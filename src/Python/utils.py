import inspect, sys, code, traceback
import json
import os
import Bio.PDB
import numpy as np
import csv
import re
import argparse

RadiiTypes = ['avgAll', 'maxCA', 'maxHeavy']
INF = float('inf')

def upsample(pts1, pts2):
    # Upsample the shape with fewer points
    pts1IsLonger = len(pts1) > len(pts2)
    N = max(len(pts1), len(pts2))

    # Use a proportion based algorithm because
    # we want to assume both shapes are roughly
    # the same length, but not exactly
    if pts1IsLonger:
        morePoints, fewerPoints = (np.copy(pts1), np.copy(pts2))
    else:
        morePoints, fewerPoints = (np.copy(pts2), np.copy(pts1))

    # N === len(morePoints)

    # Compute longer shape total length
    mpTotalLength = 0.0
    for i in xrange(1, N):
        mpTotalLength += np.linalg.norm(morePoints[i] - morePoints[i - 1])

    fpTotalLength = 0.0
    for i in xrange(1, len(fewerPoints)):
        fpTotalLength += np.linalg.norm(fewerPoints[i] - fewerPoints[i - 1])

    # Upsample fewerPoints
    upsampled = np.empty([0, 3])

    # First and last points are the same
    upsampled = np.append(upsampled, [fewerPoints[0]], axis=0)

    mpProportion = 0.0
    fpProportion = 0.0
    mpi = 1
    for i in xrange(1, len(fewerPoints)):
        baseFpPoint = fewerPoints[i - 1]
        nextFpPoint = fewerPoints[i]
        baseFpProportion = fpProportion
        fpSegment = np.linalg.norm(nextFpPoint - baseFpPoint) / fpTotalLength
        vec = nextFpPoint - baseFpPoint

        fpProportion += fpSegment
        while (mpProportion <= fpProportion and mpi < N):
            mpSegment = \
                np.linalg.norm(morePoints[mpi] - morePoints[mpi - 1]) \
                / mpTotalLength

            if (mpProportion + mpSegment) > fpProportion:
                break
            mpProportion += mpSegment

            s = (mpProportion - baseFpProportion) / fpSegment
            upsampled = np.append(upsampled, [baseFpPoint + (vec * s)], axis=0)

            mpi += 1

    # Sometimes the last node is automatically added
    if len(upsampled) < N:
        upsampled = np.append(upsampled, [fewerPoints[-1]], axis=0)

    if pts1IsLonger:
        pts2 = upsampled 
    else: 
        pts1 = upsampled

    return pts1, pts2

def readCSVPoints(csvFile):
    pts = []
    
    with open(csvFile, 'r') as file:
        pts = np.asarray([[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n') if len(l) > 0])
    
    return pts

def saveCSV(npArray, saveFile, delimiter=' '):
    with open(saveFile, 'wb') as csvFile:
        wt = csv.writer(csvFile, delimiter=delimiter)
        for row in npArray:
            wt.writerow(row)

def canConvertToFloat(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

# Credits to http://stackoverflow.com/questions/2597278/python-load-variables-in-a-dict-into-namespace
class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

def floatApproximates(a, b, error=1e-6):
    return abs(a-b) < error

def realPath(path):
    return os.path.realpath(path)

def normPath(path):
    return os.path.normpath(path)

def suffixPdb(className, fromFunction, scale, targetLen):
    return '_{}_{}_s{}_l{}.pdb'.format(
            className, 
            fromFunction,
            scale, 
            targetLen)

def minDistFromLine(point, linePoints, allowPerp=True):
    minDist = float('inf')
    for i in xrange(1, len(linePoints)):
        lineSeg = (linePoints[i-1], linePoints[i])

        # First determine whether point is outside line segment regime
        v = lineSeg[1] - lineSeg[0]
        w = point - lineSeg[0]

        if allowPerp:
            c1 = np.dot(w, v)
            if c1 <= 0: # before lineSeg[0]
                dist = np.linalg.norm(w)
            else:
                c2 = np.dot(v, v)
                if c2 <= c1: # after lineSeg[1]
                    dist = np.linalg.norm(point - lineSeg[1])
                else:
                    # If not outside, then calculate perpendicular distance
                    b = c1 / c2
                    pol = lineSeg[0] + b*v
                    dist = np.linalg.norm(point - pol)
        else:
            dist = min(np.linalg.norm(point - lineSeg[0]), np.linalg.norm(point - lineSeg[1]))

        if dist < minDist:
            minDist = dist

    return minDist

def die(condition, str):
    if condition:
        print str
        exit()

def checkCollision(xdb, collisionMeasure, nodes, newNode, shape):
    newCOM = xdb['pairsData'][nodes[-1]][newNode]['comB']

    # previous node PAIR (not just single node!) is inherently non-colliding
    for i in xrange(0, len(nodes) - 2):
        comDist = np.linalg.norm(shape[i] - newCOM)
        collisionDist = xdb['singlesData'][newNode]['radii'][collisionMeasure] + \
                            xdb['singlesData'][nodes[i]]['radii'][collisionMeasure]

        if comDist < collisionDist:
            return True

    return False

def makePdbFromNodes(xdb, nodes, pairsDir, saveFile=None, fRot=None, movieMode=False):

    # Load first PDB and clear its chains to host all other residues
    pairName = nodes[0] + '-' + nodes[1]
    pdbFile = pairsDir + '/' + pairName + '.pdb'
    motherPdb = readPdb(pairName, pdbFile)
    motherModel = motherPdb.child_list[0]
    motherModel.detach_child('A')
    motherModel.detach_child('B')
    motherChain = Bio.PDB.Chain.Chain('A')
    motherModel.add(motherChain)

    moviePdbs = []
    baseRId = 1

    comShape = np.empty([1, 3])
    startingPoint = np.zeros(3)

    chainLenDigits = len(str(len(nodes)))
    for i in xrange(1, len(nodes)):
        lastNode = nodes[i-1]
        currNode = nodes[i]
        # pymol load new pair
        pairName = lastNode + '-' + currNode
        rel = xdb['pairsData'][lastNode][currNode]
        # mother pdb append and transform
        pdbFile = pairsDir + '/' + pairName + '.pdb'
        pdbPair = readPdb(pairName, pdbFile)
        
        comShape = np.append(comShape, [[0,0,0]], axis=0)
        if i == len(nodes) - 1:
            comShape = np.append(comShape, [rel['comB']], axis=0)
        
        if movieMode:
            moviePdbs.append(pdbPair)
            for pdb in moviePdbs:
            	pdb.transform(np.asarray(rel['rot']), rel['tran'])
        else:
            childChain = pdbPair.child_list[0].child_dict['A']
            nResi = len(childChain.child_list)
            for j in xrange(0, nResi):
                r = childChain.child_list[j]
                nextId = baseRId + j
                r.id = (r.id[0], nextId, r.id[2]) 
                motherChain.add(r)
            baseRId += nResi
            
            if i == len(nodes) - 1:
                childChain = pdbPair.child_list[0].child_dict['B']
                nResi = len(childChain.child_list)
                for j in xrange(0, nResi):
                    r = childChain.child_list[j]
                    nextId = baseRId + j
                    r.id = (r.id[0], nextId, r.id[2]) 
                    motherChain.add(r)
                baseRId += nResi
            motherPdb.transform(np.asarray(rel['rot']), rel['tran'])


        comShape = np.dot(comShape, np.asarray(rel['rot'])) + rel['tran']

        startingPoint = np.dot(startingPoint, np.asarray(rel['rot'])) + rel['tran']
        print 'Pair[{}]:   {}---{}'.format(str(i).ljust(chainLenDigits),
            lastNode.ljust(16), currNode.rjust(16))

       
    if fRot is not None:
        if movieMode:
            motherPdb.transform(np.eye(3), -startingPoint)
            motherPdb.transform(np.asarray(fRot), np.zeros(3))
        else:
            for pdb in moviePdbs:
                pdb.transform(np.eye(3), -startingPoint)
                pdb.transform(np.asarray(fRot), np.zeros(3))

    if not movieMode:
        if saveFile is not None:
            savePdb(motherPdb, saveFile)
        return motherPdb, comShape
    else:
        if saveFile is not None:
            pdbId = 0
            saveFileDotIndex = saveFile.rfind('.')
            for pdb in moviePdbs:
                savePartFile = saveFile[0:saveFileDotIndex] + \
                    'part' + str(pdbId) + \
                    saveFile[saveFileDotIndex:]
                savePdb(pdb, savePartFile)
                pdbId = pdbId + 1
        return moviePdbs, comShape

def getXDBStat(xDB):
    # xdb = readJSON('res/xDB.json')

    pd = xDB['pairsData']
    dists = []
    for s1 in pd.keys():
        for s2 in pd[s1].keys():
            dists.append(np.linalg.norm(pd[s1][s2]['comB']))

    return np.average(dists), min(dists), max(dists)

def readJSON(filename):
    with open(filename, 'r') as openFile:
        return json.load(openFile)

def readPdb(customName, inFile, permissive=0):
    parser = Bio.PDB.PDBParser(permissive)
    structure = parser.get_structure(customName, inFile)
    return structure

def savePdb(struct, saveFile):
    io = Bio.PDB.PDBIO()
    io.set_structure(struct)
    io.save(saveFile)

def getPdbSingles(pStruct):
    pSingles = []
    for pSingle in pStruct.get_chains():
        pSingles.append(pSingle)

    nPs = len(pSingles)
    assert nPs == 2

    return pSingles

def interact(globalVars=None, localsVars=None):
    print "Entering interactive mode"
    print
    if(localsVars == None):
        localsVars = locals()
    if(globalVars == None):
        globalVars = globals()
    code.interact(local=dict(globalVars, **localsVars))

def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def pauseCode(frame=None):
    print '---------pauseCode()---------'
    if frame is None:
        # Use current frame (one above the exception wrapper)
        frame = inspect.currentframe().f_back
    
    ns = dict(frame.f_globals)
    ns.update(frame.f_locals)
    code.interact(local=ns)

def safeExec(func, *args):
    try:
        func(*args)
    except Exception as e:
        print '---------safeExec() caught exception---------'

        # Find last (failed) inner frame
        type, value, tb = sys.exc_info()
        last_frame = lambda tb=tb: last_frame(tb.tb_next) if tb.tb_next else tb
        frame = last_frame().tb_frame
        traceback.print_exc()
        pauseCode(frame)
