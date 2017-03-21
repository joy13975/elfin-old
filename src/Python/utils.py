import inspect, sys, code, traceback
import json
import os
import Bio.PDB
import numpy as np

RadiiTypes = ['avgAll', 'maxCA', 'maxHeavy']
INF = float('inf')

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
        comDist = np.linalg.norm(shape[i] - newCOM);
        collisionDist = xdb['singlesData'][newNode]['radii'][collisionMeasure] + \
                            xdb['singlesData'][nodes[i]]['radii'][collisionMeasure];

        if comDist < collisionDist:
            return True;

    return False;

def makePdbFromNodes(xdb, nodes, pairsDir, saveFile=None, fRot=None, movieMode=False):    
    # Load first "mother" pdb to host all subsequent chains
    pairName = nodes[0] + '-' + nodes[1]
    pdbFile = pairsDir + '/' + pairName + '.pdb'
    motherPdb = readPdb(pairName, pdbFile)

    motherChain = motherPdb.child_list[0].child_dict['A']
    baseRId = motherChain.child_list[-1].id[1]

    # Only for mother pdb, keep both chains
    motherChainB = motherPdb.child_list[0].child_dict['B']
    motherModel = motherPdb.child_list[0]
    motherModel.detach_child('B')
    if movieMode:
        movieChainID = 'B'
    else:
        nResi = len(motherChainB.child_list)
        for j in xrange(0, nResi):
            r = motherChainB.child_list[j]
            r.id = (r.id[0], 
                    (baseRId + 1 + j),
                    r.id[2]) 
            motherChain.add(r)
        baseRId += nResi

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

        
        childChain = pdbPair.child_list[0].child_dict['B']
        nResi = len(childChain.child_list)
        if movieMode:
            childChain.id = movieChainID
            motherModel.add(childChain)

            movieChainID = chr(ord(movieChainID) + 1)
            # pauseCode()
        else:
            for j in xrange(0, nResi):
                r = childChain.child_list[j]
                r.id = (r.id[0], 
                        (baseRId + 1 + j),
                        r.id[2]) 
                motherChain.add(r)

        baseRId += nResi
        motherPdb.transform(np.asarray(rel['rot']), rel['tran'])
        startingPoint = np.dot(startingPoint, np.asarray(rel['rot'])) + rel['tran']
        print 'Pair[{}]:   {}---{}'.format(str(i).ljust(chainLenDigits),
            lastNode.ljust(16), currNode.rjust(16))

       
    if fRot is not None:
        motherPdb.transform(np.eye(3), -startingPoint)
        motherPdb.transform(np.asarray(fRot), np.zeros(3))

    if saveFile is None:
        return motherPdb
    else:
        savePdb(motherPdb, saveFile)

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