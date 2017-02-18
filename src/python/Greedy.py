#!/usr/bin/env python

import glob, random
from utils import *
from Designer import *
import Kabsch
import time
import re
import numpy as np

def main():
    if(len(sys.argv) < 2):
        print './Greedy.py <specFile.{json|pdb}> s=scale l=len'
        exit()

    specFile = sys.argv[1]
    fileExt = specFile[specFile.rfind('.'):]

    if fileExt == '.json':
        spec = readJSON(specFile)
        targetLen = len(spec['nodes'])
    elif fileExt == '.csv':
        with open(specFile, 'r') as file:
            pts = [[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n')]

        spec = {'coms': pts}
        npts = np.asarray(pts)

        (avgD, minD, maxD) = getXDBStat()
        # Use total length/avgD as heuristic. avgD is average xDB pair distance
        targetLen = int(np.ceil(sum(np.linalg.norm(npts-np.roll(npts, 1, axis=0), axis=1)) / avgD))
    else:
        print 'Unknown spec file type: {}'.format(fileExt)
        exit()


    # default options
    scale = 1.0
    userLen = -1
    for i in xrange(2, len(sys.argv)):
        arg = sys.argv[i]
        if arg.startswith('s='):
            scale = float(arg[(arg.rfind('=')+1):])
        elif arg.startswith('l='):
            userLen = int(arg[(arg.rfind('=')+1):])
        else:
            print 'Unrecognised flag: \"{}\"'.format(arg)
            exit()

    die(scale < 0, 'Scale must be > 0.0')

    # User specified length should overwrite 
    if userLen != -1:
        die(userLen < 3, 'Length must be > 3')
        targetLen = userLen
    else:
        targetLen = int(round(targetLen * scale))

    xDB         = readJSON('res/xDB.json')

    designer    = GreedyDesigner(xDB, 'maxHeavy')
    startTime   = time.clock()
    (nodes,shape,score, fRot) = designer.design(spec, targetLen, scaleFactor=scale)
    print "{:.2f}s, score: {}".format(time.clock() - startTime, score)

    makePdbFromNodes(xDB, nodes, 'res/centered_pdb/pair',
        specFile.replace(fileExt, suffixPdb(
            designer.__class__.__name__, 
            'Main',
            scale, 
            targetLen)),
        fRot)
            
class GreedyDesigner(Designer):

    def __init__(self, xDB, collisionMeasure):
    	self.xDB                = xDB
        self.collisionMeasure   = collisionMeasure
        assert collisionMeasure in RadiiTypes

    def grow(self, lastNode, shape, newNode):
        rel = self.xDB['pairsData'][lastNode][newNode]
        shape = np.dot(shape, np.asarray(rel['rot'])) + rel['tran']
        shape = np.append(shape, [[0,0,0]], axis=0)

        return shape

    def upsampleTarget(self, target):
        # Distances avg: 38.0111371052, min: 17.6585177454, max: 50.4035189874
        avgD = 38.0111371052
        minD = 17.6585177454
        maxD = 50.4035189874

        result = []
        for i in xrange(1, len(target)):
            result.append(target[i-1])

            delta = target[i]-target[i-1]
            dist = np.linalg.norm(delta)

            if dist < maxD:
                continue

            # Use min or avg!?
            n = int(np.floor(dist / avgD))

            partDelta = delta / n
            for j in xrange(1, n):
                result.append(target[i-1] + (j * partDelta))

        result.append(target[-1])
        return np.asarray(result)

    def evalShape(self, candidate, target, allowPerp=True):
        sumScore = 0
        for point in candidate:
            sumScore += minDistFromLine(point, target, allowPerp=allowPerp)
        return sumScore

    def design(self, spec, targetLen, scaleFactor=1.0):
        print 'Greedy design: target length {}, scale {}'.format(
            targetLen, scaleFactor)
        minTotalScore = float('inf')
        minScoreShape = None
        minScoreNodes = None

        # Load target and shift it by the first CoM 
        # so we can calculate score progressively
        target = spec.get('coms', None)
        assert target is not None
        target = np.asarray(target)
        target = target - target[0]
        target *= scaleFactor

        target = self.upsampleTarget(target)

        assert len(target) >= 2

        # Try all starting pairs - we don't have information
        # to know which pair is best to start with
        pd = self.xDB['pairsData']
        candidate = np.asarray([[0, 0, 0]])

        bestStartPair = None
        bestStartScore = float('inf')
        for s1 in pd.keys():
            for s2 in pd[s1].keys(): 
                rel = pd[s1][s2]
                tmpShape = np.asarray([[0, 0, 0], rel['comB']])
                kR = Kabsch.kabsch(tmpShape, target[:2])
                tmpShape = np.dot(tmpShape, kR)
                s = self.evalShape(tmpShape, target, allowPerp=False)
                if s < bestStartScore:
                    bestStartPair = (s1, s2)
                    bestStartScore = s

        assert bestStartPair is not None
        candidate = self.grow(bestStartPair[0], [[0, 0, 0]], bestStartPair[1])
        nodes = [bestStartPair[0], bestStartPair[1]]

        # Main greedy construction loop
        totalScore = 0.0
        for i in xrange(2, targetLen):
            s1 = nodes[-1]

            bestNextNode = None
            bestNextScore = float('inf')
            atLeastOneNonColliding = False
            for s2 in pd[s1].keys(): 
                if checkCollision(self.xDB, self.collisionMeasure, nodes, s2, candidate):
                    continue

                atLeastOneNonColliding = True

                tmpShape = self.grow(s1, np.copy(candidate), s2)
                tmpShape = tmpShape - tmpShape[0]

                matchLen = min(len(tmpShape), len(target))
                kR = Kabsch.kabsch(tmpShape[:matchLen], target[:matchLen])
                tmpShape = np.dot(tmpShape, kR)

                s = self.evalShape(tmpShape, target, allowPerp=True)
                if s < bestNextScore:
                    bestNextNode = s2
                    bestNextScore = s

            if not atLeastOneNonColliding:
                print 'All next nodes collide! Cannot continue... Current length: {}'.format(i+1)
                print nodes
                print 'Available next nodes: {}'.format(pd[s2].keys())
                break

            assert bestNextNode is not None
            
            nodes.append(bestNextNode)
            candidate = self.grow(s1, candidate, bestNextNode)
            totalScore += bestNextScore

        candidate - candidate[0]
        matchLen = min(len(candidate), len(target))
        kR = Kabsch.kabsch(candidate[:matchLen]-candidate[0], target[:matchLen])
        candidate = np.dot(candidate, kR)
        return nodes, candidate, totalScore, kR

if __name__ =='__main__': safeExec(main)