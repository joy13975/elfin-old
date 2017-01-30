#!/usr/bin/env python

import glob, numpy, random
from utils import *
from Designer import *
import Kabsch
import shapely.geometry as geom

def main():
    xDB     = 'res/xDB.json'
    bmDir 	= 'bm/'

    gd = GreedyDesigner(xDB)

    # Process all benchmarks
    files = glob.glob(bmDir + '/*.json')
    nFiles = len(files)
    for i in range(0, nFiles):
    	file = files[i]
    	print 'Benchmarking GreedyDesigner on {}'.format(file)
    	gd.design(file)

class GreedyDesigner(Designer):

    def __init__(self, xDB):
    	self.xDB = xDB

    def parseTarget(self, spec):
        target = spec.get('coms', None)
        if target is not None:
            return numpy.asarray(target), len(target)

        target = spec.get('spline', None)
        if target is not None:
            raise NotImplementedError('spline determine length?')
            return numpy.asarray(target), -1

        target = spec.get('segments', None)
        if target is not None:
            raise NotImplementedError('segments determine length?')
            return numpy.asarray(target), -1

        assert target != None

    def getMostLikely(self, shapePt, target):
        s2Dist = numpy.linalg.norm(shapePt)
        tDist = numpy.linalg.norm(target, axis=1)
        return numpy.argmin(abs(tDist - s2Dist))

    def grow(self, lastNode, shape, newNode):
        rel = self.xDB['pairsData'][lastNode][newNode]
        shape = numpy.dot(shape, numpy.asarray(rel['rot'])) + rel['tran']
        shape = numpy.append(shape, [[0,0,0]], axis=0)

        return shape

    def design(self, spec):
        minTotalScore = float('inf')
        minScoreShape = None
        minScoreNodes = None

        # Shift target shape by its first CoM so we can calculate
        # score progressively
        target, targetLen = self.parseTarget(spec)
        target = target - target[0]

        assert len(target) >= 2

        # Try all starting pairs - we don't have information
        # to know which starting pair is "best" for greedy
        pd = self.xDB['pairsData']
        for s1 in pd.keys():
            for s2 in pd[s1].keys():
                nodes = [s1, s2]
                shape = self.grow(s1, [[0,0,0]], s2)

                # Find second point estimate
                mostLikelyId1 = self.getMostLikely(shape[1]-shape[0], target)

                # Find Third point  estimate
                minScore3 = float('inf')
                mostLikelyId2 = None
                s3 = None
                for s in pd[s2].keys():
                    tmpShape = self.grow(s2, numpy.copy(shape), s)
                    tmpShape = tmpShape - tmpShape[0]

                    # Get a 3D orientation estimation
                    id2 = self.getMostLikely(tmpShape[2]-tmpShape[0], target)
                    kR = Kabsch.kabsch(tmpShape[0:3], 
                        [   
                            target[0], 
                            target[mostLikelyId1], 
                            target[id2]
                        ])
                    sc = numpy.linalg.norm(numpy.dot(tmpShape, kR) - target[0:3])

                    if(sc < minScore3):
                        minScore3 = sc
                        mostLikelyId2 = id2
                        s3 = s

                assert s3 is not None and mostLikelyId2 is not None

                nodes.append(s3)
                shape = self.grow(s2, shape, s3)
                targetFixTrio = [target[0], target[mostLikelyId1], target[mostLikelyId2]]
                line = geom.LineString(target)

                # Main construction loop
                totalScore = 0
                for i in range(3, targetLen):
                    minScore = float('inf')
                    minScoreNode = None
                    lastNode = nodes[-1]
                    for s in pd[lastNode].keys():
                        tmpShape = self.grow(lastNode, numpy.copy(shape), s)
                        tmpShape = tmpShape - tmpShape[0]
    
                        kR = Kabsch.kabsch(tmpShape[0:3], targetFixTrio)
                        tmpShape = numpy.dot(tmpShape, kR)
                        point = geom.Point(tmpShape[-1])
                        sc = point.distance(line)

                        if(sc < minScore):
                            minScore = sc
                            minScoreNode = s

                    nodes.append(minScoreNode)
                    shape = self.grow(lastNode, shape, minScoreNode)
                    totalScore = totalScore + minScore

                if(totalScore < minTotalScore):
                    minTotalScore = totalScore
                    minScoreShape = shape
                    minScoreNodes = nodes

        assert minTotalScore != float('inf') and minScoreShape is not None and minScoreNodes is not None 
        return minScoreNodes, minScoreShape, minTotalScore

if __name__ =='__main__': safeExec(main)