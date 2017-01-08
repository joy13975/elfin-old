#!/usr/bin/env python

import Bio.PDB
import json
import numpy
from random import randint

try:
    from pymol import cmd
    cmd.reinitialize();
    def _cmd(funcName, *args):
        func = getattr(cmd, funcName);
        return func(*args)
except ImportError:
    print 'Could not import pymol cmd. Not running as pymol plugin...'
    def _cmd(funcName, *args):
        pass

elfinDir = '/Users/joy/src/elfin/'
import imp
utils = imp.load_source('utils', elfinDir + '/utils.py')

def main():
    dbFile      = elfinDir + 'res/xDB.json'
    outDir      = elfinDir + 'res/bmarks/'
    bg          = BenchmarkGenerator(dbFile, outDir, 'avgAll')
    nBmarks     = 1
    chainLen    = 20
    # bg.run(nBmarks, chainLen)
    utils.safeExec(bg.run, nBmarks, chainLen)

class BenchmarkGenerator:
    RadiiTypes = ['avgAll', 'maxCA', 'maxHeavy'];
    def __init__(self, dbFile, outDir, collisionMeasure):
        def makeSelf():
            with open(dbFile, 'r') as openFile:
                db = json.load(openFile)

            for k, v in db.iteritems():
                setattr(self, k, v)

            self.nonTerms = []
            for k, v in self.singlesData.iteritems():
                if(v['linkCount'] > 1):
                    self.nonTerms.append(k)

            print('DB has {} non-terminal nodes'.format(len(self.nonTerms)))

            self.outDir = outDir
            self.bmarks = []

            # Collision measure is the radius type used to check collision
            if(collisionMeasure not in self.RadiiTypes):
                raise ValueError('Possible values for collisionMeasure: {}'.format(self.RadiiTypes));

            self.collisionMeasure = collisionMeasure

        utils.safeExec(makeSelf)

    def checkCollision(self, newNode, nodes, shape):
        newCOM = self.pairsData[nodes[-1]][newNode]['comB']

        # previous node PAIR (not just single node!) is inherently non-colliding
        for i in xrange(0, len(nodes) - 2):
            comDist = numpy.linalg.norm(shape[i] - newCOM);
            collisionDist = self.singlesData[newNode]['radii'][self.collisionMeasure] + \
                                self.singlesData[nodes[i]]['radii'][self.collisionMeasure];

            if comDist < collisionDist:
                return True;

        return False;

    def chooseNextNode(self, nodes, shape):
        lastNode = nodes[-1]
        links = self.pairsData[lastNode].keys();

        # utils.interact(globals(), locals())
        collide = True
        while(collide):
            newNodeId = randint(0, len(links) - 1)
            newNode = links[newNodeId]
            collide = self.checkCollision(newNode, nodes, shape)
            if collide:
                links.remove(newNode)

            if len(links) == 0:
                print 'All possible links cause collision...'
                print 'Possible links: {}'.format(self.pairsData[lastNode].keys())
                return None;

        return newNode

    def gen(self, chainLen):
        nodes = []

        # Step 1: Pick starting single from non-terminal nodes
        nNonTerms = len(self.nonTerms)
        nodes.append(self.nonTerms[randint(0, nNonTerms-1)])

        # Shape starts from origin
        shape = numpy.zeros(shape=(1,3), dtype='float64')

        # Keep adding a next node from any node until either
        #   specified length is reached

        _cmd('hide', 'everything', 'all')
        for i in xrange(0, chainLen-1):
            lastNode = nodes[-1]
            newNode = self.chooseNextNode(nodes, shape)
            if(newNode is None):
                print('Collision stopped chain generation')
                break

            nodes.append(newNode)

            rel = self.pairsData[lastNode][newNode]
            shape = numpy.append(shape, [rel['comB']], axis=0)

            # pymol load new pair
            pairName = nodes[-2] + '-' + nodes[-1]
            _cmd('load', elfinDir + 'res/centered_pdb/pair/' + pairName + '.pdb',
                str(i) + '-' + pairName)

            # pymol rotate
            pmTrans = numpy.append(
                            numpy.append(
                                numpy.transpose(rel['rot']),
                                numpy.transpose([rel['tran']]),
                                axis=1),
                            [[0.0, 0.0, 0.0, 1.0]],
                            axis=0)

            # utils.interact(globals(), locals())
            _cmd('transform_selection', 'all', pmTrans.flatten().tolist())

            shape = numpy.dot(shape, numpy.asarray(rel['rot'])) + rel['tran']
            print 'Pair: ' + pairName

        _cmd('hide', 'everything', 'all')
        _cmd('show', 'cartoon', 'all')
        _cmd('reset')

        self.bmarks.append(shape)

    def run(self, nBmarks, chainLen):
        for i in xrange(0, nBmarks):
            print('Genereating #{}/{}'.format(i+1, nBmarks))
            self.gen(chainLen)

main()
