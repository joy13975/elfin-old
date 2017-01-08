#!/usr/bin/env python

import Bio.PDB
import json
import numpy
from random import randint
from pymol import cmd

elfinDir = "/Users/joy/src/elfin/"
import imp
utils = imp.load_source("utils", elfinDir + "/utils.py")

cmd.reinitialize()

def main():
    dbFile      = elfinDir + "res/xDB.json"
    outDir      = elfinDir + "res/bmarks/"
    bg          = BenchmarkGenerator(dbFile, outDir)
    nBmarks     = 1
    bmarkLen    = 20
    bg.run(nBmarks, bmarkLen)
    # utils.safeExec(bg.run, nBmarks, bmarkLen)

class BenchmarkGenerator:
    def __init__(self, dbFile, outDir):
        with open(dbFile, "r") as openFile:
            db = json.load(openFile)

        for k, v in db.iteritems():
            setattr(self, k, v)

        self.nonTerms = []
        for k, v in self.stat.iteritems():
            if(v > 1):
                self.nonTerms.append(k)

        print("DB has {} non-terminal nodes".format(len(self.nonTerms)))

        self.outDir = outDir
        self.bmarks = []

    def chooseNextNode(self, lastNode):
        nextId = randint(0, len(self.data[lastNode]) - 1)
        return self.data[lastNode].keys()[nextId]

    def gen(self, bmarkLen):
        nodes = []

        # Step 1: Pick starting single from non-terminal nodes
        nNonTerms = len(self.nonTerms)
        nodes.append(self.nonTerms[randint(0, nNonTerms-1)])

        # Shape starts from origin
        shape = numpy.zeros(shape=(1,3), dtype="float64")

        # Keep adding a next node from any node until either
        #   specified length is reached

        cmd.hide("everything", "all")
        for i in xrange(0, bmarkLen-1):
            lastNode = nodes[-1]
            currNode = self.chooseNextNode(lastNode)
            nodes.append(currNode)

            rel = self.data[lastNode][currNode]
            shape = numpy.append(shape, [rel["comB"]], axis=0)

            # pymol load new pair
            pairName = nodes[-2] + "-" + nodes[-1]
            cmd.load(elfinDir + "res/centered_pdb/pair/" + pairName + ".pdb",
                str(i) + "-" + pairName)

            # pymol rotate
            pmTrans = numpy.append(
                            numpy.append(
                                numpy.transpose(rel["rot"]),
                                numpy.transpose([rel["tran"]]),
                                axis=1),
                            [[0.0, 0.0, 0.0, 1.0]],
                            axis=0)

            # utils.interact(globals(), locals())
            cmd.transform_selection("all", pmTrans.flatten().tolist())

            shape = numpy.dot(shape, numpy.asarray(rel["rot"])) + rel["tran"]
            print "Pair: " + pairName

        cmd.hide("everything", "all")
        cmd.show("cartoon", "all")
        cmd.reset()

        self.bmarks.append(shape)

    def run(self, nBmarks, bmarkLen):
        for i in xrange(0, nBmarks):
            print("Genereating #{}/{}".format(i+1, nBmarks))
            self.gen(bmarkLen)

main()
