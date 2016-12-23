#!/usr/bin/env python

import Bio.PDB
import code
import os, glob
import numpy
import codecs, json
from collections import OrderedDict

def interact(localsVars):
    print "Intractive mode!"
    print
    code.interact(local=dict(globals(), **localsVars))

def main():
    pairDir     = "res/pair/"
    singleDir   = "res/single/"
    newLibDir   = "res/centered_pdb/"
    outFile     = open("res/xDB.json", "w")
    try:
        XDBGenrator(pairDir, singleDir, newLibDir, outFile).run()
    except Exception as e:
        print(e)
        interact(locals())
        raise e

def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

class XDBGenrator:

    def __init__(self,
                pairDir,
                singleDir,
                newLibDir,
                outFile,
                permissive=0):
        self.parser     = Bio.PDB.PDBParser(permissive)
        self.pairDir    = pairDir
        self.singleDir  = singleDir
        mkdir(newLibDir)
        mkdir(newLibDir + "/pair/")
        mkdir(newLibDir + "/single/")
        self.newLibDir  = newLibDir
        self.outFile    = outFile
        self.io         = Bio.PDB.PDBIO()
        self.si         = Bio.PDB.Superimposer()
        self.data       = []
        self.linkCount  = {}

    def getCOM(self, struct):
        CBpoints = []

        for a in struct.get_atoms():
            if(a.name == "CB"):
                CBpoints.append(a.get_coord().astype("float64"))

        return numpy.mean(CBpoints, axis=0)

    def moveToOrigin(self, struct, com=[]):
        if(len(com) == 0):
            com = self.getCOM(struct)

        for a in struct.get_atoms():
            a.coord = a.coord - com

    def savePDB(self, struct, filename):
        self.io.set_structure(struct)
        saveFile = self.newLibDir + "/" + filename
        self.io.save(saveFile)

    def getPair(self, pStruct):
        pSingles = []
        for pSingle in pStruct.get_chains():
            pSingles.append(pSingle)

        nPs = len(pSingles)
        if(nPs != 2):
            errStr = "Pair contains not 2 pSingles but " + str(nPs)
            raise ValueError(errStr)

        return pSingles

    def alignOneSingle(self, pStruct, sStructs, which=0):
        pSingles = self.getPair(pStruct)

        nSingles = len(pSingles)
        nSStructs = len(sStructs)
        if(nSingles != len(sStructs)):
            raise ValueError("nSingles(" +
                str(nSingles) + ") != nSStructs(" +
                str(nSStructs) + ")")

        if(which > nSingles - 1 or which < 0):
            raise ValueError("Argument \"which\"=" +
                str(which) + " is out of bound (0-" +
                str(nSingles) + ")")

        ma = []
        for a in pSingles[which].get_atoms(): ma.append(a)
        fa = []
        for a in sStructs[which].get_atoms(): fa.append(a)

        self.si.set_atoms(fa, ma)
        rot, tran = self.si.rotran
        pSingles[which].transform(rot, tran)

    def getRT(self, moving, fixed):
        ma = []
        for a in moving.get_atoms(): ma.append(a)
        fa = []
        for a in fixed.get_atoms(): fa.append(a)

        self.si.set_atoms(fa, ma)
        return self.si.rotran

    def processPDB(self, filename):
        pairName = filename.split("/")[-1].split(".")[0]
        singleNames = pairName.split("-")
        pStruct = self.parser.get_structure(pairName, filename)

        # Step 1: Center the corresponding singles and save as new
        psFilenames = [singleNames[0] + ".pdb", singleNames[1] + ".pdb"]
        sStructs = [self.parser.get_structure(singleNames[0],
                        self.singleDir + psFilenames[0]),
                    self.parser.get_structure(singleNames[1],
                        self.singleDir + psFilenames[1])]
        self.moveToOrigin(sStructs[0])
        self.moveToOrigin(sStructs[1])
        self.savePDB(sStructs[0], "/single/" + psFilenames[0])
        self.savePDB(sStructs[1], "/single/" + psFilenames[1])

        # Reload singles from saved PDBs (minimize error)
        sStructs = [self.parser.get_structure(singleNames[0],
                        self.newLibDir + "/single/" + psFilenames[0]),
                    self.parser.get_structure(singleNames[1],
                        self.newLibDir + "/single/" + psFilenames[1])]

        # Step 2: Move pair to align with first single and save as new
        self.alignOneSingle(pStruct, sStructs)
        self.savePDB(pStruct, "/pair/" + pairName + ".pdb")

        # Reload pair from saved PDB
        pStruct = self.parser.get_structure(pairName,
                    self.newLibDir + "/pair/" + pairName + ".pdb")
        pSingles = self.getPair(pStruct)

        # Step 3: Get COM of the 2 singles inside the pair
        pCOMs = [self.getCOM(pSingles[0]), self.getCOM(pSingles[1])]

        # Step 4: Get transformation of pair to the second single
        # Note: pair is already aligned to first single so
        #       there is no need for the first transformation
        #       You can check this is true by varifying that
        #           self.getRT(pSingles[0], sStructs[0])
        #       has identity rotation and zero translation. Also,
        #       pCOMs[0] is almost at origin.
        rot, tran = self.getRT(pSingles[1], sStructs[1])

        data = OrderedDict([
            ("pairName", pairName),
            ("comA",     pCOMs[0].tolist()),
            ("comB",     pCOMs[1].tolist()),
            ("rot",      rot.tolist()),
            ("tran",     tran.tolist())
            ])

        self.data.append(data)
        self.linkCount[singleNames[0]] = self.linkCount.get(singleNames[0], 0) + 1;

        # interact(locals())

    def dumpJSON(self):
        toDump = OrderedDict([
            ("stat",            self.linkCount),
            ("links",           len(self.linkCount)),
            ("complexity",      self.complexity),
            ("data",            self.data)
            ])

        json.dump(toDump,
            self.outFile,
            separators=(',', ':'),
            indent=4)

    def run(self):
        files = glob.glob(self.pairDir + '/*.pdb')
        nFiles = len(files)
        for i in range(0, nFiles):
            print "[XDBG] Processing file #{}/{}: {}".format(i+1, nFiles, files[i])
            self.processPDB(files[i])

        self.complexity = 1
        for s in self.linkCount:
            self.complexity = self.complexity * self.linkCount.get(s)

        print "[XDBG] Complexity: {}".format(self.complexity)

        self.dumpJSON()
        self.outFile.close()

if  __name__ =='__main__': main()