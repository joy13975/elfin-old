#!/usr/bin/env python

import Bio.PDB
import utils
import glob
import numpy as np
import codecs, json
from collections import OrderedDict

def main():
    pairDir         = 'res/mergedAndCleansed/pair/'
    singleDir       = 'res/mergedAndCleansed/single/'
    alignedLibDir   = 'res/aligned/'
    outFile         = 'res/xDB.json'
    xdbg = XDBGenrator(pairDir, singleDir, alignedLibDir, outFile)
    utils.safeExec(xdbg.run)

class XDBGenrator:

    def __init__(self,
                pairDir,
                singleDir,
                alignedLibDir,
                outFile,
                permissive=0):
        self.pairDir        = pairDir
        self.singleDir      = singleDir
        utils.mkdir(alignedLibDir)
        utils.mkdir(alignedLibDir     + '/pair/')
        utils.mkdir(alignedLibDir     + '/single/')
        self.alignedLibDir  = alignedLibDir
        self.outFile        = outFile
        self.si             = Bio.PDB.Superimposer()
        self.pairsData      = {}
        self.singlesData    = {}

    def getCOM(self, child, mother=None, childAtomOffset=0, motherAtomOffset=0):
        CAs = []
        for a in child.get_atoms():
            if(a.name == 'CA'):
                CAs.append(a.get_coord().astype('float64'))
        com = np.mean(CAs, axis=0)

        if mother is not None:
            # This is for finding COM of a single inside a pair
            _, tran = self.getRotTrans(child, mother, 
                movingAtomOffset=childAtomOffset, fixedAtomOffset=motherAtomOffset)

            com += tran

        return com

    def moveToOrigin(self, pdb):
        com = self.getCOM(pdb)

        # No rotation - just move to centre
        pdb.transform([[1,0,0],[0,1,0],[0,0,1]], -com)

    def align(self, moving, fixed,movingAtomOffset=0, fixedAtomOffset=0):
        rot, tran = self.getRotTrans(moving, fixed,
            movingAtomOffset=movingAtomOffset, fixedAtomOffset=fixedAtomOffset)
        moving.transform(rot, tran)

    def getRotTrans(self, moving, fixed, movingAtomOffset=0, fixedAtomOffset=0):
        # First push the generators to desired locations
        maGen = moving.get_atoms()
        for i in xrange(0, movingAtomOffset):
            try:
                maGen.next()
            except StopIteration:
                die(True, 'Moving PDB Atom Offset too large')    

        faGen = fixed.get_atoms()
        for i in xrange(0, fixedAtomOffset):
            try:
                faGen.next()
            except StopIteration:
                die(True, 'Fixed PDB Atom Offset too large')    

        # Then fill in the arrays until either of the generators run out
        ma = []
        fa = []
        while True:
            try:
                maNext = maGen.next()
                faNext = faGen.next()
            except StopIteration:
                break

            ma.append(maNext)
            fa.append(faNext)

        self.si.set_atoms(fa, ma)

        # Import note:
        # The rotation from BioPython seems be the
        # second dot operand instead of the 
        # conventional first dot operand!
        #
        # This means instead of R*v + T, the actual
        # transform is done with v'*R + T
        #
        # This has important to understand why I did
        # the rotation maths this way in the C++ GA
        return self.si.rotran

    def getRadii(self, pose):
        # Warning: this function assumes pose is centered!

        natoms = 0;
        rgSum = 0;
        maxCA = 0;

        nHeavy = 0;
        maxHeavy = 0;
        for a in pose.get_atoms():
            dist = np.linalg.norm(
                a.get_coord().astype('float64'));

            rgSum += dist;

            if(a.name =='CA'):
                maxCA = max(maxCA, dist);

            if(a.element != 'H'):
                maxHeavy = max(maxHeavy, dist);
                nHeavy = nHeavy + 1;

            natoms = natoms + 1;

        rg = rgSum / natoms;
        return OrderedDict([
            ('avgAll', rg),
            ('maxCA', maxCA),
            ('maxHeavy', maxHeavy)
        ]);

    def getAtomCount(self, pdb):
        i = 0
        for a in pdb.get_atoms():
            i += 1
        return i

    def processPDB(self, filename):
        # Step 0: Load pair and single structures
        pairName = filename.split('/')[-1].split('.')[0].replace('_mc_0001', '')
        pairPdb = utils.readPdb(pairName, filename)

        singleNames = pairName.split('-')
        psFilenames = [
            singleNames[0] + '_mc_0001.pdb', 
            singleNames[1] + '_mc_0001.pdb'
        ]
        singlePdbs = [utils.readPdb(singleNames[0], self.singleDir + psFilenames[0]),
                    utils.readPdb(singleNames[1], self.singleDir + psFilenames[1])]

        # Step 1: Center the corresponding singles
        self.moveToOrigin(singlePdbs[0])
        self.moveToOrigin(singlePdbs[1])

        # Step 2: Move pair to align with first single
        # Note: this aligns pair by superimposing pair[0] with singlePdbs[0]
        self.align(pairPdb, singlePdbs[0])

        # Step 3: Get COM of the second single inside the pair
        singlePdbAtomCounts = [
            self.getAtomCount(singlePdbs[0]), 
            self.getAtomCount(singlePdbs[1])
        ]
        com2 = self.getCOM(singlePdbs[1], pairPdb, motherAtomOffset=singlePdbAtomCounts[0])

        # Step 4: Get measures for collision:
        #           1. Avg dist to com (gyradius aka RG)
        #           2. Max dist from CA to com
        #           3. Max dist from any heavy stom (not H) to COM
        sRads = [self.getRadii(singlePdbs[0]),
                    self.getRadii(singlePdbs[1])]

        # Step 5: Get transformation of pair to the second single
        # Note: pair is already aligned to first single so
        #       there is no need for the first transformation
        #       You can check this is true by varifying that
        #           self.getRotTrans(pairPdb, singlePdbs[0])
        #       has identity rotation and zero translation. Also,
        #       com2 should be at the origin.
        rot, tran = self.getRotTrans(pairPdb, singlePdbs[1], 
            movingAtomOffset=singlePdbAtomCounts[0])

        # Step 6: Save the aligned molecules once
        # Note: here the PDB format adds some slight floating point error
        utils.savePdb(singlePdbs[0], self.alignedLibDir + '/single/' + singleNames[0])
        utils.savePdb(singlePdbs[1], self.alignedLibDir + '/single/' + singleNames[1])
        utils.savePdb(pairPdb, self.alignedLibDir + '/pair/' + pairName + '.pdb')

        # comA is aligned to centered singlePdbs[0] so should be at origin
        data = OrderedDict([
            ('comB',  com2.tolist()),
            ('rot',   rot.tolist()),
            ('tran',  tran.tolist())
            ])

        entry = self.pairsData.get(singleNames[0], None)
        if(entry == None):
            self.pairsData[singleNames[0]] = {}
            entry = self.pairsData.get(singleNames[0], None)

        entry[singleNames[1]] = data

        singleDataA = self.singlesData.get(singleNames[0],
             OrderedDict([
                ('linkCount', 0),
                ('radii', sRads[0])
                ]));
        singleDataA['linkCount'] = singleDataA['linkCount'] + 1;
        self.singlesData[singleNames[0]] = singleDataA;

        if(singleNames[1] != singleNames[0]):
            singleDataB = self.singlesData.get(singleNames[1],
                 OrderedDict([
                    ('linkCount', 0),
                    ('radii', sRads[0])
                    ]));
            self.singlesData[singleNames[1]] = singleDataB;

        # interact(globals(), locals())

    def dumpJSON(self):
        toDump = OrderedDict([
            ('singlesData',  self.singlesData),
            ('complexity',  self.complexity),
            ('pairsData',   self.pairsData)
            ])

        json.dump(toDump,
            open(self.outFile, 'w'),
            separators=(',', ':'),
            ensure_ascii=False,
            indent=4)

    def run(self):
        # _mc stands for (chain) Merged and Cleansed
        # _0001 means it is minimized by Rosetta
        files = glob.glob(self.pairDir + '/*_mc_0001.pdb')
        nFiles = len(files)
        for i in range(0, nFiles):
            print '[XDBG] Processing file #{}/{}: {}'.format(i+1, nFiles, files[i])
            self.processPDB(files[i])

        self.complexity = 1
        for s in self.singlesData:
            self.complexity = self.complexity * self.singlesData.get(s)['linkCount']

        print '[XDBG] Complexity: {}'.format(self.complexity)

        self.dumpJSON()

if __name__ =='__main__': main()