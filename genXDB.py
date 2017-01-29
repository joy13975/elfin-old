#!/usr/bin/env python

import Bio.PDB
from utils import *
import glob
import numpy
import codecs, json
from collections import OrderedDict

def main():
    pairDir     = 'res/pair/'
    singleDir   = 'res/single/'
    newLibDir   = 'res/centered_pdb/'
    outFile     = 'res/xDB.json'
    xdbg = XDBGenrator(pairDir, singleDir, newLibDir, outFile)
    safeExec(xdbg.run)

class XDBGenrator:

    def __init__(self,
                pairDir,
                singleDir,
                newLibDir,
                outFile,
                permissive=0):
        self.pairDir        = pairDir
        self.singleDir      = singleDir
        mkdir(newLibDir)
        mkdir(newLibDir     + '/pair/')
        mkdir(newLibDir     + '/single/')
        self.newLibDir      = newLibDir
        self.outFile        = outFile
        self.si             = Bio.PDB.Superimposer()
        self.pairsData      = {}
        self.singlesData    = {}

    def getCOM(self, struct):
        CAs = [];

        for a in struct.get_atoms():
            if(a.name == 'CA'):
                # I noticed some floating point error with float32, so use double
                CAs.append(a.get_coord().astype('float64'))

        return numpy.mean(CAs, axis=0)

    def moveToOrigin(self, struct, com=[]):
        if(len(com) == 0):
            com = self.getCOM(struct)

        # No rotation - just move to centre
        struct.transform([[1,0,0],[0,1,0],[0,0,1]], -com)

    def alignToSingle(self, pStruct, sStructs, which=0):
        pSingles = utils.getPdbSingles(pStruct)

        nSingles = len(pSingles)
        nSStructs = len(sStructs)
        if(nSingles != len(sStructs)):
            raise ValueError('nSingles(' +
                str(nSingles) + ') != nSStructs(' +
                str(nSStructs) + ')')

        if(which > nSingles - 1 or which < 0):
            raise ValueError('Argument \'which\'=' +
                str(which) + ' is out of bound (0-' +
                str(nSingles) + ')')

        ma = []
        for a in pSingles[which].get_atoms(): ma.append(a)
        fa = []
        for a in sStructs[which].get_atoms(): fa.append(a)

        self.si.set_atoms(fa, ma)
        rot, tran = self.si.rotran
        pStruct.transform(rot, tran)

    def getRotTrans(self, moving, fixed):
        ma = []
        for a in moving.get_atoms(): ma.append(a)
        fa = []
        for a in fixed.get_atoms(): fa.append(a)

        self.si.set_atoms(fa, ma)
        return self.si.rotran

    def getRadii(self, pose):
        # Warning: this function assumes pose is centered!

        natoms = 0;
        rgSum = 0;
        maxCA = 0;

        nHeavy = 0;
        maxHeavy = 0;
        for a in pose.get_atoms():
            dist = numpy.linalg.norm(
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

    def processPDB(self, filename):
        # Step 0: Load pair and single structures
        pairName = filename.split('/')[-1].split('.')[0]
        pair = utils.readPdb(pairName, filename)
        singlesInPair = utils.getPdbSingles(pair)

        singleNames = pairName.split('-')
        psFilenames = [singleNames[0] + '.pdb', singleNames[1] + '.pdb']
        singles = [utils.readPdb(singleNames[0], self.singleDir + psFilenames[0]),
                    utils.readPdb(singleNames[1], self.singleDir + psFilenames[1])]

        # Step 1: Center the corresponding singles
        self.moveToOrigin(singles[0])
        self.moveToOrigin(singles[1])

        # Step 2: Move pair to align with first single
        # Note: this aligns pair by superimposing pair[0] with singles[0]
        self.alignToSingle(pair, singles)

        # Step 3: Get COM of the 2 singles inside the pair
        sComs = [self.getCOM(singlesInPair[0]), self.getCOM(singlesInPair[1])]

        # Step 4: Get measures for collision:
        #           1. Avg dist to com (gyradius aka RG)
        #           2. Max dist from CA to com
        #           3. Max dist from any heavy stom (not H) to COM
        sRads = [self.getRadii(singles[0]),
                    self.getRadii(singles[1])]

        # Step 5: Get transformation of pair to the second single
        # Note: pair is already aligned to first single so
        #       there is no need for the first transformation
        #       You can check this is true by varifying that
        #           self.getRotTrans(singlesInPair[0], singles[0])
        #       has identity rotation and zero translation. Also,
        #       sComs[0] should be at the origin.
        rot, tran = self.getRotTrans(singlesInPair[1], singles[1])

        # Step 6: Save the centred molecules once
        # Note: here the PDB format adds some slight floating point error
        utils.savePDB(singles[0], self.newLibDir + '/single/' + psFilenames[0])
        utils.savePDB(singles[1], self.newLibDir + '/single/' + psFilenames[1])
        utils.savePDB(pair, self.newLibDir + '/pair/' + pairName + '.pdb')

        # comA is aligned to centered singles[0] so should be at origin
        data = OrderedDict([
            ('comB',  sComs[1].tolist()),
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
        files = glob.glob(self.pairDir + '/*.pdb')
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