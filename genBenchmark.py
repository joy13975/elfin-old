#!/usr/bin/env python

import Bio.PDB, json
import numpy, random, string, math, codecs
from collections import OrderedDict
from time import gmtime, strftime
import sys

haveCmd = False
try:
    from pymol import cmd
    cmd.reinitialize()
    haveCmd = True
except ImportError:
    print 'Could not import pymol cmd. Not running as pymol plugin...'

elfinDir = '/Users/joy/src/elfin/'
import imp
utils = imp.load_source('utils', elfinDir + '/utils.py')

def main():
    dbFile              = elfinDir + 'res/xDB.json'
    outDir              = elfinDir + './bm/'
    pdbPairsDir         = 'res/centered_pdb/pair/'
    bg                  = BenchmarkGenerator(dbFile, pdbPairsDir, outDir, 'maxHeavy')
    nBmarks             = 1
    if(len(sys.argv) > 1):
        try:
            nBmarks = int(sys.argv[1])
        except Exception as e:
            print e
            print 'Failed to convert argument {} to integer (number of shapes to generate)'.format(sys.argv[1])
            exit(1)

    chainLen            = 20
    maxRetries          = -1
    # bg.run(nBmarks, chainLen)
    utils.safeExec(bg.run, nBmarks, chainLen, maxRetries)

class BenchmarkGenerator:
    RadiiTypes = ['avgAll', 'maxCA', 'maxHeavy'];
    def __init__(self,
                dbFile,
                pdbPairsDir,
                outDir,
                collisionMeasure):
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

            self.pdbPairsDir = pdbPairsDir
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

        collide = True
        while(collide):
            newNodeId = random.randint(0, len(links) - 1)
            newNode = links[newNodeId]
            collide = self.checkCollision(newNode, nodes, shape)
            if collide:
                links.remove(newNode)

            if len(links) == 0:
                print 'Stopped because all links lead to collision'
                print 'Available links: {}'.format(
                    [str(k) for k in self.pairsData[lastNode].keys()])
                raise UserWarning('Genereation could not continue')

        return newNode

    def gen(self, chainLen):
        nodes = []

        # Step 1: Pick starting single from non-terminal nodes
        nNonTerms = len(self.nonTerms)
        nodes.append(self.nonTerms[random.randint(0, nNonTerms-1)])

        # Shape (array of CoMs) starts from origin
        coms = numpy.zeros(shape=(1,3), dtype='float64')
        pairs = []

        # Main structure generation loop
        # Keep adding a next node from any node until either
        #   specified length is reached
        for i in xrange(0, chainLen):
            lastNode = nodes[i]
            newNode = self.chooseNextNode(nodes, coms)

            nodes.append(newNode)

            rel = self.pairsData[lastNode][newNode]
            coms = numpy.append(coms, [rel['comB']], axis=0)
            coms = numpy.dot(coms, numpy.asarray(rel['rot'])) + rel['tran']

        # Move display/print/postprocess to after construction succeeded
        # Makes generation faster

        # Load first "mother" pdb to host all subsequent chains
        pairName = nodes[0] + '-' + nodes[1]
        pdbFile = elfinDir + self.pdbPairsDir + '/' + pairName + '.pdb'
        motherPdb = utils.readPdb(pairName, pdbFile)

        # Only for mother pdb, remove B instead of A
        chainId = 'A'
        motherPdb.child_list[0].detach_child('B')

        chainLenDigits = int(math.ceil(math.log(chainLen, 10)))
        for i in xrange(1, len(nodes)):
            lastNode = nodes[i-1]
            currNode = nodes[i]

            # pymol load new pair
            pairName = lastNode + '-' + currNode
            pairs.append(pairName)

            rel = self.pairsData[lastNode][currNode]

            # mother pdb append and transform
            pdbFile = elfinDir + self.pdbPairsDir + '/' + pairName + '.pdb'
            pdbPair = utils.readPdb(pairName, pdbFile)

            chainId = chr(ord(chainId) + 1)
            pdbPair.child_list[0].child_dict['B'].id = chainId
            motherPdb.child_list[0].add(pdbPair.child_list[0].child_dict['B'])

            motherPdb.transform(numpy.asarray(rel['rot']), rel['tran'])

            print 'Pair[{}]:   {}---{}'.format(str(i).ljust(chainLenDigits),
                lastNode.ljust(16), currNode.rjust(16))

        if haveCmd:
            tmpFile = './elfin.tmp'
            utils.savePdb(motherPdb, tmpFile)

            cmd.load(tmpFile, str(i) + '-' + pairName)
            cmd.hide('everything', 'all')
            cmd.show('cartoon', 'all')
            cmd.reset()
            cmd.util.cbc()

        self.bmarks.append({
            'pdb': motherPdb,
            'data': OrderedDict([
            ('nodes', nodes),
            ('coms', coms.tolist()),
            ('pairs', pairs)
            ])
        })

    def run(self, nBmarks, chainLen, maxRetries):
        count = 0
        retries = 0

        # Outer retry loop
        while count < nBmarks:
            print 'Attempt #{}'.format(retries)
            try:
                for i in xrange(count, nBmarks):
                    print('Genereating #{}/{}'.format(i+1, nBmarks))
                    self.gen(chainLen)
                    count = count + 1

                    # if all fine then retry is reset
                    retries = 0

            except UserWarning as uw:
                print 'Warning: {}'.format(uw)
                if maxRetries is not -1 and retries >= maxRetries:
                    utils.pauseCode()
                    print 'Warning: Maximum retries reached - stopping'
                    break
                else:
                    retries = retries + 1

        # When done trying, dump metadata and PDB
        bmNameLen = 7
        for bm in self.bmarks:
            #bmName = strftime("%m%d_%H%M%S", gmtime())
            bmName = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(bmNameLen))
            outFile = self.outDir + '/' + bmName

            json.dump(bm['data'],
                open(outFile + '.json', 'w'),
                separators=(',', ':'),
                ensure_ascii=False,
                indent=4)

            utils.savePdb(bm['pdb'], outFile + '.pdb')

main()
