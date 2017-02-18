#!/usr/bin/env python

import glob
from utils import *

def main():
    if len(sys.argv) < 2:
        print './Synth.py <specFile.{json|pdb}> <optional: outFile>'
        exit()

    specFile = sys.argv[1]
    specExt = specFile[specFile.rfind('.'):]

    if specExt == '.json':
        spec = readJSON(specFile)
        targetLen = len(spec['nodes'])
    elif specExt == '.csv':
        with open(specFile, 'r') as file:
            pts = [[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n')]

        spec = {'coms': pts}
        npts = np.asarray(pts)

        (avgD, minD, maxD) = getXDBStat()
        # Use total length/avgD as heuristic. avgD is average xDB pair distance
        targetLen = int(np.ceil(sum(np.linalg.norm(npts-np.roll(npts, 1, axis=0), axis=1)) / avgD))
    else:
        print 'Unknown spec file type: {}'.format(specExt)
        exit()

    scale = 1.0
    outFile = (sys.argv[2] if len(sys.argv) > 2 else specFile)
    outExt = outFile[outFile.rfind('.'):]
    outFile = outFile.replace(outExt, suffixPdb(
            'Synth', 
            'Main',
            scale, 
            targetLen))

    xDB = readJSON('res/xDB.json')
    makePdbFromNodes(xDB, spec['nodes'], 'res/centered_pdb/pair',
        outFile)
            

if __name__ == '__main__':
    main()