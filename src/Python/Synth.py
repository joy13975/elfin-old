#!/usr/bin/env python

import glob
from utils import *

def main():
    if len(sys.argv) < 2:
        print './Synth.py <nodes.json> <optional: outFile>'
        exit()

    specFile = sys.argv[1]
    specExt = specFile[specFile.rfind('.'):]

    if specExt == '.json':
        spec = readJSON(specFile)
        targetLen = len(spec['nodes'])
    else:
        print 'Unknown spec file type: {}'.format(specExt)
        exit()

    scale = 1.0
    outFile = (sys.argv[2] if len(sys.argv) > 2 else specFile)

    dotIndex = outFile.rfind('.')
    if dotIndex == -1:
    	outExt = ''
    	outFile = outFile + suffixPdb(
		            'Synth', 
		            'Main',
		            scale, 
		            targetLen)
    else:
		outExt = outFile[dotIndex:]
		outFile = outFile.replace(outExt, suffixPdb(
		        'Synth', 
		        'Main',
		        scale, 
		        targetLen))

    xDB = readJSON('res/xDB.json')
    makePdbFromNodes(xDB, spec['nodes'], 'res/centered_pdb/pair',
        outFile, movieMode=True)
            

if __name__ == '__main__':
    main()