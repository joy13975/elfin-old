#!/usr/bin/env python

import glob
from Designer import *
import Greedy
from utils import *

def main():
    xDB     = readJSON('res/xDB.json')
    bmDir 	= 'bm/pc10'

    designers = []
    designers.append(Greedy.GreedyDesigner(xDB))
    # MCDesigner
    # ...

    # Process all benchmarks
    for jsonFile in glob.glob(bmDir + '/*.json'):
    	spec = readJSON(jsonFile)

    	for designer in designers:
	    	print 'Benchmarking {} on {}'.format(designer.__class__.__name__, jsonFile)
    		(nodes,shape,score) = designer.design(spec)
    		if nodes == spec['nodes']:
    			print 'Pass: score {}'.format(score)
    		else:
    			print 'Failed: score {}'.format(score)
    			pauseCode()

if __name__ =='__main__': safeExec(main)