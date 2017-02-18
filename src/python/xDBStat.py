#!/usr/bin/env python
from utils import *

def main():
	(avgD, minD, maxD) = getXDBStat()
	print 'Distances avg: {}, min: {}, max: {}'.format(avgD, minD, maxD)

if __name__ =='__main__': safeExec(main)