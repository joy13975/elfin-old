#!/usr/bin/env python
import glob

files = glob.glob('scores/*.comp.sc')
nFiles = len(files)
rmsds = []

for i in range(0, nFiles):
    with open(files[i], 'r') as file:
		line = file.read().split('\n')[1]
		rmsdStr = line.split(' ')[-2]
		print '{} RMSD: {}'.format(files[i], rmsdStr)
		rmsds.append(float(rmsdStr))

maxRmsd = max(rmsds)
print 'Average: {}, Min: {}, Max: {}'.format(sum(rmsds)/nFiles, min(rmsds), maxRmsd)

if(maxRmsd > 5.0):
	print 'WARNING: One or more molecules exceed 5A RMSD!'