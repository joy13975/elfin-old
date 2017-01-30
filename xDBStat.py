#!/usr/bin/env python
from utils import *
import numpy as np

xdb = readJSON('res/xDB.json')

pd = xdb['pairsData']
dists = []
for s1 in pd.keys():
	for s2 in pd[s1].keys():
		dists.append(np.linalg.norm(pd[s1][s2]['comB']))


print 'Distances avg: {}, min: {}, max: {}'.format(np.average(dists), min(dists), max(dists))