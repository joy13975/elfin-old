#!/usr/bin/env python

import utils
import json
import argparse
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def main():
	ap = argparse.ArgumentParser(description='Compute the number of combinations for a given MMC protein length');
	ap.add_argument('--xdbFile', default='./res/xDB.json')
	ap.add_argument('--length', type=int, default=21)

	globals().update(vars(ap.parse_args()))

	print xdbFile
	with open(xdbFile, 'r') as file:
		xdb = json.load(file)
		pd = xdb['pairsData']

		singleNames = xdb['singlesData'].keys()
		dim = len(singleNames)

	dim = len(singleNames)
	adjMat = np.zeros([dim, dim])

	for pdk in pd.keys():
		i1 = singleNames.index(pdk)
		for pdkk in pd[pdk].keys():
			i2 = singleNames.index(pdkk)
			adjMat[i1][i2] = 1

	for l in xrange(1, length):
		nCombs = np.sum(np.linalg.matrix_power(adjMat, l))

		print 'L={}, NC={}'.format(l+1, nCombs)

if __name__ == '__main__':
	main()