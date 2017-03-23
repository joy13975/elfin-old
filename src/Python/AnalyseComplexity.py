#!/usr/bin/env python

import utils
import json
import argparse
import numpy as np
from decimal import Decimal

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

	mmcYs = []
	for l in xrange(1, length):
		nCombs = Decimal(np.sum(np.linalg.matrix_power(adjMat, l)))

		mmcYs.append(nCombs)
		# print 'L={}, NC={}'.format(l+1, nCombs)

	Xs = np.asarray(range(2, length + 1))
	mmcYs = np.asarray(mmcYs)
	aaYs = np.power(Decimal(20 * 80), Xs) # A typical repeat module is ~80 AA

	fig, ax1 = plt.subplots()

	ax1.set_xlabel('Design Length/AA')
	# ax2 = ax1.twinx()

	ax1.plot(Xs, aaYs, label='AA')
	ax1.set_ylabel('No. Combs (log scale)')
	ax1.set_yscale('log')

	ax1.plot(Xs, mmcYs, label='MMC')

	xTickIds = np.arange(3, len(Xs) + 1, 5)
	xTickIds = np.insert(xTickIds, 0, 0)
	plt.xticks(xTickIds+2, Xs[(xTickIds)])

	plt.legend()
	plt.show()

if __name__ == '__main__':
	main()