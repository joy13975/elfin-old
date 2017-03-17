#!/usr/bin/env python

import sys
import utils
import json
from collections import OrderedDict

def main():
	outputDir = './gridSearchConfigs/'

	if len(sys.argv) < 2:
		print 'Using DEFAULT output directory: {}'.format(outputDir)
		utils.mkdir(outputDir)

		# print 'Usage: ./GenGAConfigs.py <outputDir>'
		# exit()
	else:
		outputDir = sys.argv[1] + '/'
		print 'Using output directory: {}'.format(outputDir)
		utils.mkdir(outputDir)

	# Define the grid
	chromoLenDevs 		= [0.1, 0.2, 0.3]
	gaPopSizes 			= [int(100e3)] 					
	gaIters 			= [int(1e3)]						
	gaSurviveRates 		= [0.005, 0.01, 0.02]
	gaCrossRates 		= [0.3, 0.5, 0.7]
	gaPointMutateRates 	= []
	gaLimbMutateRates  	= []

	# Create 3 configs - for 3 different benchmarks shapes
	# using the same config
	bmNames = ['6vjex8d', '9y8hxgo', 'j0m06n4']

	# PM and LM rates depend on Cross rates
	pmRatios = (0.25, 0.5, 0.75)
	for cr in gaCrossRates:
		# Each remaining portion after cross rates
		# generate 3 ratios of RM and LM
		rem = 1 - cr
		for pmRatio in pmRatios:
			(pm, lm) = (pmRatio * rem, (1 - pmRatio) * rem)
			gaPointMutateRates.append(pm)
			gaLimbMutateRates.append(lm)

	nRuns = len(chromoLenDevs) * len(gaPopSizes) * len(gaIters) * \
		len(gaSurviveRates) * len(gaCrossRates) * len(gaPointMutateRates) * \
		len(gaLimbMutateRates) * len(bmNames)
	print 'Total runs needed: {}'.format(nRuns)

	# Write all combinations of GA parameters to output
	configId = 0
	for cld in chromoLenDevs:
		for gps in gaPopSizes:
			for gi in gaIters:
				for gsr in gaSurviveRates:
					for gcr in gaCrossRates:
						for (gpmr, glmr) in zip(gaPointMutateRates, gaLimbMutateRates):
							for bmName in bmNames:
								outputName = 'gs_{}_{}'.format(configId, bmName)
								bmOutputDir = './gs_out/{}/'.format(outputName)
								utils.mkdir(bmOutputDir)

								configJson = OrderedDict([
									('inputFile', './bm/l10/{}.json'.format(bmName)),
									('xDBFile', './res/xDB.json'),
									('outputDir', bmOutputDir),
									('randSeed', '0x600d1337'),

									('chromoLenDev', cld),
									('gaPopSize', gps),
									('gaIters', gi),
									('gaSurviveRate', gsr),
									('gaCrossRate', gcr),
									('gaPointMutateRate', gpmr),
									('gaLimbMutateRate', glmr),
								])

								json.dump(configJson,
									open(outputDir + '{}.json'.format(outputName), 'w'),
									separators=(',', ':'),
									ensure_ascii=False,
									indent=4)

							configId = configId + 1

	print 'Max GS config ID: {}'.format(configId - 1)

if __name__ == '__main__':
	main()
