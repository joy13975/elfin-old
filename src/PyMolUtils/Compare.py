
import numpy as np
from pymol import cmd

# Not sure how to just figure out where elfin is located
# So we need to load our library this way
elfinDir = '/Users/joy/src/elfin/'
elfinPyLibDir = elfinDir + '/src/python/'
elfinMovieDir = elfinDir + '/movieOutput/'
import imp
utils = imp.load_source('utils', elfinPyLibDir + '/utils.py')
Kabsch = imp.load_source('Kabsch', elfinPyLibDir + '/Kabsch.py')

cmd.run('/Users/joy/src/elfin/src/PyMolUtils/LineUtils.py');

def compare_csv(specCSV, solCSV):
    specPts = utils.readCSVPoints(specCSV)
    solPts = utils.readCSVPoints(solCSV)

    # Centre both pts
    centredSpec = specPts - np.mean(specPts, axis=0)
    centredSol = solPts - np.mean(solPts, axis=0)

    # Draw specification
    draw_pts(centredSpec, color=[0.7,0,0])

    # Equalise sample points
    specEs, solEs = utils.upsample(centredSpec, centredSol)

    draw_pts(specEs, color=[0.5,0.5,0])

    # Find Kabsch rotation for solution -> spec
    R = Kabsch.kabsch(solEs, specEs)

    solEsR = np.dot(solEs, R)

    draw_pts(solEsR, color=[0,0.5,0.7])

    cmd.reset()
    cmd.set("depth_cue", 0)


cmd.extend("compare_csv", compare_csv)