import code
import os
from collections import namedtuple
from numpy import mean

import pyrosetta as pyros
from pyrosetta import rosetta as ros
pyros.init()

def caCoords(pose, rosFormat = False):
    if rosFormat:
        coords = ros.utility.vector1_numeric_xyzVector_double_t()
    else:
        coords = [];

    for pos in range(1, pose.size() + 1):
        if rosFormat:
            coords.append(pose.residue(pos).xyz("CA"))
        else:
            xyz = pose.residue(pos).xyz("CA")
            coords.append(xyzVecToList(xyz))

    return coords;

def getCenterOfMass(pose):
    caXYZs = caCoords(pose)

    return mean(caXYZs, axis=0)

def xyzVecToList(xyz):
    return [xyz.at(0), xyz.at(1), xyz.at(2)];

def getDockingRotation(child, mother):
    # Child moves and aligns to mother
    mothercoords = caCoords(mother, True)
    childCoords = caCoords(child, True)

    R = ros.numeric.xyzMatrix_double_t();
    motherToOrigin = ros.numeric.xyzVector_double_t();
    childToOrigin = ros.numeric.xyzVector_double_t();

    ros.protocols.toolbox.superposition_transform(mothercoords,
        childCoords,
        R,
        motherToOrigin,
        childToOrigin);

    rot = [R.xx(), R.yx(), R.zx(),
    R.xy(), R.yy(), R.zy(),
    R.xz(), R.yz(), R.zz()];

    cto = xyzVecToList(childToOrigin);
    mto = [x * -1 for x in xyzVecToList(motherToOrigin)];

    return cto,rot,mto

def floatListStr(floats, precision):
    fmt = "{:." + str(precision) + "f}"
    return ", ".join(fmt.format(f) for f in floats)

def processPDB(xDB, file, import_opts):
    print "Processing " + file

    #pair is a 1-based pose vector
    in_pose = pyros.pose_from_file(file)

    pair = in_pose.split_by_chain()

    com1 = getCenterOfMass(pair[1])
    com2 = getCenterOfMass(pair[2])

    # Compute translation
    trans = com2 - com1;
    # code.interact(local=locals())

    # Load individual monomers for alignment
    _, filename = os.path.split(file)
    names = filename[0:filename.find(".")].split("-")
    file1 = "res/single/" + names[0] + ".pdb"
    m1 = pyros.pose_from_file(file1)
    cto1,rot1,mto1 = getDockingRotation(pair[1], m1)
    # code.interact(local=locals())

    file2 = "res/single/" + names[1] + ".pdb"
    m2 = pyros.pose_from_file(file2)
    cto2,rot2,mto2 = getDockingRotation(pair[2], m2)
    # code.interact(local=locals())

    # xDB.write("".format())

    # pair_name, trans, single_name_1, com1, cto1, rot1, mto1, single_name_2, com2, cto2, rot2, mto2
    pair_name = filename.replace(".pdb", "");
    single_name_1 = names[0];
    single_name_2 = names[1];
    prec = 12;

    strParts = [pair_name,
    floatListStr(trans, prec),
    single_name_1, floatListStr(com1, prec),
    floatListStr(cto1, prec), floatListStr(rot1, prec), floatListStr(mto1, prec),
    single_name_2, floatListStr(com2, prec),
    floatListStr(cto2, prec), floatListStr(rot2, prec), floatListStr(mto2, prec)];

    outStr = ", ".join("{}".format(s) for s in strParts) + "\n";
    # code.interact(local=locals())
    xDB.write(outStr);

libDirs = ["res/pair/"]
import_opts=ros.core.import_pose.ImportPoseOptions();
# import_opts.set_pack_missing_sidechains(False);

xDB = open("res/xDB.csv", "w");

for libDir in libDirs:
    for file in os.listdir(libDir):
        if file.lower().endswith(".pdb"):
            processPDB(xDB, libDir + file, import_opts)
