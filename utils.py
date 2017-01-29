#!/usr/bin/env python

import traceback, sys, code
import os
import inspect
import Bio.PDB

def readPdb(customName, inFile, permissive=0):
    parser = Bio.PDB.PDBParser(permissive)
    structure = parser.get_structure(customName, inFile)
    return structure

def savePdb(struct, saveFile):
    io = Bio.PDB.PDBIO()
    io.set_structure(struct)
    io.save(saveFile)

def getPdbSingles(pStruct):
    pSingles = []
    for pSingle in pStruct.get_chains():
        pSingles.append(pSingle)

    nPs = len(pSingles)
    assert nPs == 2

    return pSingles

def interact(globalVars=None, localsVars=None):
    print "Entering interactive mode"
    print
    if(localsVars == None):
        localsVars = locals()
    if(globalVars == None):
        globalVars = globals()
    code.interact(local=dict(globalVars, **localsVars))

def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def pauseCode(failFrame=None):
    print '---------pauseCode()---------'
    type, value, tb = sys.exc_info()
    last_frame = lambda tb=tb: last_frame(tb.tb_next) if tb.tb_next else tb
    frame = last_frame().tb_frame
    ns = dict(frame.f_globals)
    ns.update(frame.f_locals)

    if failFrame is None:
        ns.update(inspect.currentframe().f_back.f_locals)
    else:
        ns.update(failFrame.f_locals)

    code.interact(local=ns)

def safeExec(func, *args):
    try:
        func(*args)
    except Exception as e:
        print '---------safeExec() caught exception---------'
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pauseCode(inspect.currentframe().f_back)