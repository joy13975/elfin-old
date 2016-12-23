#!/usr/bin/env python

import code
import os

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

def safeExec(func, *args):
    try:
        func(*args)
    except Exception as e:
        print(e)
        interact(locals())
        raise e