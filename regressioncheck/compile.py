#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import combinations 
from execute_cmd import execute_cmd

class BasedirNotFoundException(Exception) :
    def __init__(self, cwd):
        self.cwd = cwd
    def __str__(self):
        return "No basedir found. Started searching for 'CMakeLists.txt' in '%s'" % self.cwd 

class CMakeFailedException(Exception) :
    def __init__(self, build_directory, error):
        self.build_directory = build_directory
        self.error = error
    def __str__(self):
        return "cmake command faild in directory '%s'. Error Message:\n %s" % (self.build_directory, self.error)

class MakeFailedException(Exception) :
    def __init__(self, build_directory, error):
        self.build_directory = build_directory
        self.error = error
    def __str__(self):
        return "make command faild in directory '%s'. Error Message:\n %s" % (self.build_directory, self.error)

def find_basedir() :
    """Search 'CMakeLists.txt' in directories above current working directory.
    The directory containing the 'CMakeLists.txt' is the 'basedir'."""
    basedir = os.getcwd()                                           # start with current working directory
    found = os.path.exists(os.path.join(basedir, "CMakeLists.txt")) # check if actual directory is the basedir
    while not found :                                               # look upwards until basedir found
        basedir = os.path.dirname(basedir)                              # basedir = basedir/..
        found = os.path.exists(os.path.join(basedir, "CMakeLists.txt")) # check if actual directory is the basedir
        if basedir == "/" : break                                       # check if root of filesystem is reached

    if not found :
        raise BasedirNotFoundException(os.getcwd())

    return basedir







