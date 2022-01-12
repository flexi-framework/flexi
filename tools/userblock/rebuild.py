#!/usr/bin/python
# -*- coding: utf-8 -*-
#************************************************************************************
#
# Author:       Matthias Sonntag, Thomas Bolemann
# Institution:  Inst. of Aero- and Gasdynamics, University of Stuttgart
# Date:         14.07.2016
#
# Description:  This script will rebuild a specific code configuration from
#               useblock data contained in an HDF5 file.
#
# Input:        1: directory to place build, 2: HDF5 file with userblock data
#
#************************************************************************************

import argparse
import subprocess
import sys
from extract_userblock import get_userblock, get_part
import os
import re
from distutils import spawn

parser = argparse.ArgumentParser(description='Rebuild code revision from userblock data' + \
                                             'contained in HDF5 state file.')
parser.add_argument('dir'  ,help='Name of empty directory where the rebuild will take place')
parser.add_argument('state',help='HDF5 state file containing userblock')

args = parser.parse_args()

if not os.path.exists(args.dir) :
   os.mkdir(args.dir)
else :
   if os.listdir(args.dir) != []: # is dir empty?
      print os.listdir(args.dir)
      print "Rebuild directory is not empty => exit!"
      exit(1)

# get svn
try:
  userblock = ""
  userblock = get_userblock(args.state,userblock)
except:
  print 'Error while extracting userblock.'
  exit(1)
git_url = get_part(userblock, "GIT URL")
git_url = git_url.strip()

# branch name
gitbranch = get_part(userblock,"GIT BRANCH").strip()
branch = gitbranch.split()[0]
print 'Branch to be rebuilt is ' + branch
print 'Branch revision ' + gitbranch.split()[1]

# checkout code
gitrevision = get_part(userblock,"GIT REFERENCE").strip()
rev = gitrevision.split()[1]
print 'Reference branch name ' + gitrevision.split()[0]
print 'Reference branch revision ' + rev

cmd = "git clone " + git_url + " " + args.dir
subprocess.call(cmd, shell=True)
cmd = "git checkout -b " + branch + " " + rev
os.chdir(args.dir)
subprocess.call(cmd, shell=True)

## apply git format-patch
#git_patch = get_part(userblock, "GIT FORMAT-PATCH")
#if git_patch :
#  f = open("format_patch", 'w')
#  f.write(git_patch)
#  f.close()
#  try:
#    subprocess.call("git am format_patch", shell=True)
#  except:
#    print 'Error while patching source code.'
#    exit(1)

# apply simple diff
git_diff = get_part(userblock, "GIT DIFF")
if git_diff :
  f = open("diff_patch", 'w')
  f.write(git_diff)
  f.close()
  try:
    subprocess.call("patch -p1 < diff_patch", shell=True)
  except:
    print 'Error while patching source code.'
    exit(1)

# write ini file
if not os.path.exists("ini"):
  os.mkdir("ini")
ini = get_part(userblock, "INIFILE")
f = open(os.path.join("ini", "parameter.ini"), 'w')
f.write(ini)
f.close()

# configure
builddir = "build"
if not os.path.exists(builddir):
  os.mkdir(builddir)
cmake = get_part(userblock, "CMAKE")
f = open(os.path.join(builddir, "config.cmake"), 'w')
f.write(cmake)
f.close()

cmakepath = spawn.find_executable("cmake")
if not cmakepath:
  print 'CMake not found, configuring not possible.'
  exit(1)

os.chdir(builddir)
try:
  p = subprocess.call(["cmake", "-C", "config.cmake" ,"../"])
except:
  print 'Error while configuring the build.'
  exit(1)

# make
try:
  p = subprocess.call(["make", "-j"])
except:
  print 'Error while compiling the code.'
  exit(1)

