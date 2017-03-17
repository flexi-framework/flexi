#!/bin/bash
#************************************************************************************
#
# Author:       Thomas Bolemann
# Institution:  Inst. of Aero- and Gasdynamics, University of Stuttgart
# Date:         07.07.2016
#
# Description:  This script will generate a userblock file to be appended to
#               executables or simulation results enabling the exact rebuilding of
#               of the simulation code, which was used to generate those results.
#               A patch to the remote Git branch of the current code is generated
#               or to the master, if the branch does not exist on the remote.
# 
#************************************************************************************

# $1: CMAKE_RUNTIME_OUTPUT_DIRECTORY
# $2: CMAKE_CURRENT_BINARY_DIR

if [ ! -d "$1" ]; then
  exit 1;
fi

# get branch name (only info)
BRANCHNAME=$(git rev-parse --abbrev-ref HEAD)
PARENTNAME=$BRANCHNAME
PARENTCOMMIT=$(git show-ref | grep "origin/$BRANCHNAME$" | cut -b -40)

if [ -z "$PARENTCOMMIT" ]; then
  LOCBRANCHNAME=$BRANCHNAME
  # recursively search for parent branch
  FOUND=0
  while [ $FOUND -eq 0 ]; do
    # get commit on server, where branch started
    COLUMN=$((    $(git show-branch | grep '^[^\[]*\*'  | head -1 | cut -d* -f1 | wc -c) - 1 )) 
    START_ROW=$(( $(git show-branch | grep -n "^[\-]*$" | cut -d: -f1) + 1 )) 
    PARENTNAME=$(   git show-branch | tail -n +$START_ROW | grep -v "^[^\[]*\[$LOCBRANCHNAME" | grep "^.\{$COLUMN\}[^ ]" | head -n1 | sed 's/.*\[\(.*\)\].*/\1/' | sed 's/[\^~].*//')
    if [ -z "$PARENTNAME" ]; then
      break
    fi
  
    PARENTCOMMIT=$(git show-ref | grep "origin/$PARENTNAME$" | cut -b -40)
    if [ -z "$PARENTCOMMIT" ]; then
      LOCBRANCHNAME=$PARENTNAME 
    else
      FOUND=1
      break
    fi
  done

  if [ $FOUND -eq 0 ]; then
    PARENTNAME='master'
    PARENTCOMMIT=$(git rev-parse origin/master)
    echo "WARNING: Could not find parent commit, creating userblock diff to master."
  fi
fi

cd $1
echo "{[( CMAKE )]}"               >  userblock.txt
cat configuration.cmake            >> userblock.txt
echo "{[( GIT BRANCH )]}"          >> userblock.txt
echo "$BRANCHNAME"                 >> userblock.txt
echo $(git rev-parse HEAD)         >> userblock.txt

# Reference is the start commit, which is either identical to
# the branch, if it exists on the remote or points to the first
# real parent in branch history available on remote.
echo "{[( GIT REFERENCE )]}"       >> userblock.txt
echo "$PARENTNAME"                 >> userblock.txt
echo $PARENTCOMMIT                 >> userblock.txt

#echo "{[( GIT FORMAT-PATCH )]}"    >> userblock.txt
## create format patch containing log info for commit changes
## PARENT should be identical to origin
#git format-patch $PARENTCOMMIT..HEAD --minimal --stdout >> $1/userblock.txt

# Also store binary changes in diff
echo "{[( GIT DIFF )]}"            >> userblock.txt
# commited changes
git diff -p $PARENTCOMMIT..HEAD    >> userblock.txt
# uncommited changes
git diff -p                        >> userblock.txt

echo "{[( GIT URL )]}"             >> userblock.txt
git config --get remote.origin.url >> userblock.txt

builddir=$2/CMakeFiles
echo "{[( flexilibF90.dir/flags.make )]}"  >> userblock.txt
cat $builddir/flexilibF90.dir/flags.make   >> userblock.txt
echo "{[( flexilib.dir/flags.make )]}"     >> userblock.txt
cat $builddir/flexilib.dir/flags.make      >> userblock.txt
echo "{[( flexi.dir/flags.make )]}"        >> userblock.txt
cat $builddir/flexi.dir/flags.make         >> userblock.txt

# Compress the userblock
tar cJf userblock.tar.xz userblock.txt

# Build the module
objcopy -I binary -O elf64-x86-64 -B i386 --redefine-sym _binary_userblock_tar_xz_start=userblock_start --redefine-sym _binary_userblock_tar_xz_end=userblock_end --redefine-sym _binary_userblock_tar_xz_size=userblock_size userblock.tar.xz userblock.o
rm userblock.tar.xz
