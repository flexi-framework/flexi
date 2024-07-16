#!/bin/bash
#************************************************************************************
#
# Author:       Thomas Bolemann
# Institution:  Inst. of Aero- and Gasdynamics, University of Stuttgart
# Date:         07.07.2016
#
# Description:  This script will read GIT info and print it to a specified file
#
#************************************************************************************

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
   #echo "WARNING: Could not find parent commit, creating userblock diff to master."
 fi
fi

echo "{[( GIT BRANCH )]}"
echo "$BRANCHNAME"
echo $(git rev-parse HEAD)

# Reference is the start commit, which is either identical to
# the branch, if it exists on the remote or points to the first
# real parent in branch history available on remote.
echo "{[( GIT REFERENCE )]}"
echo "$PARENTNAME"
echo $PARENTCOMMIT

#echo "{[( GIT FORMAT-PATCH )]}"
## create format patch containing log info for commit changes
## PARENT should be identical to origin
#git format-patch $PARENTCOMMIT..HEAD --minimal --stdout

# Also store binary changes in diff
echo "{[( GIT DIFF )]}"
# commited changes
git diff -p $PARENTCOMMIT..HEAD
# uncommited changes
git diff -p HEAD

echo "{[( GIT URL )]}"
git config --get remote.origin.url
