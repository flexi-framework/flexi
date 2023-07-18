#!/bin/bash
#************************************************************************************
#
# This script will output a hash of the current source code status.
# To this end, unstaged changes are added to the tree, the tree is hashed,
# and previously unstaged changes are reset. A comparison of the git status before
# and after the operation ensures that this caused no changes to the git repo.
#
#************************************************************************************

# save old git status for safety check after operation
PREV=$(git status --porcelain)
# get unstaged changes
UNSTAGED=($(git status --porcelain | grep '^\s[A-Z]\s*' | cut -c 4-))

# save old directory
PWD=$(pwd)
# change to git root directory
GITROOT=$(git rev-parse --show-toplevel)
cd $GITROOT

#add all unsategd changes
for i in "${UNSTAGED[@]}"
do
  git add "$i"
done

# get hash of git tree and output
TREE=$(git write-tree)
echo "$TREE"

# reset unstaged changes (go back to original version)
for i in "${UNSTAGED[@]}"
do
  git reset -q -- "$i"
done
# go back to old directory
cd $PWD

# compare git status with previous one
POST=$(git status --porcelain)
if [[ ${PREV} != ${POST} ]]
then
  >&2 echo "GIT ERROR. CHECK GIT STATUS."
  exit 1
fi
