#!/bin/bash
#************************************************************************************
#
# Author:       Stephen Copplestone
# Institution:  boltzplatz, Stuttgart
# Date:         08.07.2020
#
# Description:  This script will set the supplied git commit hash under src/commit.h
#
#************************************************************************************

# $1: path to commit.h

GITCOMMIT=$(git rev-parse HEAD)

if [ -f "$1" ]; then
  GITCOMMITFILE=$(grep "#define GIT_CURRENT_COMMIT" "$1" | cut -d "'" -f2)
  if [ "$GITCOMMIT" != "$GITCOMMITFILE" ]; then
    GITCOMMITQUOTES="'$GITCOMMIT'"
    sed -i -e 's/.*#define GIT_CURRENT_COMMIT.*/#define GIT_CURRENT_COMMIT  '"$GITCOMMITQUOTES"'/' "$1"
  fi
fi
