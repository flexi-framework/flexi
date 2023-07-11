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
# $2: CMAKE_CACHEFILE_DIR
# $3: CMAKE_CACHE_MAJOR_VERSION.CMAKE_CACHE_MINOR_VERSION.CMAKE_CACHE_PATCH_VERSION
# $4: CMAKE_CURRENT_SOURCE_DIR

if [ ! -d "$1" ]; then
  exit 1;
fi
if [ ! -d "$2" ]; then
  exit 1;
fi
if [ ! -d "$4" ]; then
  exit 1;
fi

#----------------------------------------------------------------------------------------------------
# Build ID and Cmake config
#----------------------------------------------------------------------------------------------------

cd "$1"
BUILD_ID=$(date +%Y-%m-%d_%H-%M-%S-%N)
echo "{[( CMAKE )]}"                > build_info.txt
cat configuration.cmake            >> build_info.txt

#----------------------------------------------------------------------------------------------------
# Git (outsourced to other file to be available for other applications)
#----------------------------------------------------------------------------------------------------

IS_GIT_DIR=$(git rev-parse --git-dir 2>/dev/null)

if [ $IS_GIT_DIR ]; then
  ## Code ID is a git hash of the working tree
  #CODE_ID=$(bash $4/tools/userblock/hash_tree.sh)
  # ATTENTION: Requires Git operations in current git repo!
  CODE_ID=$(git rev-parse HEAD)
  # TODO: Compute instead unique hash from git commit and diff!

  sh $4/tools/userblock/print_git_info.sh > code_info.txt
else
  CODE_ID="ERROR_COULD_NOT_GET_CODE_INFO_FROM_GIT"
  touch code_info.txt
fi


#----------------------------------------------------------------------------------------------------
# More detailed Cmake Info
#----------------------------------------------------------------------------------------------------

## change directory to cmake chache dir
#cd "$2/CMakeFiles"
## copy compile flags of the flexi(lib) to build_info
#echo "{[( libflexistatic.dir/flags.make )]}" >> $1/build_info.txt
#cat libflexistatic.dir/flags.make            >> $1/build_info.txt
#echo "{[( libflexishared.dir/flags.make )]}" >> $1/build_info.txt
#cat libflexishared.dir/flags.make            >> $1/build_info.txt
#echo "{[( flexi.dir/flags.make )]}"          >> $1/build_info.txt
#cat flexi.dir/flags.make                     >> $1/build_info.txt

## change directory to actual cmake version
#cd "$3"
## copy detection of compiler to build_info
#echo "{[( COMPILER VERSIONS )]}"           >> $1/build_info.txt
#cat CMakeFortranCompiler.cmake             >> $1/build_info.txt

#----------------------------------------------------------------------------------------------------
# Transfer to C fprintf commands
#----------------------------------------------------------------------------------------------------

cd "$1" # go back to the runtime output directory
# generate C print commands to print userblock:
#      replace \ by \\
#                    replace % by %%
#                                     replace " by \"
#                                                               prepend fprintf to line
#                                                                               append end of line
sed -e 's/\\/\\\\/g' -e 's/\%/\%\%/g' -e 's/"/\\"/g' -e 's/^/   fprintf(fp, "/' -e 's/$/\\n");/' build_info.txt > build_info_print.txt
sed -e 's/\\/\\\\/g' -e 's/\%/\%\%/g' -e 's/"/\\"/g' -e 's/^/   fprintf(fp, "/' -e 's/$/\\n");/' code_info.txt > code_info_print.txt
# copy empty source file template
cp "$4/src/output/print_userblock.c" .
# insert userblock print commands
sed -i -e '/INSERT_CODE_INFO_HERE/r code_info_print.txt' print_userblock.c
sed -i -e '/INSERT_BUILD_INFO_HERE/r build_info_print.txt' print_userblock.c
#insert build id to be returned by function
sed -i -e "s/DUMMY_INSERT_CODE_ID_HERE/$CODE_ID/g" print_userblock.c
sed -i -e "s/DUMMY_INSERT_BUILD_ID_HERE/$BUILD_ID/g" print_userblock.c
