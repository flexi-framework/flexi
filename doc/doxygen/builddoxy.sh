#!/bin/bash


EQUATION=navierstokes/
TESTCASE=default/
LIFTING=br1/

#=======================================
# DO NOT CHANGE ANYTHING BELOW THIS LINE
#=======================================

cd "$(dirname "$0")" # change to directory of this script
CURDIR=`readlink -f .`
cd $CURDIR
DOXYDIR=${CURDIR##*/}
if [ "$DOXYDIR" != "doxygen" ]; then
   echo "Security exit:"
   echo "builddoxy.sh must be placed in a directory with name 'doxygen'!"
   exit 1
fi   

function deletesrc {
   cd $CURDIR
   rm -rf src
}

trap deletesrc EXIT

#delete old files
rm -rf doxygen
cp -r ../../src .
cp -r ../../README.md .
cp -r ../../CONTRIBUTORS.md .
cp -r ../../INSTALL.md .
cp -r ../../LICENSE.md .
cp -r ../../REFERENCE.md .
cp -r ../../REGGIE.md .

# delete images of build status and gpl
sed -i '3,4d' README.md

# comment INTERFACE/END INTERFACE lines and MODULE PROCEDURE lines so doxygen will create call graphs
FILELIST=`find . -name '*.f90'`
for file in $FILELIST
do
	sed -i 's/END INTERFACE/!END INTERFACE/g' $file
	sed -i 's/INTERFACE/!INTERFACE/g' $file
	sed -i 's/MODULE PROCEDURE/!MODULE PROCEDURE/g' $file
done
# rename f90 -> F90: doxygen preprocessor is only run for these files
rename 's/\.f90$/\.F90/' $FILELIST

# delete all equation dirs and testcase dirs except the one specified above
function cleandirs {
DIRS=`ls -d */`
for dir in $DIRS; do
  if [ "$dir" != "$1" ]; then #delete all dirs except $1
    rm -r $dir
  fi
done
}
cd src/equations
cleandirs $EQUATION
cd ../testcase
cleandirs $TESTCASE
cd ../dg
cleandirs $LIFTING
cd ../..

# run doxygen
doxygen doxyconfig > doxylog.out

# comment in intefaces in html files
cd doxygen/html/
FILELIST=`find . -name '*.html'`
for file in $FILELIST
do
  sed -i 's/!END INTERFACE/END INTERFACE/g' $file
  sed -i 's/!INTERFACE/INTERFACE/g' $file
  sed -i 's/!MODULE PROCEDURE/MODULE PROCEDURE/g' $file
  sed -i 's/\.F90/\.f90/g' $file
done

