#!/bin/bash
#************************************************************************************
#
# Author:       Thomas Bolemann
# Institution:  Inst. of Aero- and Gasdynamics, University of Stuttgart
# Date:         18.04.2017
#
# Description:  This script to downloads, compiles and installs standard libraries
#               and tools and sets up the corresponding environment modules.
#
#               The script requires environment modules to be set up, as explained
#               in the Advanced Installation instructions in the Flexi user guide.
#               The installation of HDF5 and Paraview require MPI libraries to be
#               available for the selected compiler, the simplest way is to build
#               them using this script.
# 
#************************************************************************************

INSTALL_ROOT=/opt         # Directoy in which compiled programs and libraries are installed
                          # This may not be identical to or a subfolder of the module directory
INSTALL_REQUIRES_ROOT=0   # Set to 1 if if installation or module directory require root permissions.
                          # The script will ask for a sudo password during installation
NPROCS=4                  # Number of processes for compiling

# Items in the list:
# GNU:          Use GNU environment
# INTEL:        Use Intel environment
# Compile:      Compile using the selected environments
# Build module: Build modules for selected environments
# Set default:  Make the currently built module the default module

# NOTE: If "Compile" and "Build module" options are not set, the GNU/INTEL options are ignored

# Build with      GNU;    INTEL; Compile; Build module Set default
BUILD_OPENMPI=(    0        0       0          0            0       )
# Build with      GNU;    INTEL; Compile; Build module Set default
BUILD_HDF5=(       0        0       0          0            0       )
# Build with      GNU;           Compile; Build module Set default
BUILD_CMAKE=(      0        0       0          0            0       )
# Build with  Standard; NoVisit; Compile; Build module Set default
BUILD_PARAVIEW=(   0        0       0          0            0       )
# Build with      GNU;    INTEL; Compile; Build module Set default
BUILD_FFTW=(       0        0       0          0            0       )

# Path and version of HDF5
HDF5_NAME=hdf5
HDF5_VERSION=1.10.0-patch1
HDF5_DLPATH='http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.0-patch1/src/hdf5-1.10.0-patch1.tar.bz2'

# Path and version of MPI
# Note: MPI versions as of 1.10.3+ require autotools (aclocal) 1.15+ to be installed
OPENMPI_NAME=openmpi
OPENMPI_VERSION=2.1.0
OPENMPI_DLPATH='https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.0.tar.bz2'

# Path and version of CMake
CMAKE_NAME=cmake
CMAKE_VERSION=3.8.0
CMAKE_DLPATH='https://cmake.org/files/v3.8/cmake-3.8.0.tar.gz'

# Path and version of Paraview
PARAVIEW_NAME=paraview
PARAVIEW_VERSION=5.3.0
PARAVIEW_DLPATH='http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.3&type=source&os=all&downloadFile=ParaView-v5.3.0.tar.gz'

# Path and version of FFTW
FFTW_NAME=fftw
FFTW_VERSION=3.3.6-pl2
FFTW_DLPATH='http://www.fftw.org/fftw-3.3.6-pl2.tar.gz'

# Compiler specific subfolders for environment modules (i.e. compiler vendor + version )
GNU_PREFIX=gnu
INTEL_PREFIX=intel

#########################################
# Changes after this point at own risk
#########################################

# First check module paths
if [ -z "$MODULEPATH" ]; then
    echo "The variable MODULEPATH is not set, make sure the module environment is set up correctly."
    exit 1
fi  
IFS=':' read -ra TMP <<< "$MODULEPATH"
MYMODULEPATH=${TMP[0]}


ROOTDIR=$(pwd)
mkdir -p $INSTALL_ROOT

set -e
lockdir () {
  if [ $INSTALL_REQUIRES_ROOT == 1 ] ; then
    echo "Locking install directories."
    for MYLIB_NAME in "$@"
    do
      sudo -- sh -c \
          "find $INSTALL_ROOT/$MYLIB_NAME -type d -exec chmod o-w {} \; ;
           find $INSTALL_ROOT/$MYLIB_NAME -type f -exec chmod o-w {} \; ;
           chown -R root:root $INSTALL_ROOT/$MYLIB_NAME                 "
    done
    sudo -- sh -c \
        "find $MYMODULEPATH -type d -exec chmod o-w {} \; ;
         find $MYMODULEPATH -type f -exec chmod o-w {} \; ;
         chown -R root:root $MYMODULEPATH                 "
  fi
}

unlocklibdir () {
  if [ $INSTALL_REQUIRES_ROOT == 1 ] ; then
    echo "Unlocking install directories."
    if [ ${BUILD_MYLIB[2]} == 1 ]; then
      sudo -- sh -c \
          "mkdir -p $INSTALL_ROOT/$MYLIB_NAME                                 ;
           find     $INSTALL_ROOT/$MYLIB_NAME -type d -exec chmod o+rwx {} \; ;
           find     $INSTALL_ROOT/$MYLIB_NAME -type f -exec chmod o+rw  {} \; " 
    fi
  fi
}

unlockmoduledir () {
  if [ $INSTALL_REQUIRES_ROOT == 1 ] ; then
    echo "Unlocking module directories."
    if [ ${BUILD_MYLIB[3]} == 1 ]; then
      sudo -- sh -c \
          "mkdir -p $MYMODULEPATH                                 ;
           find     $MYMODULEPATH -type d -exec chmod o+rwx {} \; ;
           find     $MYMODULEPATH -type f -exec chmod o+rw  {} \; " 
    fi
  fi
}

trap "lockdir" EXIT
#trap "lockdir $HDF5_NAME $OPENMPI_NAME $CMAKE_NAME" EXIT

getversion_intel () {
  if hash ifort 2>/dev/null; then
    VER=$(ifort -v 2>&1 | cut -d" " -f3)
  else
    echo "Could not find Intel Fortran compiler."
    exit 1
  fi
  INTEL_DIR=$INTEL_PREFIX$VER
}

getversion_gnu () {
  if hash gfortran 2>/dev/null; then
    VER=$(gfortran -dumpversion)
  else
    echo "Could not find GNU Fortran compiler."
    exit 1
  fi
  GNU_DIR=$GNU_PREFIX$VER
}

check_module () {
  echo "Checking for module $1"
  modules=$(( module avail ) 2>&1 )
  if echo $modules | grep -q "$1" ; then
    echo 'Module found.'
  else
    echo "Could not find module $1. Make sure the module environment is set up correctly."
    exit 1
  fi
}

prepare_lib () {
  BUILD_DIR=$ROOTDIR/$MYLIB_NAME
  mkdir -p $BUILD_DIR

  cd $BUILD_DIR

  #download/unpack
  MYLIB_FILE=$(basename "$MYLIB_DLPATH")
  if [ ! -e $MYLIB_FILE ] ; then
    wget "$MYLIB_DLPATH"
  fi
  tar xf $MYLIB_FILE
  MYLIB_DIR=${BUILD_DIR}/$MYLIB_NAME-$MYLIB_VERSION
  if [ ! -d ${MYLIB_DIR} ] ; then
    MYLIB_DIR=$(find  . -maxdepth 1 -type d -iname "*$MYLIB_NAME*")
    if [ ! -d ${MYLIB_DIR} ] ; then
      echo "Could not find $MYLIB_DIR dir!"
      exit 1;
    fi
  fi
}

build_lib () {
  if [ ${BUILD_MYLIB[2]} != 1 ]; then
    return
  fi
  cd $BUILD_DIR

  MYLIB_INSTALL=${INSTALL_ROOT}/$MYLIB_NAME/$MYLIB_NAME-$MYLIB_VERSION/$COMPILER_DIR

  mkdir -p $COMPILER_DIR
  cp -r ${MYLIB_DIR}/* $COMPILER_DIR
  cd    $COMPILER_DIR
  if [ $MYLIB_USECMAKE == 1 ]; then
    mkdir -p build
    cd build
    cmake ../   ${MYLIB_OPTIONS}${MYLIB_INSTALL}
  else
    ./configure ${MYLIB_OPTIONS}${MYLIB_INSTALL}
  fi
  make -j ${NPROCS}

  if [ $INSTALL_REQUIRES_ROOT == 1 ] ; then
    sudo mkdir -p $MYLIB_INSTALL
    sudo chmod -R ugo+rwx $MYLIB_INSTALL
  else
         mkdir -p $MYLIB_INSTALL
  fi
  make install
  if [ $INSTALL_REQUIRES_ROOT == 1 ] ; then
    sudo chown -R root:root $MYLIB_INSTALL
  fi
}

build_module () {
  if [ ${BUILD_MYLIB[3]} != 1 ]; then
    return
  fi
  unlockmoduledir $MYLIB_NAME
  cd $ROOTDIR
  MYLIB_INSTALL=${INSTALL_ROOT}/$MYLIB_NAME/$MYLIB_NAME-$MYLIB_VERSION/$COMPILER_DIR
  MYLIB_INSTALL2=$(echo $MYLIB_INSTALL | sed -e 's/[\/&]/\\&/g')
  MYMODULE_DIR=$MYMODULEPATH/$MYLIB_NAME/$COMPILER_NAME
  mkdir -p $MYMODULE_DIR
  cat module_templates/${MYLIB_NAME} |
      sed -e 's/FIGIVERSION/'$MYLIB_VERSION'/g'  | \
      sed -e 's/FIGIPATH/'$MYLIB_INSTALL2'/g' \
      > $MYMODULE_DIR/${MYLIB_VERSION}_$COMPILER_DIR
  if [ ${BUILD_MYLIB[4]} == 1 ]; then
    echo '#%Module' > $MYMODULE_DIR/.version
    echo "set ModulesVersion \"${MYLIB_VERSION}_${COMPILER_DIR}\"" >> $MYMODULE_DIR/.version
  fi
  echo "Module for $MYLIB_NAME version $MYLIB_VERSION has been built."
}

#############################################################################
# 1. OpenMPI
#############################################################################
BUILD_MYLIB=("${BUILD_OPENMPI[@]}")
MYLIB_NAME=$OPENMPI_NAME
MYLIB_VERSION=$OPENMPI_VERSION
MYLIB_DLPATH=$OPENMPI_DLPATH
MYLIB_USECMAKE=0
if [ ${BUILD_MYLIB[2]} == 1 ]; then
  prepare_lib
  MYLIB_OPTIONS='--prefix='
  export CFLAGS=-fPIC
  export FCFLAGS=-fPIC
fi

if [ ${BUILD_MYLIB[0]} == 1 ]; then  # GNU
  COMPILER_NAME=gnu
  check_module env/$COMPILER_NAME
  module unload env
  module load env/$COMPILER_NAME
  export CC=gcc
  export FC=gfortran

  getversion_gnu
  COMPILER_DIR=$GNU_DIR

  build_lib
  build_module
fi

if [ ${BUILD_MYLIB[1]} == 1 ]; then  # INTEL
  COMPILER_NAME=intel
  check_module env/$COMPILER_NAME
  module unload env
  module load env/$COMPILER_NAME
  export CC=icc
  export FC=ifort

  getversion_intel
  COMPILER_DIR=$INTEL_DIR

  build_lib
  build_module
fi

#############################################################################
# 2. CMake
#############################################################################
BUILD_MYLIB=("${BUILD_CMAKE[@]}")
MYLIB_NAME=$CMAKE_NAME
MYLIB_VERSION=$CMAKE_VERSION
MYLIB_DLPATH=$CMAKE_DLPATH
MYLIB_USECMAKE=0
if [ ${BUILD_MYLIB[2]} == 1 ]; then
  prepare_lib
  MYLIB_OPTIONS='--prefix='
fi

if [ ${BUILD_MYLIB[0]} == 1 ]; then  # GNU
  COMPILER_NAME=gnu
  check_module env/$COMPILER_NAME
  module unload env
  module load env/$COMPILER_NAME
  export CC=gcc
  export FC=gfortran
  export CXX=g++

  getversion_gnu
  COMPILER_DIR=$GNU_DIR

  build_lib
  build_module
fi

#if [ ${BUILD_MYLIB[1]} == 1 ]; then  # INTEL
#  COMPILER_NAME=intel
#  COMPILER_DIR=$INTEL_DIR
#  check_module env/$COMPILER_NAME
#  module swap env env/$COMPILER_NAME
#  export CC=icc
#  export FC=ifort
#  export CXX=icpc
#
#  build_lib
#  build_module
#fi

#############################################################################
# 3. HDF5
#############################################################################
BUILD_MYLIB=("${BUILD_HDF5[@]}")
MYLIB_NAME=$HDF5_NAME
MYLIB_VERSION=$HDF5_VERSION
MYLIB_DLPATH=$HDF5_DLPATH
MYLIB_USECMAKE=0
if [ ${BUILD_MYLIB[2]} == 1 ]; then
  prepare_lib
  MYLIB_OPTIONS='--with-pic --enable-fortran --enable-parallel --disable-shared --prefix='
  export CC=mpicc
  export FC=mpif90
  export CFLAGS=-fPIC
  export FCFLAGS=-fPIC
fi

if [ ${BUILD_MYLIB[0]} == 1 ]; then  # GNU
  COMPILER_NAME=gnu
  check_module env/$COMPILER_NAME
  check_module openmpi/$COMPILER_NAME
  module unload env
  module load env/$COMPILER_NAME

  getversion_gnu
  COMPILER_DIR=$GNU_DIR

  build_lib
  build_module
fi

if [ ${BUILD_MYLIB[1]} == 1 ]; then  # INTEL
  COMPILER_NAME=intel
  check_module env/$COMPILER_NAME
  check_module openmpi/$COMPILER_NAME
  module unload env
  module load env/$COMPILER_NAME

  getversion_intel
  COMPILER_DIR=$INTEL_DIR

  build_lib
  build_module
fi

#############################################################################
# 4. Paraview
#############################################################################
BUILD_MYLIB=("${BUILD_PARAVIEW[@]}")
MYLIB_NAME=$PARAVIEW_NAME
MYLIB_VERSION=$PARAVIEW_VERSION
MYLIB_DLPATH=$PARAVIEW_DLPATH
MYLIB_USECMAKE=1
if [ ${BUILD_MYLIB[2]} == 1 ]; then
  prepare_lib
  cd $MYLIB_DIR
  patch -p1 < $ROOTDIR/patches/paraview_patch
fi

# always use GNU for Paraview
if [ ${BUILD_MYLIB[0]} == 1 ]; then  # default
  MYLIB_OPTIONS='-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DPARAVIEW_USE_VISITBRIDGE=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DCMAKE_INSTALL_PREFIX='
  COMPILER_NAME=standard
  check_module env/gnu
  module unload env
  module load env/gnu
  export CC=gcc
  export FC=gfortran
  export CXX=g++

  getversion_gnu
  COMPILER_DIR=${GNU_DIR}_standard

  build_lib
  build_module
fi

if [ ${BUILD_MYLIB[1]} == 1 ]; then  # no visit bridge (easier use of Flexi Paraview plugin)
  MYLIB_OPTIONS='-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DPARAVIEW_USE_VISITBRIDGE=OFF -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DCMAKE_INSTALL_PREFIX='
  COMPILER_NAME=novisit
  check_module env/gnu
  module unload env
  module load env/gnu
  export CC=gcc
  export FC=gfortran
  export CXX=g++

  getversion_gnu
  COMPILER_DIR=${GNU_DIR}_novisit

  build_lib
  build_module
fi

#if [ ${BUILD_MYLIB[1]} == 1 ]; then  # INTEL
#  COMPILER_NAME=intel
#  COMPILER_DIR=$INTEL_DIR
#  check_module env/$COMPILER_NAME
#  module swap env env/$COMPILER_NAME
#  export CC=icc
#  export FC=ifort
#  export CXX=icpc
#
#  build_lib
#  build_module
#fi

#############################################################################
# 5. FFTW
#############################################################################
BUILD_MYLIB=("${BUILD_FFTW[@]}")
MYLIB_NAME=$FFTW_NAME
MYLIB_VERSION=$FFTW_VERSION
MYLIB_DLPATH=$FFTW_DLPATH
MYLIB_USECMAKE=0
if [ ${BUILD_MYLIB[2]} == 1 ]; then
  prepare_lib
  MYLIB_OPTIONS='--prefix='
  export CC=mpicc
  export FC=mpif90
  export CFLAGS=-fPIC
  export FCFLAGS=-fPIC
fi

if [ ${BUILD_MYLIB[0]} == 1 ]; then  # GNU
  COMPILER_NAME=gnu
  check_module env/$COMPILER_NAME
  check_module openmpi/$COMPILER_NAME
  module unload env
  module load env/$COMPILER_NAME

  getversion_gnu
  COMPILER_DIR=$GNU_DIR

  build_lib
  build_module
fi

if [ ${BUILD_MYLIB[1]} == 1 ]; then  # INTEL
  COMPILER_NAME=intel
  check_module env/$COMPILER_NAME
  check_module openmpi/$COMPILER_NAME
  module unload env
  module load env/$COMPILER_NAME

  getversion_intel
  COMPILER_DIR=$INTEL_DIR

  build_lib
  build_module
fi
