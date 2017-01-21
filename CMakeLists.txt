# =========================================================================
# Flexi CMake configuration
# =========================================================================
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
# TODO: eventually require CMake 3.0 in future releases
# Advantages: Possible to split string
IF ("${CMAKE_VERSION}" VERSION_GREATER 3.0.0)
  cmake_policy (SET CMP0045 OLD)
  cmake_policy (SET CMP0046 OLD)
ENDIF()
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)


# =========================================================================
# Check whether we are on Github or at IAG
# =========================================================================
EXECUTE_PROCESS(COMMAND git ls-remote --get-url OUTPUT_VARIABLE GIT_ORIGIN)
STRING(FIND "${GIT_ORIGIN}" ".iag.uni-stuttgart.de" HOMERUN)
STRING(FIND "${GIT_ORIGIN}" "git@" USESSH)
IF (HOMERUN GREATER -1)
  IF(USESSH EQUAL 0)
    SET(FLEXI_DLPATH "git@gitlab.iag.uni-stuttgart.de:")
  ELSE()
    SET(FLEXI_DLPATH "https://gitlabext.iag.uni-stuttgart.de/")
  ENDIF()
  SET (HDF5DOWNLOAD ${FLEXI_DLPATH}libs/hdf5.git )
ELSE()
  IF(USESSH EQUAL 0)
    SET(FLEXI_DLPATH "git@github.com:flexi-framework")
  ELSE()
    SET(FLEXI_DLPATH "https://github.com/flexi-framework")
  ENDIF()
  SET (HDF5DOWNLOAD https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git )
ENDIF()
SET (HDF5TAG "hdf5-1_10_0-patch1" CACHE STRING "HDF5 version tag")

# =========================================================================
# Detect machine environements
# =========================================================================
SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "eslogin" ONHORNET)
# Set options that need to be enabled BEFORE ANY PROJECT OR LANGUAGE COMMAND
IF (ONHORNET GREATER -1)
  SET(CMAKE_C_COMPILER       cc)
  SET(CMAKE_CXX_COMPILER     CC)
  SET(CMAKE_Fortran_COMPILER ftn)
  SET(HDF5_USE_STATIC_LIBRARIES TRUE)
  MARK_AS_ADVANCED(C_PATH CXX_PATH Fortran_PATH)
  MESSAGE(STATUS "Compiling on Hazelhen")
  ADD_DEFINITIONS(-DLUSTRE)
ENDIF()

# =========================================================================
# Project definitions and output paths
# =========================================================================
PROJECT(Flexi)
ENABLE_LANGUAGE(Fortran)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Set custom install dir (needs to be done after project has been defined!)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()

# =========================================================================
# Options
# =========================================================================

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()
SUBDIRLIST(TESTCASES  ${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/)
SUBDIRLIST(EQNSYSTEMS ${CMAKE_CURRENT_SOURCE_DIR}/src/equations/)
SUBDIRLIST(LIFTINGS   ${CMAKE_CURRENT_SOURCE_DIR}/src/dg/)

OPTION(FLEXI_BUILDPOSTI "Build POSTI toolset" ON)
OPTION(FLEXI_MPI "Compile SINGLE or MPI version" ON)
OPTION(FLEXI_PARABOLIC "Enable parabolic terms of equation-system" ON)
OPTION(FLEXI_EDDYVISCOSITY "Enable eddy viscosity modells" OFF)
OPTION(FLEXI_MKL "Use MKL libraries instead of BLAS/LAPACK" OFF)
OPTION(FLEXI_PAPI "Use PAPI library to perform performance measurements (e.g. flop counts)." OFF)
OPTION(FLEXI_FV "Enable Finite Volume subcell shock capturing" OFF)
CMAKE_DEPENDENT_OPTION(FLEXI_FV_RECONSTRUCTION "Enable 2nd order reconstruction for FV" ON
                                               "FLEXI_FV" OFF)
SET(   FLEXI_BUILD_HDF5 ON  CACHE BOOL "Compile and build hdf5 library")

SET(   FLEXI_EQNSYSNAME "navierstokes" CACHE STRING
                        "Used equation system (navierstokes/linearscalaradvection)")
SET_PROPERTY(CACHE FLEXI_EQNSYSNAME  PROPERTY STRINGS ${EQNSYSTEMS})

SET(   FLEXI_TESTCASE   "default"      CACHE STRING
                        "Name of testcase which is used. Must be identical to directory name in testcase folder.")
SET_PROPERTY(CACHE FLEXI_TESTCASE PROPERTY STRINGS ${TESTCASES})

SET(FLEXI_POLYNOMIAL_DEGREE "N" CACHE STRING "N: compile with variable polynomial degree. 1,2,3,...: compile with fixed polynomial degree for performance (untested).")
ADD_DEFINITIONS(-DPP_N=${FLEXI_POLYNOMIAL_DEGREE})

##### Node type
SET(FLEXI_NODETYPE "GAUSS" CACHE STRING "Set of nodes used for representing the solution (GAUSS/GAUSS-LOBATTO).")
SET_PROPERTY(CACHE FLEXI_NODETYPE PROPERTY STRINGS "GAUSS" "GAUSS-LOBATTO")
IF(FLEXI_NODETYPE STREQUAL "GAUSS")
  ADD_DEFINITIONS(-DPP_NodeType=1)
ELSEIF(FLEXI_NODETYPE STREQUAL "GAUSS-LOBATTO")
  ADD_DEFINITIONS(-DPP_NodeType=2)
ELSE()
  MESSAGE(ERROR "No such node type known.")
ENDIF()

# Enable or disable unit tests
OPTION(FLEXI_UNITTESTS "Enable unit tests after build" ON)
MARK_AS_ADVANCED(FLEXI_UNITTESTS)

# =========================================================================
# Equation system
# =========================================================================

IF(FLEXI_EQNSYSNAME STREQUAL "navierstokes")
  ADD_DEFINITIONS(-DPP_nVar=5 -DPP_nVarPrim=6 -DEQNSYSNR=2)
  SET(WITHEOS 1)
ELSEIF(FLEXI_EQNSYSNAME STREQUAL "linearscalaradvection")
  ADD_DEFINITIONS(-DPP_nVar=1 -DPP_nVarPrim=1 -DEQNSYSNR=1)
  SET(WITHEOS 0)
ELSE()
  MESSAGE(ERROR "Unknown equation system.")
ENDIF()

IF(WITHEOS GREATER 0)
  SET(FLEXI_EOS "idealgas" CACHE STRING "Used equation of state (idealgas/table)")
  SET_PROPERTY(CACHE FLEXI_EOS PROPERTY STRINGS idealgas)
ENDIF()

IF (FLEXI_PARABOLIC)
  ADD_DEFINITIONS(-DPARABOLIC=1)

  #### Set the lifting
  SET(FLEXI_LIFTING "br1" CACHE STRING "Lifting method used for PARABOLIC")
  SET_PROPERTY(CACHE FLEXI_LIFTING PROPERTY STRINGS ${LIFTINGS})

  IF(FLEXI_LIFTING STREQUAL "br1")
    ADD_DEFINITIONS(-DPP_Lifting=1)
  ELSEIF(FLEXI_LIFTING STREQUAL "br2")
    ADD_DEFINITIONS(-DPP_Lifting=2)
  ELSE()
    MESSAGE(ERROR "No such LIFTING known.")
  ENDIF()

  #### Set the lifting
  IF(FLEXI_EQNSYSNAME STREQUAL "navierstokes")
    SET(FLEXI_VISCOSITY "constant" CACHE STRING "Method to calculate viscosity (constant/sutherland/powerlaw)")
    SET_PROPERTY(CACHE FLEXI_VISCOSITY PROPERTY STRINGS constant sutherland powerlaw)
    IF(FLEXI_VISCOSITY STREQUAL "constant")
      ADD_DEFINITIONS(-DPP_VISC=0)
    ELSEIF(FLEXI_VISCOSITY STREQUAL "sutherland")
      ADD_DEFINITIONS(-DPP_VISC=1)
    ELSEIF(FLEXI_VISCOSITY STREQUAL "powerlaw")
      ADD_DEFINITIONS(-DPP_VISC=2)
    ELSE()
      MESSAGE(ERROR "No such VISCOSITY known.")
    ENDIF()
    IF (FLEXI_EDDYVISCOSITY)
      ADD_DEFINITIONS(-DEDDYVISCOSITY)
    ENDIF()
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DPARABOLIC=0)
ENDIF()

# =========================================================================
# FV subcell shock capturing
# =========================================================================
IF(FLEXI_FV)
  ADD_DEFINITIONS(-DFV_ENABLED=1)
  IF (FLEXI_FV_RECONSTRUCTION)
    ADD_DEFINITIONS(-DFV_RECONSTRUCT=1)
  ELSE()
    ADD_DEFINITIONS(-DFV_RECONSTRUCT=0)
    IF (FLEXI_PARABOLIC)
      MESSAGE(ERROR " FLEXI_PARABOLIC is enabled, please also enabled FLEXI_FV_RECONSTRUCTION or disable FLEXI_PARABOLIC")
    ENDIF()
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DFV_ENABLED=0)
  ADD_DEFINITIONS(-DFV_RECONSTRUCT=0)
ENDIF()

# =========================================================================
# Exclude certain non-compatible combinations from compiling
# =========================================================================
IF(FLEXI_FV AND ${FLEXI_NODETYPE} STREQUAL "GAUSS-LOBATTO")
  MESSAGE(FATAL_ERROR " Gauss-Lobatto nodes not implemented for FV.")
ENDIF()

# =========================================================================
# USERBLOCK + Preproc_flags
# =========================================================================

# A function to get all user defined variables with a specified prefix
function (getListOfVarsStartingWith _prefix _varResult)
    GET_CMAKE_PROPERTY(_vars CACHE_VARIABLES)
    STRING(REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    SET(${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

# A function to get all user defined variables with a specified prefix
SET(configuration ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake)
FILE(WRITE ${configuration} "")
GETLISTOFVARSSTARTINGWITH("FLEXI_" flexiVars)
FOREACH (_var IN LISTS flexiVars)
  GET_PROPERTY(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
  SET(boolian (${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF")) 
  IF(${boolian})
    FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
  ELSE()
    IF(${_var})
      STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
    ENDIF()
    FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
  ENDIF()
ENDFOREACH()



# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release Profile."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_DEFINITIONS("-DDEBUG")
ENDIF()

# =========================================================================
# Location of binary and filenames
# =========================================================================

# append relative filename-macro for __FILENAME__ in Stamp of abort function (see flexi.h)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# add basedir for location of corresponding userblock-file
ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

# =========================================================================
# COMPILER FLAGS
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  #environment specific variables
  IF(ONHORNET GREATER -1)
    SET(FLEXI_INSTRUCTION "-march=core-avx2")
  ELSE()
    SET(FLEXI_INSTRUCTION "-march=native")
  ENDIF()
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -fbounds-check -finit-real=snan -fbacktrace  -Wall")
  # add flags only for compiling not linking!
  SET (FLEXI_COMPILE_FLAGS "-xf95-cpp-input -fPIC")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  #environment specific variables
  IF(ONHORNET GREATER -1)
     SET(FLEXI_INSTRUCTION "-xCORE-AVX2")
  ELSE()
     SET(FLEXI_INSTRUCTION "-xHost")
  ENDIF()
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}    -O3 ${FLEXI_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -p -O3 ${FLEXI_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
  # add flags only for compiling not linking!
  SET (FLEXI_COMPILE_FLAGS "-fpp -allow nofpp_comments -assume bscc")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -f free -s real64 -em -DCRAY")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -p . -rm")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -h profile_generate -p . -rm")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -eD -rm")
  # add flags only for compiling not linking!
  SET (FLEXI_COMPILE_FLAGS "-F")

ELSE()
   MESSAGE(ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# Add the libraries
# =========================================================================

IF(FLEXI_MPI)
  SET(FLEXI_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID}-MPI)
ELSE()
  SET(FLEXI_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID})
ENDIF()
MARK_AS_ADVANCED(FLEXI_EXTERNAL_LIB_DIR)

# =========================================================================
# HDF5 library
# =========================================================================

IF(NOT FLEXI_BUILD_HDF5)
  FIND_PACKAGE(HDF5 COMPONENTS C Fortran QUIET) #if built with configure script
  IF (HDF5_FOUND)
    SET(FLEXI_HDF5_CMAKE FALSE)
  ELSE()
    FIND_PACKAGE(HDF5 COMPONENTS C Fortran static NO_MODULE QUIET PATH_SUFFIXES share/cmake )
    IF (HDF5_FOUND)
      SET(FLEXI_HDF5_CMAKE TRUE)
      MESSAGE (STATUS "HDF5 C libs:${HDF5_FOUND} static:${HDF5_static_C_FOUND} and shared:${HDF5_shared_C_FOUND}")
      MESSAGE (STATUS "HDF5 Fortran libs: static:${HDF5_static_Fortran_FOUND} and shared:${HDF5_shared_Fortran_FOUND}")
      MESSAGE (STATUS "h5diff found:${HDF5_DIFF_EXECUTABLE}")
    ENDIF()
  ENDIF()
ELSE()

  SET(FLEXI_HDF5_DLDIR ${FLEXI_EXTERNAL_LIB_DIR}/HDF5)
  SET(FLEXI_HDF5_DIR  ${FLEXI_HDF5_DLDIR}/build)
  IF (NOT EXISTS "${FLEXI_HDF5_DIR}/lib/libhdf5.a")
    IF(FLEXI_MPI)
      SET(FLEXI_HDF5PARALLEL --enable-parallel)
      SET(FLEXI_HDF5FC ${MPI_Fortran_COMPILER})
      SET(FLEXI_HDF5CC ${MPI_C_COMPILER})
    ELSE()
      UNSET(FLEXI_HDF5PARALLEL)
      SET(FLEXI_HDF5FC ${CMAKE_Fortran_COMPILER})
      SET(FLEXI_HDF5CC ${CMAKE_C_COMPILER} )
    ENDIF()
    #NOTE: old 1.8.X versions of HDF5 require additional flag --enable-fortran2003
    ExternalProject_Add(HDF5
      GIT_REPOSITORY ${HDF5DOWNLOAD}
      GIT_TAG ${HDF5TAG}
      PREFIX ${FLEXI_HDF5_DIR}
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND FC=${FLEXI_HDF5FC} CC=${FLEXI_HDF5CC} ${FLEXI_HDF5_DIR}/src/HDF5/configure --prefix=${FLEXI_HDF5_DIR} --enable-fortran ${FLEXI_HDF5PARALLEL} --disable-shared --libdir=${FLEXI_HDF5_DIR}/lib
      BUILD_COMMAND ${MAKE}
    )
    # CMake HDF5 is fast but not that reliable
    #ExternalProject_Add(HDF5
    #  GIT_REPOSITORY ${HDF5DOWNLOAD}
    #  GIT_TAG ${HDF5TAG}
    #  PREFIX ${FLEXI_HDF5_DIR}
    #  UPDATE_COMMAND ""
    #  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${FLEXI_HDF5_DIR} -DCMAKE_C_COMPILER=${FLEXI_HDF5CC} -DCMAKE_Fortran_COMPILER=${FLEXI_HDF5FC} -DBUILD-TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_PARALLEL=ON
    #  BUILD_COMMAND ${MAKE}
    #)
  ENDIF()

  SET(HDF5_hdf5_LIBRARY_RELEASE         ${FLEXI_HDF5_DIR}/lib/libhdf5.a)
  SET(HDF5_hdf5_fortran_LIBRARY_RELEASE ${FLEXI_HDF5_DIR}/lib/libhdf5_fortran.a)
  SET(HDF5_C_INCLUDE_DIR                ${FLEXI_HDF5_DIR}/include)
  SET(HDF5_Fortran_INCLUDE_DIR          ${FLEXI_HDF5_DIR}/include)
  UNSET(HDF5_LIBRARIES)
  UNSET(HDF5_INCLUDE_DIR_FORTRAN)
  UNSET(HDF5_INCLUDE_DIR)
  UNSET(HDF5_DIFF_EXECUTABLE)
  SET(HDF5_DIFF_EXECUTABLE              ${FLEXI_HDF5_DIR}/bin/h5diff)
  LIST(APPEND HDF5_INCLUDE_DIR_FORTRAN ${HDF5_Fortran_INCLUDE_DIR} ${HDF5_C_INCLUDE_DIR})
  LIST(APPEND HDF5_INCLUDE_DIR  ${HDF5_C_INCLUDE_DIR})
  FIND_LIBRARY(HDF5_z_LIBRARY_RELEASE z)
  LIST(APPEND HDF5_LIBRARIES ${HDF5_hdf5_fortran_LIBRARY_RELEASE} ${HDF5_hdf5_LIBRARY_RELEASE} ${HDF5_z_LIBRARY_RELEASE} -ldl)
ENDIF()

IF(FLEXI_HDF5_CMAKE)
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR_FORTRAN} ${HDF5_INCLUDE_DIR})
  LIST(APPEND linkedlibs ${HDF5_FORTRAN_STATIC_LIBRARY} ${HDF5_LIBRARIES})
ELSE()
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR_FORTRAN} ${HDF5_INCLUDE_DIR})
  LIST(INSERT HDF5_LIBRARIES 0 ${HDF5_hdf5_fortran_LIBRARY_RELEASE})
  LIST(REMOVE_DUPLICATES HDF5_LIBRARIES)
  MESSAGE(STATUS "Found HDF5 Libs: " ${HDF5_LIBRARIES})
  LIST(APPEND linkedlibs ${HDF5_LIBRARIES} )
ENDIF()

# =========================================================================
# LAPACK/MKL library
# =========================================================================
# use Lapack/Blas for GNU
MESSAGE(STATUS "Using BLAS/Lapack library")
IF (FLEXI_MKL)
  SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
ENDIF()
FIND_PACKAGE(LAPACK QUIET)
IF (LAPACK_FOUND)
  LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
ENDIF()

# =========================================================================
# MPI
# =========================================================================
IF(FLEXI_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF (NOT MPI_Fortran_NO_INTERROGATE)
    FOREACH(DIR ${MPI_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    FOREACH(DIR ${MPI_Fortran_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    LIST(APPEND linkedlibs ${MPI_Fortran_LIBRARIES})
  ENDIF()
  MARK_AS_ADVANCED(MPI_LIBRARY MPI_EXTRA_LIBRARY) # don't show these ones
  MESSAGE(STATUS "Building Flexi with MPI")
  MESSAGE(STATUS "MPI Compiler: " ${MPI_Fortran_COMPILER})
  ADD_DEFINITIONS(-DUSE_MPI=1)
ELSE()
  ADD_DEFINITIONS(-DUSE_MPI=0)
ENDIF()

# ========================================================================= 
# PAPI library 
# ========================================================================= 
IF(FLEXI_PAPI)
  FIND_PACKAGE(PAPI REQUIRED)
  ADD_DEFINITIONS(-DPAPI) 
  LIST(APPEND linkedlibs ${PAPI_LIBRARIES})
  INCLUDE_DIRECTORIES(${PAPI_INCLUDE_DIRS})
  MESSAGE(STATUS "Building Flexi with PAPI benchmark support.") 
ENDIF()

# =========================================================================
# Flexi sources & targets
# =========================================================================

function(add_exec target)
  ADD_EXECUTABLE(${target} ${ARGN})
  TARGET_LINK_LIBRARIES(${target} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.o flexilib userblocklib ${linkedlibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}) 
  SET_TARGET_PROPERTIES(${target} PROPERTIES LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
endfunction()

# include flexi.h 
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src) 

# collect all F90 files needed for flexilib
FILE(GLOB_RECURSE flexiF90 ./src/interpolation/*.f90
                           ./src/readintools/*.f90
                           ./src/mortar/*.f90
                           ./src/mesh/*.f90
                           ./src/indicator/*.f90
                           ./src/mpi/*.f90
                           ./src/globals/*.f90
                           ./src/analyze/*.f90
                           ./src/recordpoints/*.f90
                           ./src/restart/*.f90
                           ./src/output/*.f90
                           ./src/sponge/*.f90
                           ./src/filter/*.f90
                           ./src/io_hdf5/*.f90
                           ./src/timedisc/*.f90
                           ./src/testcase/${FLEXI_TESTCASE}/*.f90
                           ./unitTests/generateunittestreferencedata.f90)

FILE(GLOB eqnsysF90        ./src/equations/${FLEXI_EQNSYSNAME}/*.f90)
IF (WITHEOS GREATER 0)
  FILE(GLOB eosF90         ./src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/*.f90)
  LIST(APPEND eqnsysF90 ${eosF90})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/) 
ENDIF()
IF (FLEXI_EDDYVISCOSITY)
  FILE(GLOB_RECURSE eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/*.f90)
  LIST(APPEND eqnsysF90 ${eddyviscF90})
ENDIF()

FILE(GLOB dgF90            ./src/dg/*.f90)
IF (FLEXI_PARABOLIC)
  FILE(GLOB liftingF90     ./src/dg/${FLEXI_LIFTING}/*.f90)
  LIST(APPEND dgF90 ${liftingF90})
ENDIF()
IF (FLEXI_FV)
  FILE(GLOB fvF90 ./src/fv/*.f90)
  LIST(APPEND flexiF90 ${fvF90})
ENDIF()


# define libs
ADD_LIBRARY(flexilibF90 OBJECT ${flexiF90} ${dgF90} ${eqnsysF90})
SET_TARGET_PROPERTIES(flexilibF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}) 
ADD_LIBRARY(flexilib STATIC $<TARGET_OBJECTS:flexilibF90>)
ADD_DEPENDENCIES(flexilibF90 HDF5)

ADD_LIBRARY(userblocklib STATIC ./src/output/read_userblock.c)
SET_TARGET_PROPERTIES(userblocklib PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(userblocklib PROPERTIES COMPILE_FLAGS "-fPIC")

add_exec(flexi ./src/flexi.f90)
ADD_DEPENDENCIES(flexi flexilib userblocklib HDF5)

# special compiler flags to avoid warnings for empty functions of testcase.f90 and overwrite of intrinsic abort in globals.f90
IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/${FLEXI_TESTCASE}/testcase.f90 PROPERTIES COMPILE_FLAGS "-Wno-unused-dummy-argument")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals.f90 PROPERTIES COMPILE_FLAGS "-Wno-intrinsic-shadow")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/${FLEXI_TESTCASE}/testcase.f90 PROPERTIES COMPILE_FLAGS "-warn nounused")
ENDIF()

# -------------------------------------------------------------------------
# echo COMPILE_DEFINITIONS
# -------------------------------------------------------------------------
GET_DIRECTORY_PROPERTY(comp_defs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
SET(comp_defs_str "COMPILE_DEFINITIONS = ")
FOREACH (def ${comp_defs})
   SET(comp_defs_str "${comp_defs_str} -D${def}")
ENDFOREACH()
ADD_CUSTOM_TARGET(preproc_defines COMMAND echo ${comp_defs_str})

# =========================================================================
# Add FLEXI2VTK tool
# =========================================================================
add_exec(flexi2vtk ./src/flexi2vtk.f90)
ADD_DEPENDENCIES(flexi2vtk flexilib HDF5)

# =========================================================================
# Posti
# =========================================================================
IF(FLEXI_BUILDPOSTI)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/CMakeLists.txt)
ENDIF()

# =========================================================================
# RegressionCheck tool && Examples
# =========================================================================
add_exec(regressioncheck ./regressioncheck/regressioncheck.f90 
                         ./regressioncheck/regressioncheck_vars.f90 
                         ./regressioncheck/regressioncheck_build.f90 
                         ./regressioncheck/regressioncheck_compare.f90 
                         ./regressioncheck/regressioncheck_run.f90 
                         ./regressioncheck/regressioncheck_tools.f90 )
ADD_DEFINITIONS(-DH5DIFF="${HDF5_DIFF_EXECUTABLE}")
SET(FLEXI_REGRESSIONCHECK_EXAMPLES "${CMAKE_CURRENT_SOURCE_DIR}/regressioncheck/examples")

# =========================================================================
# Userblock
# =========================================================================
ADD_CUSTOM_COMMAND(TARGET flexilib PRE_BUILD COMMAND 
  sh ${CMAKE_CURRENT_SOURCE_DIR}/tools/userblock/generateuserblock.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# =========================================================================
# Install
# =========================================================================
INSTALL(FILES    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt DESTINATION bin)

# =========================================================================
# Ctags
# =========================================================================
FIND_PROGRAM(CTAGS_PATH ctags)
IF(CTAGS_PATH)
  SET(EQNS_WO_ACTUAL ${EQNSYSTEMS})
  LIST(REMOVE_ITEM EQNS_WO_ACTUAL ${FLEXI_EQNSYSNAME})
  SET(EXLUDES "")
  FOREACH(item ${EQNS_WO_ACTUAL})
     SET(EXLUDES "--exclude=src/equations/${item}/*" ${EXLUDES})
  ENDFOREACH()
  SET(LIFTINGS_WO_ACTUAL ${LIFTINGS})
  LIST(REMOVE_ITEM LIFTINGS_WO_ACTUAL ${FLEXI_LIFTING})
  FOREACH(item ${LIFTINGS_WO_ACTUAL})
     SET(EXLUDES "--exclude=src/dg/${item}/*" ${EXLUDES})
  ENDFOREACH()
  SET(TESTCASES_WO_ACTUAL ${TESTCASES})
  LIST(REMOVE_ITEM TESTCASES_WO_ACTUAL ${FLEXI_TESTCASE})
  FOREACH(item ${TESTCASES_WO_ACTUAL})
     SET(EXLUDES "--exclude=src/testcase/${item}/*" ${EXLUDES})
  ENDFOREACH()
  ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ctags -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} src && ctags -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} --append posti)
  ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
  MESSAGE(STATUS "Found ctags: ${CTAGS_PATH}")
  ADD_DEPENDENCIES(flexi tags)
ENDIF(CTAGS_PATH)

# =========================================================================
# Tests
# =========================================================================

# Deactivate tests on hornet since no programs are allowed to run on frontend
IF(ONHORNET GREATER -1)
  SET(FLEXI_UNITTESTS OFF CACHE BOOL "Enable unit tests after build" FORCE)
ENDIF()
IF(FLEXI_UNITTESTS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/unitTests/CMakeLists.txt)
ENDIF()

