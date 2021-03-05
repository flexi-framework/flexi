# =========================================================================
# FLEXI Options
# =========================================================================

# Build directory list
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
SUBDIRLIST(EQNSYSTEMS ${CMAKE_CURRENT_SOURCE_DIR}/src/equations/)
SUBDIRLIST(LIFTINGS   ${CMAKE_CURRENT_SOURCE_DIR}/src/dg/)
SUBDIRLIST(TESTCASES  ${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/)

# =========================================================================
# Node type
# =========================================================================
SET(FLEXI_NODETYPE "GAUSS" CACHE STRING "Set of nodes used for representing the solution (GAUSS/GAUSS-LOBATTO).")
SET_PROPERTY(CACHE FLEXI_NODETYPE PROPERTY STRINGS "GAUSS" "GAUSS-LOBATTO")
IF(FLEXI_NODETYPE STREQUAL "GAUSS")
  ADD_DEFINITIONS(-DPP_NodeType=1)
ELSEIF(FLEXI_NODETYPE STREQUAL "GAUSS-LOBATTO")
  ADD_DEFINITIONS(-DPP_NodeType=2)
ELSE()
  MESSAGE(SEND_ERROR "No such node type known.")
ENDIF()

# =========================================================================
# SplitDG
# =========================================================================
OPTION(FLEXI_SPLIT_DG "Use split DG formulation" OFF)
IF(FLEXI_SPLIT_DG)
  ADD_DEFINITIONS(-DSPLIT_DG)
  IF(FLEXI_NODETYPE STREQUAL "GAUSS")
    MESSAGE(SEND_ERROR "FLEXI_SPLIT_DG only works on FLEXI_NODETYPE==GAUSS-LOBATTO points.")
  ENDIF()
  IF(NOT FLEXI_EQNSYSNAME STREQUAL "navierstokes")
    MESSAGE(SEND_ERROR "FLEXI_SPLIT_DG only works for FLEXI_EQNSYSNAME==navierstokes.")
  ENDIF()
ENDIF()

# =========================================================================
# Equation system
# =========================================================================
SET(FLEXI_EQNSYSNAME "navierstokes" CACHE STRING
                     "Used equation system (navierstokes/linearscalaradvection/rans_sa)")
SET_PROPERTY(CACHE FLEXI_EQNSYSNAME  PROPERTY STRINGS ${EQNSYSTEMS})
IF(   FLEXI_EQNSYSNAME STREQUAL "navierstokes")
  ADD_DEFINITIONS(-DPP_nVar=5 -DPP_nVarPrim=6 -DEQNSYSNR=2)
  SET(WITHEOS 1)
ELSEIF(FLEXI_EQNSYSNAME STREQUAL "linearscalaradvection")
  ADD_DEFINITIONS(-DPP_nVar=1 -DPP_nVarPrim=1 -DEQNSYSNR=1)
  SET(WITHEOS 0)
ELSEIF(FLEXI_EQNSYSNAME STREQUAL "rans_sa")
  ADD_DEFINITIONS(-DPP_nVar=6 -DPP_nVarPrim=7 -DEQNSYSNR=3)
  SET(WITHEOS 1)
ELSE()
  MESSAGE(SEND_ERROR "Unknown equation system.")
ENDIF()

IF(WITHEOS GREATER 0)
  SET(FLEXI_EOS "idealgas" CACHE STRING "Used equation of state (idealgas/table)")
  SET_PROPERTY(CACHE FLEXI_EOS PROPERTY STRINGS idealgas)
  MARK_AS_ADVANCED(FORCE FLEXI_EOS)
ENDIF()

# =========================================================================
# Polynomial degree
# =========================================================================
SET(FLEXI_POLYNOMIAL_DEGREE "N" CACHE STRING "N: compile with variable polynomial degree. 1,2,3,...: compile with fixed polynomial degree for performance (untested).")
ADD_DEFINITIONS(-DPP_N=${FLEXI_POLYNOMIAL_DEGREE})

# =========================================================================
# Parabolic terms
# =========================================================================
OPTION(FLEXI_PARABOLIC "Enable parabolic terms of equation-system" ON)
IF (FLEXI_PARABOLIC)
  ADD_DEFINITIONS(-DPARABOLIC=1)

  # Set the lifting
  LIST(APPEND FLEXI_LIFTING_LIST  "br1"
                                  "br2")
  SET(FLEXI_LIFTING "br1" CACHE STRING "Lifting method used for PARABOLIC")
  SET_PROPERTY(CACHE FLEXI_LIFTING PROPERTY STRINGS ${FLEXI_LIFTING_LIST})

  IF(    FLEXI_LIFTING STREQUAL "br1")
    ADD_DEFINITIONS(-DPP_Lifting=1)
  ELSEIF(FLEXI_LIFTING STREQUAL "br2")
    ADD_DEFINITIONS(-DPP_Lifting=2)
  ELSE()
    MESSAGE(SEND_ERROR "No such LIFTING known.")
  ENDIF()

  # Set the law for viscosity
  IF((FLEXI_EQNSYSNAME STREQUAL "navierstokes") OR (FLEXI_EQNSYSNAME STREQUAL "rans_sa"))
    SET(FLEXI_VISCOSITY "constant" CACHE STRING "Method to calculate viscosity (constant/sutherland/powerlaw)")
    SET_PROPERTY(CACHE FLEXI_VISCOSITY PROPERTY STRINGS constant sutherland powerlaw)
    IF(    FLEXI_VISCOSITY STREQUAL "constant")
      ADD_DEFINITIONS(-DPP_VISC=0)
    ELSEIF(FLEXI_VISCOSITY STREQUAL "sutherland")
      ADD_DEFINITIONS(-DPP_VISC=1)
    ELSEIF(FLEXI_VISCOSITY STREQUAL "powerlaw")
      ADD_DEFINITIONS(-DPP_VISC=2)
    ELSE()
      MESSAGE(SEND_ERROR "No such VISCOSITY known.")
    ENDIF()
    IF (FLEXI_EDDYVISCOSITY)
      IF(FLEXI_2D)
        MESSAGE(SEND_ERROR " FLEXI_EDDYVISCOSITY is enabled, please disable FLEXI_2D")
      ENDIF()
    ENDIF()
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DPARABOLIC=0)
ENDIF()

# =========================================================================
# Eddy viscosity models
# =========================================================================
CMAKE_DEPENDENT_OPTION(FLEXI_EDDYVISCOSITY "Enable eddy viscosity modells" OFF
                                           "FLEXI_PARABOLIC"               OFF)
IF (FLEXI_EDDYVISCOSITY)
  ADD_DEFINITIONS(-DEDDYVISCOSITY=1)
ELSE()
  ADD_DEFINITIONS(-DEDDYVISCOSITY=0)
ENDIF()

# =========================================================================
# Dimension related preproc variables
# =========================================================================
OPTION(FLEXI_2D "Enable two dimensional computation" OFF)
IF(FLEXI_2D)
  ADD_DEFINITIONS(-DPP_dim=2)
ELSE()
  ADD_DEFINITIONS(-DPP_dim=3)
ENDIF()

# =========================================================================
# FV subcell shock capturing
# =========================================================================
OPTION(FLEXI_FV "Enable Finite Volume subcell shock capturing" OFF)
CMAKE_DEPENDENT_OPTION(FLEXI_FV_RECONSTRUCTION "Enable 2nd order reconstruction for FV" ON
                                               "FLEXI_FV" OFF)
IF(FLEXI_FV)
  ADD_DEFINITIONS(-DFV_ENABLED=1)
  IF (FLEXI_FV_RECONSTRUCTION)
    ADD_DEFINITIONS(-DFV_RECONSTRUCT=1)
  ELSE()
    ADD_DEFINITIONS(-DFV_RECONSTRUCT=0)
    IF (FLEXI_PARABOLIC)
      MESSAGE(SEND_ERROR " FLEXI_PARABOLIC is enabled, please also enabled FLEXI_FV_RECONSTRUCTION or disable FLEXI_PARABOLIC")
    ENDIF()
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DFV_ENABLED=0)
  ADD_DEFINITIONS(-DFV_RECONSTRUCT=0)
ENDIF()

# =========================================================================
# Test case
# =========================================================================
SET(   FLEXI_TESTCASE   "default"      CACHE STRING
                        "Name of testcase which is used. Must be identical to directory name in testcase folder.")
SET_PROPERTY(CACHE FLEXI_TESTCASE PROPERTY STRINGS ${TESTCASES})

# =========================================================================
# Unit test at end of compilation
# =========================================================================
OPTION(FLEXI_UNITTESTS "Enable unit tests after build" ON)
MARK_AS_ADVANCED(FORCE FLEXI_UNITTESTS)
IF(FLEXI_UNITTESTS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/unitTests/CMakeLists.txt)
ENDIF()

# =========================================================================
# Flexi sources & targets
# =========================================================================
FUNCTION(add_exec target)
  ADD_EXECUTABLE(       ${target} ${ARGN})
  TARGET_LINK_LIBRARIES(${target} libflexishared)
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                             LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
ENDFUNCTION()

FUNCTION(add_lib_shared target)
  ADD_LIBRARY(          ${target} SHARED ${ARGN})
  TARGET_LINK_LIBRARIES(${target} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.o libflexistatic userblocklib stacksizelib ${linkedlibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES PREFIX ""
                                             COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                             LINKER_LANGUAGE Fortran)
  INSTALL(TARGETS ${target} LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                                                GROUP_READ GROUP_EXECUTE
                                                WORLD_READ WORLD_EXECUTE
                            DESTINATION lib)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
ENDFUNCTION()

# include flexi.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)

# collect all F90 files needed for flexilib, except the extra files below
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
                           ./src/flexilib.f90
                           ./unitTests/unittest.f90
                           ./unitTests/unittest_vars.f90)

# The following files/directories where excluded from the above global recursive search,
# since some of the subdirectories should only be included if the respective feature is enabled.

# >> collect F90 files of the equation system (non-recursive) in a new list 'eqnsysF90'
FILE(GLOB eqnsysF90        ./src/equations/${FLEXI_EQNSYSNAME}/*.f90)
# append F90 files of the equation of state to 'eqnsysF90'-list
IF (WITHEOS GREATER 0)
  FILE(GLOB eosF90         ./src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/*.f90)
  LIST(APPEND eqnsysF90 ${eosF90})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/)
ELSE()
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/equations/${FLEXI_EQNSYSNAME}/) 
ENDIF()

# >> append F90 files of the eddy viscosity to 'eqnsysF90'-list
IF (FLEXI_EDDYVISCOSITY)
  FILE(GLOB eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/default/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/smagorinsky/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/sigmaModel/*.f90)
  LIST(APPEND eqnsysF90 ${eddyviscF90})
ENDIF()

# >> append 'eqnsysF90'-list to 'flexiF90'-list
LIST(APPEND flexiF90 ${eqnsysF90})

# >> collect F90 files of 'dg' directory (non-recursive) in a new list 'dgF90'
FILE(GLOB dgF90            ./src/dg/*.f90)

# >> append F90 files of lifting to 'dgF90'-list
IF (FLEXI_PARABOLIC)
  FILE(GLOB liftingF90     ./src/dg/lifting/*.f90)
  LIST(APPEND dgF90 ${liftingF90})
ENDIF()

# >> append 'dgF90'-list to 'flexiF90'-list
LIST(APPEND flexiF90 ${dgF90})

# >> collect F90 files of 'fv' directory (non-recursive) and append to 'flexiF90' list
IF (FLEXI_FV)
  FILE(GLOB fvF90 ./src/fv/*.f90)
  LIST(APPEND flexiF90 ${fvF90})
ENDIF()

# define libs
ADD_LIBRARY(          libflexistaticF90 OBJECT ${flexiF90})
SET_TARGET_PROPERTIES(libflexistaticF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(          libflexistatic STATIC $<TARGET_OBJECTS:libflexistaticF90>)
SET_TARGET_PROPERTIES(libflexistatic PROPERTIES PREFIX ""
                                     OUTPUT_NAME "libflexif90")
IF(SELFBUILTEXTERNALS)
  ADD_DEPENDENCIES(libflexistaticF90 ${SELFBUILTEXTERNALS})
ENDIF()

ADD_LIBRARY(userblocklib STATIC ./src/output/read_userblock.c)
ADD_LIBRARY(stacksizelib STATIC ./src/globals/stacksize.c)
SET_TARGET_PROPERTIES(stacksizelib userblocklib PROPERTIES LINKER_LANGUAGE C
                                                           COMPILE_FLAGS "-fPIC")

# we link the statically built libs
ADD_LIB_SHARED(       libflexishared ./src/flexilib.f90)
SET_TARGET_PROPERTIES(libflexishared PROPERTIES OUTPUT_NAME "libflexi")
ADD_DEPENDENCIES(     libflexishared libflexistatic userblocklib stacksizelib ${INTERNALLIBS})
ADD_EXEC(        flexi ./src/flexi.f90)
ADD_DEPENDENCIES(flexi libflexistatic userblocklib stacksizelib ${INTERNALLIBS})

# special compiler flags to avoid warnings for empty functions of testcase.f90 and overwrite of intrinsic abort in globals.f90
IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/${FLEXI_TESTCASE}/testcase.f90 PROPERTIES COMPILE_FLAGS "-Wno-unused-dummy-argument")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals.f90 PROPERTIES COMPILE_FLAGS "-Wno-intrinsic-shadow")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/${FLEXI_TESTCASE}/testcase.f90 PROPERTIES COMPILE_FLAGS "-warn nounused")
ENDIF()

# =========================================================================
# Ctags
# =========================================================================
MARK_AS_ADVANCED(FORCE CTAGS_PATH)
FIND_PROGRAM(CTAGS_PATH ctags-exuberant)
IF(CTAGS_PATH)
  SET(             EQNS_WO_ACTUAL ${EQNSYSTEMS})
  LIST(REMOVE_ITEM EQNS_WO_ACTUAL ${FLEXI_EQNSYSNAME})
  SET(EXLUDES "")
  FOREACH(item ${EQNS_WO_ACTUAL})
     SET(EXLUDES "--exclude=src/equations/${item}/*" ${EXLUDES})
  ENDFOREACH()

  IF(FLEXI_PARABOLIC)
    SET(             LIFTINGS_WO_ACTUAL ${LIFTINGS})
    LIST(REMOVE_ITEM LIFTINGS_WO_ACTUAL ${FLEXI_LIFTING})
    FOREACH(item ${LIFTINGS_WO_ACTUAL})
       SET(EXLUDES "--exclude=src/dg/${item}/*" ${EXLUDES})
    ENDFOREACH()
  ENDIF()

  SET(             TESTCASES_WO_ACTUAL ${TESTCASES})
  LIST(REMOVE_ITEM TESTCASES_WO_ACTUAL ${FLEXI_TESTCASE})
  FOREACH(item ${TESTCASES_WO_ACTUAL})
     SET(EXLUDES "--exclude=src/testcase/${item}/*" ${EXLUDES})
  ENDFOREACH()

  ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ${CTAGS_PATH} -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} src && ${CTAGS_PATH} -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} --append posti)
  ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
  MESSAGE(STATUS "Found ctags: ${CTAGS_PATH}")
  ADD_DEPENDENCIES(flexi tags)
ENDIF(CTAGS_PATH)
