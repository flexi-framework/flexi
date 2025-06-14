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
  ADD_COMPILE_DEFINITIONS(PP_NodeType=1)
ELSEIF(FLEXI_NODETYPE STREQUAL "GAUSS-LOBATTO")
  ADD_COMPILE_DEFINITIONS(PP_NodeType=2)
ELSE()
  MESSAGE(SEND_ERROR "No such node type known.")
ENDIF()

# =========================================================================
# Equation system
# =========================================================================
SET(FLEXI_EQNSYSNAME "navierstokes" CACHE STRING
                     "Used equation system (navierstokes/linearscalaradvection/rans_sa)")
SET_PROPERTY(CACHE FLEXI_EQNSYSNAME  PROPERTY STRINGS ${EQNSYSTEMS})
IF(   FLEXI_EQNSYSNAME STREQUAL "navierstokes")
  ADD_COMPILE_DEFINITIONS(PP_nVar=5 PP_nVarPrim=6 EQNSYSNR=2)
  SET(WITHEOS 1)
ELSEIF(FLEXI_EQNSYSNAME STREQUAL "linearscalaradvection")
  ADD_COMPILE_DEFINITIONS(PP_nVar=1 PP_nVarPrim=1 EQNSYSNR=1)
  SET(WITHEOS 0)
ELSEIF(FLEXI_EQNSYSNAME STREQUAL "rans_sa")
  ADD_COMPILE_DEFINITIONS(PP_nVar=6 PP_nVarPrim=7 EQNSYSNR=3)
  SET(WITHEOS 1)
ELSEIF(FLEXI_EQNSYSNAME STREQUAL "burgers")
  ADD_DEFINITIONS(-DPP_nVar=3 -DPP_nVarPrim=3 -DEQNSYSNR=4)
  SET(WITHEOS 0)
ELSE()
  MESSAGE(SEND_ERROR "Unknown equation system.")
ENDIF()

IF(WITHEOS GREATER 0)
  SET(FLEXI_EOS "idealgas" CACHE STRING "Used equation of state (idealgas/table)")
  SET_PROPERTY(CACHE FLEXI_EOS PROPERTY STRINGS idealgas)
  MARK_AS_ADVANCED(FORCE FLEXI_EOS)
ENDIF()

# =========================================================================
# SplitDG
# =========================================================================
OPTION(FLEXI_SPLIT_DG "Use split DG formulation" OFF)
IF(FLEXI_SPLIT_DG)
  ADD_COMPILE_DEFINITIONS(SPLIT_DG)
  IF((FLEXI_EQNSYSNAME STREQUAL "navierstokes") AND (FLEXI_NODETYPE STREQUAL "GAUSS"))
    ADD_COMPILE_DEFINITIONS(PP_EntropyVars=1)
  ELSE()
    ADD_COMPILE_DEFINITIONS(PP_EntropyVars=0)
  ENDIF()
  IF(NOT FLEXI_EQNSYSNAME STREQUAL "navierstokes")
    MESSAGE(SEND_ERROR "FLEXI_SPLIT_DG only works for FLEXI_EQNSYSNAME==navierstokes.")
  ENDIF()
ELSE()
  ADD_COMPILE_DEFINITIONS(PP_EntropyVars=0)
ENDIF()

# =========================================================================
# PPLimiter
# =========================================================================
OPTION(FLEXI_PP_LIMITER "Enable Positivity Preserving Limiter" OFF)
IF(FLEXI_PP_LIMITER)
  ADD_COMPILE_DEFINITIONS(PP_LIMITER)
  IF(FLEXI_NODETYPE STREQUAL "GAUSS")
    MESSAGE(WARNING "Positivity Preserving Limiter only works on volume interpolation points. Positivity on faces is not guaranteed with FLEXI_NODETYPE==GAUSS." )
  ENDIF()
  IF(NOT FLEXI_EQNSYSNAME STREQUAL "navierstokes")
    MESSAGE(SEND_ERROR "FLEXI_PP_LIMITER only works for FLEXI_EQNSYSNAME==navierstokes.")
  ENDIF()
ENDIF()

# =========================================================================
# ExactMM
# =========================================================================
# exact mass matrix
OPTION(FLEXI_EXACT_MASSMATRIX "Build and invert exact mass Matrix" OFF)
IF(FLEXI_EXACT_MASSMATRIX)
  ADD_COMPILE_DEFINITIONS(EXACT_MM)
  IF(FLEXI_NODETYPE STREQUAL "GAUSS")
    MESSAGE(SEND_ERROR "FLEXI_EXACT_MASSMATRIX only works on FLEXI_NODETYPE==GAUSS-LOBATTO points.")
  ENDIF()
ENDIF()

# =========================================================================
# Polynomial degree
# =========================================================================
SET(FLEXI_POLYNOMIAL_DEGREE "N" CACHE STRING "N: compile with variable polynomial degree. 1,2,3,...: compile with fixed polynomial degree for performance (untested).")
ADD_COMPILE_DEFINITIONS(PP_N=${FLEXI_POLYNOMIAL_DEGREE})

# =========================================================================
# Parabolic terms
# =========================================================================
OPTION(FLEXI_PARABOLIC "Enable parabolic terms of equation-system" ON)
IF (FLEXI_PARABOLIC)
  ADD_COMPILE_DEFINITIONS(PARABOLIC=1)

  # Set the lifting
  LIST(APPEND FLEXI_LIFTING_LIST  "br1"
                                  "br2")
  SET(FLEXI_LIFTING "br1" CACHE STRING "Lifting method used for PARABOLIC")
  SET_PROPERTY(CACHE FLEXI_LIFTING PROPERTY STRINGS ${FLEXI_LIFTING_LIST})

  IF(    FLEXI_LIFTING STREQUAL "br1")
    ADD_COMPILE_DEFINITIONS(PP_Lifting=1)
  ELSEIF(FLEXI_LIFTING STREQUAL "br2")
    ADD_COMPILE_DEFINITIONS(PP_Lifting=2)
  ELSE()
    MESSAGE(SEND_ERROR "No such LIFTING known.")
  ENDIF()

  # Set the law for viscosity
  IF((FLEXI_EQNSYSNAME STREQUAL "navierstokes") OR (FLEXI_EQNSYSNAME STREQUAL "rans_sa"))
    SET(FLEXI_VISCOSITY "constant" CACHE STRING "Method to calculate viscosity (constant/sutherland/powerlaw)")
    SET_PROPERTY(CACHE FLEXI_VISCOSITY PROPERTY STRINGS constant sutherland powerlaw)
    IF(    FLEXI_VISCOSITY STREQUAL "constant")
      ADD_COMPILE_DEFINITIONS(PP_VISC=0)
    ELSEIF(FLEXI_VISCOSITY STREQUAL "sutherland")
      ADD_COMPILE_DEFINITIONS(PP_VISC=1)
    ELSEIF(FLEXI_VISCOSITY STREQUAL "powerlaw")
      ADD_COMPILE_DEFINITIONS(PP_VISC=2)
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
  ADD_COMPILE_DEFINITIONS(PARABOLIC=0)
ENDIF()

# =========================================================================
# Eddy viscosity models
# =========================================================================
CMAKE_DEPENDENT_OPTION(FLEXI_EDDYVISCOSITY "Enable eddy viscosity models" OFF
                                           "FLEXI_PARABOLIC"              OFF)
IF (FLEXI_EDDYVISCOSITY)
  ADD_COMPILE_DEFINITIONS(EDDYVISCOSITY=1)
ELSE()
  ADD_COMPILE_DEFINITIONS(EDDYVISCOSITY=0)
ENDIF()

# =========================================================================
# Dimension related preproc variables
# =========================================================================
OPTION(FLEXI_2D "Enable two dimensional computation" OFF)
IF(FLEXI_2D)
  ADD_COMPILE_DEFINITIONS(PP_dim=2)
ELSE()
  ADD_COMPILE_DEFINITIONS(PP_dim=3)
ENDIF()

# =========================================================================
# FV subcell shock capturing
# SWITCH: switching
# BLEND/SBLEND: element-wise/DOF-wise blending
# =========================================================================
SET(FLEXI_FV "OFF" CACHE STRING "Enable Finite Volume subcell shock capturing")
SET_PROPERTY(CACHE FLEXI_FV PROPERTY STRINGS OFF SWITCH BLEND SBLEND)

CMAKE_DEPENDENT_OPTION(FLEXI_FV_RECONSTRUCTION "Enable 2nd order reconstruction for FV" ON
                                               "FLEXI_FV" ON) # Everything except "OFF" is logically "ON"

IF(NOT (FLEXI_FV STREQUAL "OFF"))
  IF(   FLEXI_FV STREQUAL "SWITCH")
    ADD_COMPILE_DEFINITIONS(FV_ENABLED=1)
  ELSEIF(FLEXI_FV STREQUAL "BLEND")
    IF(NOT FLEXI_EQNSYSNAME STREQUAL "navierstokes")
      MESSAGE(SEND_ERROR "FLEXI_FV=BLEND only works for FLEXI_EQNSYSNAME==navierstokes.")
    ELSEIF(FLEXI_NODETYPE STREQUAL "GAUSS")
      MESSAGE(SEND_ERROR "FLEXI_FV=BLEND only works on FLEXI_NODETYPE==GAUSS-LOBATTO points.")
    ENDIF()
    ADD_COMPILE_DEFINITIONS(FV_ENABLED=2)
  ELSEIF(FLEXI_FV STREQUAL "SBLEND")
    IF(NOT FLEXI_EQNSYSNAME STREQUAL "navierstokes")
      MESSAGE(SEND_ERROR "FLEXI_FV=SBLEND only works for FLEXI_EQNSYSNAME==navierstokes.")
    ELSEIF(FLEXI_NODETYPE STREQUAL "GAUSS")
      MESSAGE(SEND_ERROR "FLEXI_FV=SBLEND only works on FLEXI_NODETYPE==GAUSS-LOBATTO points.")
    ENDIF()
    ADD_COMPILE_DEFINITIONS(FV_ENABLED=3)
  ELSE()
    MESSAGE(SEND_ERROR "No such FV shock capturing known.")
  ENDIF()

  IF (FLEXI_FV_RECONSTRUCTION)
    ADD_COMPILE_DEFINITIONS(FV_RECONSTRUCT=1)
  ELSE()
    ADD_COMPILE_DEFINITIONS(FV_RECONSTRUCT=0)
    IF (FLEXI_PARABOLIC AND (NOT FLEXI_FV STREQUAL "BLEND"))
      MESSAGE(SEND_ERROR " FLEXI_PARABOLIC is enabled, please also enabled FLEXI_FV_RECONSTRUCTION or disable FLEXI_PARABOLIC")
    ENDIF()
  ENDIF()
ELSE()
  ADD_COMPILE_DEFINITIONS(FV_ENABLED=0)
  ADD_COMPILE_DEFINITIONS(FV_RECONSTRUCT=0)
ENDIF()

# =========================================================================
# Reduced lifting in performance mode
# =========================================================================
CMAKE_DEPENDENT_OPTION(FLEXI_PERFORMANCE_OPTLIFT "Lift only needed variables" OFF
                                                 "FLEXI_PERFORMANCE;NOT POSTI" OFF)
IF(FLEXI_PERFORMANCE_OPTLIFT)
  ADD_COMPILE_DEFINITIONS(PP_OPTLIFT=1)
ELSE()
  ADD_COMPILE_DEFINITIONS(PP_OPTLIFT=0)
ENDIF()

# =========================================================================
# Preconditioner for implicit
# =========================================================================
OPTION(FLEXI_USE_PRECONDITIONER "Use preconditioner for implicit" OFF)
IF(FLEXI_USE_PRECONDITIONER)
  IF (NOT (FLEXI_FV STREQUAL "OFF"))
    IF (FLEXI_FV_RECONSTRUCTION AND FLEXI_PARABOLIC)
      MESSAGE(SEND_ERROR "FLEXI_USE_PRECONDITIONER not supported for FLEXI_PARABOLIC and FLEXI_FV_RECONSTRUCTION enabled!")
    ELSE()
      ADD_COMPILE_DEFINITIONS(USE_PRECOND=1)
    ENDIF()
  ELSE()
    ADD_COMPILE_DEFINITIONS(USE_PRECOND=1)
  ENDIF()
ELSE()
  ADD_COMPILE_DEFINITIONS(USE_PRECOND=0)
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
# Check is unit tests can be activated
IF(BUILD_TYPE_LC MATCHES "sanitize" AND LIBS_USE_MPI OR BUILD_TYPE_LC MATCHES "profile")
  SET(ACTIVATE_UNIT_TEST       OFF)
  # Store if we force-disabled the unit tests
  SET(UNIT_TEST_FORCE_DISABLED ON  CACHE BOOL "Force disable unit tests" FORCE)
  # Check if unit tests were previously activated
  IF(FLEXI_UNITTESTS)
    IF (BUILD_TYPE_LC MATCHES "profile")
      MESSAGE(WARNING "Unit tests automatically disabled (FLEXI_UNITTESTS=OFF) with\n CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}\nas the tests will always fail due to memory errors returned from the profiling library.")
    ELSEIF (LIBS_MPI_NAME MATCHES "OpenMPI")
      MESSAGE(WARNING "Unit tests automatically disabled (FLEXI_UNITTESTS=OFF) with\n CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}\n LIBS_USE_MPI = ${LIBS_USE_MPI}\nas the tests will always fail due to memory errors returned from the MPI library.")
    ENDIF()
    SET(FLEXI_UNITTESTS OFF)
  ENDIF()
ELSE()
  SET(ACTIVATE_UNIT_TEST ON)
ENDIF()
MARK_AS_ADVANCED(FORCE UNIT_TEST_FORCE_DISABLED)

# Enable unit testing only when ACTIVATE_UNIT_TEST=ON
CMAKE_DEPENDENT_OPTION(FLEXI_UNITTESTS "Enable unit tests after build" ON "ACTIVATE_UNIT_TEST" OFF)
# Re-enable if previously force disabled
IF (ACTIVATE_UNIT_TEST AND UNIT_TEST_FORCE_DISABLED)
  SET(FLEXI_UNITTESTS ON)
  SET(UNIT_TEST_FORCE_DISABLED OFF CACHE BOOL "Force disable unit tests" FORCE)
ENDIF()

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
                                             COMPILE_OPTIONS "-fPIC"
                                             CMAKE_POSITION_INDEPENDENT_CODE ON
                                             LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
  # Strip executable library when building in release mode, except on macOS
  # IF(NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
  IF(NOT CMAKE_HOST_APPLE)
    ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND $<$<CONFIG:release>:${CMAKE_STRIP}> ARGS --strip-unneeded $<TARGET_FILE:${target}>)
  ENDIF()
ENDFUNCTION()

FUNCTION(add_lib_shared target)
  ADD_LIBRARY(          ${target} SHARED ${ARGN})
  TARGET_LINK_LIBRARIES(${target} libflexistatic userblocklib stacksizelib ${linkedlibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES PREFIX ""
                                             COMPILE_FLAGS   ${FLEXI_COMPILE_FLAGS}
                                             COMPILE_OPTIONS "-fPIC"
                                             CMAKE_POSITION_INDEPENDENT_CODE ON
                                             LINKER_LANGUAGE Fortran)
  INSTALL(TARGETS ${target} LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                                                GROUP_READ GROUP_EXECUTE
                                                WORLD_READ WORLD_EXECUTE
                            DESTINATION lib)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
  # Strip shared library when building in release mode, except on macOS
  # IF(NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
  IF(NOT CMAKE_HOST_APPLE)
    ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND $<$<CONFIG:release>:${CMAKE_STRIP}> ARGS --strip-unneeded $<TARGET_FILE:${target}>)
  ENDIF()
ENDFUNCTION()

# include flexi.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)

# collect all F90 files needed for flexilib, except the extra files below
FILE(GLOB_RECURSE flexiF90 ./src/interpolation/*.f90
                           ./src/readintools/*.f90
                           ./src/mortar/*.f90
                           ./src/mesh/*.f90
                           ./src/mpi/*.f90
                           ./src/globals/*.f90
                           ./src/analyze/*.f90
                           ./src/recordpoints/*.f90
                           ./src/restart/*.f90
                           ./src/output/*.f90
                           ./src/baseflow/*.f90
                           ./src/sponge/*.f90
                           ./src/avg2d/*.f90
                           ./src/filter/*.f90
                           ./src/io_hdf5/*.f90
                           ./src/timedisc/*.f90
                           ./src/implicit/*.f90
                           ./src/areas/*.f90
                           ./src/testcase/${FLEXI_TESTCASE}/*.f90
                           # ./src/flexilib.f90
                           ./unitTests/unittest.f90
                           ./unitTests/unittest_vars.f90)

# The following files/directories where excluded from the above global recursive search,
# since some of the subdirectories should only be included if the respective feature is enabled.

# >> collect F90 files of the equation system (non-recursive) in a new list 'eqnsysF90'
FILE(GLOB eqnsysF90        ./src/equations/${FLEXI_EQNSYSNAME}/*.f90)
# append F90 files of the equation of state to 'eqnsysF90'-list
IF (WITHEOS GREATER 0)
  IF(POSTI)
    FILE(GLOB eosF90         ./src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/*.f90
              eosF90         ./src/equations/${FLEXI_EQNSYSNAME}/posti/*.f90)
  ELSE()
    FILE(GLOB eosF90         ./src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/*.f90)
  ENDIF()
  LIST(APPEND eqnsysF90 ${eosF90})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/equations/${FLEXI_EQNSYSNAME}/${FLEXI_EOS}/)
ELSE()
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src/equations/${FLEXI_EQNSYSNAME}/)
ENDIF()

# >> remove F90 files of the splitflux from 'eqnsysF90'-list (if file exists for selected equation system)
FILE(GLOB splitfluxF90 ./src/equations/${FLEXI_EQNSYSNAME}/splitflux.f90)
IF (NOT FLEXI_SPLIT_DG AND EXISTS "${splitfluxF90}")
    LIST(REMOVE_ITEM eqnsysF90 ${splitfluxF90})
ENDIF()

# >> append F90 files of the eddy viscosity to 'eqnsysF90'-list
IF (FLEXI_EDDYVISCOSITY)
  FILE(GLOB eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/default/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/smagorinsky/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/dynSmagorinsky/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/vreman/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/sigmaModel/*.f90
            eddyviscF90 ./src/equations/${FLEXI_EQNSYSNAME}/eddyVisc/wale/*.f90)
  LIST(APPEND eqnsysF90 ${eddyviscF90})
ENDIF()

# >> append 'eqnsysF90'-list to 'flexiF90'-list
LIST(REMOVE_DUPLICATES eqnsysF90)
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
  FILE(GLOB fvF90 ./src/indicator/*.f90
                  ./src/fv/*.f90)
  LIST(APPEND flexiF90 ${fvF90})
ENDIF()

IF (FLEXI_USE_PRECONDITIONER)
  FILE(GLOB precondF90     ./src/precond/*.f90
                           ./src/precond/jac_ex/*.f90
                           ./src/precond/jac_fd/*.f90)
  LIST(APPEND flexiF90 ${precondF90})
ENDIF()

# define libs
ADD_LIBRARY(          libflexistaticF90 OBJECT ${flexiF90})
SET_TARGET_PROPERTIES(libflexistaticF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                                   COMPILE_OPTIONS "-fPIC"
                                                   CMAKE_POSITION_INDEPENDENT_CODE ON)
ADD_LIBRARY(          libflexistatic    STATIC $<TARGET_OBJECTS:libflexistaticF90>)
SET_TARGET_PROPERTIES(libflexistatic    PROPERTIES PREFIX ""
                                                   COMPILE_OPTIONS "-fPIC"
                                                   CMAKE_POSITION_INDEPENDENT_CODE ON
                                                   OUTPUT_NAME "libflexif90")
# Make sure that the script is executed even though cmake is not explicitly executed after the commit is changed in git
ADD_DEPENDENCIES(libflexistaticF90 UpdateGitCommitHash)
IF(SELFBUILTEXTERNALS)
  ADD_DEPENDENCIES(libflexistaticF90 ${SELFBUILTEXTERNALS})
ENDIF()

# Userblock
#EXECUTE_PROCESS(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/output/print_userblock.c ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/output/print_userblock.c ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COPYONLY)
ADD_CUSTOM_TARGET(generateUserblock COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/tools/userblock/generate_userblock.sh
                                               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
                                               ${CMAKE_BINARY_DIR}
                                              "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}"
                                               ${CMAKE_CURRENT_SOURCE_DIR})
ADD_LIBRARY(userblocklib STATIC ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/print_userblock.c)
ADD_DEPENDENCIES(userblocklib generateUserblock)

ADD_LIBRARY(stacksizelib STATIC ./src/globals/stacksize.c)
SET_TARGET_PROPERTIES(stacksizelib userblocklib PROPERTIES LANGUAGE C
                                                           LINKER_LANGUAGE C
                                                           COMPILE_OPTIONS "-fPIC"
                                                           CMAKE_POSITION_INDEPENDENT_CODE ON)

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
