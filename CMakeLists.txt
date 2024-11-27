# =========================================================================
# Flexi information
# =========================================================================
SET(PROJECT_NAME "Flexi")
SET(PROJECT_DESC "A high order discontinuous Galerkin framework for hyperbolic–parabolic conservation laws")
SET(PROJECT_VER  "22.03")
SET(PROJECT_URL  "https://www.flexi-project.org")

# =========================================================================
# CMake configuration
# =========================================================================
MESSAGE(STATUS "This is ${PROJECT_NAME} (v${PROJECT_VER})")
MESSAGE(STATUS "Generating with [cmake] (v${CMAKE_VERSION})")
# NOTE: use minimum CMake version required
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0)
# Policies (see https://cmake.org/cmake/help/v3.13/manual/cmake-policies.7.html)
  SET(GITSHALLOW "GIT_SHALLOW ON")
IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()
# Required to link visulib against visuReader
IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13.0")
  CMAKE_POLICY(SET CMP0079 NEW)
ENDIF()

# =========================================================================
# Load CMake modules
# =========================================================================
# Directory for custom CMake source files
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(CheckOrigin)

# Load CMake internal extensions
INCLUDE(CMakeDependentOption)
INCLUDE(ExternalProject)

# =========================================================================
# Performance mode
# =========================================================================
OPTION(FLEXI_PERFORMANCE "Enable performance optimizations" OFF)

# =========================================================================
# Machine environment
# =========================================================================

INCLUDE(SetLinker)
# Languages are specified in CMakeListMachine.txt as Hawk requires them
# after specifiying CMAKE_COMPILER
INCLUDE(SetMachine)

PROJECT(Flexi
        VERSION      ${PROJECT_VER}
        DESCRIPTION  ${PROJECT_DESC}
        HOMEPAGE_URL ${PROJECT_URL})

INCLUDE(SetCompiler)
# Load CMake internal extensions
# > must be loaded after target architecture is known
INCLUDE(GNUInstallDirs)

# =========================================================================
# Store the current commit information
# =========================================================================
EXECUTE_PROCESS(COMMAND git rev-parse HEAD OUTPUT_VARIABLE GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} RESULT_VARIABLE ERR_NO_GIT)
IF(${ERR_NO_GIT})
  MESSAGE("No git information found. Continue without git.")
ELSE()
  SET(commit ${CMAKE_CURRENT_SOURCE_DIR}/src/commit.h)
  MESSAGE(STATUS "Current git commit ${GIT_COMMIT} will be written to ${commit}")
  FILE(WRITE ${commit} "! Define the current commit hash. The default must remain empty, i.e., ''. Do not commit the changed file!\n#define GIT_CURRENT_COMMIT  ''\n")
  EXECUTE_PROCESS(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/cmake/setCommitHash.sh ${commit} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()
# Make sure that the script is executed even though cmake is not explicitly executed after the commit is changed in git
ADD_CUSTOM_TARGET(UpdateGitCommitHash COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/cmake/setCommitHash.sh ${commit})

# =========================================================================
# Check IPO support:
# =========================================================================
# we need to have languages enabled and compilers defined for this
OPTION(FLEXI_PERFORMANCE_IPO "Enable Interprocedural Optimization for build (if available)." ON)
MARK_AS_ADVANCED(FORCE FLEXI_PERFORMANCE_IPO)
IF(FLEXI_PERFORMANCE_IPO AND NOT(${CMAKE_VERSION} VERSION_LESS "3.9.0"))
  CMAKE_POLICY(SET CMP0069 NEW)
  INCLUDE(CheckIPOSupported)
  CHECK_IPO_SUPPORTED(RESULT HASIPO OUTPUT error)
ELSE()
  SET(HASIPO FALSE)
ENDIF()

# =========================================================================
# Output paths
# =========================================================================
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Set custom install dir (needs to be done after project has been defined!)
# MARK_AS_ADVANCED(FORCE CMAKE_INSTALL_PREFIX)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()

# =========================================================================
# Build type
# =========================================================================
INCLUDE(SetBuildType)

# =========================================================================
# Location of binary and filenames
# =========================================================================
# append relative filename-macro for __FILENAME__ in Stamp of abort function (see flexi.h)
# subst replacement is GNU make specific
IF("${CMAKE_GENERATOR}" MATCHES "Unix Makefiles")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")
ENDIF()

# add basedir for location of corresponding userblock-file
ADD_COMPILE_DEFINITIONS("BASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

# =========================================================================
# echo COMPILE_DEFINITIONS
# =========================================================================
GET_DIRECTORY_PROPERTY(comp_defs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
SET(comp_defs_str "COMPILE_DEFINITIONS = ")
FOREACH (def ${comp_defs})
   SET(comp_defs_str "${comp_defs_str} -D${def}")
ENDFOREACH()
ADD_CUSTOM_TARGET(preproc_defines COMMAND echo ${comp_defs_str})

# =========================================================================
# USERBLOCK + PREPROC_FLAGS
# =========================================================================
# A function to get all user defined variables with a specified prefix
FUNCTION (GETLISTOFVARSSTARTINGWITH _PREFIX _VARRESULT)
    GET_CMAKE_PROPERTY(_VARS CACHE_VARIABLES)
    STRING(REGEX MATCHALL "(^|;)${_PREFIX}[A-Za-z0-9_]*" _MATCHEDVARS "${_VARS}")
    SET(${_VARRESULT} ${_MATCHEDVARS} PARENT_SCOPE)
ENDFUNCTION()

# A function to get all user defined variables with a specified prefix
SET(configuration ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake)
FILE(WRITE ${configuration} "")
GETLISTOFVARSSTARTINGWITH("CMAKE_BUILD_TYPE" flexiVars1)
GETLISTOFVARSSTARTINGWITH("FLEXI_" flexiVars2)
GETLISTOFVARSSTARTINGWITH("LIBS_"  flexiVars3)
GETLISTOFVARSSTARTINGWITH("POSTI_" flexiVars4)
FOREACH (_var IN LISTS flexiVars1 flexiVars2 flexiVars3 flexiVars4)
  GET_PROPERTY(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
  # Skip empty variables
  IF(NOT ${${_var}} STREQUAL "")
    SET(boolian (${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF"))
    IF(${boolian})
      FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
    ELSE()
      IF(${_var})
        STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
      ENDIF()
      FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
    ENDIF()
  ENDIF()
ENDFOREACH()

# =========================================================================
# ADD LIBRARIES
# =========================================================================
INCLUDE(SetLibraries)

# =========================================================================
# FLEXI
# =========================================================================
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/src/CMakeLists.txt)

# =========================================================================
# POSTI
# =========================================================================
OPTION(POSTI "Build POSTI toolset" OFF)
IF(POSTI)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/CMakeLists.txt)
ELSE()
  GETLISTOFVARSSTARTINGWITH("POSTI_" postiVars)
  FOREACH (_var IN LISTS postiVars)
    UNSET(${_var} CACHE)
  ENDFOREACH()
ENDIF()

# =========================================================================
# Install
# =========================================================================
INSTALL(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake DESTINATION bin)
