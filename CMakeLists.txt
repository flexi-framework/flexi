# =========================================================================
# Flexi CMake configuration
# =========================================================================
MESSAGE(STATUS "Generating with [cmake] (${CMAKE_VERSION})")
# NOTE: use minimum CMake version required
CMAKE_MINIMUM_REQUIRED(VERSION 3.6.3)
# Policies (see https://cmake.org/cmake/help/v3.13/manual/cmake-policies.7.html)
# VERSION_GREATER_EQUAL is incompatible with CMAKE versions < 3.7
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.6.0"))
  SET(GITSHALLOW "GIT_SHALLOW ON")
ENDIF()
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.12.0"))
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()
# Required to link visulib against visuReader
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.13.0"))
  CMAKE_POLICY(SET CMP0079 NEW)
ENDIF()

# =========================================================================
# Git configuration
# =========================================================================
# Check where the code originates
EXECUTE_PROCESS(COMMAND git ls-remote --get-url OUTPUT_VARIABLE GIT_ORIGIN)

# Setup git hooks
SET(PRECOMMIT_FILE ".githooks/pre-commit")

# Perform checks only if origin points to IAG, other origins can't commit
IF("${GIT_ORIGIN}" MATCHES ".iag.uni-stuttgart.de")
  # Check if the pre-commit hooks exits
  IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git/hooks/pre-commit)
    # Create otherwise
    EXECUTE_PROCESS(COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/.git/hooks)
    EXECUTE_PROCESS(COMMAND ln -s ${CMAKE_CURRENT_SOURCE_DIR}/${PRECOMMIT_FILE} ${CMAKE_CURRENT_SOURCE_DIR}/.git/hooks/pre-commit)
  ELSE()
    # Check if the hook is the correct symlink and warn otherwise
    EXECUTE_PROCESS(COMMAND readlink ${CMAKE_CURRENT_SOURCE_DIR}/.git/hooks/pre-commit OUTPUT_VARIABLE PRECOMMIT_LINK OUTPUT_STRIP_TRAILING_WHITESPACE)
    IF (NOT ${PRECOMMIT_LINK} MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/${PRECOMMIT_FILE}")
      MESSAGE (WARNING "Custom git pre-commit hook detected. Please ensure to call ${PRECOMMIT_FILE} manually.")
    ENDIF()
  ENDIF()

  # Check if the hook actually gets loaded
  EXECUTE_PROCESS(COMMAND git config --get core.hooksPath OUTPUT_VARIABLE HOOKSPATH OUTPUT_STRIP_TRAILING_WHITESPACE)
  IF (DEFINED HOOKSPATH  AND NOT "${HOOKSPATH}" STREQUAL "" AND NOT "${HOOKSPATH}" STREQUAL ".git/hooks")
    # STRING(ASCII 27 ESCAPE)
    # MESSAGE (STATUS "${ESCAPE}[34mCustom hooks path detected. Please ensure to call ${PRECOMMIT_FILE} manually.${ESCAPE}[0m")
    MESSAGE (WARNING "Custom git hooks path detected. Please ensure to call ${PRECOMMIT_FILE} manually.")
  ENDIF()
ENDIF()

# =========================================================================
# Project definition
# =========================================================================
# PROJECT(Flexi Fortran C CXX)

# Folder for custom CMake source files
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)

# =========================================================================
# Performance mode
# =========================================================================
OPTION(FLEXI_PERFORMANCE "Enable performance optimizations" OFF)

# =========================================================================
# Check lld / gold support
# =========================================================================
# EXECUTE_PROCESS(COMMAND ld.lld --version COMMAND grep "^LLD" COMMAND grep -Eo "[0-9]+\.[0-9]+\.[0-9]+" OUTPUT_VARIABLE LLD_VERSION)
# # lld should be faster than gold
# IF (DEFINED LLD_VERSION  AND NOT "${LLD_VERSION}" STREQUAL "")
#   STRING(STRIP "${LLD_VERSION}" LLD_VERSION)
#   MESSAGE(STATUS "Setting linker to [lld] (v${LLD_VERSION})")
#   # Shift responsibility of driving the final stages of compilation from collect2 to gold via the linker plugin
#   # More information at: https://gcc.gnu.org/wiki/LinkTimeOptimization
#   IF(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
#     ADD_LINK_OPTIONS("-fuse-ld=lld")
#     # Make it abundantly clear we want to use lld
#     FIND_PROGRAM(LLD_LINKER NAMES ${_CMAKE_TOOLCHAIN_PREFIX}ld.lld${_CMAKE_TOOLCHAIN_SUFFIX} HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
#     SET (CMAKE_LINKER "${LLD_LINKER}" CACHE FILEPATH "" FORCE)
#     MARK_AS_ADVANCED(FORCE LLD_LINKER)
#   ELSE()
#     SET (CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    -fuse-ld=lld")
#     SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=lld")
#     SET (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -fuse-ld=lld")
#     SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fuse-ld=lld")
#     # Make it abundantly clear we want to use lld
#     FIND_PROGRAM(LLD_LINKER NAMES ${_CMAKE_TOOLCHAIN_PREFIX}ld.lld${_CMAKE_TOOLCHAIN_SUFFIX} HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
#     SET (CMAKE_LINKER "${LLD_LINKER}" CACHE FILEPATH "" FORCE)
#     MARK_AS_ADVANCED(FORCE LLD_LINKER)
#   ENDIF()
# # gold should be faster than GNU ld
# ELSE()
  EXECUTE_PROCESS(COMMAND ld.gold --version COMMAND grep "^GNU gold" COMMAND sed "s/^.* //g" OUTPUT_VARIABLE GNU_GOLD_VERSION)
  IF (DEFINED GNU_GOLD_VERSION  AND NOT "${GNU_GOLD_VERSION}" STREQUAL "")
    STRING(STRIP "${GNU_GOLD_VERSION}" GNU_GOLD_VERSION)
    MESSAGE(STATUS "Linking with [gold] (v${GNU_GOLD_VERSION})")
    # Shift responsibility of driving the final stages of compilation from collect2 to gold via the linker plugin
    # More information at: https://gcc.gnu.org/wiki/LinkTimeOptimization
    IF(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
      ADD_LINK_OPTIONS("-fuse-ld=gold")
      # Make it abundantly clear we want to use gold
      FIND_PROGRAM(CMAKE_GOLD_LINKER NAMES ${_CMAKE_TOOLCHAIN_PREFIX}ld.gold${_CMAKE_TOOLCHAIN_SUFFIX} HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
      SET (CMAKE_LINKER "${CMAKE_GOLD_LINKER}" CACHE FILEPATH "" FORCE)
      MARK_AS_ADVANCED(FORCE CMAKE_GOLD_LINKER)
    ELSE()
      SET (CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    -fuse-ld=gold")
      SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold")
      # This currently breaks ar (binutils archiver)
      # SET (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -fuse-ld=gold")
      SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fuse-ld=gold")
      # Make it abundantly clear we want to use gold
      FIND_PROGRAM(CMAKE_GOLD_LINKER NAMES ${_CMAKE_TOOLCHAIN_PREFIX}ld.gold${_CMAKE_TOOLCHAIN_SUFFIX} HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
      SET (CMAKE_LINKER "${CMAKE_GOLD_LINKER}" CACHE FILEPATH "" FORCE)
      MARK_AS_ADVANCED(FORCE CMAKE_GOLD_LINKER)
    ENDIF()
# Found neither lld nor gold, output GNU ld
  ELSE()
    EXECUTE_PROCESS(COMMAND ld --version OUTPUT_VARIABLE GNU_LD_STRING)
    # Check if we actually got GNU ld
    IF (GNU_LD_STRING MATCHES "^GNU ld")
      EXECUTE_PROCESS(COMMAND ld --version COMMAND grep "^GNU ld" COMMAND sed "s/^.* //g" OUTPUT_VARIABLE GNU_LD_VERSION)
      STRING(STRIP "${GNU_LD_VERSION}" GNU_LD_VERSION)
      MESSAGE(STATUS "Linking with [ld] (v${GNU_LD_VERSION})")
    # ... or if we are running mold
    ELSEIF (GNU_LD_STRING MATCHES "^mold")
      EXECUTE_PROCESS(COMMAND ld --version COMMAND grep "^mold" COMMAND grep -Eo "[0-9]+\.[0-9]+\.[0-9]+ " OUTPUT_VARIABLE GNU_MOLD_VERSION)
      STRING(STRIP "${GNU_MOLD_VERSION}" GNU_MOLD_VERSION)
      MESSAGE(STATUS "Linking with [mold] (v${GNU_MOLD_VERSION})")
    ENDIF()
  ENDIF()
# ENDIF()

# =========================================================================
# Machine environment
# =========================================================================
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMakeListsMachine.txt)

# Hawk requires this after specifiying CMAKE_COMPILER
PROJECT(Flexi)

# =========================================================================
# Store the current commit information
# =========================================================================
EXECUTE_PROCESS(COMMAND git rev-parse HEAD OUTPUT_VARIABLE GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE ERR_NO_GIT)
IF(${ERR_NO_GIT})
  MESSAGE("No git information found. Continue without git.")
ELSE()
  SET(commit ${CMAKE_CURRENT_SOURCE_DIR}/src/commit.h)
  MESSAGE(STATUS "Current git commit ${GIT_COMMIT} will be written to ${commit}")
  FILE(WRITE ${commit} "! Define the current commit hash. The default must remain empty, i.e., ''. Do not commit the changed file!\n#define GIT_CURRENT_COMMIT  ''\n")
  EXECUTE_PROCESS(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/cmake/setCommitHash.sh ${commit})
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
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug, Release, Profile, Sanitize (only GNU))." FORCE)
   IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile Sanitize)
   ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
   ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
   ENDIF()
ENDIF (NOT CMAKE_BUILD_TYPE)

STRING(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE_LC)
IF (BUILD_TYPE_LC MATCHES "debug" OR BUILD_TYPE_LC MATCHES "sanitize")
  ADD_DEFINITIONS("-DDEBUG")
ELSE()
  IF (HASIPO)
    # enable IPO globally (IPO branding: Intel => IPO, GNU => LTO, PGI => IPA)
    IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      # Check the GCC wrapper for the complete toolchain
      FIND_PROGRAM(CMAKE_GCC_AR     NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-ar${_CMAKE_TOOLCHAIN_SUFFIX} HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
      FIND_PROGRAM(CMAKE_GCC_NM     NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-nm                           HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
      FIND_PROGRAM(CMAKE_GCC_RANLIB NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-ranlib                       HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
      MARK_AS_ADVANCED(FORCE CMAKE_GCC_AR)
      MARK_AS_ADVANCED(FORCE CMAKE_GCC_NM)
      MARK_AS_ADVANCED(FORCE CMAKE_GCC_RANLIB)
      # Do not use the standard CMake LTO option for GNU (-flto -fno-fat-lto-objects), as it does not allow speed-up during linking
      IF( CMAKE_GCC_AR AND CMAKE_GCC_NM AND CMAKE_GCC_RANLIB )
        MESSAGE(STATUS "Found GCC binutils wrappers for LTO. Enabling LTO linker plugin.")
        # Do not use the standard CMake LTO option for GNU (-flto -fno-fat-lto-objects), as it does not allow speed-up during linking
        SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
        # Static libraries require either fat LTO objects (increases compilation time) or the use of linker plugins (per default enabled); the jobserver option reduces linking time
        # More information at: https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html#Optimize-Options
        IF (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
          # GCC10 introduced the -flto=auto option, https://gcc.gnu.org/legacy-ml/gcc-patches/2019-07/msg01488.html
          SET(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -flto=auto -fuse-linker-plugin")
        ELSE()
          # Otherwise we must rely on hardcoded GCC jobserver. If it is not available, we will get warnings in the console but it still works, albeit in serial mode
          SET(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -flto=jobserver -fuse-linker-plugin")
        ENDIF()
        # Use the GCC wrapper for the complete toolchain to provide the path to the linker plugin (this might be the problem with using the CMAKE IPO)
        SET(CMAKE_AR     "${CMAKE_GCC_AR}"     CACHE FILEPATH "" FORCE)
        SET(CMAKE_NM     "${CMAKE_GCC_NM}"     CACHE FILEPATH "" FORCE)
        SET(CMAKE_RANLIB "${CMAKE_GCC_RANLIB}" CACHE FILEPATH "" FORCE)
        MARK_AS_ADVANCED(FORCE CMAKE_AR)
        MARK_AS_ADVANCED(FORCE CMAKE_NM)
        MARK_AS_ADVANCED(FORCE CMAKE_RANLIB)
      ELSE()
        MESSAGE(WARNING "GCC indicates LTO support, but binutils wrappers could not be found. Disabling LTO linker plugin." )
        SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)  # enable IPO globally
      ENDIF()
    ELSE()
      SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)  # enable IPO globally
    ENDIF()
  ENDIF()
ENDIF()
MESSAGE(STATUS "Compiling with [${BUILD_TYPE_LC}] build type")

# =========================================================================
# Location of binary and filenames
# =========================================================================
# append relative filename-macro for __FILENAME__ in Stamp of abort function (see flexi.h)
# subst replacement is GNU make specific
IF("${CMAKE_GENERATOR}" MATCHES "Unix Makefiles")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")
ENDIF()

# add basedir for location of corresponding userblock-file
ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

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
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMakeListsLib.txt)

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
