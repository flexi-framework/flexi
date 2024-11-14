# =========================================================================
# Detect machine environments
# =========================================================================
CMAKE_HOST_SYSTEM_INFORMATION(RESULT CMAKE_FQDN_HOST QUERY FQDN)
MARK_AS_ADVANCED(FORCE CMAKE_FQDN_HOST)
MARK_AS_ADVANCED(FORCE CMAKE_HOSTNAME)
SITE_NAME(CMAKE_HOSTNAME)

# =========================================================================
# CMake generator settings
# =========================================================================
SET(USED_CMAKE_GENERATOR "${CMAKE_GENERATOR}" CACHE STRING "Expose CMAKE_GENERATOR (cannot be changed here)" FORCE)
IF("${CMAKE_GENERATOR}" MATCHES "Ninja")
  # CMake introduced the CMAKE_COLOR_DIAGNOSTICS flag with 3.24.0, https://gitlab.kitware.com/cmake/cmake/-/merge_requests/6990
  IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.24.0"))
    SET(CMAKE_COLOR_DIAGNOSTICS ON CACHE INTERNAL "Flag if CMake should attempt to color output")
  ELSE()
    SET(NINJA_COLOR_DIAGNOSTICS "-fdiagnostics-color=always" CACHE INTERNAL "Flag if Ninja should attempt to color output")
  ENDIF()
ENDIF()
MESSAGE(STATUS "Generating for [${CMAKE_GENERATOR}] build system")

# =========================================================================
# Some clusters requires setting the compilers by hand and invoking
# ENABLE_LANGUAGE afterwards, which is required for
# CMAKE_Fortran_COMPILER_ID that is used below
# > This block must be called before ENABLE_LANGUAGE
# =========================================================================
# HLRS HAWK
IF (CMAKE_FQDN_HOST MATCHES "hawk\.hww\.hlrs\.de$")
  SET(CMAKE_C_COMPILER       mpicc)
  SET(CMAKE_CXX_COMPILER     mpicxx)
  SET(CMAKE_Fortran_COMPILER mpif90) # mpif08 wrapper seems to have issue
# SuperMUC
# ELSEIF(CMAKE_FQDN_HOST MATCHES "sng\.lrz\.de$"
# LUMI
ELSEIF(CMAKE_FQDN_HOST MATCHES "\.can$")
  SET(CMAKE_C_COMPILER       cc)
  SET(CMAKE_CXX_COMPILER     CC)
  SET(CMAKE_Fortran_COMPILER ftn)
# IAG Prandtl
ELSEIF(CMAKE_FQDN_HOST MATCHES "^(prandtl|grafik.*)\.iag\.uni\-stuttgart\.de")
  SET(CMAKE_C_COMPILER       gcc)
  SET(CMAKE_CXX_COMPILER     c++)
  SET(CMAKE_Fortran_COMPILER gfortran)
# IAG Grafik01/Grafik02
ELSEIF (CMAKE_FQDN_HOST MATCHES "^ila(head.*|cfd.*)\.ila.uni\-stuttgart\.de")
  SET(CMAKE_C_COMPILER       mpicc)
  SET(CMAKE_CXX_COMPILER     mpicxx)
  SET(CMAKE_Fortran_COMPILER mpif90) # mpif08 wrapper seems to have issue
ELSEIF (CMAKE_FQDN_HOST MATCHES "^(xenon.*|argon.*)\.ila.uni\-stuttgart\.de")
  SET(CMAKE_C_COMPILER       mpicc)
  SET(CMAKE_CXX_COMPILER     mpicxx)
  SET(CMAKE_Fortran_COMPILER mpif90) # mpif08 wrapper seems to have issue
ENDIF()

# =========================================================================
# Score-P instrumentation infrastructure
# > This option must be called before ENABLE_LANGUAGE, thus is only available
# > through -DMACHINE_USE_SCOREP=ON
# =========================================================================
IF (MACHINE_USE_SCOREP)
  FIND_PROGRAM(SCOREP_C_COMPILER scorep-${CMAKE_C_COMPILER})
  MARK_AS_ADVANCED(FORCE SCOREP_C_COMPILER )
  IF (SCOREP_C_COMPILER MATCHES "NOTFOUND")
    MESSAGE (FATAL_ERROR "Score-P not available in PATH. Did you load the module?")
  ENDIF()

  # Set default build type to profile
  IF (NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE Profile CACHE STRING "Choose the type of build, options are: Release RelWithDebInfo Profile Debug Sanitize." FORCE)
  ENDIF (NOT CMAKE_BUILD_TYPE)
  IF (CMAKE_BUILD_TYPE MATCHES "Release")
    MESSAGE (WARNING "Score-P requires debug compile flags which are not available with BUILD_TYPE='Release'")
  ENDIF()

  SET(CMAKE_C_COMPILER       "scorep-${CMAKE_C_COMPILER}")
  SET(CMAKE_CXX_COMPILER     "scorep-${CMAKE_CXX_COMPILER}")
  SET(CMAKE_Fortran_COMPILER "scorep-${CMAKE_Fortran_COMPILER}")
ENDIF()
