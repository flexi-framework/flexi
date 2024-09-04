# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Release, RelWithDebInfo, Profile, Debug, Sanitize (only GNU))." FORCE)
   IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo Profile Sanitize)
   ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo Profile)
   ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo Profile)
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

