# =========================================================================
# After settings specific compilers, enable named languages for cmake
# =========================================================================
ENABLE_LANGUAGE(Fortran C CXX)
MARK_AS_ADVANCED(FORCE C_PATH CXX_PATH Fortran_PATH)

# =========================================================================
# Set machine-specific definitions and settings
# =========================================================================
# SuperMUC
IF (CMAKE_FQDN_HOST MATCHES "sng\.lrz\.de$")
  MESSAGE(STATUS "Compiling on SuperMUC")
  # Overwrite compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    SET(FLEXI_INSTRUCTION "-march=skylake-avx512 -mtune=skylake-avx512" CACHE STRING "Compiler optimization options")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(FLEXI_INSTRUCTION "-xSKYLAKE-AVX512" CACHE STRING "Compiler optimization options")
    # Explicitely enable usage of AVX512 registers
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopt-zmm-usage=high")
  ENDIF()
  # Set LUSTRE definition to account for filesystem and MPI implementation
  ADD_COMPILE_DEFINITIONS(LUSTRE)

# LUMI
ELSEIF(CMAKE_FQDN_HOST MATCHES "\.can$")
  MESSAGE(STATUS "Compiling on LUMI")
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    SET(FLEXI_INSTRUCTION "-march=znver3 -mtune=znver3" CACHE STRING "Compiler optimization options")
  ELSE()
    MESSAGE(FATAL_ERROR "LUMI currently only supported using the GNU Compiler Collection (GCC). Please load/swap the following modules: LUMI PrgEnv-gnu cray-hdf5-parallel")
  ENDIF()

# IAG Prandtl
ELSEIF(CMAKE_FQDN_HOST MATCHES "^(prandtl|grafik.*)\.iag\.uni\-stuttgart\.de")
  MESSAGE(STATUS "Compiling on ${CMAKE_HOSTNAME}")
  SET(FLEXI_INSTRUCTION "-march=native -mtune=native" CACHE STRING "Compiler optimization options")
  # Set LUSTRE definition to account for filesystem
  ADD_COMPILE_DEFINITIONS(LUSTRE)

# macOS
# ELSEIF(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
ELSEIF(CMAKE_HOST_APPLE)
  MESSAGE(STATUS "Compiling on generic ${CMAKE_HOST_SYSTEM_NAME} ${CMAKE_HOST_SYSTEM_VERSION} machine [${CMAKE_HOSTNAME}]")
  # Skip -mtune on macOS
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang" OR CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    SET(FLEXI_INSTRUCTION "-march=native" CACHE STRING "Compiler optimization options")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(FLEXI_INSTRUCTION "-xHost" CACHE STRING "Compiler optimization options")
  ENDIF()

# Generic machine
ELSE()
  MESSAGE(STATUS "Compiling on generic ${CMAKE_HOST_SYSTEM_NAME} ${CMAKE_HOST_SYSTEM_VERSION} machine [${CMAKE_HOSTNAME}]")
  # Set compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang" OR CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    SET(FLEXI_INSTRUCTION "-march=native -mtune=native" CACHE STRING "Compiler optimization options")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(FLEXI_INSTRUCTION "-xHost" CACHE STRING "Compiler optimization options")
  ENDIF()
ENDIF()

MESSAGE(STATUS "Compiling with [${CMAKE_Fortran_COMPILER_ID}] (v${CMAKE_Fortran_COMPILER_VERSION}) fortran compiler using [${FLEXI_INSTRUCTION}] instruction")

# =========================================================================
# COMPILER FLAGS
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# CMake can always request position independent code
# SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

# GNU Compiler Collection (GCC)
IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # set Flags (disable lto type warnings due to false positives with MATMUL, which is a known bug)
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -finit-real=snan -finit-integer=snan -Wno-lto-type-mismatch -Wno-missing-include-dirs -lstdc++ -DGNU")
  ENDIF()
  # initialize all variables as signalling NaNs to force the user to correctly initialize these data types
  SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS}     -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS} -g  -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays -ffpe-trap=invalid,zero,overflow -fbacktrace")
  SET (CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS} -pg -O3 ${FLEXI_INSTRUCTION} -finline-functions -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS} -g  -Og -ggdb3 -ffpe-trap=invalid,zero,overflow -fbounds-check -fbacktrace -Wall")
  SET (CMAKE_Fortran_FLAGS_SANITIZE       "${CMAKE_Fortran_FLAGS} -g  -Og -ggdb3 -ffpe-trap=invalid,zero,overflow,denorm -fbounds-check -fbacktrace -Wall -fsanitize=address,undefined,leak -fno-omit-frame-pointer -Wc-binding-type -Wuninitialized -pedantic")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (FLEXI_COMPILE_FLAGS "-xf95-cpp-input")
  ELSE()
    # Trailing white space required in case variable is unset!
    SET (FLEXI_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} ")
  ENDIF()

# AMD Optimized LLVM/CLANG
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS              "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -std=f2008 -lstdc++ -DFLANG")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS}     -O3 ${FLEXI_INSTRUCTION} -finline-functions ")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS}     -O3 ${FLEXI_INSTRUCTION} -finline-functions -ffpe-trap=invalid,zero,overflow -fbacktrace")
  SET (CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS} -pg -O3 ${FLEXI_INSTRUCTION} -finline-functions ")
  SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS} -g  -O0 -ggdb3 -ffpe-trap=invalid,zero,overflow -fbounds-check -finit-real=snan -fbacktrace -Wall")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (FLEXI_COMPILE_FLAGS "-xf95-cpp-input")
  ELSE()
    # Trailing white space required in case variable is unset!
    SET (FLEXI_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} ")
  ENDIF()

# Intel Compiler
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS              "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -lstdc++ -DINTEL")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS}    -O3 ${FLEXI_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS}    -O3 ${FLEXI_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div -fpe0 -traceback")
  SET (CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS} -p -O3 ${FLEXI_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (FLEXI_COMPILE_FLAGS "-fpp -allow nofpp_comments -assume bscc")
  ELSE()
    SET (FLEXI_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} -allow nofpp_comments -assume bscc")
  ENDIF()

# Cray Compiler
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS              "${CMAKE_Fortran_FLAGS} -ffree -s real64 -s integer64 -em -lstdc++ -hfp0 -DCRAY")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS} -s integer32  -O2 -hfp3 -p . -rm")
  SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS} -s integer32  -O2 -hfp3 -p . -rm -eD")
  SET (CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS} -s integer32  -O2 -hfp3 -h profile_generate -p . -rm")
  SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS} -s integer32 -g -O0 -eD -rm")
  # add flags only for compiling not linking!
  SET (FLEXI_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} -F")
ELSE()
  MESSAGE(SEND_ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# Profile-Guided Optimization (PGO)
# =========================================================================
CMAKE_DEPENDENT_OPTION(FLEXI_PERFORMANCE_PGO "Enable profile-guided optimization (Only GNU Compiler supported)" OFF
                                             "FLEXI_PERFORMANCE" OFF)
IF (FLEXI_PERFORMANCE_PGO)
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    SET(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}       -fprofile-use")
    SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_REWITHDEBINFO} -fprofile-use")
    SET(CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS_PROFILE}       -fprofile-generate")
  ELSE()
    MESSAGE(SEND_ERROR "Profile-guided optimization (PGO) currently only supported for GNU compiler. Either set FLEXI_PERFORMANCE_PGO=OFF or use the GNU compiler." )
  ENDIF()
ENDIF()

# =========================================================================
# SAVE CURRENT COMPILER FLAGS TO THE CACHE
# =========================================================================
MARK_AS_ADVANCED(FORCE FLEXI_INSTRUCTION)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_RELEASE)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_PROFILE)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_DEBUG)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_SANITIZE)
SET(CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS}"                CACHE STRING "Default compiler flags"        FORCE)
SET(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}"        CACHE STRING "Release compiler flags"        FORCE)
SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" CACHE STRING "RelWithDebInfo compiler flags" FORCE)
SET(CMAKE_Fortran_FLAGS_PROFILE        "${CMAKE_Fortran_FLAGS_PROFILE}"        CACHE STRING "Profile compiler flags"        FORCE)
SET(CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}"          CACHE STRING "Debug compiler flags"          FORCE)
SET(CMAKE_Fortran_FLAGS_SANITIZE       "${CMAKE_Fortran_FLAGS_SANITIZE}"       CACHE STRING "Sanitize compiler flags"       FORCE)
