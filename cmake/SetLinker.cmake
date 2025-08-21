# =========================================================================
# Check/set the linker
# > The following code only works with CMake 3.29+ and defaults to ld otherwise
# =========================================================================
IF(${CMAKE_VERSION} VERSION_LESS "3.29.0")
  RETURN()
ENDIF()

# Try to read from environment variable
IF(DEFINED ENV{CMAKE_LINKER_TYPE} AND NOT "$ENV{CMAKE_LINKER_TYPE}" STREQUAL "")
  SET(CMAKE_LINKER_TYPE  "$ENV{CMAKE_LINKER_TYPE}")
ENDIF()

# Default is already set
IF (CMAKE_LINKER_TYPE)
  # Use the standard linker, typically GNU ld
  IF ("${CMAKE_LINKER_TYPE}" STREQUAL "DEFAULT" OR "${CMAKE_LINKER_TYPE}" STREQUAL "SYSTEM")
    EXECUTE_PROCESS(COMMAND ld       --version OUTPUT_VARIABLE GNU_LD_STRING)
    # # Check if we really got GNU ld
    IF (GNU_LD_STRING MATCHES "^GNU ld")
      EXECUTE_PROCESS(COMMAND ld     --version COMMAND grep "^GNU ld" COMMAND sed "s/^.* //g"                    OUTPUT_VARIABLE GNU_LD_VERSION)
      STRING(STRIP "${GNU_LD_VERSION}" GNU_LD_VERSION)
      MESSAGE(STATUS "Linking with [ld] (v${GNU_LD_VERSION})")
    # ... or if we are running mold
    ELSEIF (GNU_LD_STRING MATCHES "^mold")
      EXECUTE_PROCESS(COMMAND ld     --version COMMAND grep "^mold"   COMMAND grep -Eo "[0-9]+\.[0-9]+\.[0-9]+ " OUTPUT_VARIABLE GNU_MOLD_VERSION)
      STRING(STRIP "${GNU_MOLD_VERSION}" GNU_MOLD_VERSION)
      MESSAGE(STATUS "Linking with [mold] (v${GNU_MOLD_VERSION})")
    ENDIF()
  # Use the LLVM linker
  ELSEIF("${CMAKE_LINKER_TYPE}" STREQUAL "LLD")
    EXECUTE_PROCESS(COMMAND lld-link --version COMMAND grep "^LLD"    COMMAND grep -Eo "[0-9]+\.[0-9]+\.[0-9]+"  OUTPUT_VARIABLE LLVM_LLD_VERSION)
    STRING(STRIP "${LLVM_LLD_VERSION}" LLVM_LLD_VERSION)
    MESSAGE(STATUS "Linking with [lld] (v${LLVM_LLD_VERSION})")
  # Use the mold linker
  ELSEIF("${CMAKE_LINKER_TYPE}" STREQUAL "MOLD")
    EXECUTE_PROCESS(COMMAND mold     --version COMMAND grep "^mold"   COMMAND grep -Eo "[0-9]+\.[0-9]+\.[0-9]+ " OUTPUT_VARIABLE GNU_MOLD_VERSION)
    STRING(STRIP "${GNU_MOLD_VERSION}" GNU_MOLD_VERSION)
    MESSAGE(STATUS "Linking with [mold] (v${GNU_MOLD_VERSION})")
  ENDIF()
# linker type is not set
ELSE()
  # Try to use GOLD, which should be faster than LD (and doesn't cause segfaults)
  EXECUTE_PROCESS(COMMAND ld.gold --version COMMAND grep "^GNU gold" COMMAND sed "s/^.* //g"                     OUTPUT_VARIABLE GNU_GOLD_VERSION)
  IF (DEFINED GNU_GOLD_VERSION  AND NOT "${GNU_GOLD_VERSION}" STREQUAL "")
    SET(CMAKE_LINKER_TYPE "GOLD")
    STRING(STRIP "${GNU_GOLD_VERSION}" GNU_GOLD_VERSION)
    MESSAGE(STATUS "Linking with [gold] (v${GNU_GOLD_VERSION})")
  ENDIF()
ENDIF()
