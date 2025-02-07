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
# # > Gold is deprecated as of GNU Binutils 2.44
# # > https://lists.gnu.org/archive/html/info-gnu/2025-02/msg00001.html
# ELSE()
#   EXECUTE_PROCESS(COMMAND ld.gold --version COMMAND grep "^GNU gold" COMMAND sed "s/^.* //g" OUTPUT_VARIABLE GNU_GOLD_VERSION)
#   IF (DEFINED GNU_GOLD_VERSION  AND NOT "${GNU_GOLD_VERSION}" STREQUAL "")
#     STRING(STRIP "${GNU_GOLD_VERSION}" GNU_GOLD_VERSION)
#     MESSAGE(STATUS "Linking with [gold] (v${GNU_GOLD_VERSION})")
#     # Shift responsibility of driving the final stages of compilation from collect2 to gold via the linker plugin
#     # More information at: https://gcc.gnu.org/wiki/LinkTimeOptimization
#     IF(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
#       ADD_LINK_OPTIONS("-fuse-ld=gold")
#       # Make it abundantly clear we want to use gold
#       FIND_PROGRAM(CMAKE_GOLD_LINKER NAMES ${_CMAKE_TOOLCHAIN_PREFIX}ld.gold${_CMAKE_TOOLCHAIN_SUFFIX} HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
#       SET (CMAKE_LINKER "${CMAKE_GOLD_LINKER}" CACHE FILEPATH "" FORCE)
#       MARK_AS_ADVANCED(FORCE CMAKE_GOLD_LINKER)
#     ELSE()
#       SET (CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    -fuse-ld=gold")
#       SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold")
#       # This currently breaks ar (binutils archiver)
#       # SET (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -fuse-ld=gold")
#       SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fuse-ld=gold")
#       # Make it abundantly clear we want to use gold
#       FIND_PROGRAM(CMAKE_GOLD_LINKER NAMES ${_CMAKE_TOOLCHAIN_PREFIX}ld.gold${_CMAKE_TOOLCHAIN_SUFFIX} HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
#       SET (CMAKE_LINKER "${CMAKE_GOLD_LINKER}" CACHE FILEPATH "" FORCE)
#       MARK_AS_ADVANCED(FORCE CMAKE_GOLD_LINKER)
#     ENDIF()
# Found neither lld nor gold, output GNU ld
#   ELSE()
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
#   ENDIF()
# ENDIF()
