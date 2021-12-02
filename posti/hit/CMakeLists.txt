# =========================================================================
# HIT Tools
# =========================================================================

# =========================================================================
# Build shared hit library
# =========================================================================
# collect F90 files for the common HIT modules
FILE(GLOB_RECURSE hitF90 ./posti/hit/shared/*.f90)
ADD_LIBRARY(hitlib STATIC ${hitF90})
SET_TARGET_PROPERTIES(hitlib PROPERTIES OUTPUT_NAME "hitlib")
SET_TARGET_PROPERTIES(hitlib PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
TARGET_LINK_LIBRARIES(hitlib fftwlib libflexishared)
ADD_DEPENDENCIES(hitlib fftwlib libflexishared)

IF(POSTI_HIT_INIT)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/init_hit/CMakeLists.txt)
ENDIF()

IF(POSTI_HIT_ANALYZE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/analyze_hit/CMakeLists.txt)
ENDIF()

IF(POSTI_HIT_FILTER)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/filter_hit/CMakeLists.txt)
ENDIF()
