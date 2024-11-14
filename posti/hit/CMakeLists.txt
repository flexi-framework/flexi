# =========================================================================
# HIT Tools
# =========================================================================

# =========================================================================
# Build shared hit fft library
# =========================================================================
# collect F90 files for the common HIT modules
FILE(GLOB_RECURSE hitF90 ./posti/hit/fft/*.f90)
ADD_LIBRARY(hit_fftlib STATIC ${hitF90})
SET_TARGET_PROPERTIES(hit_fftlib PROPERTIES OUTPUT_NAME "hit_fftlib")
SET_TARGET_PROPERTIES(hit_fftlib PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                            COMPILE_OPTIONS "-fPIC"
                                            CMAKE_POSITION_INDEPENDENT_CODE ON)
TARGET_LINK_LIBRARIES(hit_fftlib fftwlib libflexishared)
ADD_DEPENDENCIES(hit_fftlib fftwlib libflexishared)

# Add single HIT tools
IF(POSTI_HIT_INIT)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/hit_init/CMakeLists.txt)
ENDIF()

IF(POSTI_HIT_ANALYZE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/hit_analyze/CMakeLists.txt)
ENDIF()

IF(POSTI_HIT_FILTER)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/hit_filter/CMakeLists.txt)
ENDIF()
