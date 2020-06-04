# =========================================================================
# INIT HIT
# =========================================================================
SET(TOOL_NAME "filter_hit")
SET(TOOL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/posti/filter_hit)

MESSAGE(STATUS "Posti: using filter_hit tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/posti_filter_hit.f90
                          ${TOOL_DIR}/filter_hit.f90
                          ${TOOL_DIR}/filter_hit_vars.f90
                          ${TOOL_DIR}/fft.f90
                          ${TOOL_DIR}/../common/fftw/fftw.f90)

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared ${SELFBUILTEXTERNALS})
TARGET_INCLUDE_DIRECTORIES(${TOOL_NAME}libF90 PRIVATE ${FFTW_INCLUDE_DIR})

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/filter_hit.f90)
IF (ENABLE_OPENMP)
  TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared ${FFTW_LIBRARIES} ${FFTW_OMP_LIBRARIES})
ELSE()
  TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared ${FFTW_LIBRARIES})
ENDIF()
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared ${SELFBUILTEXTERNALS})
