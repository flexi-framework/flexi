# =========================================================================
# INIT HIT
# =========================================================================
SET(TOOL_NAME "init_hit")
SET(TOOL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/init_hit)

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/posti_init_hit.f90
                          ${TOOL_DIR}/init_hit.f90
                          ${TOOL_DIR}/init_hit_vars.f90
                          ${TOOL_DIR}/../shared/fft.f90
                          ${TOOL_DIR}/../shared/fft_vars.f90

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
TARGET_LINK_LIBRARIES(${TOOL_NAME}lib fftwlib)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared fftwlib)

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/posti_init_hit.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared ${SELFBUILTEXTERNALS})
