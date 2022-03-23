# =========================================================================
# SWAPMESH
# =========================================================================
SET(TOOL_NAME "swapmesh")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/swapmesh")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/swapmesh.f90
                          ${TOOL_DIR}/swapmesh_vars.f90
                          ${TOOL_DIR}/parametriccoords.f90
                          ${TOOL_DIR}/interpolatesolution.f90)

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
TARGET_LINK_LIBRARIES(${TOOL_NAME}lib libposti)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared libposti)

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/posti_swapmesh.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
