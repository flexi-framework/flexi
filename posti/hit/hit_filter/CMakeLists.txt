# =========================================================================
# HIT FILTER
# =========================================================================
SET(TOOL_NAME "hit_filter")
SET(TOOL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/hit_filter)
SET(TOOL_MAIN "posti_hit_filter.f90")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/posti_hit_filter.f90
                          ${TOOL_DIR}/hit_filter.f90
                          ${TOOL_DIR}/hit_filter_vars.f90)

# Build an object library from the Fortran source files
ADD_LIBRARY(               ${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
TARGET_INCLUDE_DIRECTORIES(${TOOL_NAME}libF90 PRIVATE ${TOOL_DIR})
SET_TARGET_PROPERTIES(     ${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                                         COMPILE_OPTIONS "-fPIC"
                                                         CMAKE_POSITION_INDEPENDENT_CODE ON)

# Build a static library from the object library
ADD_LIBRARY(          ${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
TARGET_LINK_LIBRARIES(${TOOL_NAME}lib hit_fftlib)
ADD_DEPENDENCIES(     ${TOOL_NAME}libF90 libflexishared hit_fftlib)

# When using Ninja, CMake tries to parallelize the build resulting in missing module files. Add a dummy target to ensure synchronization
ADD_CUSTOM_TARGET(target_${TOOL_NAME} COMMENT "Synchronizing build for ${TOOL_NAME} tool")
ADD_DEPENDENCIES( target_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES( ${TOOL_NAME}libF90 target_${TOOL_NAME})

# Standalone Posti executable
ADD_EXEC(                  posti_${TOOL_NAME} ${TOOL_DIR}/${TOOL_MAIN})
TARGET_INCLUDE_DIRECTORIES(posti_${TOOL_NAME} PRIVATE ${TOOL_DIR})
TARGET_LINK_LIBRARIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(          posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared ${SELFBUILTEXTERNALS})
