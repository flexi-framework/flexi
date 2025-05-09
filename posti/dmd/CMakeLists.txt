# =========================================================================
# DMD tool
# =========================================================================
SET(TOOL_NAME "dmd")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/dmd")
SET(TOOL_MAIN "posti_dmd.f90")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/dmd.f90
                          ${TOOL_DIR}/dmd_vars.f90
                          ${TOOL_DIR}/navierstokes/equationDmd.f90
                          ${TOOL_DIR}/navierstokes/equationDmd_vars.f90)

# Build an object library from the Fortran source files
ADD_LIBRARY(               ${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
TARGET_INCLUDE_DIRECTORIES(${TOOL_NAME}libF90 PRIVATE ${TOOL_DIR})
SET_TARGET_PROPERTIES(     ${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                                         COMPILE_OPTIONS "-fPIC"
                                                         CMAKE_POSITION_INDEPENDENT_CODE ON)

# Build a static library from the object library
ADD_LIBRARY(          ${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(     ${TOOL_NAME}libF90 libflexishared)

# When using Ninja, CMake tries to parallelize the build resulting in missing module files. Add a dummy target to ensure synchronization
ADD_CUSTOM_TARGET(target_${TOOL_NAME} COMMENT "Synchronizing build for ${TOOL_NAME} tool")
ADD_DEPENDENCIES( target_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES( ${TOOL_NAME}libF90 target_${TOOL_NAME})

# Standalone Posti executable
ADD_EXEC(                  posti_${TOOL_NAME} ${TOOL_DIR}/${TOOL_MAIN})
TARGET_INCLUDE_DIRECTORIES(posti_${TOOL_NAME} PRIVATE ${TOOL_DIR})
TARGET_LINK_LIBRARIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(          posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
