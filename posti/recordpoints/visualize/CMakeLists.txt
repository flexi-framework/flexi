# =========================================================================
# VISUALIZE RECORDPOINTS
# =========================================================================
SET(TOOL_NAME "visualizerecordpoints")
SET(TOOL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/visualize)
SET(TOOL_MAIN "visualizerecordpoints.f90")

MESSAGE(STATUS "Posti: using visualize record points tool.")

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/parameters.f90
                          ${TOOL_DIR}/fft/*.f90
                          ${TOOL_DIR}/filterrp.f90
                          ${TOOL_DIR}/ensemble/*.f90
                          ${TOOL_DIR}/rpdata/*.f90
                          ${TOOL_DIR}/rpset/*.f90
                          ${TOOL_DIR}/interpolation/*.f90
                          ${TOOL_DIR}/navierstokes/*.f90
                          ${TOOL_DIR}/turbulence/*.f90
                          ${TOOL_DIR}/output/*.f90)

# Build an object library from the Fortran source files
ADD_LIBRARY(               ${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
TARGET_INCLUDE_DIRECTORIES(${TOOL_NAME}libF90 PRIVATE ${TOOL_DIR})
SET_TARGET_PROPERTIES(     ${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                                         COMPILE_OPTIONS "-fPIC"
                                                         CMAKE_POSITION_INDEPENDENT_CODE ON)

# Build a static library from the object library
ADD_LIBRARY(          ${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
TARGET_LINK_LIBRARIES(${TOOL_NAME}lib libposti fftwlib)
ADD_DEPENDENCIES(     ${TOOL_NAME}libF90 libflexishared libposti fftwlib)

# When using Ninja, CMake tries to parallelize the build resulting in missing module files. Add a dummy target to ensure synchronization
ADD_CUSTOM_TARGET(target_${TOOL_NAME} COMMENT "Synchronizing build for ${TOOL_NAME} tool")
ADD_DEPENDENCIES( target_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES( ${TOOL_NAME}libF90 target_${TOOL_NAME})

# Standalone Posti executable
ADD_EXEC(                  posti_${TOOL_NAME} ${TOOL_DIR}/${TOOL_MAIN})
TARGET_INCLUDE_DIRECTORIES(posti_${TOOL_NAME} PRIVATE ${TOOL_DIR})
TARGET_LINK_LIBRARIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(          posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
