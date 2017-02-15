# =========================================================================
# PREPARE RECORDPOINTS
# =========================================================================
SET(TOOL_NAME "preparerecordpoints")
SET(TOOL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/prepare)

MESSAGE(STATUS "Posti: using prepare record points tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR}) 

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/parameters.f90 
                          ${TOOL_DIR}/parametriccoords.f90 
                          ${TOOL_DIR}/rpset.f90 
                          ${TOOL_DIR}/rpset_vars.f90 
                          ${TOOL_DIR}/blprojection.f90
                          ${TOOL_DIR}/spline.f90
                          ${TOOL_DIR}/hdf5_outputrp.f90
                          ${TOOL_DIR}/visurp.f90 )

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}) 
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 flexilibF90 HDF5)

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/preparerecordpoints.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs})
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs})
