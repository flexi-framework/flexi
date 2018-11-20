# =========================================================================
# MergeTimeAverage
# =========================================================================
SET(TOOL_NAME "dmd")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/dmd")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR}) 

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/dmd.f90
                          ${TOOL_DIR}/dmd_vars.f90
                          ${TOOL_DIR}/navierstokes/equationDmd.f90
                          ${TOOL_DIR}/navierstokes/equationDmd_vars.f90)
                          #${TOOL_DIR}/parametriccoords.f90
                          #${TOOL_DIR}/interpolatesolution.f90)

#SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
#ADD_DEPENDENCIES(${TOOL_NAME}libF90 flexilibF90 HDF5)

#ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
#SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
#ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
#ADD_DEPENDENCIES(${TOOL_NAME}libF90 flexilibF90 HDF5)
ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 flexilibF90 HDF5)

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/posti_dmd.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs} ${LAPACK_LIBRARIES})
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs} )
SET_TARGET_PROPERTIES(posti_${TOOL_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

