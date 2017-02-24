# =========================================================================
# VISU
# =========================================================================
SET(TOOL_NAME "visu")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/visu")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR}) 

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/avg2d.f90
                          ${TOOL_DIR}/calc.f90
                          ${TOOL_DIR}/CMakeLists.txt
                          ${TOOL_DIR}/converttovisu.f90
                          ${TOOL_DIR}/mappings.f90
                          ${TOOL_DIR}/paraviewReader/
                          ${TOOL_DIR}/pluginTypes_visu.h
                          ${TOOL_DIR}/plugin_visu.h
                          ${TOOL_DIR}/posti_visu.f90
                          ${TOOL_DIR}/readstate.f90
                          ${TOOL_DIR}/visu_Cwrapper.f90
                          ${TOOL_DIR}/visu.f90
                          ${TOOL_DIR}/visumesh.f90
                          ${TOOL_DIR}/visu_vars.f90)


ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}) 
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 flexilibF90 HDF5)

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/posti_visu.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs})
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs})

# Posti visu plugin for Paraview
CMAKE_DEPENDENT_OPTION(POSTI_USE_PARAVIEW "Build the Paraview visualization plugin for Posti." OFF "POSTI_BUILD_VISU" OFF)
IF (POSTI_USE_PARAVIEW)
  FIND_PACKAGE(ParaView)
ENDIF()
IF (POSTI_USE_PARAVIEW AND ParaView_FOUND)
   INCLUDE(${PARAVIEW_USE_FILE})
   MESSAGE(STATUS ${PARAVIEW_USE_FILE})

   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   INCLUDE_DIRECTORIES(${TOOL_DIR}/paraviewReader/)
   ADD_PARAVIEW_PLUGIN(visu3DReader "1.0"
     SERVER_MANAGER_XML     ${TOOL_DIR}/paraviewReader/visuReader.xml
     SERVER_MANAGER_SOURCES ${TOOL_DIR}/paraviewReader/visuReader.cxx
     REQUIRED_ON_SERVER
     DEPENDS               "${TOOL_DIR}/paraviewReader/visuReader.h"
     )

   ADD_DEPENDENCIES(     visu3DReader ${TOOL_NAME}lib flexilib userblocklib_dummy HDF5)
   TARGET_LINK_LIBRARIES(visu3DReader PUBLIC ${TOOL_NAME}lib flexilib ${linkedlibs})
ENDIF()
