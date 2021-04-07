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
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared)

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/posti_visu.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)

# Posti visu plugin for Paraview
CMAKE_DEPENDENT_OPTION(POSTI_VISU_PARAVIEW "Build the Paraview visualization plugin for Posti." OFF "POSTI_VISU" OFF)
IF (POSTI_VISU_PARAVIEW)
  FIND_PACKAGE(ParaView)
ENDIF()
IF (POSTI_VISU_PARAVIEW AND ParaView_FOUND)
   SET(POSTI_VISU_PARAVIEW_SUFFIX "" CACHE STRING "Pass custom suffix to plugin to be able to load several plugins simultaneously")
   IF(POSTI_VISU_PARAVIEW_SUFFIX STREQUAL "")
     SET(TOOL_DIR_PV ${TOOL_DIR}/paraviewReader)
   ELSE()
     SET(TOOL_DIR_PV ${CMAKE_CURRENT_BINARY_DIR})
     FILE(WRITE ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.h   "")
     FILE(WRITE ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.cxx "")
     FILE(WRITE ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.xml "")
   ENDIF()
   INCLUDE(${PARAVIEW_USE_FILE})
   MESSAGE(STATUS ${PARAVIEW_USE_FILE})

   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   INCLUDE_DIRECTORIES(${TOOL_DIR}/paraviewReader/)
   ADD_PARAVIEW_PLUGIN(visuReader${POSTI_VISU_PARAVIEW_SUFFIX} "1.0"
      SERVER_MANAGER_XML     ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.xml
      SERVER_MANAGER_SOURCES ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.cxx
     REQUIRED_ON_SERVER
     DEPENDS                "${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.h" 
     )

  IF(NOT POSTI_VISU_PARAVIEW_SUFFIX STREQUAL "")
     SET(EXTENSIONS "h" "xml" "cxx")
     FOREACH(EX ${EXTENSIONS})
       SET(OLDFILE ${TOOL_DIR}/paraviewReader/visuReader.${EX})
       SET(NEWFILE ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.${EX})
       ADD_CUSTOM_TARGET(
         ${EX}_COPY
         COMMAND ${CMAKE_COMMAND} -E copy ${OLDFILE} ${NEWFILE}
         COMMAND sed -i 's/visuReader/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}/g' ${NEWFILE}
         COMMAND sed -i 's/2 Posti Visu/2 Posti Visu ${POSTI_VISU_PARAVIEW_SUFFIX}/g' ${NEWFILE}
         )
      ADD_DEPENDENCIES(visuReader${POSTI_VISU_PARAVIEW_SUFFIX} ${EX}_COPY)
     ENDFOREACH()
   ENDIF()

   ADD_DEPENDENCIES(     visuReader${POSTI_VISU_PARAVIEW_SUFFIX} ${TOOL_NAME}lib libflexishared)
   TARGET_LINK_LIBRARIES(visuReader${POSTI_VISU_PARAVIEW_SUFFIX} PUBLIC ${TOOL_NAME}lib libflexishared)


   mark_as_advanced([CLEAR|FORCE] ParaView_DIR)
   mark_as_advanced([CLEAR|FORCE] QT_IMPORTS_DIR)
   mark_as_advanced([CLEAR|FORCE] QT_QMAKE_EXECUTABLE)
ENDIF()
