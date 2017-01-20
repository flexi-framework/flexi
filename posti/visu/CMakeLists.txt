# =========================================================================
# Visu
# =========================================================================
# collect all F90 files needed for visulib
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/posti/visu/) 
FILE(GLOB visuF90 ./posti/visu/calc.f90
                  ./posti/visu/converttovisu.f90
                  ./posti/visu/mappings.f90
                  ./posti/visu/readstate.f90
                  ./posti/visu/visu3D.f90
                  ./posti/visu/visu3D_vars.f90
                  ./posti/visu/visumesh.f90)
ADD_LIBRARY(visulibF90 OBJECT ${visuF90})
SET_TARGET_PROPERTIES(visulibF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}) 
ADD_LIBRARY(visulib STATIC $<TARGET_OBJECTS:visulibF90>)
ADD_DEPENDENCIES(visulibF90 flexilibF90 HDF5)

# Standalone Posti executable
add_exec(posti_visu3D ./posti/visu/posti_visu3D.f90)
TARGET_LINK_LIBRARIES(posti_visu3D visulib flexilib userblocklib_dummy ${linkedlibs})
ADD_DEPENDENCIES(     posti_visu3D visulib flexilib userblocklib_dummy ${linkedlibs})

# Posti visu plugin for Paraview
CMAKE_DEPENDENT_OPTION(POSTI_USE_PARAVIEW "Build the Paraview visualization plugin for Posti." OFF "POSTI_BUILD_VISU" OFF)
IF (POSTI_USE_PARAVIEW)
  FIND_PACKAGE(ParaView)
ENDIF()
IF (POSTI_USE_PARAVIEW AND ParaView_FOUND)
   INCLUDE(${PARAVIEW_USE_FILE})
   MESSAGE(STATUS ${PARAVIEW_USE_FILE})

   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/posti/visu/paraviewReader/)
   ADD_PARAVIEW_PLUGIN(visu3DReader "1.0"
      SERVER_MANAGER_XML ${CMAKE_CURRENT_SOURCE_DIR}/posti/visu/paraviewReader/visu3DReader.xml
      SERVER_MANAGER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/posti/visu/paraviewReader/visu3DReader.cxx
      REQUIRED_ON_SERVER
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/posti/visu/paraviewReader/visu3DReader.h"
      )

   ADD_DEPENDENCIES(visu3DReader visulib flexilib userblocklib_dummy HDF5)
   TARGET_LINK_LIBRARIES(visu3DReader PUBLIC visulib flexilib ${linkedlibs})
ENDIF()
