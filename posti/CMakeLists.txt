# =========================================================================
# POSTI
# =========================================================================

# use dummy userblock, as userblocks are often not required for Posti tools
ADD_LIBRARY(userblocklib_dummy SHARED ./src/output/read_userblock_dummy.c)
SET_TARGET_PROPERTIES(userblocklib_dummy PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(userblocklib_dummy PROPERTIES COMPILE_FLAGS "-fPIC")

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
CMAKE_DEPENDENT_OPTION(POSTI_USE_PARAVIEW "Build the Paraview visualization plugin for Posti." OFF "FLEXI_BUILDPOSTI" OFF)
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

# =========================================================================
# FFTW (FFT tool)
# =========================================================================
SET(FFTWDOWNLOAD   ${FLEXI_DLPATH}libs/FFTW.git)
SET(POSTI_FFTWVERSION "3.3.6-pl1" CACHE STRING "FFTW version tag")
SET(POSTI_FFTWDIR ${FLEXI_EXTERNAL_LIB_DIR}/fft)
FIND_PACKAGE(FFTW QUIET)
IF (NOT FFTW_FOUND)
  OPTION(POSTI_BUILDFFTW "Compile and build fftw library" ON)
ELSE()
  OPTION(POSTI_BUILDFFTW "Compile and build fftw library" OFF)
ENDIF()
IF (POSTI_BUILDFFTW)
  IF (NOT EXISTS ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
    ExternalProject_Add(fftw
      GIT_REPOSITORY ${FFTWDOWNLOAD}
      PREFIX ${POSTI_FFTWDIR}
      CONFIGURE_COMMAND tar -xzf fftw-${POSTI_FFTWVERSION}.tar.gz &&
                        FC=${CMAKE_Fortran_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_F77_COMPILER}
                        ${POSTI_FFTWDIR}/src/fftw/fftw-${POSTI_FFTWVERSION}/configure
                        --prefix=${POSTI_FFTWDIR}/build
      BUILD_IN_SOURCE 1
    )
  ENDIF()
  INCLUDE_DIRECTORIES(${POSTI_FFTWDIR}/build/include)
  SET(FFTW_LIBRARY ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
  SET(FFTW_LIBRARIES ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
ELSE()
  INCLUDE_DIRECTORIES (${FFTW_INCLUDE_DIR})
  MESSAGE(STATUS "FFTW include dir: " ${FFTW_INCLUDE_DIR})
ENDIF()

# =========================================================================
# Channel FFT
# =========================================================================
FILE(GLOB channel_fftF90 ./posti/channel_fft/fft_vars.f90
                         ./posti/channel_fft/fft.f90
			 ./posti/channel_fft/fftw/fftw.f90)

ADD_LIBRARY(channel_fftlibF90 OBJECT ${channel_fftF90})
SET_TARGET_PROPERTIES(channel_fftlibF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(channel_fftlib STATIC $<TARGET_OBJECTS:channel_fftlibF90>)
ADD_DEPENDENCIES(channel_fftlibF90 flexilibF90 HDF5 fftw)

# add channel fft executable
add_exec(posti_channel_fft ./posti/channel_fft/channel_fft.f90)
TARGET_LINK_LIBRARIES(posti_channel_fft channel_fftlib flexilib userblocklib_dummy ${FFTW_LIBRARIES} ${linkedlibs})
ADD_DEPENDENCIES(     posti_channel_fft channel_fftlib flexilib userblocklib_dummy FFTW ${linkedlibs})
