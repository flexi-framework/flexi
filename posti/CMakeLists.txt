# =========================================================================
# POSTI
# =========================================================================

CMAKE_DEPENDENT_OPTION(POSTI_VISU                  "Build the visu tool"                                ON  "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_MERGETIMEAVERAGES     "Build the time-averaging tool"                      OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_DMD                   "Build the DMD tool"                                 OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_RP_PREPARERE          "Build the Record Point preprocessing tool"          OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_RP_VISUALIZE          "Build the Record Point visualization tool"          OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_RP_EVALUATE           "Build the Record Point aposteriori evaluation tool" OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_SWAPMESH              "Build the Swapmesh tool"                            OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_CHANNEL_FFT           "Build the channel FFT tool"                         OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_WALLDISTANCE          "Build the walldistance tool"                        OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_TO3D                  "Build the 3D extrusion tool"                        OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_AVG2D                 "Build the 2D average tool"                          OFF "POSTI" OFF)


IF(POSTI_RP_VISUALIZE OR POSTI_CHANNEL_FFT)
  # =========================================================================
  # FFTW (FFT tool)
  # =========================================================================
  SET (FFTWDOWNLOAD ${LIBS_DLPATH}libs/FFTW.git )
  SET(POSTI_FFTWTAG "fftw-3.3.9" CACHE STRING "FFTW version tag")
  SET(POSTI_FFTWDIR ${LIBS_EXTERNAL_LIB_DIR}/fftw)
  # Try to find the library in the share directory or in a global, user-defined $FFTW_DIR
  FIND_PATH(FFTW_INCLUDES fftw3.h         HINTS $ENV{FFTW_DIR}/include ${LIBS_EXTERNAL_LIB_DIR}/fftw/build/include)
  FIND_LIBRARY(FFTW_LIBRARIES NAMES fftw3 HINTS $ENV{FFTW_DIR}/lib     ${LIBS_EXTERNAL_LIB_DIR}/fftw/build/lib    )
  IF(FFTW_INCLUDES STREQUAL "FFTW_INCLUDES-NOTFOUND" OR FFTW_LIBRARIES STREQUAL "FFTW_LIBRARIES-NOTFOUND")
    MESSAGE(STATUS "FFTW not found, building library")
      ExternalProject_Add(fftw
        GIT_REPOSITORY ${FFTWDOWNLOAD}
        GIT_TAG ${POSTI_FFTWTAG}
        ${${GITSHALLOW}}
        PREFIX  ${POSTI_FFTWDIR}
	UPDATE_COMMAND ""
        CONFIGURE_COMMAND tar xzf ${POSTI_FFTWTAG}.tar.gz && FC=${CMAKE_Fortran_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_F77_COMPILER} ${POSTI_FFTWDIR}/src/fftw/${POSTI_FFTWTAG}/configure --prefix=${POSTI_FFTWDIR}/build
        BUILD_IN_SOURCE 1
      )
    SET(FFTW_INCLUDE_DIR ${POSTI_FFTWDIR}/build/include)
    SET(FFTW_LIBRARY     ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
    SET(FFTW_LIBRARIES   ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
    SET_PROPERTY(DIRECTORY PROPERTY CLEAN_NO_CUSTOM 1)
    LIST(APPEND SELFBUILTEXTERNALS fftw)
  ELSE()
    MESSAGE(STATUS "FFTW includes found: " ${FFTW_INCLUDES})
    SET(FFTW_INCLUDE_DIR ${FFTW_INCLUDES})
    SET(FFTW_LIBRARY     ${FFTW_LIBRARIES})
  ENDIF()

ENDIF()


IF(POSTI_VISU)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/visu/CMakeLists.txt)
ENDIF()

IF(POSTI_MERGETIMEAVERAGES)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/mergetimeaverages/CMakeLists.txt)
ENDIF()

IF(POSTI_DMD)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/dmd/CMakeLists.txt)
ENDIF()

IF(POSTI_RP_PREPARERE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/prepare/CMakeLists.txt)
ENDIF()

IF(POSTI_RP_VISUALIZERE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/visualize/CMakeLists.txt)
ENDIF()

IF(POSTI_RP_EVALUATE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/evaluate/CMakeLists.txt)
ENDIF()

IF(POSTI_SWAPMESH)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/swapmesh/CMakeLists.txt)
ENDIF()

IF(POSTI_CHANNEL_FFT)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/channel_fft/CMakeLists.txt)
ENDIF()

IF(POSTI_WALLDISTANCE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/walldistance/CMakeLists.txt)
ENDIF()

IF(POSTI_TO3D)
	INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/to3D/CMakeLists.txt)
ENDIF()

IF(POSTI_AVG2D)
	INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/avg2D/CMakeLists.txt)
ENDIF()
