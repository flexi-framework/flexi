# =========================================================================
# POSTI
# =========================================================================

CMAKE_DEPENDENT_OPTION(POSTI_BUILD_VISU                  "Build the visu tool"                                ON  "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_MERGETIMEAVERAGES     "Build the time-averaging tool"                      ON  "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_DMD                   "Build the DMD tool"                                 OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_PREPARERECORDPOINTS   "Build the Record point preprocessing tool"          ON  "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_VISUALIZERECORDPOINTS "Build the Record Point visualization tool"          OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_EVALUATERECORDPOINTS  "Build the Record Point aposteriori evaluation tool" ON  "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_SWAPMESH              "Build the Swapmesh tool"                            OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_CHANNEL_FFT           "Build the channel FFT tool"                         OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_WALLDISTANCE          "Build the walldistance tool"                        OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_TO3D                  "Build the 3D extrusion tool"                        OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_AVG2D                 "Build the 2D average tool"                          OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_INIT_HIT              "Build the init hit tool"                            OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_FILTER_HIT            "Build the filter hit tool"                          OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(ENABLE_OPENMP                     "Enable OpenMP in posti tools"                       OFF "FLEXI_BUILDPOSTI" OFF)

IF(ENABLE_OPENMP)
  IF ("${CMAKE_VERSION}" VERSION_LESS 3.1.0)
    MESSAGE(WARNING "For finding OpenMP Fortran flags at least CMake version 3.1.0 is required. Please specify flags manually or use newer CMake version.")
  ENDIF()
  FIND_PACKAGE(OpenMP REQUIRED)
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG}   ${OpenMP_Fortran_FLAGS}")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
  SET (CMAKE_CXX_FLAGS_DEBUG       "${CMAKE_CXX_FLAGS_DEBUG}       ${OpenMP_CXX_FLAGS}")
  SET (CMAKE_CXX_FLAGS_RELEASE     "${CMAKE_CXX_FLAGS_RELEASE}     ${OpenMP_CXX_FLAGS}")
  SET (CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS}      ${OpenMP_EXE_LINKER_FLAGS}")
  ADD_DEFINITIONS(-DUSE_OPENMP=1)
ELSE()
  ADD_DEFINITIONS(-DUSE_OPENMP=0)
ENDIF()


SET(   POSTI_WHICH_EQNSYS "navierstokes" CACHE STRING
                        "Used equation system (navierstokes)")
SET(   POSTI_WHICH_EOS    "idealgas" CACHE STRING
                        "Used equation of state (idealgas)")

IF(POSTI_BUILD_VISUALIZERECORDPOINTS OR POSTI_BUILD_CHANNEL_FFT OR POSTI_BUILD_INIT_HIT OR POSTI_BUILD_FILTER_HIT)
  # =========================================================================
  # FFTW (FFT tool)
  # =========================================================================
  # TODO: Upload new version to FlexiDL
  IF (HOMERUN GREATER -1)
    SET (FFTWDOWNLOAD ${FLEXI_DLPATH}libs/FFTW.git )
  ELSE()
    SET (FFTWDOWNLOAD ${FLEXI_DLPATH}libs/FFTW.git )
    #SET (FFTWDOWNLOAD https://github.com/FFTW/fftw3.git ) ! Github repo does not contain all data, use fftw.org instead
    #SET (FFTWDOWNLOAD http://fftw.org/ ) fftw-3.3.6-pl1.tar.gz
  ENDIF()
  SET(POSTI_FFTWTAG "fftw-3.3.6-pl1" CACHE STRING "FFTW version tag")
  SET(POSTI_FFTWDIR ${FLEXI_EXTERNAL_LIB_DIR}/fftw)
  # Try to find the library in the share directory or in a global, user-defined $FFTW_DIR
  FIND_PATH(FFTW_INCLUDES fftw3.h         HINTS $ENV{FFTW_DIR}/include ${FLEXI_EXTERNAL_LIB_DIR}/fftw/build/include)
  FIND_LIBRARY(FFTW_LIBRARIES NAMES fftw3 HINTS $ENV{FFTW_DIR}/lib     ${FLEXI_EXTERNAL_LIB_DIR}/fftw/build/lib    )
  IF (ENABLE_OPENMP)
    FIND_LIBRARY(FFTW_OMP_LIBRARIES NAMES fftw3_omp HINTS $ENV{FFTW_DIR}/lib     ${FLEXI_EXTERNAL_LIB_DIR}/fftw/build/lib    )
  ENDIF()
  IF(FFTW_INCLUDES STREQUAL "FFTW_INCLUDES-NOTFOUND" OR FFTW_LIBRARIES STREQUAL "FFTW_LIBRARIES-NOTFOUND" OR FFTW_OMP_LIBRARIES STREQUAL "FFTW_OMP_LIBRARIES-NOTFOUND")
    MESSAGE(STATUS "FFTW not found, building library")
      ExternalProject_Add(fftw
        GIT_REPOSITORY ${FFTWDOWNLOAD}
        GIT_TAG ${POSTI_FFTWTAG}
        ${${GITSHALLOW}}
        PREFIX  ${POSTI_FFTWDIR}
	UPDATE_COMMAND ""
        CONFIGURE_COMMAND tar xzf ${POSTI_FFTWTAG}.tar.gz && FC=${CMAKE_Fortran_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_F77_COMPILER} ${POSTI_FFTWDIR}/src/fftw/${POSTI_FFTWTAG}/configure --prefix=${POSTI_FFTWDIR}/build --enable-openmp BUILD_IN_SOURCE 1
      )
    SET(FFTW_INCLUDE_DIR ${POSTI_FFTWDIR}/build/include)
    SET(FFTW_LIBRARY     ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
    SET(FFTW_LIBRARIES   ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
    IF (ENABLE_OPENMP)
      SET(FFTW_OMP_LIBRARY   ${POSTI_FFTWDIR}/build/lib/libfftw3_omp.a)
      SET(FFTW_OMP_LIBRARIES ${POSTI_FFTWDIR}/build/lib/libfftw3_omp.a)
    ENDIF()
    SET_PROPERTY(DIRECTORY PROPERTY CLEAN_NO_CUSTOM 1)
    LIST(APPEND SELFBUILTEXTERNALS fftw)
  ELSE()
    MESSAGE(STATUS "FFTW includes found: " ${FFTW_INCLUDES})
    SET(FFTW_INCLUDE_DIR ${FFTW_INCLUDES})
    SET(FFTW_LIBRARY     ${FFTW_LIBRARIES})
    IF (ENABLE_OPENMP)
      SET(FFTW_OMP_LIBRARY   ${POSTI_FFTWDIR}/build/lib/libfftw3_omp.a)
    ENDIF()
  ENDIF()

ENDIF()


IF(POSTI_BUILD_VISU)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/visu/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_MERGETIMEAVERAGES)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/mergetimeaverages/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_DMD)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/dmd/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_PREPARERECORDPOINTS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/prepare/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_VISUALIZERECORDPOINTS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/visualize/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_EVALUATERECORDPOINTS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/evaluate/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_SWAPMESH)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/swapmesh/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_CHANNEL_FFT)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/channel_fft/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_WALLDISTANCE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/walldistance/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_TO3D)
	INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/to3D/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_AVG2D)
	INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/avg2D/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_INIT_HIT)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/init_hit/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_FILTER_HIT)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/filter_hit/CMakeLists.txt)
ENDIF()
