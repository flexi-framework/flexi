# =========================================================================
# POSTI
# =========================================================================
CMAKE_DEPENDENT_OPTION(POSTI_VISU                  "Build the visu tool"                                ON  "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_MERGETIMEAVERAGES     "Build the time-averaging tool"                      OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_DMD                   "Build the DMD tool"                                 OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_RP_PREPARE            "Build the Record Point preprocessing tool"          OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_RP_VISUALIZE          "Build the Record Point visualization tool"          OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_RP_EVALUATE           "Build the Record Point aposteriori evaluation tool" OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_SWAPMESH              "Build the Swapmesh tool"                            OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_CHANNEL_FFT           "Build the channel FFT tool"                         OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_WALLDISTANCE          "Build the walldistance tool"                        OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_TO3D                  "Build the 3D extrusion tool"                        OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_TODG                  "Build the DG conversion tool"                       OFF "POSTI;FLEXI_FV" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_AVG2D                 "Build the 2D average tool"                          OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_HIT_INIT              "Build the HIT initialization tool"                  OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_HIT_FILTER            "Build the HIT spectral filter tool"                 OFF "POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_HIT_ANALYZE           "Build the HIT analyze tool"                         OFF "POSTI" OFF)

# =========================================================================
# Build shared posti library
# =========================================================================
# collect F90 files for the common Posti modules
FILE(GLOB_RECURSE postiF90 ./posti/common/vtkOutput/*.f90
                           ./posti/common/newton/*.f90
                           ./posti/common/spline/*.f90)
ADD_LIBRARY(libposti STATIC ${postiF90})
SET_TARGET_PROPERTIES(libposti PROPERTIES OUTPUT_NAME "libposti")
SET_TARGET_PROPERTIES(libposti PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
TARGET_LINK_LIBRARIES(libposti libflexishared)
ADD_DEPENDENCIES(libposti libflexishared)

# =========================================================================
# FFTW library
# =========================================================================
IF(POSTI_RP_VISUALIZE OR POSTI_CHANNEL_FFT OR POSTI_HIT_INIT OR POSTI_HIT_FILTER OR POSTI_HIT_ANALYZE)
  SET (FFTWDOWNLOAD ${LIBS_DLPATH}libs/FFTW.git )
  SET(POSTI_FFTWTAG "fftw-3.3.9" CACHE STRING "FFTW version tag")
  IF (LIBS_USE_OPENMP)
    SET(POSTI_FFTWDIR ${LIBS_EXTERNAL_LIB_DIR}/fftw_omp)
  ELSE()
    SET(POSTI_FFTWDIR ${LIBS_EXTERNAL_LIB_DIR}/fftw)
  ENDIF()
  # Try to find the library in the share directory or in a global, user-defined $FFTW_DIR
  FIND_PATH(FFTW_INCLUDES fftw3.h       HINTS $ENV{FFTW_DIR}/include ${LIBS_EXTERNAL_LIB_DIR}/fftw/build/include)
  FIND_LIBRARY(FFTW_LIBRARY NAMES fftw3 HINTS $ENV{FFTW_DIR}/lib     ${LIBS_EXTERNAL_LIB_DIR}/fftw/build/lib    )
  IF (LIBS_USE_OPENMP)
    FIND_LIBRARY(FFTW_OMP_LIBRARY NAMES fftw3_omp HINTS $ENV{FFTW_DIR}/lib ${LIBS_EXTERNAL_LIB_DIR}/fftw/build/lib)
  ENDIF()
  IF(FFTW_INCLUDES STREQUAL "FFTW_INCLUDES-NOTFOUND" OR  FFTW_LIBRARY STREQUAL "FFTW_LIBRARY-NOTFOUND"
		                                     OR (FFTW_OMP_LIBRARY STREQUAL "FFTW_OMP_LIBRARY-NOTFOUND" AND LIBS_USE_OPENMP))
    SET(FFTW_FLAGS tar xzf ${POSTI_FFTWTAG}.tar.gz && FC=${CMAKE_Fortran_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_F77_COMPILER} ${POSTI_FFTWDIR}/src/fftw/${POSTI_FFTWTAG}/configure --prefix=${POSTI_FFTWDIR}/build)
    IF (LIBS_USE_OPENMP)
      LIST(APPEND FFTW_FLAGS --enable-openmp)
    ENDIF()
    MESSAGE(STATUS "FFTW not found, building library")
      ExternalProject_Add(fftw
        GIT_REPOSITORY ${FFTWDOWNLOAD}
        GIT_TAG ${POSTI_FFTWTAG}
        ${${GITSHALLOW}}
        PREFIX  ${POSTI_FFTWDIR}
        UPDATE_COMMAND ""
	CONFIGURE_COMMAND ${FFTW_FLAGS}
        BUILD_IN_SOURCE 1
      )
    SET(FFTW_INCLUDE_DIR ${POSTI_FFTWDIR}/build/include)
    SET(FFTW_LIBRARY     ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
    SET(FFTW_OMP_LIBRARY ${POSTI_FFTWDIR}/build/lib/libfftw3_omp.a)
    SET_PROPERTY(DIRECTORY PROPERTY CLEAN_NO_CUSTOM 1)
    LIST(APPEND SELFBUILTEXTERNALS fftw)
  ELSE()
    MESSAGE(STATUS "FFTW includes found: " ${FFTW_INCLUDES})
    SET(FFTW_INCLUDE_DIR ${FFTW_INCLUDES})
  ENDIF()

  SET(FFTW_LIBRARIES   ${FFTW_LIBRARY})
  IF (LIBS_USE_OPENMP)
    LIST(APPEND FFTW_LIBRARIES ${FFTW_OMP_LIBRARY})
  ENDIF()

  # Build static FFTW library to be included in each POSTI tool
  ADD_LIBRARY(fftwlib STATIC ${CMAKE_CURRENT_SOURCE_DIR}/posti/common/fftw/fftw.f90)
  SET_TARGET_PROPERTIES(fftwlib PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
  TARGET_LINK_LIBRARIES(fftwlib ${FFTW_LIBRARIES})
  TARGET_INCLUDE_DIRECTORIES(fftwlib PRIVATE ${FFTW_INCLUDE_DIR})
  IF(SELFBUILTEXTERNALS)
    ADD_DEPENDENCIES(fftwlib ${SELFBUILTEXTERNALS})
  ENDIF()
  MARK_AS_ADVANCED(FORCE FFTW_INCLUDES)
  MARK_AS_ADVANCED(FORCE FFTW_LIBRARY)
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

IF(POSTI_RP_PREPARE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/prepare/CMakeLists.txt)
ENDIF()

IF(POSTI_RP_VISUALIZE)
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

IF(POSTI_TODG)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/toDG/CMakeLists.txt)
ENDIF()

IF(POSTI_AVG2D)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/avg2D/CMakeLists.txt)
ENDIF()

IF(POSTI_HIT_INIT OR POSTI_HIT_FILTER OR POSTI_HIT_ANALYZE)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/hit/CMakeLists.txt)
ENDIF()
