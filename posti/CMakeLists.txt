# =========================================================================
# POSTI
# =========================================================================

# use dummy userblock, as userblocks are often not required for Posti tools
ADD_LIBRARY(userblocklib_dummy SHARED ./src/output/read_userblock_dummy.c)
SET_TARGET_PROPERTIES(userblocklib_dummy PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(userblocklib_dummy PROPERTIES COMPILE_FLAGS "-fPIC")

CMAKE_DEPENDENT_OPTION(POSTI_BUILD_VISU                  "Build the visu tool"                        ON  "FLEXI_BUILDPOSTI" OFF)
#CMAKE_DEPENDENT_OPTION(POSTI_BUILD_MERGETIMEAVERAGE      "Build the time-averaging tool"              ON  "FLEXI_BUILDPOSTI" OFF)
#CMAKE_DEPENDENT_OPTION(POSTI_BUILD_FLUCTUATIONS          "Build the fluctuations tool"                ON  "FLEXI_BUILDPOSTI" OFF)
#CMAKE_DEPENDENT_OPTION(POSTI_BUILD_SWAPMESH              "Build the swapmesh tool"                    ON  "FLEXI_BUILDPOSTI" OFF)
#CMAKE_DEPENDENT_OPTION(POSTI_BUILD_WALLDISTANCE          "Build the walldistance tool"                OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_PREPARERECORDPOINTS   "Build the Record point preprocessing tool"          ON  "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_VISUALIZERECORDPOINTS "Build the Record Point post-processing tool"        OFF "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_EVALUATERECORDPOINTS  "Build the Record Point aposteriori evaluation tool" ON  "FLEXI_BUILDPOSTI" OFF)


SET(   POSTI_WHICH_EQNSYS "navierstokes" CACHE STRING
                        "Used equation system (navierstokes)")
SET(   POSTI_WHICH_EOS    "idealgas" CACHE STRING
                        "Used equation of state (idealgas)")

IF(POSTI_BUILD_VISUALIZERECORDPOINTS)
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
        GIT_TAG ${POSTI_FFTWTAG}
        PREFIX  ${POSTI_FFTWDIR}
        CONFIGURE_COMMAND tar xzf ${POSTI_FFTWTAG}.tar.gz && FC=${CMAKE_Fortran_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_F77_COMPILER} ${POSTI_FFTWDIR}/src/fftw/${POSTI_FFTWTAG}/configure --prefix=${POSTI_FFTWDIR}/build
        BUILD_IN_SOURCE 1
      )
    ENDIF()
    SET(FFTW_INCLUDE_DIR ${POSTI_FFTWDIR}/build/include)
    SET(FFTW_LIBRARY     ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
    SET(FFTW_LIBRARIES   ${POSTI_FFTWDIR}/build/lib/libfftw3.a)
  ELSE()
    MESSAGE(STATUS "FFTW include dir: " ${FFTW_INCLUDE_DIR})
  ENDIF()
ENDIF()


IF(POSTI_BUILD_VISU)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/posti/visu/CMakeLists.txt)
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


