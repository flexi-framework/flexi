# =========================================================================
# VISU
# =========================================================================
SET(TOOL_NAME "visu")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/visu")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/avg2d/avg2d.f90
                          ${TOOL_DIR}/io_hdf5/visu_hdf5_output.f90
                          ${TOOL_DIR}/mesh/mesh.f90
                          ${TOOL_DIR}/calc.f90
                          ${TOOL_DIR}/CMakeLists.txt
                          ${TOOL_DIR}/converttovisu.f90
                          ${TOOL_DIR}/mappings.f90
                          ${TOOL_DIR}/paraviewReader/Reader/
                          ${TOOL_DIR}/pluginTypes_visu.h
                          ${TOOL_DIR}/plugin_visu.h
                          ${TOOL_DIR}/posti_visu.f90
                          ${TOOL_DIR}/visu.f90
                          ${TOOL_DIR}/visu_Cwrapper.f90
                          ${TOOL_DIR}/visu_init.f90
                          ${TOOL_DIR}/visu_mesh.f90
                          ${TOOL_DIR}/visu_state.f90
                          ${TOOL_DIR}/visu_vars.f90)

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}
                                                    COMPILE_OPTIONS "-fPIC"
                                                    CMAKE_POSITION_INDEPENDENT_CODE ON)
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared)

# When using Ninja, CMake tries to parallelize the build resulting in missing module files. Add a dummy target to ensure synchronization
ADD_CUSTOM_TARGET(${TOOL_NAME}target)
ADD_DEPENDENCIES(${TOOL_NAME}target libflexishared)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 ${TOOL_NAME}target)

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/posti_visu.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib libflexishared)

# Posti visu plugin for Paraview
CMAKE_DEPENDENT_OPTION(POSTI_VISU_PARAVIEW "Build the Paraview visualization plugin for Posti." OFF "POSTI_VISU" OFF)
IF (POSTI_VISU_PARAVIEW)
  # OSPRay causes developer warnings because of duplicate 'optimized' flag. CMake fixes itself, so suppress it for the time being
  # Unfortunately, we have to suppress them before actually finding Paraview since the find itself causes the warnings
  MARK_AS_ADVANCED(FORCE CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
  SET(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE BOOL "Suppress warnings that are meant for the author of the CMakeLists.txt files" FORCE)

  # Required to build visuReader with ParaView >= 5.7
  INCLUDE(GNUInstallDirs)

  # If ParaView is compiled with Ospray, it must be found first
  FIND_PACKAGE(TBB QUIET)
  IF(TBB_FOUND)
    MESSAGE(STATUS "Found TBB version ${TBB_VERSION}")

    # ParaView cannot find the correct path to TBB, so give it some help
    FIND_PATH (TBB_ROOT_DIR
      NAMES include/tbb/tbb.h
      PATHS ENV TBBROOT
            ENV TBB40_INSTALL_DIR
            ENV TBB30_INSTALL_DIR
            ENV TBB22_INSTALL_DIR
            ENV TBB21_INSTALL_DIR
            ENV TBB_ROOT_DIR
      DOC "TBB root directory")

    SET(TBB_ROOT ${TBB_ROOT_DIR})
  ENDIF()

  # However, if we locate OSPRay before ParaView, ParaView complains because it wants to add all of it to its path
  # FIND_PACKAGE(ospray QUIET)
  # IF(ospray_FOUND)
  #   TARGET_LINK_LIBRARIES(${client_target} ospray::ospray)
  #   MESSAGE(STATUS "Found Ospray version ${OSPRAY_VERSION}")
  # ENDIF()

  FIND_PACKAGE(ParaView)
ENDIF()

IF (POSTI_VISU_PARAVIEW AND ParaView_FOUND)
   MESSAGE(STATUS "Found ParaView version ${ParaView_VERSION}")
   SET(POSTI_VISU_PARAVIEW_SUFFIX "" CACHE STRING "Pass custom suffix to plugin to be able to load several plugins simultaneously")

  # High-order elements only supported with ParaView 5.5+
  # IF(${ParaView_VERSION} VERSION_LESS "5.5")
  #   ADD_COMPILE_DEFINITIONS(PV_HIGHORDER=0)
  # ELSE()
  #   ADD_COMPILE_DEFINITIONS(PV_HIGHORDER=1)
  # ENDIF()

   # For Paraview >= 5.7, we need to stash the files one directory deeper
   IF(${ParaView_VERSION} VERSION_LESS "5.7")
     IF(POSTI_VISU_PARAVIEW_SUFFIX STREQUAL "")
        SET(TOOL_DIR_PV ${TOOL_DIR}/paraviewReader/Reader)
     ELSE()
       SET(TOOL_DIR_PV ${CMAKE_CURRENT_BINARY_DIR})
       FILE(WRITE ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.h   "")
       FILE(WRITE ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.cxx "")
       FILE(WRITE ${TOOL_DIR_PV}/visuReader${POSTI_VISU_PARAVIEW_SUFFIX}.xml "")
     ENDIF()
   ELSE()
     SET(TOOL_DIR_PV ${TOOL_DIR}/paraviewReader)
   ENDIF()

   # Only required for ParaView <= 5.6
   IF(${ParaView_VERSION} VERSION_LESS "5.7")
   INCLUDE(${PARAVIEW_USE_FILE})
   MESSAGE(STATUS ${PARAVIEW_USE_FILE})
   ELSE()
     SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
     SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
     SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
     SET(BUILD_SHARED_LIBS ON)
   ENDIF()

   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   INCLUDE_DIRECTORIES(${TOOL_DIR}/paraviewReader/Reader/)

   # Changed beginning form Paraview 5.7
   IF(${ParaView_VERSION} VERSION_LESS "5.7")
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

   ELSE()
     PARAVIEW_PLUGIN_SCAN(
       PLUGIN_FILES "${TOOL_DIR_PV}/visuReader.plugin"
       PROVIDES_PLUGINS plugins
       REQUIRES_MODULES  required_modules
       ENABLE_BY_DEFAULT ON)

     PARAVIEW_PLUGIN_BUILD(
       PLUGINS ${plugins})
   ENDIF()

   ADD_DEPENDENCIES(     visuReader${POSTI_VISU_PARAVIEW_SUFFIX} ${TOOL_NAME}lib libflexishared)
   TARGET_LINK_LIBRARIES(visuReader${POSTI_VISU_PARAVIEW_SUFFIX} PUBLIC ${TOOL_NAME}lib libflexishared)

  # Enable VTKm optimizations
  SET(VTKm_Vectorization "native" CACHE STRING "Level of compiler vectorization support" FORCE)

   # ParaView really adds a lot of CMake variables
   # > Adios2
   MARK_AS_ADVANCED(FORCE ADIOS2_DIR)
   # > BLAS
   MARK_AS_ADVANCED(FORCE libLAS_DIR)
   # > Boost
   MARK_AS_ADVANCED(FORCE Boost_FILESYSTEM_LIBRARY_RELEASE)
   MARK_AS_ADVANCED(FORCE Boost_INCLUDE_DIR)
   MARK_AS_ADVANCED(FORCE Boost_IOSTREAMS_LIBRARY_RELEASE)
   MARK_AS_ADVANCED(FORCE Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE)
   MARK_AS_ADVANCED(FORCE Boost_SYSTEM_LIBRARY_RELEASE)
   MARK_AS_ADVANCED(FORCE Boost_THREAD_LIBRARY_RELEASE)
   # > CLI11
   MARK_AS_ADVANCED(FORCE CLI11_DIR)
   # > Catalyst
   MARK_AS_ADVANCED(FORCE catalyst_DIR)
   # > Embree
   MARK_AS_ADVANCED(FORCE embree_DIR)
   # > FastFloat
   MARK_AS_ADVANCED(FORCE FastFloat_DIR)
   # > FreeType
   MARK_AS_ADVANCED(FORCE freetype_DIR)
   # > HDF5
   MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_dl)
   MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_hdf5)
   MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_hdf5_hl)
   MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_m)
   MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_sz)
   MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_z)
   # > ICU
   MARK_AS_ADVANCED(FORCE ICU_DERB_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENBRK_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENCCODE_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENCFU_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENCMN_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENCNVAL_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENDICT_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENNORM2_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENRB_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_GENSPREP_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_ICU-CONFIG_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_ICUINFO_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_ICUPKG_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_INCLUDE_DIR)
   MARK_AS_ADVANCED(FORCE ICU_MAKECONV_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_MAKEFILE_INC)
   MARK_AS_ADVANCED(FORCE ICU_PKGDATA_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE ICU_PKGDATA_INC)
   MARK_AS_ADVANCED(FORCE ICU_UCONV_EXECUTABLE)
   # > JSON
   MARK_AS_ADVANCED(FORCE nlohmann_json_DIR)
   # > LZMA
   MARK_AS_ADVANCED(FORCE LZMA_INCLUDE_DIR)
   MARK_AS_ADVANCED(FORCE LZMA_LIBRARY)
   # > netCDF
   MARK_AS_ADVANCED(FORCE netCDF_DIR)
   # > ParaView
   MARK_AS_ADVANCED(FORCE ParaView_DIR)
   MARK_AS_ADVANCED(FORCE PARAVIEW_PLUGIN_ENABLE_visuReader)
   # > PDAL
   MARK_AS_ADVANCED(FORCE PDAL_DIR)
   # > pugixml
   MARK_AS_ADVANCED(FORCE pugixml_DIR)
   # > OpenGL
   MARK_AS_ADVANCED(FORCE OPENGL_GLES2_INCLUDE_DIR)
   MARK_AS_ADVANCED(FORCE OPENGL_GLES3_INCLUDE_DIR)
   # > OpenVKL
   MARK_AS_ADVANCED(FORCE openvkl_DIR)
   # > OSPRay
   MARK_AS_ADVANCED(FORCE ospray_DIR)
   # > OIDN
   MARK_AS_ADVANCED(FORCE OpenImageDenoise_DIR)
   # > Qt
   MARK_AS_ADVANCED(FORCE QT_IMPORTS_DIR)
   MARK_AS_ADVANCED(FORCE QT_QMAKE_EXECUTABLE)
   MARK_AS_ADVANCED(FORCE Qt5_DIR)
   MARK_AS_ADVANCED(FORCE Qt5Core_DIR)
   MARK_AS_ADVANCED(FORCE Qt5Gui_DIR)
   MARK_AS_ADVANCED(FORCE Qt5Help_DIR)
   MARK_AS_ADVANCED(FORCE Qt5Network_DIR)
   MARK_AS_ADVANCED(FORCE Qt5OpenGL_DIR)
   MARK_AS_ADVANCED(FORCE Qt5Sql_DIR)
   MARK_AS_ADVANCED(FORCE Qt5Svg_DIR)
   MARK_AS_ADVANCED(FORCE Qt5Widgets_DIR)
   MARK_AS_ADVANCED(FORCE Qt5X11Extras_DIR)
   # > rkcommon
   MARK_AS_ADVANCED(FORCE rkcommon_DIR)
   # > TBB
   MARK_AS_ADVANCED(FORCE TBB_DIR)
   MARK_AS_ADVANCED(FORCE TBB_ROOT)
   MARK_AS_ADVANCED(FORCE TBB_ROOT_DIR)
   # > TIFF
   MARK_AS_ADVANCED(FORCE tiff_DIR)
   MARK_AS_ADVANCED(FORCE Tiff_DIR)
   # > VTK
   MARK_AS_ADVANCED(FORCE VTK_DIR)
   MARK_AS_ADVANCED(FORCE VTK_MPI_NUMPROCS)
   MARK_AS_ADVANCED(FORCE VTKm_DIR)
   MARK_AS_ADVANCED(FORCE VTKm_Vectorization)
   # > Verdict
   MARK_AS_ADVANCED(FORCE Verdict_DIR)
   # > XCB
   MARK_AS_ADVANCED(FORCE X11_xcb_xkb_INCLUDE_PATH)
ENDIF()
