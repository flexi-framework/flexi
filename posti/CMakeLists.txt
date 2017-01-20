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
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_PREPARERECORDPOINTS   "Build the record point preprocessing tool"  ON  "FLEXI_BUILDPOSTI" OFF)
#CMAKE_DEPENDENT_OPTION(POSTI_BUILD_VISUALIZERECORDPOINTS "Build the Record Point Postprocessing tool" ON  "FLEXI_BUILDPOSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_EVALUATERECORDPOINTS  "Build the Record Point Postprocessing tool" ON  "FLEXI_BUILDPOSTI" OFF)

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

