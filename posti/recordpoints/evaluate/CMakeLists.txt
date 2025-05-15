# =========================================================================
# EVALUATE RECORDPOINTS
# =========================================================================
SET(TOOL_NAME "evaluaterecordpoints")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/evaluate")
SET(TOOL_MAIN "evaluaterecordpoints.f90")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")

# Standalone Posti executable
ADD_EXEC(                  posti_${TOOL_NAME} ${TOOL_DIR}/${TOOL_MAIN})
TARGET_INCLUDE_DIRECTORIES(posti_${TOOL_NAME} PRIVATE ${TOOL_DIR})
TARGET_LINK_LIBRARIES(     posti_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES(          posti_${TOOL_NAME} libflexishared)

# When using Ninja, CMake tries to parallelize the build resulting in missing module files. Add a dummy target to ensure synchronization
ADD_CUSTOM_TARGET(target_${TOOL_NAME} COMMENT "Synchronizing build for ${TOOL_NAME} tool")
ADD_DEPENDENCIES( target_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES(  posti_${TOOL_NAME} target_${TOOL_NAME})
