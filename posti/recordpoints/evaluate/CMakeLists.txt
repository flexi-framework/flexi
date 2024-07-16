# =========================================================================
# EVALUATE RECORDPOINTS
# =========================================================================
SET(TOOL_NAME "evaluaterecordpoints")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/evaluate")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

# Standalone Posti executable
ADD_EXEC(             posti_${TOOL_NAME} ${TOOL_DIR}/evaluaterecordpoints.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} libflexishared)

# When using Ninja, CMake tries to parallelize the build resulting in missing module files. Add a dummy target to ensure synchronization
ADD_CUSTOM_TARGET(${TOOL_NAME}target)
ADD_DEPENDENCIES(${TOOL_NAME}target libflexishared)
ADD_DEPENDENCIES(posti_${TOOL_NAME} ${TOOL_NAME}target)
