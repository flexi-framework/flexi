# =========================================================================
# EVALUATE RECORDPOINTS
# =========================================================================
SET(TOOL_NAME "evaluaterecordpoints")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/evaluate")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/evaluaterecordpoints.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} libflexishared)
