# =========================================================================
# toDG
# =========================================================================
SET(TOOL_NAME "toDG")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/toDG")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

# Standalone Posti executable
ADD_EXEC(             posti_${TOOL_NAME} ${TOOL_DIR}/toDG.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} libflexishared)
