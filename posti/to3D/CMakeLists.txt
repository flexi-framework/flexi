# =========================================================================
# to3D
# =========================================================================
SET(TOOL_NAME "to3D")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/to3D")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

# Standalone Posti executable
ADD_EXEC(             posti_${TOOL_NAME} ${TOOL_DIR}/posti_to3D.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} libflexishared)
