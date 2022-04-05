# =========================================================================
# MergeTimeAverage
# =========================================================================
SET(TOOL_NAME "mergetimeaverages")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/mergetimeaverages")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/timeavg.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} libflexishared)

# =========================================================================
# CalcFluctuations
# =========================================================================
SET(TOOL_NAME "calcfluctuations")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/calcfluc.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} libflexishared)
ADD_DEPENDENCIES(     posti_${TOOL_NAME} libflexishared)
