# =========================================================================
# CHANNEL_FFT
# =========================================================================
SET(TOOL_NAME "init_hit")
SET(TOOL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/posti/init_hit)

MESSAGE(STATUS "Posti: using init_hit tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/init_hit.f90
                          ${TOOL_DIR}/init_hit_vars.f90
                          ${TOOL_DIR}/fft.f90
                          ${TOOL_DIR}/../common/fftw/fftw.f90)

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 flexilibF90 HDF5 fftw)
TARGET_INCLUDE_DIRECTORIES(${TOOL_NAME}libF90 PRIVATE ${FFTW_INCLUDE_DIR})

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/init_hit.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${FFTW_LIBRARIES} ${linkedlibs})
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy FFTW ${linkedlibs})
