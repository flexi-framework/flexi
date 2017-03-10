# =========================================================================
# PREPARE RECORDPOINTS
# =========================================================================
SET(TOOL_NAME "visualizerecordpoints")
SET(TOOL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/posti/recordpoints/visualize)

MESSAGE(STATUS "Posti: using visualize record points tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/parameters.f90
                          ${TOOL_DIR}/fft/*.f90
                          ${TOOL_DIR}/filterrp.f90
                          ${TOOL_DIR}/rpdata/*.f90
                          ${TOOL_DIR}/rpset/*.f90
                          ${TOOL_DIR}/interpolation/*.f90
                          ${TOOL_DIR}/turbulence/*.f90
                          ${TOOL_DIR}/varnamemappings/*.f90
                          ${TOOL_DIR}/output/*.f90
                          ${TOOL_DIR}/../prepare/spline.f90
                          ${TOOL_DIR}/../../common/fftw/fftw.f90
                          ${TOOL_DIR}/${POSTI_WHICH_EQNSYS}/*.f90 )

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 flexilibF90 HDF5 fftw)
TARGET_INCLUDE_DIRECTORIES(${TOOL_NAME}libF90 PRIVATE ${FFTW_INCLUDE_DIR})

#check whether these two lines are correct (why link externals to lib AND exec??)
SET_TARGET_PROPERTIES(${TOOL_NAME}lib PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(${TOOL_NAME}lib ${POSTI_HDF5_LIBRARY} ${FFTW_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# Standalone Posti executable
add_exec(             posti_${TOOL_NAME} ${TOOL_DIR}/visualizerecordpoints.f90)
TARGET_LINK_LIBRARIES(posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs} ${LAPACK_LIBRARIES} ${CGNS_LIBRARY} ${FFTW_LIBRARIES})
ADD_DEPENDENCIES(     posti_${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs} HDF5 fftw)
SET_TARGET_PROPERTIES(posti_${TOOL_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

