cmake_minimum_required (VERSION 2.8.11)
project (DUMMY_FLEXI)
ENABLE_LANGUAGE(Fortran)

OPTION(FLEXI_MPI "Compile SINGLE or MPI version" ON)
OPTION(FLEXI_FV "Enable Finite Volume subcell shock capturing" OFF)

SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")

SET (FLEXI_COMPILE_FLAGS "-xf95-cpp-input")

IF(FLEXI_MPI)
  ADD_DEFINITIONS(-DUSE_MPI=1)                                                                                                       
ELSE()                                                                                                                               
  ADD_DEFINITIONS(-DUSE_MPI=0)                                                                                                       
ENDIF()

IF(FLEXI_FV)
  ADD_DEFINITIONS(-DFV_ENABLED=1)
ELSE()
  ADD_DEFINITIONS(-DFV_ENABLED=0)
ENDIF()

IF(FLEXI_FAIL)
  ADD_DEFINITIONS(-DFAIL=1)
ELSE()
  ADD_DEFINITIONS(-DFAIL=0)
ENDIF()

add_executable (dummy dummy.f90)
SET_TARGET_PROPERTIES(dummy PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS})

