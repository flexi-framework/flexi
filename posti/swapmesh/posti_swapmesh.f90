!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!===================================================================================================================================
!> Tool to non-conservatively interpolate the solution from one mesh to another.
!> General process is as follows:
!>   * Read in the mesh coordinates of the old and the new mesh and store as CL points.
!>   * Search for the parametric coordinates and the elemID in the old mesh of all interpolation points of the new mesh.
!>     The interpolation points are CL points of a polynomial degree NInter which is usually higher that the new solution
!>     polynomial degree NNew.
!>     For the search, a Newton algorithm is employed. To get a good starting point for the Newton algorithm, the old mesh
!>     is supersampled on a polynomial degree NSuper.
!>     The found reference coordinates are marked as invalid if the reference coordinates are outside of [-1,1] more than
!>     maxTol and the code will abort if they are outside more than abortTol.
!>   * The old state files will then be interpolated to the new interpolation nodes and transformed to the new polynomial degree
!>     NNew. If some points are marked as invalid and a reference solution has been specified, this reference solution will be
!>     used at those points.
!===================================================================================================================================
PROGRAM swapMesh
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_SwapMesh_Vars
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_SwapMesh,                ONLY: InitSwapmesh,ReadOldStateFile,WriteNewStateFile,FinalizeSwapMesh
USE MOD_InterpolateSolution,     ONLY: InterpolateSolution
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg
INTEGER                            :: limit,nTotalNew,nTotalOld
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

! Define parameters needed
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()

CALL prms%SetSection("swapMesh")
CALL prms%CreateStringOption(   "MeshFileOld"        , "Old mesh file (if different than the one found in the state file)")
CALL prms%CreateStringOption(   "MeshFileNew"        , "New mesh file")
CALL prms%CreateLogicalOption(  "useCurvedsOld"      , "Controls usage of high-order information in old mesh. Turn off to discard "//&
                                                       "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption(  "useCurvedsNew"      , "Controls usage of high-order information in new mesh. Turn off to discard "//&
                                                       "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateIntOption(      "NInter"             , "Polynomial degree used for interpolation on new mesh (should be equal or  "//&
                                                       "higher than NNew) - the state will be interpolated to this degree and then "//& 
                                                       "projected down to NNew")
CALL prms%CreateIntOption(      "NNew"               , "Polynomial degree used in new state files")
CALL prms%CreateIntOption(      "NSuper"             , "Polynomial degree used for supersampling on the old mesh, used to get an "//&
                                                       "initial guess for Newton's method - should be higher than NGeo of old mesh")
CALL prms%CreateRealOption(     "maxTolerance"       , "Tolerance used to mark points as invalid if outside of reference element "//&
                                                       "more than maxTolerance",'5.e-2')
CALL prms%CreateLogicalOption(  "printTroublemakers" , "Turn output of not-found points on or off",'.TRUE.')
CALL prms%CreateRealArrayOption("RefState"           , "If a RefState is defined, this state will be used at points that are "// &
                                                        "marked as invalid - without a RefState, the program will abort in this case")
CALL prms%CreateRealOption(     "abortTolerance"     , "Tolerance used to decide if the program should abort if no "// &
                                                       "RefState is given")
CALL prms%CreateLogicalOption(  "ExtrudeTo3D"        , "Perform an extrusion of a one-layer mesh to the 3D version",'.FALSE.')
CALL prms%CreateIntOption(      "ExtrudeK"           , "Layer which is used in extrusion")

! Parse parameters
! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.2).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: swapmesh prm-file statefile [statefiles]')
END IF
CALL prms%read_options(Args(1))


SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(A)') &
"     _______.____    __    ____  ___      .______   .___  ___.  _______     _______. __    __  "
SWRITE(UNIT_stdOut,'(A)') &
"    /       |\\   \\  /  \\  /   / /   \\     |   _  \\  |   \\/   | |   ____|   /       ||  |  |  | "
SWRITE(UNIT_stdOut,'(A)') &
"   |   (----` \\   \\/    \\/   / /  ^  \\    |  |_)  | |  \\  /  | |  |__     |   (----`|  |__|  | "
SWRITE(UNIT_stdOut,'(A)') &
"    \\   \\      \\            / /  /_\\  \\   |   ___/  |  |\\/|  | |   __|     \\   \\    |   __   | "
SWRITE(UNIT_stdOut,'(A)') &
".----)   |      \\    /\\    / /  _____  \\  |  |      |  |  |  | |  |____.----)   |   |  |  |  | "
SWRITE(UNIT_stdOut,'(A)') &
"|_______/        \\__/  \\__/ /__/     \\__\\ | _|      |__|  |__| |_______|_______/    |__|  |__| "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')

! Initialization
CALL InitIOHDF5()
CALL InitSwapmesh()
#if USE_MPI
CALL InitMPIvars()
#endif

#if USE_MPI
nTotalNew=REAL(nVar_State*(NNew+1)**3*nElemsNew)
nTotalOld=REAL(nVar_State*(NState+1)**3*nElemsOld)
!limit=(2**31-1)/8.
limit=2**28-1/8. ! max. 32 bit integer / 8
IF((nTotalNew.GT.limit).OR.(nTotalNew.GT.limit))THEN
  WRITE(UNIT_StdOut,'(A,F13.0,A)')' New or old state file size is too big! Total array size may not exceed', limit, ' entries!'
  WRITE(UNIT_StdOut,'(A)')' Lower number of elements or NNew! Alternative: compile swapmesh without MPI'
  STOP
END IF
#endif

! Evaluate solution at new solution nodes
DO iArg=2,nArgs
  ! Check if a .h5 file has been given to the swapmesh tool
  IF (.NOT.(STRICMP(GetFileExtension(Args(iArg)),'h5'))) THEN
    CALL CollectiveStop(__STAMP__,'ERROR - Must specify .h5 files!')
  END IF

  ! Read in the old state
  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' READING STATE FILE ',iArg-1,' of ',nArgs-1,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(Args(iArg)),'" )'
  SWRITE(UNIT_stdOut,'(132("="))')
  CALL ReadOldStateFile(Args(iArg))


  SWRITE(UNIT_stdOut,'(A)') ' EVALUATING SOLUTION ON NEW MESH ...'
  CALL InterpolateSolution()

  SWRITE(UNIT_stdOut,'(A)') ' WRITING NEW SOLUTION ...'
  CALL WriteNewStateFile()
END DO

CALL FinalizeSwapMesh()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' SWAPMESH TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM swapMesh
