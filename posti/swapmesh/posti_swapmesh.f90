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
!> Tool to non-conservatively interpolate the solution from one mesh to another
!===================================================================================================================================
PROGRAM swapMesh
! MODULES
USE MOD_Globals
USE MOD_SwapMesh_Vars
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_SwapMesh,                ONLY: InitSwapmesh,ReadOldStateFile,WriteNewStateFile,FinalizeSwapMesh
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,FinalizeInterpolation
USE MOD_InterpolateSolution,     ONLY: InterpolateSolution
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg
!===================================================================================================================================
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()
! check if parameter file is given
IF ((nArgs.LT.2).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: swapmesh prm-file statefile [statefiles]')
END IF

! Define parameters needed
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()

CALL prms%SetSection("swapMesh")
CALL prms%CreateStringOption(   "MeshFileOld"        , "Old mesh file (different than the one found in the state file)")
CALL prms%CreateStringOption(   "MeshFileNew"        , "New mesh file")
CALL prms%CreateLogicalOption(  "useCurvedsOld"      , "Controls usage of high-order information in old mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption(  "useCurvedsNew"      , "Controls usage of high-order information in new mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateIntOption(      "NInter"             , "Polynomial degree used for interpolation (should be equal or higher than "//&
                                                       "N)",'-1')
CALL prms%CreateIntOption(      "NNew"               , "Polynomial degree used for solution",'-1')
CALL prms%CreateIntOption(      "NSuper"             , "Polynomial degree used for supersampling",'-1')
CALL prms%CreateRealOption(     "maxTolerance"       , "Tolerance used for coarse point search",'5.e-2')
CALL prms%CreateLogicalOption(  "printTroublemakers" , "Turn output of not-found points on or off",'.TRUE.')
CALL prms%CreateRealArrayOption("RefState"           , "If a RefState is defined, this state will be used at points that are "// &
                                                       "not found - without a RefState, the program will abort")
CALL prms%CreateRealOption(     "abortTolerance"     , "Tolerance used to decide if the program should abort if no "// &
                                                       "RefState is given")

! Parse parameters
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

!#ifdef MPI
!nTotal=REAL(nVar_HDF5*(NNew+1)**3*nElemsNew)
!!limit=(2**31-1)/8.
!limit=2**28-1/8. ! max. 32 bit integer / 8
!IF(nTotal.GT.limit)THEN
  !WRITE(UNIT_StdOut,'(A,F13.0,A)')' New state file size is too big! Total array size may not exceed', limit, ' entries!'
  !WRITE(UNIT_StdOut,'(A)')' Lower number of elements or NNew! Alternative: compile swapmesh without MPI'
  !STOP
!END IF
!#endif

! Evaluate solution at new GP
DO iArg=2,nArgs
  IF (.NOT.(STRICMP(GetFileExtension(Args(iArg)),'h5'))) THEN
    CALL CollectiveStop(__STAMP__,'ERROR - Must specify .h5 files!')
  END IF

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' READING STATE FILE ',iArg-1,' of ',nArgs-1,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(Args(iArg)),'" )'
  SWRITE(UNIT_stdOut,'(132("="))')

  ! Read in the old state
  CALL ReadOldStateFile(Args(iArg))
  SWRITE(UNIT_stdOut,'(A)') ' EVALUATING SOLUTION ON NEW MESH ...'

  ! Interpolate solution to new mesh
  CALL InterpolateSolution()

  ! Write new solution
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
