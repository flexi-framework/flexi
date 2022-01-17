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
!> Tool used to pre-calculate the distance to the closest solid wall for each volume gauss point when needed e.g. for the 
!> RANS SA turbulence model.
!> General process is as follows:
!>   * Read in of mesh (global, only single execution)
!>   * For each volume gauss points: 
!>       * Coarse search. Super-sampling of all solid wall surfaces to find the closest surface.
!>       * Minimization of the square of the distance function using a simple conjugate gradient algorithm on that specific surface.
!>         The coarse search is used as a starting point.
!===================================================================================================================================
PROGRAM walldistance
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_walldistance_Vars
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_Walldistance,            ONLY: InitWalldistance,FinalizeWalldistance,CalcWalldistance
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,FinalizeInterpolation
#if FV_ENABLED
USE MOD_FV,                ONLY:DefineParametersFV,InitFV,FinalizeFV
USE MOD_FV_Basis,          ONLY:InitFV_Basis,FinalizeFV_Basis
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

! Define parameters needed
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersMesh()

CALL prms%SetSection("walldistance")
CALL prms%CreateIntOption("NSuper" , "Polynomial degree used for supersampling on the surface in the coarse search.")
CALL prms%CreateIntOption("NVisu" ,  "Polynomial degree used for visualization.")
CALL prms%CreateLogicalOption("DebugVisu" ,  "Visualize the walldistance in .vtu format.")
CALL prms%CreateLogicalOption("includeTrip" ,"Switch on to project a trip point to the surface.",'.FALSE.')
CALL prms%CreateRealArrayOption("TripX"     ,"2D coordinates of trip point.")

! Parse parameters
! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.NE.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: walldistance prm-file')
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
CALL InitInterpolation()
#if FV_ENABLED
CALL InitFV_Basis()
#endif
#if USE_MPI
CALL InitMPIvars()
#endif
CALL InitMesh(meshMode=2)
CALL InitWalldistance()

CALL CalcWalldistance()

CALL FinalizeWalldistance()
CALL FinalizeMesh()
#if FV_ENABLED
CALL FinalizeFV_Basis()
#endif
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL Abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' WALLDISTANCE TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM walldistance
