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
!> This tool will take a time-sequence of equi-spaced state files
!> and perform global dynamic mode decomposition
!===================================================================================================================================
PROGRAM DMD
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DMD_Vars
USE MOD_DMD,                     ONLY: DefineParametersDMD,InitDMD,performDMD,WriteDmdStateFile,FinalizeDMD
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Output,                  ONLY: DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_EOS                     ,ONLY:DefineParametersEOS
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: limit,nTotalNew
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
CALL DefineParametersEOS()
CALL DefineParametersMesh()
CALL DefineParametersDMD()

! Parse parameters
! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.2).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: dmd prm-file statefile [statefiles]')
END IF
CALL prms%read_options(Args(1))

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(A)') &
"      ___           ___           ___        "
SWRITE(UNIT_stdOut,'(A)') &
"     /\  \         /\__\         /\  \       "
SWRITE(UNIT_stdOut,'(A)') &
"    /::\  \       /::|  |       /::\  \      "
SWRITE(UNIT_stdOut,'(A)') &
"   /:/\:\  \     /:|:|  |      /:/\:\  \     "
SWRITE(UNIT_stdOut,'(A)') &
"  /:/  \:\__\   /:/|:|__|__   /:/  \:\__\    "
SWRITE(UNIT_stdOut,'(A)') &
" /:/__/ \:|__| /:/ |::::\__\ /:/__/ \:|__|   "
SWRITE(UNIT_stdOut,'(A)') &
" \:\  \ /:/  / \/__/~~/:/  / \:\  \ /:/  /   "
SWRITE(UNIT_stdOut,'(A)') &
"  \:\  /:/  /        /:/  /   \:\  /:/  /    "
SWRITE(UNIT_stdOut,'(A)') &
"   \:\/:/  /        /:/  /     \:\/:/  /     "
SWRITE(UNIT_stdOut,'(A)') &
"    \::/__/        /:/  /       \::/__/      "
SWRITE(UNIT_stdOut,'(A)') &
"     ~~            \/__/         ~~          "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')

! Initialization
CALL InitIOHDF5()
CALL InitInterpolation()

CALL InitDMD()
! Init interpolation on new polynomial degree, will set PP_N to NNew
#if USE_MPI
CALL InitMPIvars()
#endif

#if USE_MPI
nTotalNew=REAL(nVar_State*(N_State+1)**3*nElems_State)
!limit=(2**31-1)/8.
limit=(2**28-1)/8. ! max. 32 bit integer / 8
IF((nTotalNew.GT.limit))THEN
  WRITE(UNIT_stdOut,'(A,F13.0,A)')' Resulting file size is too big! Total array size may not exceed', limit, ' entries!'
  WRITE(UNIT_stdOut,'(A)')' Compile dmd without MPI'
  STOP
END IF
#endif

!===================================================================================================================================
CALL performDMD()
CALL WriteDmdStateFile()
!===================================================================================================================================

CALL FinalizeDMD()
CALL FinalizeMesh()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL Abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' DMD TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM DMD
