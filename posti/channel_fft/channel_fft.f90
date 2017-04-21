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
!> This tool will read in state files (single or multiple at once) from channel flow simulations and will perform 
!> a FFT (using external libraries) to generate spectra as well as mean profiles for fluctuations and the mean velocity.
!> If several state files are given, an average over all of them is calculated.
!> The mesh needs to be ijk sorted. 
!> To perform the FFT, the mesh and the solution will be interpolated to a equidistant FFT grid.
!===================================================================================================================================
PROGRAM channel_fft
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Readin,             ONLY: ReadIJKSorting
USE MOD_Mesh_Vars,               ONLY: nElems,OffsetElem
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5,File_ID
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_HDF5_Input,              ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,ReadArray
USE MOD_Interpolation_Vars,      ONLY: NodeType
USE MOD_DG_Vars,                 ONLY: U
USE MOD_FFT,                     ONLY: InitFFT,PerformFFT,FFTOutput,FinalizeFFT,PrimStateAtFFTCoords
USE MOD_FFT_Vars,                ONLY: ProjectName,Time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg
INTEGER                            :: nElems_State,nVar_State,N_State
CHARACTER(LEN=255)                 :: MeshFile_state,NodeType_State
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()
! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: channel_fft [prm-file] statefile [statefiles]')
END IF

! Define parameters needed
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersMesh()

CALL prms%SetSection("channelFFT")
CALL prms%CreateIntOption( "NCalc",  "Polynomial degree to perform DFFT on.")
CALL prms%CreateRealOption("Re_tau", "Reynolds number based on friction velocity and channel half height.")

! Parse parameters
CALL prms%read_options(Args(1))

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(A)') &
"   _____ _    _          _   _ _   _ ______ _        ______ ______ _______ "
SWRITE(UNIT_stdOut,'(A)') &
"  / ____| |  | |   /\   | \ | | \ | |  ____| |      |  ____|  ____|__   __|"
SWRITE(UNIT_stdOut,'(A)') &
" | |    | |__| |  /  \  |  \| |  \| | |__  | |      | |__  | |__     | |   "
SWRITE(UNIT_stdOut,'(A)') &
" | |    |  __  | / /\ \ | . ` | . ` |  __| | |      |  __| |  __|    | |   "
SWRITE(UNIT_stdOut,'(A)') &
" | |____| |  | |/ ____ \| |\  | |\  | |____| |____  | |    | |       | |   "
SWRITE(UNIT_stdOut,'(A)') &
"  \_____|_|  |_/_/    \_\_| \_|_| \_|______|______| |_|    |_|       |_|   "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')


! Initialization
CALL InitIOHDF5()

! Open the first statefile to read necessary attributes
CALL OpenDataFile(Args(2),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=MeshFile_state)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',        1,RealScalar=Time)
CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State)
CALL CloseDataFile()
IF (.NOT.STRICMP(NodeType_State,NodeType)) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Node type of statefile not equal to node type FLEXI is compiled with')
END IF

CALL InitInterpolation(N_State)
CALL InitMesh(meshMode=1,MeshFile_IN=MeshFile_state)
CALL ReadIJKSorting()
#if USE_MPI
CALL InitMPIvars()
#endif
CALL InitFFT()

! Loop over all statefiles
DO iArg=2,nArgs
  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iArg-1,' of ',nArgs-1,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(Args(iArg)),'" )'
  SWRITE(UNIT_stdOut,'(132("="))')

  ! Get Solution
  CALL OpenDataFile(Args(iArg),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nElems/),OffsetElem,5,RealArray=U)
  CALL CloseDataFile()

  ! Interpolate solution to FFT grid
  ! FFT grid is globally equidistant
  CALL PrimStateAtFFTCoords()

  ! Perform actual FFT
  CALL PerformFFT()
END DO

! Do output of results
CALL FFTOutput()

! Finalize 
CALL FinalizeFFT()
CALL FinalizeMesh()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' Channel_FFT TOOL FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM channel_fft

