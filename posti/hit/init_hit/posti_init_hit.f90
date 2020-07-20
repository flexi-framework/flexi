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
!> This tool generates a state file which serves as an initial condition to compute the HIT test case.
!> It generates a random velocity field based on a specified turbulent kinetic energy distribution.
!> The phase of each velocity mode is chosen random.
!===================================================================================================================================
PROGRAM posti_init_hit
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Init_HIT_Vars
USE MOD_Init_HIT
USE MOD_ReadInTools
USE MOD_Commandline_Arguments
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Mesh_Vars,               ONLY: nElems
USE MOD_FFT,                     ONLY: InitFFT,FinalizeFFT,ComputeFFT_R2C,ComputeFFT_C2R,EvalFourierAtDGCoords
USE MOD_FFT_Vars,                ONLY: N_FFT,Endw,II
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Output,                  ONLY: DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
USE MOD_HDF5_Output,             ONLY: WriteState
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_MPI,                     ONLY: InitMPI,DefineParametersMPI
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE FFTW3
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
  'This tool does currently not work with more than 1 Proc!')

CALL ParseCommandlineArguments()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
    " ||==========================================||                                                           "
SWRITE(UNIT_stdOut,'(A)') &
    " || Generate Spec Data for FLEXI Hit Restart ||                                                           "
SWRITE(UNIT_stdOut,'(A)') &
    " ||==========================================||                                                           "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')


! Define Parameters
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
!================================
CALL prms%SetSection("initHIT")
CALL prms%CreateStringOption( "MeshFile"   , "Desired mesh file for initial hit solution.")
CALL prms%CreateIntOption(    "N_FFT"      , "Polynomial degree to perform DFFT on")
CALL prms%CreateIntOption(    "InitSpec"   , "Initial energy spectrum (1) Rogallo, (2) blaisdell, (3) Chasnov,&
                                             & (4) Inf interial range .")

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: init_hit [prm-file]')
END IF
! Parse parameters
CALL prms%read_options(Args(1))
ParameterFile = Args(1)

! Readin Parameters
N_FFT    = GETINT('N_FFT')
InitSpec = GETINT('InitSpec')
MeshFile = GETSTR('MeshFile')

CALL InitIOHDF5()
CALL InitInterpolation()
CALL InitOutput()
#if USE_MPI
CALL InitMPIvars()
#endif

CALL InitMesh(meshMode=2,MeshFile_IN=MeshFile)

CALL InitFFT()
CALL Init_InitHit()

ALLOCATE(U_FFT(1:PP_nVar,1:Endw(1),1:Endw(2),1:Endw(3)))

CALL Rogallo()

ALLOCATE(U(     1:PP_nVar,0:N,0:N,0:N,nElems))
ALLOCATE(Uloc_c(1:PP_nVar,0:N,0:N,0:N))
ALLOCATE(U_k(   1:PP_nVar,1:endw(1),1:endw(2),0:N))
ALLOCATE(U_j(   1:PP_nVar,1:endw(1),0:N,0:N))

U      = 0.
U_FFT  = 0.
Uloc_c = 0.

CALL ComputeFFT_R2C(5,Uloc,U_FFT)
U_FFT=U_FFT/(N_FFT**3)

Call EvalFourierAtDGCoords(5,U_FFT,U)

! Write State-File to initialize HIT
CALL WriteState(TRIM(MeshFile),0.,0.,.FALSE.)

CALL Finalize_InitHIT()
CALL FinalizeParameters()
CALL FinalizeInterpolation()
CALL FinalizeOutput()
CALL FinalizeMesh()
CALL FinalizeFFT()
#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
CALL FinalizeMPI()
#endif

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' initHIT FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM posti_init_hit
