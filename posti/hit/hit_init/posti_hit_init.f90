!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
!> This tool generates a state file which serves as an initial condition to compute the HIT test case. It generates a random
!> velocity field based on a specified turbulent kinetic energy distribution. The phase of each velocity mode is chosen random.
!===================================================================================================================================
PROGRAM HIT_Init
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_HIT_Init
USE MOD_HIT_Init_Vars
USE MOD_ReadInTools
USE MOD_Commandline_Arguments
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_HIT_FFT,                 ONLY: InitFFT,FinalizeFFT,EvalFourierAtDGCoords
USE MOD_HIT_FFT_Vars,            ONLY: N_FFT,Endw
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars,               ONLY: nElems
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5,FinalizeIOHDF5
USE MOD_HDF5_Output,             ONLY: WriteState
USE MOD_Output,                  ONLY: DefineParametersOutput,InitOutput,FinalizeOutput
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
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(20X,A)') &
    "     ____  ____ _____ _________                      _____ ____  _____ _____ _________  "
SWRITE(UNIT_stdOut,'(20X,A)') &
    "    |_   ||   _|_   _|  _   _  |                    |_   _|_   \|_   _|_   _|  _   _  | "
SWRITE(UNIT_stdOut,'(20X,A)') &
    "      | |__| |   | | |_/ | | \_|       ______         | |   |   \ | |   | | |_/ | | \_| "
SWRITE(UNIT_stdOut,'(20X,A)') &
    "      |  __  |   | |     | |          |______|        | |   | |\ \| |   | |     | |     "
SWRITE(UNIT_stdOut,'(20X,A)') &
    "     _| |  | |_ _| |_   _| |_                        _| |_ _| |_\   |_ _| |_   _| |_    "
SWRITE(UNIT_stdOut,'(20X,A)') &
    "    |____||____|_____| |_____|                      |_____|_____|\____|_____| |_____|   "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')

! Define Parameters
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersOutput()
CALL DefineParametersMesh()

! Define Parameters HIT_Init
CALL prms%SetSection("HIT_Init")
CALL prms%CreateIntOption(    "Seed"       , "Seed for random number generator for Rogallo precedure (Only Debug)")
CALL prms%CreateIntOption(    "N_FFT"      , "Number of global interpolation points to perform DFFT on.")
CALL prms%CreateIntOption(    "InitSpec"   , "Initial energy spectrum (1) Rogallo,&
                                                                     &(2) Blaisdell,&
                                                                     &(3) Chasnov,&
                                                                     &(4) Inf intertial range,&
                                                                     &(5) Karman-Pao.")

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti_hit_init parameter.ini')
END IF
! Parse parameter file
CALL prms%read_options(Args(1))
ParameterFile = Args(1)

! Readin Parameters
N_FFT    = GETINT('N_FFT')
InitSpec = GETINT('InitSpec')
Seed     = GETINT('Seed','0')
MeshFile = GETSTR('MeshFile')

CALL InitIOHDF5()
CALL InitInterpolation()
CALL InitOutput()
#if USE_MPI
CALL InitMPIvars()
#endif
CALL InitMesh(meshMode=0,MeshFile_IN=MeshFile)
CALL InitFFT()

! Allocate the solution arrays in Fourier and DG space
ALLOCATE(U_FFT(1:PP_nVar,1:Endw(1),1:Endw(2),1:Endw(3)))
ALLOCATE(U(1:PP_nVar,0:N,0:N,0:N,nElems))

! Build random flow realization with desired energy distribution
CALL Rogallo(U_FFT)

! Evaluate solution in Fourier space at Gauss points
CALL EvalFourierAtDGCoords(5,U_FFT,U,doPrintTime=.TRUE.)

! Write State-File to initialize HIT
CALL WriteState(TRIM(MeshFile),0.,0.,.FALSE.)

! Deallocate solution arrays
DEALLOCATE(U_FFT)
DEALLOCATE(U)

CALL FinalizeParameters()
CALL FinalizeInterpolation()
CALL FinalizeOutput()
CALL FinalizeMesh()
CALL FinalizeFFT()
CALL FinalizeIOHDF5
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' HIT_INIT FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM HIT_Init
