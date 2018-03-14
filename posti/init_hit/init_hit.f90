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
PROGRAM init_hit
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Init_Hit_Vars
USE MOD_DG_Vars,                 ONLY: U              
USE MOD_Interpolation_Vars      ,ONLY: NodeType
USE MOD_FFT,                     ONLY: InitFFT,Rogallo,FinalizeFFT                     
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars,               ONLY: nElems,Elem_xGP
USE MOD_Output,                  ONLY: DefineParametersOutput,InitOutput
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_HDF5_Output,             ONLY: WriteState
USE MOD_Commandline_Arguments
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_ReadInTools
USE FFTW3
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,j,k,iw,jw,kw,iElem,q
INTEGER                            :: nVar=5
REAL                               :: wavenumber,Time
COMPLEX                            :: basis
!===================================================================================================================================

CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

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
CALL prms%CreateStringOption(  "MeshFile"        , "Desired mesh file for initial hit solution.")
CALL prms%CreateIntOption(     "N_FFT"           , "Polynomial degree to perform DFFT on")
CALL prms%CreateIntOption(     "InitSpec"        , "Initial energy spectrum (1) Rogallo, (2) blaisdell, (3) Chasnov,&
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
N_FFT         = GETINT('N_FFT')
InitSpec      = GETINT('InitSpec')
MeshFile      = GETSTR('MeshFile')

CALL InitIOHDF5()
CALL InitInterpolation()
CALL InitOutput()
#if USE_MPI
CALL InitMPIvars()
#endif

CALL InitMesh(meshMode=2,MeshFile_IN=MeshFile)

CALL InitFFT()  
CALL Rogallo()

IF (.NOT.(ALLOCATED(Uloc_c))) ALLOCATE(Uloc_c(1:nVar,0:N,0:N,0:N))
IF (.NOT.(ALLOCATED(U_k))) ALLOCATE(U_k(1:nVar,1:endw(1),1:endw(2),0:N))
IF (.NOT.(ALLOCATED(U_j))) ALLOCATE(U_j(1:nVar,1:endw(1),0:N,0:N))
IF (.NOT.(ALLOCATED(U))) ALLOCATE(U(1:nVar,0:N,0:N,0:N,nElems))

U = 0.
U_FFT = 0.
Uloc_c = 0.
DO q=1,PP_nVar! fft of old solution to fourier modes
  CALL DFFTW_PLAN_DFT_R2C_3D(plan,N_FFT,N_FFT,N_FFT,Uloc(q,:,:,:),U_FFT(q,:,:,:),FFTW_ESTIMATE)
  CALL DFFTW_Execute(plan,Uloc(q,:,:,:),U_FFT(q,:,:,:))
END DO
U_FFT=U_FFT/(N_FFT**3)

DO iElem=1,nElems
  U_k=0.
  DO k=0,N
    DO kw=1,endw(3)
      wavenumber=kw-1
      IF (kw.GE.Nc+2) wavenumber=-(2*Nc+1-kw)
      basis=EXP(II*(wavenumber)*Elem_xGP(3,0,0,k,iElem))
      DO iw=1,endw(1)
        DO jw=1,endw(2)
          U_k(:,iw,jw,k) = U_k(:,iw,jw,k)+U_FFT(:,iw,jw,kw)*basis
        END DO !k
      END DO !kw
    END DO !jw
  END DO !iw
  U_j=0.
  DO j=0,N
    DO jw=1,endw(2)
      wavenumber=jw-1
      IF (jw.GE.Nc+2) wavenumber=-(2*Nc+1-jw)
      basis=EXP(II*(wavenumber)*Elem_xGP(2,0,j,0,iElem))
      DO k=0,N
        DO iw=1,endw(1)
          U_j(:,iw,j,k) = U_j(:,iw,j,k)+U_k(:,iw,jw,k)*basis
        END DO !k
      END DO !j
    END DO !jw
  END DO !iw
  Uloc_c=0.
  U_j(:,1,:,:)=U_j(:,1,:,:)/2.
  DO i=0,N
    DO iw=1,endw(1)
      basis=EXP(II*(iw-1)*Elem_xGP(1,i,0,0,iElem))
      DO k=0,N
        DO j=0,N
          Uloc_c(:,i,j,k) = Uloc_c(:,i,j,k)+U_j(:,iw,j,k)*basis
        END DO !k
      END DO !j
    END DO !i
  END DO !iw

  U(:,:,:,:,iElem) = 2*REAL(Uloc_c)
END DO !iElem


! Write State-File to initialize HIT
Time=0.
CALL WriteState(TRIM(MeshFile),Time,Time,.FALSE.)


CALL FinalizeMesh()
CALL FinalizeFFT()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
    CALL abort(__STAMP__,'MPI finalize error',iError,999.)
#endif

WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' initHIT FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM init_hit





