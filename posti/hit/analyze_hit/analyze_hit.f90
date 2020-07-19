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

MODULE MOD_Analyze_Hit
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitAnalyze
  MODULE PROCEDURE InitAnalyze
END INTERFACE

INTERFACE FinalizeAnalyze
  MODULE PROCEDURE FinalizeAnalyze
END INTERFACE

INTERFACE AnalyzeTGV
  MODULE PROCEDURE AnalyzeTGV
END INTERFACE

INTERFACE ReadOldStateFile
  MODULE PROCEDURE ReadOldStateFile
END INTERFACE

PUBLIC::InitAnalyze, FinalizeAnalyze, AnalyzeTGV, ReadOldStateFile
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
! Get some parameters needed by equation modules and initialize variables for turbulence statistics
!===================================================================================================================================
SUBROUTINE InitAnalyze()
! MODULES
USE MOD_Globals
USE MOD_ANALYZE_HIT_Vars
USE MOD_FFT_Vars,                ONLY: N_FFT
USE MOD_Readintools
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: connected
!===================================================================================================================================
IF(AnalyzeInitIsDone)THEN
  SWRITE(*,*) "InitAnalyze not ready to be called or already called."
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ANALYZE...'

! Find free file unit for analysis data
FileUnit_HIT=55
INQUIRE(UNIT=FileUnit_HIT, OPENED=connected)
IF (Connected) THEN
  DO
    FileUnit_HIT=Fileunit_HIT+1
    INQUIRE(UNIT=FileUnit_HIT, OPENED=connected)
    IF (.NOT. connected) EXIT
  END DO
END IF

FileName_HIT=TRIM(ProjectName_HDF5)//'_Dissrate.dat'

! delete old file and open file in append-mode
OPEN(FileUnit_HIT,FILE=Filename_HIT,STATUS="REPLACE")
WRITE(FileUnit_HIT,*)'TITLE     = "Dissipation Rate and other turbulent data "'
WRITE(FileUnit_HIT,'(a)')'VARIABLES = "Time" "Ekin" "EkinWave" "Dissipation" "KolmogorovLength" "KolmogorovLength*K"&
  & "TaylorMicroScale" "TaylorMicroScale*K" "Int_Length" "Int_Length*K" "U_RMS" "Re_lambda"'
CLOSE(FILEUnit_HIT)

kmax=NINT(sqrt(REAL(((N_FFT*N_FFT)*3))))+1

ALLOCATE(E_k(0:kmax))
ALLOCATE(numDiss_k(0:kmax))
ALLOCATE(MuSGS_k(0:kmax))
ALLOCATE(T_k(0:kmax))
ALLOCATE(eps_k(0:kmax))
ALLOCATE(eMean_k(0:kmax))

AnalyzeInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitAnalyze

!===================================================================================================================================
! Deallocate Variables and Finalize Analyze
!===================================================================================================================================
SUBROUTINE FinalizeAnalyze()
! MODULES
USE MOD_Globals
USE MOD_ANALYZE_HIT_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(AnalyzeInitIsDone)THEN
  AnalyzeInitIsDone=.FALSE.
END IF

SDEALLOCATE(E_k)
SDEALLOCATE(numDiss_k)
SDEALLOCATE(MuSGS_k)
SDEALLOCATE(T_k)
SDEALLOCATE(eps_k)
SDEALLOCATE(eMean_k)

SWRITE(UNIT_stdOut,'(A)')' FINALIZE ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE FinalizeAnalyze

!===================================================================================================================================
! Main Routine for Analysis of turbulence
!===================================================================================================================================
SUBROUTINE AnalyzeTGV(Time,nVar_In,U_In)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE FFTW3
USE MOD_ANALYZE_HIT_Vars
USE MOD_FFT,                ONLY: ComputeFFT_R2C, ComputeFFT_C2R, Interpolate_DG2FFT
USE MOD_FFT_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: Time
INTEGER,INTENT(IN)    :: nVar_In
REAL,INTENT(INOUT)    :: U_in(1:nVar_In,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< elementwise DG solution from state file
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: U_Global(1:nVar_In,1:N_FFT  ,1:N_FFT  ,1:N_FFT  )        ! Real global DG solution
COMPLEX :: U_FFT(   1:nVar_In,1:Endw(1),1:Endw(2),1:Endw(3))        ! Complex FFT solution
REAL    :: v(1:5,1:N_FFT,1:N_FFT,1:N_FFT)
REAL    :: E_K_loc(0:kmax),dummysum,Lambda, L_int, ETA,eps,Ekinglob
REAL    :: IntEps, IntInt, Lambda_K, L_Int_K, ETA_K, urms ,Re_lambda
LOGICAL :: connected
! Transfer function
INTEGER :: intLim
! Pi and kmax
REAL    :: TwoPi
REAL    :: Pi
!===================================================================================================================================
TwoPi  = 2.*PP_Pi
Pi     = PP_Pi

Ekin=0.
Ekinglob=0.
U_FFT=0.
U_Global = 0.

! Convert Cons to Prims
DO iElem=1,nElems_HDF5
  DO i=0,N_HDF5;  DO j=0,N_HDF5;  DO k=0,N_HDF5
    U_in(2:4,i,j,k,iElem) = U_in(2:4,i,j,k,iElem)/U_in(1,i,j,k,iElem)
  END DO;END DO;END DO
END DO

! 1. Interpolate DG solution to equidistant points
CALL Interpolate_DG2FFT(NodeType_HDF5,nVar_HDF5,U_in,U_Global)

! 2. Apply complex Fourier-Transform on solution from state file
CALL ComputeFFT_R2C(nVar_HDF5,U_Global,U_FFT)

U_FFT=U_FFT/REAL(N_FFT**3)

! 3. Fourier cutoff filter (usefull for DNS to LES resolution for example)
!Due to conjugate symmetry in x (here i) direction only half of the modes are stored
!the others have the property u(k_x) = u*(-k_x), meaning that the sum we compute here
!is only half the sum of all modes!

IF (N_Filter.GT.-1) THEN
  DO i=startw(1),endw(1)
    DO j=startw(2),endw(2)
      DO k=startw(3),endw(3)
        IF (ABS(localk(4,i,j,k)).GT.N_Filter) U_FFT(:,i,j,k)=0.
      END DO
    END DO
  END DO
ELSE
  DO i=startw(1),endw(1)
    DO j=startw(2),endw(2)
      DO k=startw(3),endw(3)
        IF (ABS(localk(4,i,j,k)).GT.Nc) U_FFT(:,i,j,k)=0.
      END DO
    END DO
  END DO
END IF

! 4. Compute kinetic energy per wavelength
E_k_loc=0.
CALL Compute_Spectrum(3,U_FFT(2:4,:,:,:),E_k_loc)
v=0.

! 5. Apply inverse Fourier-Transform back into physical space
CALL ComputeFFT_C2R(nVar_HDF5,U_FFT,U_Global)

! 6. Compute other turbulence statistics
Ekin=0.
DO k=startw(3),endw(3)
  DO j=startw(2),endw(2)
    DO i=startw(1),endw(1)
      Ekin=Ekin+SUM(U_Global(2:4,i,j,k)*U_Global(2:4,i,j,k))
    END DO
  END DO
END DO
Ekin=0.5*Ekin/REAL(N_VISU**3)

E_k_loc(0) = 0.5*E_k_loc(0)
FileUnit_EK=155
INQUIRE(UNIT=FileUnit_EK, OPENED=connected)

IF (Connected) THEN
  DO
    FileUnit_EK=Fileunit_EK+1
    INQUIRE(UNIT=FileUnit_EK, OPENED=connected)
    IF (.NOT. connected) EXIT
  END DO
END IF

FileName_EK=TIMESTAMP(TRIM(ProjectName_HDF5)//'_EnergySpectrum',time)
FileName_EK=TRIM(Filename_EK)//'.dat'
OPEN(FileUnit_Ek,FILE=Filename_EK,STATUS="REPLACE")
WRITE(FileUnit_EK,*)'TITLE     = "Energy Spectrum "'

WRITE(FileUnit_EK,'(a)')'VARIABLES = "Wavenumber k" "E(k)"'

IntEps=1E-16
IntInt=0.
dummysum=0.
intLim = MIN(Nyq,N_Filter)

DO k=1,Nyq ! IntLim!-1
  dummysum=dummysum+E_k_loc(k)
  ! Get integrand for mean dissipation rate
  IntEps=IntEps+0.5*(E_k_loc(k)*k**2+E_k_loc(k+1)*(k+1)**2)
  ! Get integrand for integral scale
  IntInt=IntInt+0.5*(E_k_loc(k)/(real(k)+1E-16)+E_k_loc(k+1)/(real(k+1)+1E-16))
  WRITE(FileUnit_EK,'(I5.5,1(E20.12,X))')k,E_k_loc(k)
END DO

CLOSE(FILEUnit_EK)

EPS=IntEps*2.*mu0                         ! Dissipation
ETA=mu0**(3./4.)/(IntEps*2.*mu0)**(1./4.) ! KolmogorovLength
ETA_K=TwoPi/ETA                           ! KolmogorovLength*K
Lambda=(5.*dummysum/IntEps)**0.5          ! TaylorMicroScale
Lambda_K=2.*Pi/Lambda                     ! TaylorMicroScale*K
L_int=3.*Pi/4.*IntInt/dummysum            ! Int_Length
L_int_K=2*Pi/L_int                        ! Int_Length*K
urms=(2./3.*dummysum)**0.5                ! U_RMS
Re_lambda=urms*lambda/mu0                 ! Re_lambda
OPEN(FileUnit_HIT,FILE=Filename_HIT,STATUS="UNKNOWN",POSITION="APPEND")
WRITE(FileUnit_HIT,'(12(E20.12,X))')Time,Ekin,dummysum,EPS,ETA,eta_k,lambda,lambda_k,L_int,L_int_k,urms,Re_lambda

CLOSE(FILEUnit_HIT)
WRITE(UNIT_StdOut,*)'-------------------------------------------------------------------------------------------'
WRITE(Unit_StdOut,'(A50)')'Some output regarding turbulence statistics'
WRITE(UNIT_StdOut,'(A14,E20.10)')'  Epsilon   = ',EPS
WRITE(UNIT_StdOut,'(A14,E20.10,A13,E20.10)')'  ETA       = ',ETA,   '  ETA_K    = ',ETA_K
WRITE(UNIT_StdOut,'(A14,E20.10,A13,E20.10)')'  Lambda    = ',Lambda,'  Lambda_K = ',Lambda_K
WRITE(UNIT_StdOut,'(A14,E20.10,A13,E20.10)')'  L_int     = ',L_int, '  L_int_K  = ',L_int_K
WRITE(UNIT_StdOut,'(A14,E20.10)')'  urms      = ',urms
WRITE(UNIT_StdOut,'(A14,E20.10)')'  Re_lambda = ',Re_lambda
WRITE(UNIT_StdOut,*)'-------------------------------------------------------------------------------------------'

END SUBROUTINE AnalyzeTGV

!===================================================================================================================================
!> Compute energy per wavelength
!===================================================================================================================================
SUBROUTINE Compute_Spectrum(nVar_In,U_In,E_k)
! MODULES
USE MOD_Globals
USE MOD_ANALYZE_HIT_Vars,  ONLY: N_Filter, kmax
USE MOD_FFT_Vars,          ONLY: startw,endw,localk,Nc
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
COMPLEX,INTENT(IN)  :: U_In(nVar_In,startw(1):endw(1),startw(2):endw(2),startw(3):endw(3))
INTEGER,INTENT(IN)  :: nVar_In
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT),DIMENSION(0:kmax) :: E_k
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER ::i,j,k,k_eff
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
E_k = 0.
IF (N_Filter.GT.-1) THEN
  DO k=startw(3),endw(3)
    DO j=startw(2),endw(2)
      DO i=startw(1),endw(1)
        IF (localk(4,i,j,k).GT.N_Filter) CYCLE
        k_eff = NINT(localk(4,i,j,k))
        E_k(k_eff) =   E_k(k_eff)+&
            SUM(U_In(:,i,j,k)*conjg(U_In(:,i,j,k)))
      END DO
    END DO
  END DO
ELSE
  DO k=startw(3),endw(3)
    DO j=startw(2),endw(2)
      DO i=startw(1),endw(1)
        IF (localk(4,i,j,k).GT.Nc) CYCLE
        k_eff = NINT(localk(4,i,j,k))
        E_k(k_eff) =   E_k(k_eff)+&
            SUM(U_In(:,i,j,k)*conjg(U_In(:,i,j,k)))
      END DO
    END DO
  END DO
END IF

END SUBROUTINE Compute_Spectrum

!===================================================================================================================================
!> Open a state file, read the old state and store the information later needed to write a new state.
!===================================================================================================================================
SUBROUTINE ReadOldStateFile(StateFile)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_DG_Vars,          ONLY: U
USE MOD_HDF5_Input,       ONLY: OpenDataFile,CloseDataFile,ISVALIDHDF5FILE
USE MOD_HDF5_Input,       ONLY: ReadArray,ReadAttribute,GetDataProps,GetDataSize,DataSetExists
USE MOD_IO_HDF5,          ONLY: File_ID
USE MOD_Analyze_Hit_Vars, ONLY: nVar_HDF5,N_HDF5,nElems_HDF5
USE MOD_Analyze_Hit_Vars, ONLY: Time_HDF5,NodeType_HDF5,ProjectName_HDF5
USE MOD_ReadInTools,      ONLY: ExtractParameterFile
USE MOD_Output_Vars,      ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Output,           ONLY: insert_userblock
USE MOD_Mesh_Vars,        ONLY: MeshFile
USE MOD_StringTools,      ONLY: STRICMP,GetFileExtension
USE ISO_C_BINDING,        ONLY: C_NULL_CHAR
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)      :: StateFile !< State file to be read
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: userblockFound
CHARACTER(LEN=255)               :: prmfile=".parameter.ini"
INTEGER                          :: iElem
!===================================================================================================================================
SWRITE(*,*) "READING SOLUTION FROM STATE FILE """,TRIM(StateFile), """"

! Get start index of file extension to check if it is a h5 file
IF (.NOT.STRICMP(GetFileExtension(StateFile), 'h5')) &
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid file extension of input state file!')

! Check if state file is a valid state
IF(.NOT.ISVALIDHDF5FILE(StateFile)) &
  CALL CollectiveStop(__STAMP__,'ERROR - Not a valid state file!')

! Open the data file
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Get data size of solution array
CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5)

! Allocate solution array in correct size
ALLOCATE(U( 1:nVar_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5,1:nElems_HDF5))

! Read the DG solution and store in UNew
CALL ReadArray('DG_Solution',5,&
               (/nVar_HDF5,N_HDF5+1,N_HDF5+1,N_HDF5+1,nElems_HDF5/),0,5,RealArray=U)

! Check for FV-Subcells
DO iElem=1,nElems_HDF5
  IF (FV_Elems(iElem).EQ.1) THEN ! FV Element
    SWRITE(*,*)'FV Subcells detected in Element', iElem
    CALL CollectiveStop(__STAMP__,'ERROR - Programm cannot handle FV Subcells!')
  END IF
END DO

! Read the attributes from file
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_HDF5)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName_HDF5)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)

! Extract parameter file from userblock (if found)
CALL ExtractParameterFile(StateFile,TRIM(prmfile),userblockFound)
! prepare userblock file
CALL insert_userblock(TRIM(UserBlockTmpFile)//C_NULL_CHAR,TRIM(prmfile)//C_NULL_CHAR)
INQUIRE(FILE=TRIM(UserBlockTmpFile),SIZE=userblock_total_len)

! Close the data file
CALL CloseDataFile()

SWRITE(*,*) "READING SOLUTION DONE!"
END SUBROUTINE ReadOldStateFile


END MODULE MOD_Analyze_Hit


