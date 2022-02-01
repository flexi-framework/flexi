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

MODULE MOD_HIT_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE ComputeSpectrum
  MODULE PROCEDURE ComputeSpectrum
END INTERFACE

INTERFACE WriteSpectrum
  MODULE PROCEDURE WriteSpectrum
END INTERFACE

INTERFACE WriteTurbulenceData
  MODULE PROCEDURE WriteTurbulenceData
END INTERFACE

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE AnalyzeTGV
  MODULE PROCEDURE AnalyzeTGV
END INTERFACE

INTERFACE ReadOldStateFile
  MODULE PROCEDURE ReadOldStateFile
END INTERFACE

PUBLIC:: AnalyzeTGV, ReadOldStateFile
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
! Main Routine for Analysis of turbulence
!===================================================================================================================================
SUBROUTINE AnalyzeTGV(U_In)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HIT_Analyze_Vars,   ONLY: N_Filter,nElems_HDF5,N_HDF5,NodeType_HDF5,nVar_HDF5
USE MOD_HIT_FFT,            ONLY: ComputeFFT_R2C,ComputeFFT_C2R,Interpolate_DG2FFT
USE MOD_HIT_FFT_Vars,       ONLY: N_FFT,Endw,kmax,Nc,localk
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(nVar_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5,1:nElems_HDF5) !< elementwise DG solution from state file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: E_K(0:kmax)
REAL    :: U_Global(1:PP_nVar,1:N_FFT  ,1:N_FFT  ,1:N_FFT  )  ! Real global DG solution
COMPLEX :: U_FFT(   1:PP_nVar,1:Endw(1),1:Endw(2),1:Endw(3))  ! Complex FFT solution
!===================================================================================================================================
! 1. Convert Cons to Prims
DO iElem=1,nElems_HDF5
  DO i=0,N_HDF5;  DO j=0,N_HDF5;  DO k=0,N_HDF5
    U_in(2:4,i,j,k,iElem) = U_in(2:4,i,j,k,iElem)/U_in(1,i,j,k,iElem)
  END DO;END DO;END DO
END DO

! 2. Interpolate DG solution to equidistant points
CALL Interpolate_DG2FFT(NodeType_HDF5,nVar_HDF5,U_in,U_Global)

! 3. Apply complex Fourier-Transform on solution from state file
CALL ComputeFFT_R2C(nVar_HDF5,U_Global,U_FFT)

! 4. Fourier cutoff filter (usefull for DNS to LES resolution for example)
IF (N_Filter.GT.-1) THEN
  DO k=1,Endw(3); DO j=1,Endw(2); DO i=1,Endw(1)
    IF(localk(4,i,j,k).GT.N_Filter) U_FFT(:,i,j,k) = 0.
  END DO; END DO; END DO
ELSE ! Nyquist filter
  DO k=1,Endw(3); DO j=1,Endw(2); DO i=1,Endw(1)
    IF(localk(4,i,j,k).GT.Nc) U_FFT(:,i,j,k) = 0.
  END DO; END DO; END DO
END IF

! 5. Transform filtered solution back into physical space
CALL ComputeFFT_C2R(nVar_HDF5,U_FFT,U_Global)

! 6. Compute kinetic energy per wavelength and write to file
CALL ComputeSpectrum(3,U_FFT(2:4,:,:,:),E_k)
CALL WriteSpectrum(E_k)

! 7. Compute and write additional turbulence output to file
CALL WriteTurbulenceData(E_k,U_Global)

END SUBROUTINE AnalyzeTGV


!===================================================================================================================================
!> Compute energy per wavelength
!===================================================================================================================================
PPURE SUBROUTINE ComputeSpectrum(nVar_In,U_In,E_k)
! MODULES
USE MOD_HIT_Analyze_Vars,  ONLY: N_Filter
USE MOD_HIT_FFT_Vars,      ONLY: endw,localk,Nc,kmax
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: nVar_In
COMPLEX,INTENT(IN)  :: U_In(nVar_In,1:endw(1),1:endw(2),1:endw(3))
REAL,INTENT(OUT)    :: E_k(0:kmax)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k
INTEGER :: N_max
!===================================================================================================================================
N_max = MERGE(N_Filter,Nc,(N_Filter.GT.-1))

E_k = 0.
DO k=1,endw(3); DO j=1,endw(2); DO i=1,endw(1)
  IF (localk(4,i,j,k).GT.N_max) CYCLE
  E_k(localk(4,i,j,k)) = E_k(localk(4,i,j,k)) + REAL(SUM(U_In(:,i,j,k)*CONJG(U_In(:,i,j,k))))
END DO; END DO; END DO

END SUBROUTINE ComputeSpectrum


!===================================================================================================================================
!> Writes energy spectrum over wave lengths to file
!===================================================================================================================================
SUBROUTINE WriteSpectrum(E_k)
! MODULES
USE MOD_Globals
USE MOD_HIT_Analyze_Vars,  ONLY: N_Filter,ProjectName_HDF5,Time_HDF5
USE MOD_HIT_FFT_Vars,      ONLY: Nc,kmax
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: E_k(0:kmax)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: connected
INTEGER            :: k,N_max
INTEGER            :: FileUnit
CHARACTER(LEN=255) :: Filename
!===================================================================================================================================
FileUnit = 155
INQUIRE(UNIT=FileUnit, OPENED=connected)
IF (Connected) THEN
  DO
    FileUnit=FileUnit+1
    INQUIRE(UNIT=FileUnit, OPENED=connected)
    IF (.NOT. connected) EXIT
  END DO
END IF

! Specify file name
FileName=TIMESTAMP(TRIM(ProjectName_HDF5)//'_EnergySpectrum',Time_HDF5)
FileName=TRIM(FileName)//'.dat'

! Open file and write header
OPEN(FileUnit,FILE=FileName,STATUS="REPLACE")
WRITE(FileUnit,'(a)') 'TITLE     = "Energy Spectrum "'
WRITE(FileUnit,'(a)') 'VARIABLES = "Wavenumber k" "E(k)"'

! Get maximum wave number and write to file
N_max = MERGE(N_Filter,Nc,(N_Filter.GT.-1))
DO k=1,N_max
  WRITE(FileUnit,'(I5.5,1(E20.12,1X))') k,E_k(k)
END DO

CLOSE(FileUnit)

END SUBROUTINE WriteSpectrum


!===================================================================================================================================
!> Computes and writes turbulence statistics like dissipation rate, the Kolmogorov length and the Reynolds number
!===================================================================================================================================
SUBROUTINE WriteTurbulenceData(E_k,U_In)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HIT_Analyze_Vars,  ONLY: ProjectName_HDF5,Time_HDF5,nVar_HDF5
USE MOD_HIT_Analyze_Vars,  ONLY: N_Filter,mu0
USE MOD_HIT_FFT_Vars,      ONLY: N_FFT,NCalc,Nc,kmax,Endw
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: E_k(0:kmax)                               !< Energy spectra over the wavenumbers
REAL,INTENT(IN)  :: U_in(1:nVar_HDF5,1:N_FFT,1:N_FFT,1:N_FFT) !< Global DG solution on visu nodes
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Lambda,  L_int,  Eta
REAL               :: Lambda_K,L_Int_K,Eta_K
REAL               :: IntE_k, IntEps, IntInt
REAL               :: U_rms,Re_lambda,eps,Ekin
INTEGER            :: i,j,k,N_max
INTEGER            :: FileUnit
LOGICAL            :: connected
CHARACTER(LEN=255) :: FileName
!===================================================================================================================================
! Compute integral kinetic energy in physical space
! Due to conjugate symmetry in x (here i) direction only half of the modes are stored
! the others have the property u(k_x) = u*(-k_x), meaning that the sum we compute here
! is only half the sum of all modes!
Ekin=0.
DO k=1,Endw(3); DO j=1,Endw(2); DO i=1,Endw(1)
  Ekin=Ekin+SUM(U_In(2:4,i,j,k)*U_In(2:4,i,j,k))
END DO; END DO; END DO
Ekin=0.5*Ekin/REAL(NCalc**3)

! Compute some characteristic turbulence quantities
IntE_k = 0.
IntEps = 0.
IntInt = 0.

N_max = MERGE(N_Filter,Nc,(N_Filter.GT.-1))
DO k=1,N_max
  ! Get integral energy
  IntE_k = IntE_k + E_k(k)
  ! Get integrand for mean dissipation rate
  IntEps = IntEps + 0.5*(E_k(k)*k**2+E_k(k+1)*(k+1)**2)         ! trapezoidal rule
  ! Get integrand for integral scale
  IntInt = IntInt + 0.5*( E_k(k)/REAL(k) + E_k(k+1)/REAL(k+1) ) ! trapezoidal rule
END DO

! Compute turbulence statistics
U_rms     = (2./3.*IntE_k)**0.5             ! U_RMS
Lambda    = (5.*IntE_k/IntEps)**0.5         ! TaylorMicroScale
Lambda_K  = 2.*PP_PI/Lambda                 ! TaylorMicroScale*K
L_int     = 3.*PP_PI/4.*IntInt/IntE_k       ! Int_Length
L_int_K   = 2*PP_PI/L_int                   ! Int_Length*K
#if PARABOLIC
Eps       = IntEps*2.*mu0                   ! Dissipation
Eta       = SQRT(mu0)/(IntEps*2.)**(1./4.)  ! KolmogorovLength
Eta_K     = 2.*PP_PI/Eta                    ! KolmogorovLength*K
Re_lambda = U_rms*lambda/mu0                ! Re_lambda
#else
Eps       = 0.       ! Dissipation
Eta       = 0.       ! KolmogorovLength
Eta_K     = HUGE(1.) ! KolmogorovLength*K
Re_lambda = HUGE(1.) ! Re_lambda
#endif

! Find free file unit for write process
FileUnit=55
INQUIRE(UNIT=FileUnit, OPENED=connected)
IF (Connected) THEN
  DO
    FileUnit=FileUnit+1
    INQUIRE(UNIT=FileUnit, OPENED=connected)
    IF (.NOT. connected) EXIT
  END DO
END IF

! Specify name of file
FileName=TRIM(ProjectName_HDF5)//'_Dissrate.dat'

! Check if file exists to either...
IF (FileExists(FileName)) THEN
  ! ... open either in append mode
  OPEN(FileUnit,FILE=Filename,STATUS="UNKNOWN",POSITION="APPEND")
ELSE
  ! ... or create new file with header
  OPEN( FileUnit,FILE=Filename,STATUS="REPLACE")
  WRITE(FileUnit,'(A)') 'TITLE     = "Dissipation Rate and other turbulent data "'
  WRITE(FileUnit,'(A)') 'VARIABLES = "Time" "Ekin" "EkinWave" "Dissipation" "KolmogorovLength" "KolmogorovLength*K"&
    & "TaylorMicroScale" "TaylorMicroScale*K" "Int_Length" "Int_Length*K" "U_RMS" "Re_lambda"'
END IF

! Write and close file
WRITE(FileUnit,'(12(E20.12,1X))') Time_hdf5,Ekin,IntE_k,Eps,Eta,Eta_k,lambda,lambda_k,L_int,L_int_k,U_rms,Re_lambda
CLOSE(FILEUnit)

! Also dump data to stdout
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A50)')'Turbulence Statistics'
WRITE(UNIT_stdOut,'(A14,E20.10)')           '  Epsilon   = ',Eps
WRITE(UNIT_stdOut,'(A14,E20.10,A13,E20.10)')'  E_kin     = ',Ekin,  '  E_kin_K  = ',IntE_k
WRITE(UNIT_stdOut,'(A14,E20.10,A13,E20.10)')'  Eta       = ',Eta,   '  Eta_K    = ',Eta_K
WRITE(UNIT_stdOut,'(A14,E20.10,A13,E20.10)')'  Lambda    = ',Lambda,'  Lambda_K = ',Lambda_K
WRITE(UNIT_stdOut,'(A14,E20.10,A13,E20.10)')'  L_Int     = ',L_int, '  L_Int_K  = ',L_int_K
WRITE(UNIT_stdOut,'(A14,E20.10)')           '  U_RMS     = ',U_RMS
WRITE(UNIT_stdOut,'(A14,E20.10)')           '  Re_Lambda = ',Re_Lambda
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE WriteTurbulenceData


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
USE MOD_HIT_Analyze_Vars, ONLY: nVar_HDF5,N_HDF5,nElems_HDF5
USE MOD_HIT_Analyze_Vars, ONLY: Time_HDF5,NodeType_HDF5,ProjectName_HDF5
USE MOD_ReadInTools,      ONLY: ExtractParameterFile
USE MOD_Output_Vars,      ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Output,           ONLY: insert_userblock
USE MOD_Mesh_Vars,        ONLY: MeshFile
USE MOD_StringTools,      ONLY: STRICMP,GetFileExtension
USE ISO_C_BINDING,        ONLY: C_NULL_CHAR
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)      :: StateFile !< State file to be read
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: userblockFound
CHARACTER(LEN=255)               :: prmfile=".parameter.ini"
!===================================================================================================================================
SWRITE(UNIT_stdOut,('(3A)')) "READING SOLUTION FROM STATE FILE """,TRIM(StateFile), """"

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

SWRITE(UNIT_stdOut,'(A)') "READING SOLUTION DONE!"
END SUBROUTINE ReadOldStateFile

END MODULE MOD_HIT_Analyze
