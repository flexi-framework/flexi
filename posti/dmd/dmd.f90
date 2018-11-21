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
MODULE MOD_DMD
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitDMD
  MODULE PROCEDURE InitDMD
END INTERFACE

INTERFACE DefineParametersDMD
  MODULE PROCEDURE DefineParametersDMD
END INTERFACE
  
INTERFACE performDMD
  MODULE PROCEDURE performDMD
END INTERFACE

INTERFACE WriteDMDStateFile
  MODULE PROCEDURE WriteDMDStateFile
END INTERFACE

INTERFACE FinalizeDMD
  MODULE PROCEDURE FinalizeDMD
END INTERFACE

PUBLIC::InitDMD,DefineParametersDMD,performDMD,WriteDMDStateFile,FinalizeDMD
CONTAINS

SUBROUTINE DefineParametersDMD() 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools,             ONLY: prms
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension,INTTOSTR
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL prms%SetSection("DMD")
CALL prms%CreateStringOption( "VarNameDMD"         , "Name of variable to perform DMD.")
CALL prms%CreateRealOption(   "SvdThreshold"       , "Define relative lower bound of singular values.")
CALL prms%CreateIntOption(    "nModes"             , "Number of Modes to visualize.")

END SUBROUTINE DefineParametersDMD 

SUBROUTINE InitDMD() 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Commandline_Arguments
USE MOD_Readintools             ,ONLY: prms,GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY,CountOption
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension,INTTOSTR
USE MOD_DMD_Vars,                ONLY: K,nFiles,nVar_State,N_State,nElems_State,NodeType_State,nDoFs,MeshFile_state
USE MOD_DMD_Vars,                ONLY: dt,nModes,SvdThreshold,VarNameDMD,VarNames_State,Time_State,TimeEnd_State
USE MOD_IO_HDF5,                 ONLY: File_ID
USE MOD_Output_Vars,             ONLY: ProjectName,NOut
USE MOD_HDF5_Input,              ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,ReadArray
USE MOD_Mesh_Vars,               ONLY: nElems,OffsetElem,nGlobalElems,MeshFile
USE MOD_EquationDMD,             ONLY: InitEquationDMD,CalcEquationDMD
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iFile,iVar,offset
REAL,ALLOCATABLE                 :: Utmp(:,:,:,:,:)
REAL                             :: time
!===================================================================================================================================
IF(nArgs.LT.3) CALL CollectiveStop(__STAMP__,'At least two files required for dmd!')
nFiles=nArgs-1

VarNameDMD   = GETSTR ('VarNameDMD'  ,'Density')
SvdThreshold = GETREAL('SvdThreshold','0.0')
nModes       = GETINT ('nModes'      ,'9999')

IF (nModes .GT. (nFiles-1)) THEN
  nModes = nFiles - 1  
END IF ! nModes .GT. (nFiles-1)

! Open the first statefile to read necessary attributes
CALL OpenDataFile(Args(2),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=MeshFile_state)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',    1,RealScalar=starttime)
CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State)
ALLOCATE(VarNames_State(nVar_State))
CALL ReadAttribute(File_ID,'VarNames',nVar_State,StrArray=VarNames_State)
CALL CloseDataFile()
! Set everything for output
Time_State = starttime
MeshFile = MeshFile_state
NOut = N_State
nGlobalElems = nElems_State
nElems       = nElems_State
OffsetElem   = 0 ! OffsetElem is 0 since the tool only works on singel

nDoFs = (N_State+1)**3 * nElems_State

CALL InitEquationDMD()

ALLOCATE(K    (nDoFs,nFiles))
ALLOCATE(Utmp (nVar_State   ,0:N_State,0:N_State,0:N_State,nElems_State))

DO iFile = 2, nFiles+1
  WRITE(UNIT_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iFile,' of ',nFiles,' FILES.'
  WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(Args(iFile)),'" )'
  WRITE(UNIT_stdOut,'(132("="))')

  ! Read Statefiles
  CALL OpenDataFile(Args(iFile),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  offset=0
  CALL ReadAttribute(File_ID,'Time',    1,RealScalar=time)
  CALL ReadArray('DG_Solution',5,&
                (/nVar_State,N_State+1,N_State+1,N_State+1,nElems_State/),0,5,RealArray=Utmp)
  CALL CloseDataFile()

  IF (iFile .EQ. 3) THEN
    dt = (time - starttime) / (iFile-2)
  END IF
  IF (iFile .GT. 3) THEN
    IF ( ABS((time - starttime) / (iFile-2) - dt) .GT. 1E-10) THEN
      CALL CollectiveStop(__STAMP__,'Given file are not an equispaced time series')
    END IF
  END IF ! iFile .GT. 2

  CALL CalcEquationDMD(Utmp,K(:,iFile-1))
END DO ! iFile = 1,nFiles 
TimeEnd_State = time
DEALLOCATE(Utmp)
END SUBROUTINE InitDMD 

SUBROUTINE performDMD() 
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_DMD_Vars,                ONLY: K,nDoFs,nFiles,Phi,lambda,freq,dt,alpha,sigmaSort
USE MOD_Mathtools,               ONLY: INVERSE
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! External Subroutines from LAPACK
EXTERNAL         DGESDD
EXTERNAL         ZGEEV
EXTERNAL         ZGESV
! LOCAL VARIABLES
INTEGER                          :: i,M,N,LDA,LDU,LDVT,LWMAX,INFO,LWORK
INTEGER,ALLOCATABLE              :: IWORK(:),IPIV(:)
DOUBLE PRECISION,ALLOCATABLE     :: WORK(:)
DOUBLE PRECISION,ALLOCATABLE     :: SigmaInv(:,:),SigmaMat(:,:),USVD(:,:),WSVD(:,:),SigmaSVD(:)              
COMPLEX,ALLOCATABLE              :: WORKC(:)
REAL,ALLOCATABLE                 :: RWORK(:),alphaSort(:,:)
DOUBLE PRECISION,ALLOCATABLE     :: KTmp(:,:),RGlobMat(:,:),Rglob,Rsnap1,Rsnap2
COMPLEX,ALLOCATABLE              :: STilde(:,:),STildeWork(:,:),eigSTilde(:),VL(:,:),VR(:,:)                
COMPLEX,ALLOCATABLE              :: Vand(:,:),qTmp(:,:),q(:)
COMPLEX,ALLOCATABLE              :: P(:,:),Ptmp(:,:),PhiTmp(:,:),alphaTmp(:)
REAL                             :: pi = acos(-1.)
!===================================================================================================================================
LWMAX = nDoFs+1000
M = nDoFs
N = nFiles-1
LDA  = M
LDU  = M 
LDVT = N
LWORK = -1
ALLOCATE(WORK(LWMAX))
ALLOCATE(USVD(M,N))
ALLOCATE(WSVD(N,N))
ALLOCATE(SigmaSVD(min(N,M)))

ALLOCATE(KTmp(nDoFs,nFiles-1))
ALLOCATE(IWORK(8*N))
KTmp = K(:,1:N)
CALL DGESDD( 'S', M, N, KTmp, LDA, SigmaSVD, USVD, LDU, WSVD, LDVT,WORK, LWORK, IWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
CALL DGESDD( 'S', M, N, KTmp, LDA, SigmaSVD, USVD, LDU, WSVD, LDVT,WORK, LWORK, IWORK,INFO )
DEALLOCATE(WORK)
DEALLOCATE(KTmp)
DEALLOCATE(IWORK)

ALLOCATE(SigmaInv(N,N))
ALLOCATE(SigmaMat(N,N))
ALLOCATE(STilde(N,N))

SigmaInv=0.
SigmaMat=0.
DO i = 1, N
  SigmaInv(i,i)=1./SigmaSVD(i)
  SigmaMat(i,i)=SigmaSVD(i)
END DO 

STilde=dcmplx(MATMUL(MATMUL(MATMUL(TRANSPOSE(USVD),K(:,2:N+1)),TRANSPOSE(WSVD)),SigmaInv))
DEALLOCATE(SigmaInv)

LWMAX = 2*N
LWORK = -1
ALLOCATE(eigSTilde(N))
ALLOCATE(VL(N,N))
ALLOCATE(VR(N,N))
ALLOCATE(WORKC(LWMAX))
ALLOCATE(RWORK(LWMAX))
ALLOCATE(STildeWork(N,N))
STildeWork = STilde
CALL ZGEEV('V','V',N,STildeWork,N,eigSTilde,VL,N,VR,N,WORKC,LWORK,RWORK,INFO)
LWORK=MIN( LWMAX, INT( WORKC( 1 ) ) )
CALL ZGEEV('V','V',N,STildeWork,N,eigSTilde,VL,N,VR,N,WORKC,LWORK,RWORK,INFO)
DEALLOCATE(VL)              
DEALLOCATE(STildeWork)
DEALLOCATE(WORKC)
DEALLOCATE(RWORK)


ALLOCATE(Vand(N,N))
ALLOCATE(P(N,N))
ALLOCATE(qTmp(N,N))
ALLOCATE(alpha(N))
Vand = 0
DO i = 1, N
  Vand(:,i)=eigSTilde(:)**(i-1)
END DO ! i = 1, N
P = MATMUL(CONJG(TRANSPOSE(VR)),VR) * CONJG(MATMUL(VAND,CONJG(TRANSPOSE(VAND))))

qTmp = MATMUL(MATMUL(MATMUL(Vand,dcmplx(TRANSPOSE(WSVD))),dcmplx(SigmaMat)),VR)

DO i = 1, N
  alpha(i) = CONJG(qTmp(i,i))
END DO ! i = 1, N

!Compute alpha=P/alpha
ALLOCATE(IPIV(N))
CALL ZGESV(N,1,P,N,IPIV,alpha,N,INFO)
DEALLOCATE(IPIV)
DEALLOCATE(Vand)
DEALLOCATE(qTmp)

ALLOCATE(alphaSort(N,2))
ALLOCATE(sigmaSort(N))
DO i = 1, N
  alphaSort(i,1) = ABS(alpha(i))
  alphaSort(i,2) = i
END DO ! i = 1, N

CALL quicksort(alphaSort,1,N,1,2)

ALLOCATE(alphaTmp(N))
alphaTmp = alpha
DO i = 1, N
  alpha(N+1-i) = alphaTmp(alphaSort(i,2))
  sigmaSort(N+1-i) = eigSTilde(alphaSort(i,2))
END DO ! i = 1, N
DEALLOCATE(alphaTmp)

ALLOCATE(PhiTmp(M,N))
ALLOCATE(Phi(M,N))
PhiTmp = MATMUL(USVD,VR)
DO i = 1, N
  Phi(:,N+1-i) = PhiTmp(:,alphaSort(i,2))
END DO ! i = 1, N
DEALLOCATE(PhiTmp)

ALLOCATE(lambda(N))
ALLOCATE(freq(N))
lambda = log(sigmaSort)/dt
freq   = aimag(lambda)/(2*PI)


ALLOCATE(RGlobMat(M,N))
RglobMat = K(:,2:N+1)-MATMUL(MATMUL(MATMUL(USVD,STilde),SigmaMat),WSVD)
LDA  = M
LDU  = M 
LDVT=N
LWMAX = nDoFs+1000
ALLOCATE(WORK(LWMAX))
ALLOCATE(IWORK(8*N))
LWORK = -1
!Calculate norm2 of matrix, simalar to matlab matrix norm
CALL DGESDD( 'S', M, N, RglobMat, LDA, SigmaSVD, USVD, LDU, WSVD, LDVT,WORK, LWORK, IWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
CALL DGESDD( 'S', M, N, RglobMat, LDA, SigmaSVD, USVD, LDU, WSVD, LDVT,WORK, LWORK, IWORK,INFO )
DEALLOCATE(WORK)
DEALLOCATE(IWORK)
DEALLOCATE(RGlobMat)

Rglob = MAXVAL(SigmaSVD)

Rsnap1=norm2(K(:,1)-REAL(MATMUL(Phi,alpha)))/norm2(K(:,1))
Rsnap2=norm2(K(:,N+1)-REAL(MATMUL(Phi,(alpha*sigmaSort**(N)))))/norm2(K(:,N+1))

WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A,F8.6)') 'Global Residual Norm:    ',RGlob
WRITE(UNIT_stdOut,'(A,F8.6)') 'Projection error K(1):   ',Rsnap1
WRITE(UNIT_stdOut,'(A,F8.6)') 'Projection error K(N+1): ',Rsnap2
WRITE(UNIT_stdOut,'(132("="))')


DEALLOCATE(SigmaSVD)
DEALLOCATE(STilde)
DEALLOCATE(SigmaMat)
DEALLOCATE(USVD)
DEALLOCATE(WSVD)
DEALLOCATE(VR)              
DEALLOCATE(eigSTilde)          
END SUBROUTINE performDMD

SUBROUTINE WriteDmdStateFile() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_IO_HDF5
USE MOD_HDF5_Output,        ONLY: WriteState,WriteTimeAverage,GenerateFileSkeleton,WriteAttribute,WriteArray
USE MOD_Output,             ONLY: InitOutput
USE MOD_Output_Vars          
USE MOD_DMD_Vars,           ONLY: Phi,N_State,nElems_State,nModes,freq,alpha,lambda,sigmaSort,Time_State,VarNameDMD
USE MOD_DMD_Vars,           ONLY: dt,TimeEnd_State
USE MOD_Mesh_Vars,          ONLY: offsetElem,nElems,MeshFile
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)   :: DataSetName
INTEGER              :: nVal(4),i,j,k
CHARACTER(LEN=255)   :: FileName,iMode,varFreq
REAL                 :: PhiGlobal(1,0:N_State,0:N_State,0:N_State,nElems_State)
INTEGER              :: Fileunit_DMD
CHARACTER(LEN=255)   :: FileName_DMD
LOGICAL              :: connected
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE DMD MODES TO HDF5 FILE "',TRIM(ProjectName),'" ...'


FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DMD_'//TRIM(VarNameDMD),Time_State))//'.h5'

CALL GenerateFileSkeleton(TRIM(FileName),'DMD',1 ,N_State,(/'DUMMY_DO_NOT_VISUALIZE'/),&
     MeshFile,Time_State,Time_State,create=.TRUE.) ! dummy DG_Solution to fix Posti error !!!

CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
CALL WriteAttribute(File_ID,'DMD Start Time',1,RealScalar=Time_State)
CALL WriteAttribute(File_ID,'DMD END Time',1,RealScalar=TimeEnd_State)
CALL WriteAttribute(File_ID,'DMD dt',1,RealScalar=dt)
CALL CloseDataFile()
  
!================= Actual data output =======================!
! Write DMD
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
k = 0
DO i = 1, nModes
  k=k+1
  IF (freq(i) .LT. 0.) THEN
    k=k-1
    CYCLE
  END IF

  
  PhiGlobal(1,:,:,:,:) = real(RESHAPE( Phi(:,i), (/N_State+1,N_State+1,N_State+1,nElems_State/)))
  WRITE(iMode,'(I0.3)')k
  WRITE(varFreq,'(F8.4)')freq(i)
  ! Replace spaces with 0's
  DO j=1,LEN(TRIM(varFreq))
    IF(varFreq(j:j).EQ.' ') varFreq(j:j)='0'
  END DO
  DataSetName = 'Mode_'//TRIM(iMode)//'_Real_'//TRIM(varFreq)
  CALL WriteArray(        DataSetName=DataSetName, rank=5,&
                          nValGlobal=(/1,N_State+1,N_State+1,N_State+1,nElems_State/),&
                          nVal=      (/1,N_State+1,N_State+1,N_State+1,nElems_State/),&
                          offset=    (/0,      0,     0,     0,     offsetElem/),&
                          collective=.TRUE.,RealArray=PhiGlobal(1:1,:,:,:,:))
  
  PhiGlobal(1,:,:,:,:) = aimag(RESHAPE( Phi(:,i), (/N_State+1,N_State+1,N_State+1,nElems_State/)))
  DataSetName = 'Mode_'//TRIM(iMode)//'_Img_'//TRIM(varFreq)
  CALL WriteArray(        DataSetName=DataSetName, rank=5,&
                          nValGlobal=(/1,N_State+1,N_State+1,N_State+1,nElems_State/),&
                          nVal=      (/1,N_State+1,N_State+1,N_State+1,nElems_State/),&
                          offset=    (/0,      0,     0,     0,     offsetElem/),&
                          collective=.TRUE.,RealArray=PhiGlobal(1:1,:,:,:,:))
END DO ! i = 1, nModes
CALL CloseDataFile()


!Write DMD spectrum to file
FileUnit_DMD=155
INQUIRE(UNIT=FileUnit_DMD, OPENED=connected)
IF (Connected) THEN
  DO
    FileUnit_DMD=Fileunit_DMD+1
    INQUIRE(UNIT=FileUnit_DMD, OPENED=connected)
    IF (.NOT. connected) EXIT
  END DO
END IF
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DMD_Spec',Time_State))//'.dat'
OPEN(FileUnit_DMD,FILE=Filename,STATUS="REPLACE")
WRITE(FileUnit_DMD,*)'TITLE = "DMD Spectrum"'
WRITE(FileUnit_DMD,'(a)')'VARIABLES ="realAlpha"'
WRITE(FileUnit_DMD,'(a)')'"imagAlpha"'
WRITE(FileUnit_DMD,'(a)')'"realLambda"'
WRITE(FileUnit_DMD,'(a)')'"imagLambda"'
WRITE(FileUnit_DMD,'(a)')'"realSigma"'
WRITE(FileUnit_DMD,'(a)')'"imagSigma"'
WRITE(FileUnit_DMD,'(a)')'ZONE T="ZONE1"'
WRITE(FileUnit_DMD,'(a)')' STRANDID=0, SOLUTIONTIME=0'
WRITE(FileUnit_DMD,*)' I=',nModes,', J=1, K=1, ZONETYPE=Ordered'
WRITE(FileUnit_DMD,'(a)')' DATAPACKING=POINT'
WRITE(FileUnit_DMD,'(a)')' DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
DO j=1,nModes
  WRITE(FileUnit_DMD,'(6(E20.12,X))') &
  real(alpha(j)),aimag(alpha(j)),real(lambda(j)),aimag(lambda(j)),real(sigmaSort(j)),aimag(sigmaSort(j))
END DO
CLOSE(FILEUnit_DMD)

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
!!===================================================================================================================================

END SUBROUTINE WriteDmdStateFile

SUBROUTINE FinalizeDMD()
! MODULES
USE MOD_Globals
USE MOD_DMD_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DEALLOCATE(VarNames_State)
DEALLOCATE(Phi)             
DEALLOCATE(freq)
DEALLOCATE(lambda)            
DEALLOCATE(alpha)            
DEALLOCATE(sigmaSort)            
DEALLOCATE(K)

WRITE(UNIT_stdOut,'(A)') '  DMD FINALIZED'
END SUBROUTINE FinalizeDMD

RECURSIVE SUBROUTINE quicksort(a, first, last, icolumn, columns) 
!----------------------------------------------------------------------------------------------------------------------------------!
! Quicksort algorythm: sorts a from first to last entry by the icolumn,  
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)    :: first, last, icolumn,columns
REAL,INTENT(INOUT)    :: a(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i, j
REAL    :: x,t(columns)
!===================================================================================================================================
x = a( (first+last) / 2 ,icolumn)
i = first
j = last
DO
   DO WHILE (a(i,icolumn) < x)
      i=i+1
   END DO
   DO WHILE (x < a(j,icolumn))
     j=j-1
   END DO
   IF (i >= j) EXIT
   t(:) = a(i,:);  a(i,:) = a(j,:);  a(j,:) = t(:)
   i=i+1
   j=j-1
END DO
IF (first < i-1) CALL quicksort(a(:,:), first, i-1, icolumn, columns)
IF (j+1 < last)  CALL quicksort(a(:,:), j+1,  last, icolumn, columns)

END SUBROUTINE quicksort


END MODULE MOD_DMD
