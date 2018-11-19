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

PUBLIC::InitDMD,DefineParametersDMD,performDMD,WriteDMDStateFile
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

END SUBROUTINE DefineParametersDMD 

SUBROUTINE InitDMD() 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Commandline_Arguments
USE MOD_ReadInTools,             ONLY: prms
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension,INTTOSTR
USE MOD_DMD_Vars,                ONLY: nDmdVars,K,KCons,nFiles,nVar_State,N_State,nElems_State,NodeType_State,nDoFs,MeshFile_state
USE MOD_IO_HDF5,                 ONLY: File_ID
USE MOD_Output_Vars,             ONLY: ProjectName,NOut
USE MOD_HDF5_Input,              ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,ReadArray
USE MOD_Mesh_Vars,               ONLY: nElems,OffsetElem,nGlobalElems,MeshFile
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iFile,iVar,offset
REAL,ALLOCATABLE                 :: Utmp(:,:,:,:,:)
!===================================================================================================================================
IF(nArgs.LT.3) CALL CollectiveStop(__STAMP__,'At least two files required for dmd!')
nFiles=nArgs-1

! Open the first statefile to read necessary attributes
CALL OpenDataFile(Args(2),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=MeshFile_state)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State)
CALL CloseDataFile()
! Set everything for output
MeshFile = MeshFile_state
NOut = N_State
nGlobalElems = nElems_State
nElems       = nElems_State
OffsetElem   = 0 ! OffsetElem is 0 since the tool only works on singel

nDoFs = (N_State+1)**3 * nElems_State

nDmdVars = nVar_State
ALLOCATE(K    (nDmdVars,nDoFs,nFiles))
ALLOCATE(KCons(nVar_State   ,nDoFs,nFiles))
ALLOCATE(Utmp (nVar_State   ,0:N_State,0:N_State,0:N_State,nElems_State))

DO iFile = 2, nFiles+1
  WRITE(UNIT_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iFile,' of ',nFiles,' FILES.'
  WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(Args(iFile)),'" )'
  WRITE(UNIT_stdOut,'(132("="))')

  ! Read Statefiles
  CALL OpenDataFile(Args(iFile),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  offset=0
  CALL ReadArray('DG_Solution',5,&
                (/nVar_State,N_State+1,N_State+1,N_State+1,nElems_State/),0,5,RealArray=Utmp)
  CALL CloseDataFile()
  DO iVar=1,nVar_State
    KCons(iVar,:,iFile-1) = RESHAPE( Utmp(iVar,:,:,:,:), (/nDoFs/))
  END DO ! iElem
  !WRITE(*,*) KCons(1,:,:) 
  !DO iVar=1,nDmdVars
    !! Mapping to desired dmd variables
  !END DO ! iElem
END DO ! iFile = 1,nFiles 

END SUBROUTINE InitDMD 

SUBROUTINE performDMD() 
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_DMD_Vars,                ONLY: KCons,nDoFs,nFiles,USVD,SigmaSVD,WSVD,STilde,eigSTilde,VL,VR,Phi
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
DOUBLE PRECISION,ALLOCATABLE     :: SigmaInv(:,:),SigmaMat(:,:)               
COMPLEX,ALLOCATABLE              :: WORKC(:)
REAL,ALLOCATABLE                 :: RWORK(:),alphaSort(:,:)
DOUBLE PRECISION,ALLOCATABLE     :: KTmp(:,:)
COMPLEX,ALLOCATABLE              :: Vand(:,:),qTmp(:,:),q(:),SigmaSort(:)
COMPLEX,ALLOCATABLE              :: P(:,:),Ptmp(:,:),alpha(:),PhiTmp(:,:)
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
KTmp = KCons(2,:,1:N)
CALL DGESDD( 'S', M, N, KTmp, LDA, SigmaSVD, USVD, LDU, WSVD, LDVT,WORK, LWORK, IWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
CALL DGESDD( 'S', M, N, KTmp, LDA, SigmaSVD, USVD, LDU, WSVD, LDVT,WORK, LWORK, IWORK,INFO )


ALLOCATE(SigmaInv(N,N))
ALLOCATE(SigmaMat(N,N))
ALLOCATE(STilde(N,N))

SigmaInv=0.
SigmaMat=0.
DO i = 1, N
  SigmaInv(i,i)=1./SigmaSVD(i)
  SigmaMat(i,i)=SigmaSVD(i)
END DO 

STilde=dcmplx(MATMUL(MATMUL(MATMUL(TRANSPOSE(USVD),KCONS(2,:,2:N+1)),TRANSPOSE(WSVD)),SigmaInv))

LWMAX = 2*N
LWORK = -1
ALLOCATE(eigSTilde(N))
ALLOCATE(VL(N,N))
ALLOCATE(VR(N,N))
ALLOCATE(WORKC(LWMAX))
ALLOCATE(RWORK(LWMAX))
CALL ZGEEV('V','V',N,STilde,N,eigSTilde,VL,N,VR,N,WORKC,LWORK,RWORK,INFO)
LWORK=MIN( LWMAX, INT( WORKC( 1 ) ) )
CALL ZGEEV('V','V',N,STilde,N,eigSTilde,VL,N,VR,N,WORKC,LWORK,RWORK,INFO)
!VR(:,1)=VR(:,1)*(-1)


ALLOCATE(Vand(N,N))
ALLOCATE(P(N,N))
ALLOCATE(qTmp(N,N))
ALLOCATE(q(N))
ALLOCATE(alpha(N))
ALLOCATE(IPIV(N))
Vand = 0
DO i = 1, N
  Vand(:,i)=eigSTilde(:)**(i-1)
END DO ! i = 1, N
P = MATMUL(CONJG(TRANSPOSE(VR)),VR) * CONJG(MATMUL(VAND,CONJG(TRANSPOSE(VAND))))

qTmp = MATMUL(MATMUL(MATMUL(Vand,dcmplx(TRANSPOSE(WSVD))),dcmplx(SigmaMat)),VR)
DO i = 1, N
  q(i) = CONJG(qTmp(i,i))
END DO ! i = 1, N

alpha = q
!Compute alpha=P/q
ALLOCATE(Ptmp(N,N))
Ptmp = P
CALL ZGESV(N,1,Ptmp,N,IPIV,alpha,N,INFO)

ALLOCATE(alphaSort(N,2))
ALLOCATE(SigmaSort(N))
DO i = 1, N
  alphaSort(i,1) = ABS(alpha(i))
  alphaSort(i,2) = i
END DO ! i = 1, N

CALL quicksort(alphaSort,1,N,1,2)

DO i = 1, N
  alpha(N+1-i) = alpha(alphaSort(i,2))
  SigmaSort(N+1-i) = eigSTilde(alphaSort(i,2))
END DO ! i = 1, N
ALLOCATE(PhiTmp(M,N))
ALLOCATE(Phi(M,N))
PhiTmp = MATMUL(USVD,VR)
DO i = 1, N
  Phi(:,N+1-i) = PhiTmp(:,alphaSort(i,2))
END DO ! i = 1, N
print*,'U',USVD(1,1)
print*,'U',USVD(1,2)
print*,'U',USVD(2,1)
!print*,'VR',VR(2,1)






DEALLOCATE(WORK)
DEALLOCATE(USVD)
DEALLOCATE(WSVD)
DEALLOCATE(SigmaSVD)
DEALLOCATE(STilde)
DEALLOCATE(KCons)

END SUBROUTINE performDMD

SUBROUTINE WriteDmdStateFile() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_IO_HDF5
USE MOD_HDF5_Output,        ONLY: WriteState,WriteTimeAverage,GenerateFileSkeleton,WriteAttribute,WriteArray
USE MOD_Output,             ONLY: InitOutput
USE MOD_Output_Vars          
USE MOD_DMD_Vars,           ONLY: Phi,N_State,nElems_State
USE MOD_Mesh_Vars,          ONLY: offsetElem,nElems,MeshFile
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)   :: VarNames(6)
INTEGER              :: nVal(4)
CHARACTER(LEN=255)   :: FileName
REAL                 :: PhiGlobal(1,0:N_State,0:N_State,0:N_State,nElems_State),Time_State
!===================================================================================================================================
Time_State=0.0
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE DMD MODES TO HDF5 FILE "',TRIM(ProjectName),'" ...'
VarNames( 1)= 'Mode1'  

PhiGlobal(1,:,:,:,:) = RESHAPE( Phi(:,1), (/N_State+1,N_State+1,N_State+1,nElems_State/))
print*,'Phi',PhiGlobal(1,:,:,:,1)

!CALL InitOutput()
!nVal(1)=SIZE(VarNames)
!nVal(2)=N_State+1
!nVal(3)=N_State+1
!nVal(4)=N_State+1

FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DMD',Time_State))//'.h5'

CALL GenerateFileSkeleton(TRIM(FileName),'DMD',1 ,N_State,(/'DUMMY_DO_NOT_VISUALIZE'/),&
     MeshFile,Time_State,Time_State,create=.TRUE.) ! dummy DG_Solution to fix Posti error !!!

CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
CALL WriteAttribute(File_ID,'DMD Time',1,RealScalar=Time_State)
CALL CloseDataFile()
  
!================= Actual data output =======================!
! Write DMD
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
CALL WriteArray(        DataSetName='Mode1', rank=5,&
                        nValGlobal=(/1,N_State+1,N_State+1,N_State+1,nElems_State/),&
                        nVal=      (/1,N_State+1,N_State+1,N_State+1,nElems_State/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=PhiGlobal(1:1,:,:,:,:))
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
!!===================================================================================================================================

END SUBROUTINE WriteDmdStateFile

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
