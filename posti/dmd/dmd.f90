!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz 
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
CALL prms%CreateStringOption( "VarNameDMD"         , "Name of variable to perform DMD.", multiple=.TRUE.)
CALL prms%CreateRealOption(   "SvdThreshold"       , "Define relative lower bound of singular values.")
CALL prms%CreateIntOption(    "nModes"             , "Number of Modes to visualize.")
CALL prms%CreateLogicalOption("sortFreq"           , "Sort modes by Frequency.")
CALL prms%CreateLogicalOption("PlotSingleMode"     , "Decide if a single mode is plotted.")
CALL prms%CreateRealOption(   "ModeFreq"           , "Specify the mode frequency.")
CALL prms%CreateLogicalOption("UseBaseflow"        , "Subtract the Baseflow-File from every single Snapshot",'.FALSE.')
CALL prms%CreateStringOption( "Baseflow"           , "Name of the Baseflow-File")
CALL prms%CreateLogicalOption("use2D"              , "Calculate DMD on 2D-Data",'.FALSE.')
CALL prms%CreateLogicalOption('output2D'           , "Set true to activate hdf5 data output with flat third dimension.",'.FALSE.')

END SUBROUTINE DefineParametersDMD 

SUBROUTINE InitDMD() 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Commandline_Arguments
USE MOD_Readintools             ,ONLY: GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY,CountOption
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension,INTTOSTR
USE MOD_DMD_Vars,                ONLY: K,nFiles,nVar_State,N_State,nElems_State,NodeType_State,nDoFs,MeshFile_state
USE MOD_DMD_Vars,                ONLY: dt,nModes,SvdThreshold,VarNameDMD,VarNames_State,Time_State,TimeEnd_State,sortFreq
USE MOD_DMD_Vars,                ONLY: ModeFreq,PlotSingleMode,useBaseflow,Baseflow,VarNames_TimeAvg,nVarDMD,use2D,N_StateZ
USE MOD_DMD_Vars,                ONLY: N_StateZ,N_StateZ_out
USE MOD_IO_HDF5,                 ONLY: File_ID,HSize,output2D
USE MOD_Output_Vars,             ONLY: ProjectName,NOut
USE MOD_HDF5_Input,              ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,ReadArray,GetDataSize
USE MOD_Mesh_Vars,               ONLY: OffsetElem,nGlobalElems,MeshFile
USE MOD_EquationDMD,             ONLY: InitEquationDMD,CalcEquationDMD
USE MOD_Interpolation_Vars ,     ONLY: NodeType
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iFile,iVar,offset,nDim,i,j
REAL,ALLOCATABLE                 :: Utmp(:,:,:,:,:)
REAL                             :: time
DOUBLE PRECISION,ALLOCATABLE     :: Ktmp(:)
CHARACTER(LEN=255)               :: FileType
INTEGER,ALLOCATABLE              :: VarSortTimeAvg(:)
!===================================================================================================================================
IF(nArgs.LT.3) CALL CollectiveStop(__STAMP__,'At least two files required for dmd!')
nFiles=nArgs-1

nVarDMD=CountOption("VarNameDMD")
ALLOCATE(VarNameDMD(nVarDMD))
DO iVar=1,nVarDMD
  VarNameDMD(iVAr)   = GETSTR ('VarNameDMD'  ,'Density')
END DO
SvdThreshold = GETREAL('SvdThreshold','0.0')
nModes       = GETINT ('nModes'      ,'9999')
sortFreq     = GETLOGICAL ('sortFreq','False')
PlotSingleMode = GETLOGICAL ('PlotSingleMode','False')
IF (PlotSingleMode) ModeFreq = GETREAL('ModeFreq','0.0')
useBaseflow  = GETLOGICAL('UseBaseflow','.FALSE.')
IF (useBaseflow) Baseflow    = GETSTR    ('Baseflow','none')
use2D        = GETLOGICAL('use2D','.FALSE.')
output2D     = GETLOGICAL('output2D','.FALSE.')

IF (nModes .GT. (nFiles-1)) THEN
  nModes = nFiles - 1  
END IF ! nModes .GT. (nFiles-1)

! Open the first statefile to read necessary attributes
CALL OpenDataFile(Args(2),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=MeshFile_state)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',    1,RealScalar=starttime)
CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State)
IF ( NodeType .NE. NodeType_State  ) THEN
  CALL CollectiveStop(__STAMP__,'Node Type of given state File is not applicable to compiled version')
END IF
IF (.NOT. use2D) THEN
  N_StateZ     = N_State
ELSE
  N_StateZ     = 0
END IF

IF (output2D) THEN
  N_StateZ_out = 0
ELSE
  N_StateZ_out = N_State
END IF

ALLOCATE(VarNames_State(nVar_State))
CALL ReadAttribute(File_ID,'VarNames',nVar_State,StrArray=VarNames_State)
CALL GetDataSize(File_ID,'DG_Solution',nDim,HSize)
IF (HSize(4) .EQ. 1) THEN
  use2D=.TRUE.
  N_StateZ=0
  WRITE(UNIT_stdOut,'(A)') 'Provided states contain 2D-data, continue with two dimensional DMD'
END IF
CALL CloseDataFile()
! Set everything for output
Time_State = starttime
MeshFile = MeshFile_state
NOut = N_State
nGlobalElems = nElems_State
OffsetElem   = 0 ! OffsetElem is 0 since the tool only works on singel

IF (.NOT. use2D) THEN
  nDoFs = (N_State+1)**3 * nElems_State
ELSE
  nDoFs = (N_State+1)**2 * nElems_State
END IF


CALL InitEquationDMD()

ALLOCATE(K    (nDoFs*nVarDMD,nFiles))
ALLOCATE(Utmp (nVar_State   ,0:N_State,0:N_State,0:N_StateZ,nElems_State))


DO iFile = 2, nFiles+1
  WRITE(UNIT_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iFile-1,' of ',nFiles,' FILES.'
  WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(Args(iFile)),'" )'
  WRITE(UNIT_stdOut,'(132("="))')

  ! Read Statefiles
  CALL OpenDataFile(Args(iFile),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  offset=0
  CALL ReadAttribute(File_ID,'Time',    1,RealScalar=time)
  CALL ReadArray('DG_Solution',5,&
                (/nVar_State,N_State+1,N_State+1,N_StateZ+1,nElems_State/),0,5,RealArray=Utmp)


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

IF (useBaseflow) THEN
  CALL OpenDataFile(Baseflow,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadAttribute(File_ID,'File_Type',1,StrScalar =FileType)
  SELECT CASE(TRIM(FileType))
  CASE('State')
    ALLOCATE(Ktmp(nDoFs*nVarDMD))
    CALL ReadAttribute(File_ID,'Time',    1,RealScalar=time)
    CALL ReadArray('DG_Solution',5,&
                  (/nVar_State,N_State+1,N_State+1,N_StateZ+1,nElems_State/),0,5,RealArray=Utmp)
    
    CALL CloseDataFile()
    
    CALL CalcEquationDMD(Utmp,Ktmp(:))
    DO iFile = 2,nFiles+1
      K(:,iFile-1) = K(:,iFile-1) - Ktmp
    END DO
    DEALLOCATE(ktmp)
  CASE('TimeAvg')
    ALLOCATE(VarSortTimeAvg(nVarDMD))
    CALL GetDataSize(File_ID,'Mean',nDim,HSize)
    ALLOCATE(VarNames_TimeAvg(INT(HSize(1))))
    CALL ReadAttribute(File_ID,'VarNames_Mean',INT(HSize(1)),StrArray=VarNames_TimeAvg)
    DO j=1,nVarDMD
      DO i=1,INT(HSize(1))
        IF (TRIM(VarNameDMD(j)) .EQ. TRIM(VarNames_TimeAvg(i))) THEN
          VarSortTimeAvg(i)=j
          EXIT
        ELSE IF(i .EQ. INT(HSize(1))) THEN 
          CALL Abort(__STAMP__,'Required DMD-Var not provided by TimeAvg-File')
        END IF 
      END DO
    END DO

    ALLOCATE(Ktmp(nDoFs*nVarDMD))
    DEALLOCATE(Utmp)
    ALLOCATE(Utmp (HSize(1),0:N_State,0:N_State,0:N_StateZ,nElems_State))
    
    CALL ReadArray('Mean',5,&
                  (/INT(HSize(1)),N_State+1,N_State+1,N_StateZ+1,nElems_State/),0,5,RealArray=Utmp)
    CALL CloseDataFile()

    ktmp(:) = RESHAPE(Utmp(VarSortTimeAvg,:,:,:,:), (/nDoFs*nVarDMD/))
    
    DO iFile = 2,nFiles+1
      K(:,iFile-1) = K(:,iFile-1) - Ktmp
    END DO
      
    DEALLOCATE(ktmp)
  END SELECT

END IF

DEALLOCATE(Utmp)

END SUBROUTINE InitDMD 

SUBROUTINE performDMD() 
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_DMD_Vars,                ONLY: K,nDoFs,nFiles,Phi,lambda,freq,dt,alpha,sigmaSort,sortFreq,SvdThreshold,nModes,nVarDMD
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
INTEGER                          :: i,M,N,LDA,LDU,LDVT,LWMAX,INFO,LWORK,nFilter
INTEGER,ALLOCATABLE              :: IWORK(:),IPIV(:),filter(:)
DOUBLE PRECISION,ALLOCATABLE     :: WORK(:)
DOUBLE PRECISION,ALLOCATABLE     :: SigmaInv(:,:),SigmaMat(:,:),USVD(:,:),WSVD(:,:),SigmaSVD(:)
DOUBLE PRECISION,ALLOCATABLE     :: USVDTmp(:,:),WSVDTmp(:,:),SigmaSVDTmp(:)
COMPLEX,ALLOCATABLE              :: WORKC(:)
REAL,ALLOCATABLE                 :: RWORK(:),SortVar(:,:)
DOUBLE PRECISION,ALLOCATABLE     :: KTmp(:,:),RGlobMat(:,:),Rglob,Rsnap1,Rsnap2
COMPLEX,ALLOCATABLE              :: STilde(:,:),STildeWork(:,:),eigSTilde(:),VL(:,:),VR(:,:)                
COMPLEX,ALLOCATABLE              :: Vand(:,:),qTmp(:,:)
COMPLEX,ALLOCATABLE              :: P(:,:),PhiTmp(:,:),alphaTmp(:)
REAL                             :: pi = acos(-1.)
!===================================================================================================================================
M = nDoFs*nVarDMD
N = nFiles-1
LDA  = M
LDU  = M 
LDVT = N
LWORK = -1
LWMAX = 3*MIN(M,N)*MIN(M,N) + MAX(MAX(M,N),4*MIN(M,N)*MIN(M,N)+4*MIN(M,N))
ALLOCATE(WORK(LWMAX))
ALLOCATE(USVD(M,N))
ALLOCATE(WSVD(N,N))
ALLOCATE(SigmaSVD(min(N,M)))

ALLOCATE(KTmp(nDoFs*nVarDMD,nFiles-1))
ALLOCATE(IWORK(8*N))
KTmp = K(:,1:N)
CALL DGESDD( 'S', M, N, KTmp, LDA, SigmaSVD, USVD, LDU, WSVD, LDVT,WORK, LWORK, IWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
CALL DGESDD( 'S', M, N, KTmp, LDA, SigmaSVD, USVD, LDU, WSVD, LDVT,WORK, LWORK, IWORK,INFO )

DEALLOCATE(WORK)
DEALLOCATE(KTmp)
DEALLOCATE(IWORK)
IF( SvdThreshold .GT. 0.) THEN
  nFilter = 0
  ALLOCATE(filter(N))
  SvdThreshold=MAXVAL(SigmaSVD(2:N))*SvdThreshold
  DO i = 1, N
    IF (SigmaSVD(i) .GT. SvdThreshold) THEN
      nFilter = nFilter+1
      filter(nFilter) = i
    END IF
  END DO ! i = 1, N

  ALLOCATE(WSVDTmp(nFilter,nFiles-1))
  ALLOCATE(USVDTmp(M,nFilter))
  ALLOCATE(SigmaSVDTmp(min(nFilter,M)))
  DO i = 1, nFilter
    USVDTmp(:,i)=USVD(:,filter(i))
    WSVDTmp(i,:)=WSVD(filter(i),:)
    SigmaSVDTmp(i)=SigmaSVD(filter(i))
  END DO ! i = 1, nFilter
  DEALLOCATE(USVD)
  DEALLOCATE(WSVD)
  DEALLOCATE(SigmaSVD)
  ALLOCATE(WSVD(nFilter,N))
  
  N=nFilter
  ALLOCATE(USVD(M,N))
  ALLOCATE(SigmaSVD(min(N,M)))
  USVD    = USVDTmp
  WSVD    = WSVDTmp
  SigmaSVD= SigmaSVDTmp
  DEALLOCATE(USVDTmp)
  DEALLOCATE(WSVDTmp)
  DEALLOCATE(SigmaSVDTmp)
  DEALLOCATE(filter)
  IF (nModes .GT. nFilter) THEN
    nModes = nFilter  
  END IF 
END IF

ALLOCATE(SigmaInv(N,N))
ALLOCATE(SigmaMat(N,N))
ALLOCATE(STilde(N,N))

SigmaInv=0.
SigmaMat=0.
DO i = 1, N
  SigmaInv(i,i)=1./SigmaSVD(i)
  SigmaMat(i,i)=SigmaSVD(i)
END DO 

STilde=dcmplx(MATMUL(MATMUL(MATMUL(TRANSPOSE(USVD),K(:,2:nFiles)),TRANSPOSE(WSVD)),SigmaInv))
DEALLOCATE(SigmaInv)

LWMAX = 3*MIN(M,N)*MIN(M,N) + MAX(MAX(M,N),4*MIN(M,N)*MIN(M,N)+4*MIN(M,N))
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


ALLOCATE(Vand(N,nFiles-1))
ALLOCATE(P(N,N))
ALLOCATE(qTmp(N,N))
ALLOCATE(alpha(N))
Vand = 0
DO i = 1, nFiles-1
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

ALLOCATE(SortVar(N,2))
ALLOCATE(sigmaSort(N))
ALLOCATE(freq(N))
freq   = aimag(log(eigSTilde)/dt)/(2*PI)
IF (sortFreq) THEN
  DO i = 1, N
    SortVar(i,1) = 1./freq(i)
    SortVar(i,2) = i
  END DO ! i = 1, N
ELSE
  DO i = 1, N
    SortVar(i,1) = ABS(alpha(i))
    SortVar(i,2) = i
  END DO ! i = 1, N
END IF ! sortFreq
CALL quicksort(SortVar,1,N,1,2)

ALLOCATE(alphaTmp(N))
alphaTmp = alpha
DO i = 1, N
  alpha(N+1-i) = alphaTmp(INT(SortVar(i,2)))
  sigmaSort(N+1-i) = eigSTilde(INT(SortVar(i,2)))
END DO ! i = 1, N
DEALLOCATE(alphaTmp)

ALLOCATE(PhiTmp(M,N))
ALLOCATE(Phi(M,N))
PhiTmp = MATMUL(USVD,VR)
DO i = 1, N
  Phi(:,N+1-i) = PhiTmp(:,INT(SortVar(i,2)))
END DO ! i = 1, N
DEALLOCATE(PhiTmp)

ALLOCATE(lambda(N))
lambda = log(sigmaSort)/dt
freq   = aimag(lambda)/(2*PI)


ALLOCATE(RGlobMat(M,nFiles-1))
RglobMat = K(:,2:nFiles)-DBLE(MATMUL(MATMUL(MATMUL(USVD,STilde),SigmaMat),WSVD))
LDA  = M
LDU  = M 
LDVT = N
LWMAX = nDoFs*nVarDMD*5
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
Rsnap2=norm2(K(:,nFiles)-REAL(MATMUL(Phi,(alpha*sigmaSort**(N)))))/norm2(K(:,nFiles))

WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A,F12.6)') 'Global Residual Norm:    ',RGlob
WRITE(UNIT_stdOut,'(A,F12.6)') 'Projection error K(1):   ',Rsnap1
WRITE(UNIT_stdOut,'(A,F12.6)') 'Projection error K(nFiles): ',Rsnap2
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
USE MOD_DMD_Vars,           ONLY: Phi,N_State,N_StateZ,nElems_State,nModes,freq,alpha,lambda,sigmaSort,Time_State,VarNameDMD
USE MOD_DMD_Vars,           ONLY: dt,TimeEnd_State,ModeFreq,PlotSingleMode,nVarDMD,nDofs,N_StateZ_out
USE MOD_Mesh_Vars,          ONLY: offsetElem,MeshFile
USE MOD_2D                 ,ONLY: ExpandArrayTo3D
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)   :: DataSetName
INTEGER              :: i,j,k
CHARACTER(LEN=255)   :: FileName,iMode,varFreq
REAL,TARGET          :: PhiGlobal(nVarDMD,0:N_State,0:N_State,0:N_StateZ,nElems_State)
REAL,TARGET          :: PhiGlobal_expand3D(nVarDMD,0:N_State,0:N_State,0:N_State,nElems_State)
REAL,POINTER         :: PhiGlobal_out(:,:,:,:,:)
INTEGER              :: Fileunit_DMD
LOGICAL              :: connected
COMPLEX              :: Phi_out(nDofs,nVarDMD,nModes)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE DMD MODES TO HDF5 FILE "',TRIM(ProjectName),'" ...'


FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DMD',Time_State))//'.h5'

CALL GenerateFileSkeleton(TRIM(FileName),'DMD',1 ,N_State,(/'DUMMY_DO_NOT_VISUALIZE'/),&
     MeshFile,Time_State,Time_State,create=.TRUE.) ! dummy DG_Solution to fix Posti error !!!

CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
CALL WriteAttribute(File_ID,'DMD Start Time',1,RealScalar=Time_State)
CALL WriteAttribute(File_ID,'DMD END Time',1,RealScalar=TimeEnd_State)
CALL WriteAttribute(File_ID,'DMD dt',1,RealScalar=dt)
CALL CloseDataFile()

DO i=1,nModes
  Phi_out(:,:,i)=TRANSPOSE(RESHAPE(Phi(:,i), (/nVarDMD,nDofs/)))
END DO

!================= Actual data output =======================!
! Write DMD
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
IF (PlotSingleMode) THEN
  k=1
  DO i = 1, nModes
    IF (ABS(freq(i)-ModeFreq) .LT. ABS(freq(k)-ModeFreq)) THEN
      k=i
    END IF
  END DO
  WRITE(iMode,'(I0.3)')1
  WRITE(varFreq,'(F12.5)')freq(k)
  ! Replace spaces with 0's
  DO j=1,LEN(TRIM(varFreq))
    IF(varFreq(j:j).EQ.' ') varFreq(j:j)='0'
  END DO
  DO i=1,nVarDMD
    PhiGlobal(i,:,:,:,:) = real(RESHAPE( Phi_out(:,i,k), (/N_State+1,N_State+1,N_StateZ+1,nElems_State/)))
    DataSetName = 'Mode_'//TRIM(iMode)//'_'//TRIM(varNameDMD(i))//'_Real'
    IF(output2D) THEN
      PhiGlobal_out=>PhiGlobal
    ELSE
      CALL ExpandArrayTo3D(5,nVal=(/1,N_State+1,N_State+1,N_StateZ+1,nElems_State/),index3D=4,size3D=N_State+1,arrayIn=PhiGlobal(j,:,:,:,:),arrayOut=PhiGlobal_expand3D(j,:,:,:,:))
      PhiGlobal_out=>PhiGlobal_expand3D
    END IF
    CALL WriteArray(        DataSetName=DataSetName, rank=5,&
                            nValGlobal=(/1,N_State+1,N_State+1,N_StateZ_out+1,nElems_State/),&
                            nVal=      (/1,N_State+1,N_State+1,N_StateZ_out+1,nElems_State/),&
                            offset=    (/0,      0,     0,     0,     offsetElem/),&
                            collective=.TRUE.,RealArray=PhiGlobal_out(i,:,:,:,:))

    PhiGlobal(i,:,:,:,:) = aimag(RESHAPE( Phi_out(:,i,k), (/N_State+1,N_State+1,N_StateZ+1,nElems_State/)))
    DataSetName = 'Mode_'//TRIM(iMode)//'_'//TRIM(varNameDMD(i))//'_Img'
    IF(output2D) THEN
      PhiGlobal_out=>PhiGlobal
    ELSE
      CALL ExpandArrayTo3D(5,nVal=(/1,N_State+1,N_State+1,N_StateZ+1,nElems_State/),index3D=4,size3D=N_State+1,arrayIn=PhiGlobal(j,:,:,:,:),arrayOut=PhiGlobal_expand3D(j,:,:,:,:))
      PhiGlobal_out=>PhiGlobal_expand3D
    END IF
    CALL WriteArray(        DataSetName=DataSetName, rank=5,&
                            nValGlobal=(/1,N_State+1,N_State+1,N_StateZ_out+1,nElems_State/),&
                            nVal=      (/1,N_State+1,N_State+1,N_StateZ_out+1,nElems_State/),&
                            offset=    (/0,      0,     0,     0,     offsetElem/),&
                            collective=.TRUE.,RealArray=PhiGlobal_out(i,:,:,:,:))
  END DO
ELSE
  k = 0
  DO i = 1, nModes
    k=k+1
    IF (freq(i) .LT. 0.) THEN
      k=k-1
      CYCLE
    END IF
  
    
    WRITE(iMode,'(I0.3)')k
    WRITE(varFreq,'(F12.5)')freq(i)
    ! Replace spaces with 0's
    DO j=1,LEN(TRIM(varFreq))
      IF(varFreq(j:j).EQ.' ') varFreq(j:j)='0'
    END DO
    DO j=1,nVarDMD
      PhiGlobal(j,:,:,:,:) = real(RESHAPE( Phi_out(:,j,i), (/N_State+1,N_State+1,N_StateZ+1,nElems_State/)))
      DataSetName = 'Mode_'//TRIM(iMode)//'_'//TRIM(varNameDMD(j))//'_Real'
      IF(output2D) THEN
        PhiGlobal_out=>PhiGlobal
      ELSE
        CALL ExpandArrayTo3D(5,nVal=(/1,N_State+1,N_State+1,N_StateZ+1,nElems_State/),index3D=4,size3D=N_State+1,arrayIn=PhiGlobal(j,:,:,:,:),arrayOut=PhiGlobal_expand3D(j,:,:,:,:))
        PhiGlobal_out=>PhiGlobal_expand3D
      END IF
      CALL WriteArray(        DataSetName=DataSetName, rank=5,&
                              nValGlobal=(/1,N_State+1,N_State+1,N_StateZ_out+1,nElems_State/),&
                              nVal=      (/1,N_State+1,N_State+1,N_StateZ_out+1,nElems_State/),&
                              offset=    (/0,      0,     0,     0,     offsetElem/),&
                              collective=.TRUE.,RealArray=PhiGlobal_out(1,:,:,:,:))
    
      PhiGlobal(j,:,:,:,:) = aimag(RESHAPE( Phi_out(:,j,i), (/N_State+1,N_State+1,N_StateZ+1,nElems_State/)))
      IF(output2D) THEN
        PhiGlobal_out=>PhiGlobal
      ELSE
        CALL ExpandArrayTo3D(5,nVal=(/1,N_State+1,N_State+1,N_StateZ+1,nElems_State/),index3D=4,size3D=N_State+1,arrayIn=PhiGlobal(j,:,:,:,:),arrayOut=PhiGlobal_expand3D(j,:,:,:,:))
        PhiGlobal_out=>PhiGlobal_expand3D
      END IF
      DataSetName = 'Mode_'//TRIM(iMode)//'_'//TRIM(varNameDMD(j))//'_Img'
      CALL WriteArray(        DataSetName=DataSetName, rank=5,&
                            nValGlobal=(/1,N_State+1,N_State+1,N_StateZ_out+1,nElems_State/),&
                            nVal=      (/1,N_State+1,N_State+1,N_StateZ_out+1,nElems_State/),&
                            offset=    (/0,      0,     0,     0,     offsetElem/),&
                            collective=.TRUE.,RealArray=PhiGlobal_out(1,:,:,:,:))
    END DO
  END DO ! i = 1, nModes
END IF ! PlotSingleMode
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
DEALLOCATE(VarNameDMD)
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
