!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz 
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

!==================================================================================================================================
!> Module that provides functions for computing the solutions time history at a defined set of points ("recordpoints")
!==================================================================================================================================
MODULE MOD_RecordPoints
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitRecordPoints
  MODULE PROCEDURE InitRecordPoints
END INTERFACE

INTERFACE RecordPoints
  MODULE PROCEDURE RecordPoints
END INTERFACE

INTERFACE WriteRP
  MODULE PROCEDURE WriteRP
END INTERFACE

INTERFACE FinalizeRecordPoints
  MODULE PROCEDURE FinalizeRecordPoints
END INTERFACE

PUBLIC::InitRecordPoints,RecordPoints,FinalizeRecordPoints,WriteRP
!==================================================================================================================================

PUBLIC::DefineParametersRecordPoints
CONTAINS


!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersRecordPoints()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("RecordPoints")
CALL prms%CreateLogicalOption('RP_inUse',          "Set true to compute solution history at points defined in recordpoints file.",&
                                                   '.FALSE.')
CALL prms%CreateStringOption( 'RP_DefFile',        "File containing element-local parametric recordpoint coordinates and structure.")
CALL prms%CreateIntOption(    'RP_MaxMemory',      "Maximum memory in MiB to be used for storing recordpoint state history. "//&
                                                   "If memory is exceeded before regular IO level states are written to file.",&
                                                   '100')
CALL prms%CreateIntOption(    'RP_SamplingOffset', "Multiple of timestep at which recordpoints are evaluated.",&
                                                   '1')
END SUBROUTINE DefineParametersRecordPoints


!==================================================================================================================================
!> Read RP parameters from ini file and RP definitions from HDF5
!==================================================================================================================================
SUBROUTINE InitRecordPoints(RPDefFileOpt)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools         ,ONLY: GETSTR,GETINT,GETLOGICAL,GETREAL
USE MOD_Interpolation_Vars  ,ONLY: InterpolationInitIsDone
USE MOD_RecordPoints_Vars
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: RPDefFileOpt
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: RP_maxMemory
INTEGER               :: maxRP
!==================================================================================================================================
! check if recordpoints are activated
RP_inUse=(GETLOGICAL('RP_inUse','.FALSE.') .OR. PRESENT(RPDefFileOpt))
IF(.NOT.RP_inUse) RETURN

IF((.NOT.InterpolationInitIsDone) .OR. RecordPointsInitIsDone)THEN
   CALL Abort(__STAMP__,&
     "InitRecordPoints not ready to be called or already called.")
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RECORDPOINTS...'

IF(PRESENT(RPDefFileOpt))THEN
  RPDefFile=RPDefFileOpt
ELSE
  RPDefFile=GETSTR('RP_DefFile')                        ! Filename with RP coords
END IF
CALL ReadRPList(RPDefFile) ! RP_inUse is set to FALSE by ReadRPList if no RP is on proc.
maxRP=nGlobalRP
#if USE_MPI
CALL InitRPCommunicator()
#endif /*USE_MPI*/

IF(RP_onProc)THEN
  RP_maxMemory=GETINT('RP_MaxMemory','100')           ! Max buffer (100MB)
  RP_SamplingOffset=GETINT('RP_SamplingOffset','1')   ! Sampling offset (iteration)
  maxRP=nGlobalRP
#if USE_MPI
  CALL MPI_ALLREDUCE(nRP,maxRP,1,MPI_INTEGER,MPI_MAX,RP_COMM,iError)
#endif /*USE_MPI*/
  RP_MaxBufferSize = RP_MaxMemory*131072/(maxRP*(PP_nVar+1)) != size in bytes/(real*maxRP*nVar)
  ALLOCATE(lastSample(0:PP_nVar,nRP))
  lastSample=0.
END IF

RecordPointsInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RECORDPOINTS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRecordPoints


#if USE_MPI
!==================================================================================================================================
!> Read RP parameters from ini file and RP definitions from HDF5
!==================================================================================================================================
SUBROUTINE InitRPCommunicator()
! MODULES
USE MOD_Globals
USE MOD_RecordPoints_Vars   ,ONLY: RP_onProc,myRPrank,RP_COMM,nRP_Procs
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: color,iProc
INTEGER                   :: noRPrank,RPrank
LOGICAL                   :: hasRP
!==================================================================================================================================
color=MPI_UNDEFINED
IF(RP_onProc) color=2

! create ranks for RP communicator
IF(MPIRoot) THEN
  RPrank=-1
  noRPrank=-1
  myRPRank=0
  IF(RP_onProc) THEN
    RPrank=0
  ELSE
    noRPrank=0
  END IF
  DO iProc=1,nProcessors-1
    CALL MPI_RECV(hasRP,1,MPI_LOGICAL,iProc,0,MPI_COMM_FLEXI,MPIstatus,iError)
    IF(hasRP) THEN
      RPrank=RPrank+1
      CALL MPI_SEND(RPrank,1,MPI_INTEGER,iProc,0,MPI_COMM_FLEXI,iError)
    ELSE
      noRPrank=noRPrank+1
      CALL MPI_SEND(noRPrank,1,MPI_INTEGER,iProc,0,MPI_COMM_FLEXI,iError)
    END IF
  END DO
ELSE
    CALL MPI_SEND(RP_onProc,1,MPI_LOGICAL,0,0,MPI_COMM_FLEXI,iError)
    CALL MPI_RECV(myRPrank,1,MPI_INTEGER,0,0,MPI_COMM_FLEXI,MPIstatus,iError)
END IF

! create new RP communicator for RP output
CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI, color, myRPrank, RP_COMM,iError)
IF(RP_onProc) CALL MPI_COMM_SIZE(RP_COMM, nRP_Procs,iError)
IF(myRPrank.EQ.0 .AND. RP_onProc) WRITE(*,*) 'RP COMM:',nRP_Procs,'procs'

END SUBROUTINE InitRPCommunicator
#endif /*USE_MPI*/


!==================================================================================================================================
!> Read Recordpoint coordinates from HDF5 file
!==================================================================================================================================
SUBROUTINE ReadRPList(FileString)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input
USE MOD_Mesh_Vars             ,ONLY:MeshFile,nGlobalElems
USE MOD_Mesh_Vars             ,ONLY:OffsetElem
USE MOD_Mesh_Vars             ,ONLY:nElems
USE MOD_RecordPoints_Vars     ,ONLY:RP_onProc,L_xi_RP,L_eta_RP,L_zeta_RP
USE MOD_RecordPoints_Vars     ,ONLY:OffsetRP,RP_ElemID,nRP,nGlobalRP,offsetRP
#if FV_ENABLED
USE MOD_RecordPoints_Vars     ,ONLY:FV_RP_ijk
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileString
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: MeshFile_RPList
INTEGER                       :: nGlobalElems_RPList
INTEGER                       :: iElem,iRP,iRP_glob
INTEGER                       :: OffsetRPArray(2,nElems)
REAL,ALLOCATABLE              :: xi_RP(:,:)
!==================================================================================================================================
IF(MPIRoot)THEN
  IF(.NOT.FILEEXISTS(FileString))  CALL abort(__STAMP__, &
          'RPList from data file "'//TRIM(FileString)//'" does not exist')
END IF

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Read recordpoint definitions from data file "'//TRIM(FileString)//'" ...'
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! compare mesh file names
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile_RPList)
IF(TRIM(MeshFile_RPList).NE.TRIM(MeshFile)) THEN
  SWRITE(UNIT_stdOut,*) ' WARNING: MeshFileName from RPList differs from Mesh File specified in parameterfile!'
END IF

! Readin OffsetRP
CALL GetDataSize(File_ID,'OffsetRP',nDims,HSize)
CHECKSAFEINT(HSize(2),4)
nGlobalElems_RPList=INT(HSize(2),4) !global number of elements
DEALLOCATE(HSize)
IF(nGlobalElems_RPList.NE.nGlobalElems) CALL abort(__STAMP__, &
          'nGlobalElems from RPList differs from nGlobalElems from Mesh File!')

CALL ReadArray('OffsetRP',2,(/2,nElems/),OffsetElem,2,IntArray=OffsetRPArray)

! Check if local domain contains any record points
! OffsetRP: first index: 1: offset in RP list for first RP on elem,
!                        2: offset in RP list for last RP on elem
! If these offsets are equal, no RP on elem.
nRP=OffsetRPArray(2,nElems)-OffsetRPArray(1,1)
offsetRP = OffsetRPArray(1,1)
! Read in RP reference coordinates
CALL GetDataSize(File_ID,'xi_RP',nDims,HSize)
CHECKSAFEINT(HSize(2),4)
nGlobalRP=INT(HSize(2),4) !global number of RecordPoints
DEALLOCATE(HSize)
ALLOCATE(xi_RP(3,nRP))
CALL ReadArray('xi_RP',2,(/3,nRP/),offsetRP,2,RealArray=xi_RP)


IF(nRP.LT.1) THEN
  RP_onProc=.FALSE.
ELSE
  RP_onProc=.TRUE.
  ! create mapping to elements
  ALLOCATE(RP_ElemID(nRP))
  DO iRP=1,nRP
    iRP_glob=offsetRP+iRP
    DO iElem=1,nElems
      IF(iRP_glob .LE. OffsetRPArray(2,iElem) .AND. iRP_glob .GT. OffsetRPArray(1,iElem)) &
        RP_ElemID(iRP)=iElem
    END DO
  END DO
END IF
CALL CloseDataFile()

IF(RP_onProc)THEN
  ALLOCATE(L_xi_RP(0:PP_N,nRP), L_eta_RP(0:PP_N,nRP), L_zeta_RP(0:PP_N,nRP))
  CALL InitRPBasis(nRP,xi_RP,L_xi_RP,L_eta_RP,L_zeta_RP)
#if FV_ENABLED
  ALLOCATE(FV_RP_ijk(3,nRP))
  !=====================================================================================
  ! Two variants possible for FV:
  ! 1. The RP state is the nearest (reference-space) FV cells cell average:
  !    + Most simple solution
  !    - Spatial accuracy
  !    - Parameter space estimate can be wrong in case of strongly deformed meshes
  ! 2. The RP state is obtained by tri-linear interpolation from the 8 nearest
  !    (physical space) FV cell averages:
  !    + Probably gives the best quality
  !    - Requires geometry info in physical space
  !    - General implementation requires interpolation across macro-cell boundaries.
  !      Very difficult especially in an MPI setting.
  !
  ! We implement the first variant, the second may be an option if higher accuracy
  ! is desired, possibly with only element local interpolation.
  !=====================================================================================
  FV_RP_ijk=INT((xi_RP+1.)*0.5*(PP_N+1))
  FV_RP_ijk=MAX(FV_RP_ijk,0)
  FV_RP_ijk=MIN(FV_RP_ijk,PP_N)

#if PP_dim==2
  FV_RP_ijk(3,:)=0
#endif

#endif
END IF
DEALLOCATE(xi_RP)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' done.'
END SUBROUTINE ReadRPList



!==================================================================================================================================
!> Precompute Lagrange basis function values at recordpoints
!==================================================================================================================================
SUBROUTINE InitRPBasis(nRP,xi_RP,L_xi_RP,L_eta_RP,L_zeta_RP)
! MODULES
USE MOD_PreProc
USE MOD_Interpolation_Vars    ,ONLY: xGP,wBary
USE MOD_Basis                 ,ONLY: LagrangeInterpolationPolys
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nRP                      !< size of recordpointarray
REAL,INTENT(IN)               :: xi_RP(3,nRP)             !< coordinates of recordpoints in reference space
REAL,INTENT(OUT)              :: L_xi_RP(  0:PP_N,nRP)    !< Lagrange basis evaluated at recordpoints (\f$\xi\f$-direction)
REAL,INTENT(OUT)              :: L_eta_RP( 0:PP_N,nRP)    !< Lagrange basis evaluated at recordpoints (\f$\eta\f$-direction)
REAL,INTENT(OUT)              :: L_zeta_RP(0:PP_N,nRP)    !< Lagrange basis evaluated at recordpoints (\f$\zeta\f$-direction)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iRP
!==================================================================================================================================
! build local basis for Recordpoints
DO iRP=1,nRP
  CALL LagrangeInterpolationPolys(xi_RP(1,iRP),PP_N,xGP,wBary,L_xi_RP(:,iRP))
  CALL LagrangeInterpolationPolys(xi_RP(2,iRP),PP_N,xGP,wBary,L_eta_RP(:,iRP))
#if PP_dim == 3
  CALL LagrangeInterpolationPolys(xi_RP(3,iRP),PP_N,xGP,wBary,L_zeta_RP(:,iRP))
#endif
END DO

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  L_zeta_RP = 0.
END IF
#endif
END SUBROUTINE InitRPBasis


!==================================================================================================================================
!> Evaluate solution at current time t at recordpoint positions and fill output buffer
!==================================================================================================================================
SUBROUTINE RecordPoints(iter,t,forceSampling)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars          ,ONLY: U
USE MOD_Timedisc_Vars,    ONLY: dt
USE MOD_Analyze_Vars,     ONLY: WriteData_dt,tWriteData
USE MOD_RecordPoints_Vars,ONLY: RP_Data,RP_ElemID
USE MOD_RecordPoints_Vars,ONLY: RP_Buffersize,RP_MaxBuffersize,RP_SamplingOffset,iSample
USE MOD_RecordPoints_Vars,ONLY: l_xi_RP,l_eta_RP,nRP
#if PP_dim==3
USE MOD_RecordPoints_Vars,ONLY: l_zeta_RP
#endif
#if FV_ENABLED
USE MOD_RecordPoints_Vars,ONLY: FV_RP_ijk
USE MOD_FV_Vars          ,ONLY: FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)     :: iter                    !< current number of timesteps
REAL,INTENT(IN)                :: t                       !< current time t
LOGICAL,INTENT(IN)             :: forceSampling           !< force sampling (e.g. at first/last timestep of computation)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k,iRP
REAL                    :: u_RP(PP_nVar,nRP)
REAL                    :: l_eta_zeta_RP
!----------------------------------------------------------------------------------------------------------------------------------
IF(MOD(iter,RP_SamplingOffset).NE.0 .AND. .NOT. forceSampling) RETURN
IF(.NOT.ALLOCATED(RP_Data))THEN
  ! Compute required buffersize from timestep and add 20% tolerance
  ! +1 is added to ensure a minimum buffersize of 2
  RP_Buffersize = MIN(CEILING((1.2*WriteData_dt)/(dt*RP_SamplingOffset))+1,RP_MaxBuffersize)
  ALLOCATE(RP_Data(0:PP_nVar,nRP,RP_Buffersize))
END IF

! evaluate state at RP
iSample=iSample+1
U_RP=0.
DO iRP=1,nRP
#if FV_ENABLED
  IF (FV_Elems(RP_ElemID(iRP)).EQ.0)THEN ! DG
#endif
    DO k=0,PP_NZ; DO j=0,PP_N
#if PP_dim==3
      l_eta_zeta_RP=l_eta_RP(j,iRP)*l_zeta_RP(k,iRP)
#else
      l_eta_zeta_RP=l_eta_RP(j,iRP)
#endif
      DO i=0,PP_N
        U_RP(:,iRP)=U_RP(:,iRP) + U(:,i,j,k,RP_ElemID(iRP))*l_xi_RP(i,iRP)*l_eta_zeta_RP
      END DO !i
    END DO; END DO !k
#if FV_ENABLED
  ELSE                                   ! FV
    ! RP value is cell average of nearest cell
    U_RP(:,iRP)=U(:,FV_RP_ijk(1,iRP),&
                    FV_RP_ijk(2,iRP),&
                    FV_RP_ijk(3,iRP),RP_ElemID(iRP))
  END IF
#endif
END DO ! iRP
RP_Data(1:PP_nVar,:,iSample)=U_RP
RP_Data(0,        :,iSample)=t

! dataset is full, write data and reset
IF(iSample.EQ.RP_Buffersize) CALL WriteRP(tWriteData,.FALSE.)
END SUBROUTINE RecordPoints


!==================================================================================================================================
!> Writes the time history of the solution at the recordpoints to an HDF5 file
!==================================================================================================================================
SUBROUTINE WriteRP(OutputTime,resetCounters)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE HDF5
USE MOD_IO_HDF5           ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Equation_Vars     ,ONLY: StrVarNames
USE MOD_HDF5_Output       ,ONLY: WriteAttribute,WriteArray,MarkWriteSuccessfull
USE MOD_Output_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars         ,ONLY: MeshFile
#if USE_MPI
USE MOD_Recordpoints_Vars ,ONLY: RP_COMM
#endif
USE MOD_Recordpoints_Vars ,ONLY: myRPrank,lastSample
USE MOD_Recordpoints_Vars ,ONLY: RPDefFile,RP_Data,iSample,nSamples
USE MOD_Recordpoints_Vars ,ONLY: offsetRP,nRP,nGlobalRP
USE MOD_Recordpoints_Vars ,ONLY: RP_Buffersize,RP_Maxbuffersize,RP_fileExists,chunkSamples
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,   INTENT(IN)             :: OutputTime            !< time
LOGICAL,INTENT(IN)             :: resetCounters         !< flag to reset sample counters and reallocate buffers, once file is done
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileString
REAL                           :: startT,endT
!==================================================================================================================================
IF(myRPrank.EQ.0)THEN
  WRITE(UNIT_stdOut,'(a)')' WRITE RECORDPOINT DATA TO HDF5 FILE...'
  WRITE(UNIT_stdOut,'(a,I4,a,I4,a)')' RP Buffer  : ',iSample,'/',RP_Buffersize,' samples.'
  GETTIME(startT)
END IF

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_RP',OutputTime))//'.h5'
! init file or just update time
IF(myRPrank.EQ.0)THEN
  CALL OpenDataFile(Filestring,create=.NOT.RP_fileExists,single=.TRUE.,readOnly=.FALSE.)
  IF(.NOT.RP_fileExists)THEN
    ! Create dataset attributes
    CALL WriteAttribute(File_ID,'File_Type'  ,1,StrScalar=(/TRIM('RecordPoints_Data')/))
    CALL WriteAttribute(File_ID,'MeshFile'   ,1,StrScalar=(/TRIM(MeshFile)/))
    CALL WriteAttribute(File_ID,'ProjectName',1,StrScalar=(/TRIM(ProjectName)/))
    CALL WriteAttribute(File_ID,'RPDefFile'  ,1,StrScalar=(/TRIM(RPDefFile)/))
    CALL WriteAttribute(File_ID,'VarNames'   ,PP_nVar,StrArray=StrVarNames)
    CALL WriteAttribute(File_ID,'Time'       ,1,RealScalar=OutputTime)
  END IF
  CALL CloseDataFile()
END IF

#if USE_MPI
CALL MPI_BARRIER(RP_COMM,iError)
CALL OpenDataFile(Filestring,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=RP_COMM)
#else
CALL OpenDataFile(Filestring,create=.FALSE.,single=.TRUE. ,readOnly=.FALSE.)
#endif

IF(iSample.GT.0)THEN
  IF(.NOT.RP_fileExists) chunkSamples=iSample
  ! write buffer into file, we need two offset dimensions (one buffer, one processor)
  nSamples=nSamples+iSample
  CALL WriteArray(DataSetName='RP_Data', rank=3,&
                        nValGlobal=(/PP_nVar+1,nGlobalRP,nSamples/),&
                        nVal=      (/PP_nVar+1,nRP      ,iSample/),&
                        offset=    (/0        ,offsetRP ,nSamples-iSample/),&
                        resizeDim= (/.FALSE.  ,.FALSE.  ,.TRUE./),&
                        chunkSize= (/PP_nVar+1,nGlobalRP,chunkSamples/),&
                        RealArray=RP_Data(:,:,1:iSample),&
                        collective=.TRUE.)
  lastSample=RP_Data(:,:,iSample)
END IF
! Reset buffer
RP_Data=0.

iSample=0
RP_fileExists=.TRUE.
IF(resetCounters)THEN
  ! Recompute required buffersize from timestep and add 10% tolerance
  IF((nSamples.GE.RP_Buffersize).AND.(RP_Buffersize.LT.RP_Maxbuffersize))THEN
    RP_Buffersize=MIN(CEILING(1.1*nSamples)+1,RP_MaxBuffersize)
    DEALLOCATE(RP_Data)
    ALLOCATE(RP_Data(0:PP_nVar,nRP,RP_Buffersize))
  END IF
  RP_fileExists=.FALSE.
  iSample=1
  nSamples=0
  ! last sample of previous file = first sample of next file
  RP_Data(:,:,1)=lastSample
END IF
CALL CloseDataFile()

IF(myRPrank.EQ.0)THEN
  CALL MarkWriteSuccessfull(Filestring)
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' DONE  [',EndT-StartT,'s]'
END IF
END SUBROUTINE WriteRP


!==================================================================================================================================
!> Deallocate recordpoint arrays
!==================================================================================================================================
SUBROUTINE FinalizeRecordPoints()
! MODULES
USE MOD_RecordPoints_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(RP_Data)
SDEALLOCATE(RP_ElemID)
SDEALLOCATE(L_xi_RP)
SDEALLOCATE(L_eta_RP)
SDEALLOCATE(L_zeta_RP)
SDEALLOCATE(lastSample)
#if FV_ENABLED
SDEALLOCATE(FV_RP_ijk)
#endif
RecordPointsInitIsDone = .FALSE.
END SUBROUTINE FinalizeRecordPoints


END MODULE MOD_RecordPoints
