!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
#include "eos.h"

!==================================================================================================================================
!> \brief Contains the routines that set up communicators and control non-blocking communication
!==================================================================================================================================
MODULE MOD_MPI
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: DefineParametersMPI
PUBLIC:: InitMPI
#if USE_MPI
PUBLIC:: InitMPIvars
PUBLIC:: StartReceiveMPIData
PUBLIC:: StartSendMPIData
#if FV_ENABLED
PUBLIC:: StartExchange_FV_Elems
#endif /*FV_ENABLED*/
#if FV_ENABLED == 2
PUBLIC:: StartExchange_FV_alpha
#endif /*FV_ENABLED == 2*/
PUBLIC:: FinishExchangeMPIData
PUBLIC:: FinalizeMPI
#endif /*USE_MPI*/
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersMPI()
! MODULES
USE MOD_ReadInTools,              ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("MPI")
CALL prms%CreateIntOption('GroupSize', "Define size of MPI subgroups, used to e.g. perform grouped IO, where group master\n"//&
                                       "collects and outputs data.",&
                                       '0')
END SUBROUTINE DefineParametersMPI


!==================================================================================================================================
!> Basic MPI initialization. Calls initialization routine of the MPI library and sets myRank, nProcessors and MPIRoot. If the code
!> is not compiled with MPI, InitMPI sets standard values for these variables.
!==================================================================================================================================
SUBROUTINE InitMPI(            &
#if USE_MPI
                   mpi_comm_IN &
#endif /*USE_MPI*/
                  )
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
#if USE_MPI
TYPE(MPI_Comm),INTENT(IN),OPTIONAL      :: mpi_comm_IN !< MPI communicator
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: initDone,foundAttr
INTEGER :: color
INTEGER(KIND=MPI_ADDRESS_KIND) :: myApp
!==================================================================================================================================
IF (PRESENT(mpi_comm_IN)) THEN
  MPI_COMM_FLEXI = mpi_comm_IN
ELSE
  CALL MPI_INIT(iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_INIT',iError)
  CALL MPI_INITIALIZED(initDone,iError)
  IF(.NOT.initDone) CALL MPI_INIT(iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_INITIALIZED',iError)

  ! Get number of own app if multiple apps have been launched in single mpirun command
  CALL MPI_COMM_GET_ATTR(MPI_COMM_WORLD,MPI_APPNUM,myApp,foundAttr,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_GET_ATTR',iError)
  IF (foundAttr) THEN
    ! Split communicator to obtain own MPI_COMM_FLEXI per executable (explicit cast, since API requires INT().)
    color = MAX(INT(myApp), 0)
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,0,MPI_COMM_FLEXI,iError)
  ELSE
    ! Duplicate communicator instead of just copying it. Creates a clean copy with all the cached information intact
    CALL MPI_COMM_DUP(MPI_COMM_WORLD,MPI_COMM_FLEXI,iError)
  END IF
END IF

CALL MPI_COMM_RANK(MPI_COMM_FLEXI, myRank     , iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_RANK',iError)
CALL MPI_COMM_SIZE(MPI_COMM_FLEXI, nProcessors, iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SIZE',iError)
MPIRoot=(myRank .EQ. 0)
#else  /*USE_MPI*/
myRank      = 0
myLocalRank = 0
nProcessors = 1
MPIRoot     =.TRUE.
MPILocalRoot=.TRUE.
#endif  /*USE_MPI*/

! At this point the initialization is not completed. We first have to create a new MPI communicator.
END SUBROUTINE InitMPI


#if USE_MPI
!==================================================================================================================================
!> Initialize derived MPI variables used for communication
!==================================================================================================================================
SUBROUTINE InitMPIVars()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation_Vars,      ONLY: InterpolationInitIsDone
USE MOD_MPI_Vars
USE MOD_ReadinTools,             ONLY: GETINT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: color,groupsize
!==================================================================================================================================
IF(.NOT.InterpolationInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,'InitMPITypes called before InitInterpolation')
END IF

ALLOCATE(MPIRequest_U(nNbProcs,2)    )
ALLOCATE(MPIRequest_Flux(nNbProcs,2) )
MPIRequest_U      = MPI_REQUEST_NULL
MPIRequest_Flux   = MPI_REQUEST_NULL
#if FV_ENABLED
ALLOCATE(MPIRequest_FV_Elems(nNbProcs,2) )
ALLOCATE(MPIRequest_FV_gradU(nNbProcs,2) )
MPIRequest_FV_Elems = MPI_REQUEST_NULL
MPIRequest_FV_gradU = MPI_REQUEST_NULL
#if ((FV_ENABLED == 2) && (PP_NodeType == 1))
ALLOCATE(MPIRequest_FV_U(nNbProcs,2)    )
ALLOCATE(MPIRequest_FV_Flux(nNbProcs,2) )
MPIRequest_FV_U     = MPI_REQUEST_NULL
MPIRequest_FV_Flux  = MPI_REQUEST_NULL
#endif
#if FV_RECONSTRUCT
ALLOCATE(MPIRequest_Rec_MS(nNbProcs,2))
ALLOCATE(MPIRequest_Rec_SM(nNbProcs,2))
MPIRequest_Rec_MS = MPI_REQUEST_NULL
MPIRequest_Rec_SM = MPI_REQUEST_NULL
#endif
#endif
#if EDDYVISCOSITY
ALLOCATE(MPIRequest_SGS(nNbProcs,2) )
MPIRequest_SGS     = MPI_REQUEST_NULL
#endif

#if PARABOLIC
ALLOCATE(MPIRequest_gradU(nNbProcs,3,2))
MPIRequest_gradU = MPI_REQUEST_NULL
#endif /*PARABOLIC*/

#if EDDYVISCOSITY
DataSizeSideSGS= (PP_N+1)*(PP_NZ+1)
#endif
DataSizeSide             =PP_nVar*(PP_N+1)*(PP_NZ+1)
DataSizeSidePrim         =PP_nVarPrim*(PP_N+1)*(PP_NZ+1)
DataSizeSideGrad         =PP_nVarLifting*(PP_N+1)*(PP_NZ+1)
DataSizeSideGradParabolic=PP_nVarLifting*(PP_N+1)*(PP_NZ+1)*3

! split communicator into smaller groups (e.g. for local nodes)
GroupSize=GETINT('GroupSize')
IF(GroupSize.LT.1)THEN ! group procs by node
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,myRank,0,MPI_COMM_NODE,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SPLIT',iError)
ELSE ! use groupsize
  color=myRank/GroupSize
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,color,0,MPI_COMM_NODE,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SPLIT',iError)
END IF
CALL MPI_COMM_RANK(MPI_COMM_NODE,myLocalRank,iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_RANK',iError)
CALL MPI_COMM_SIZE(MPI_COMM_NODE,nLocalProcs,iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SIZE',iError)
MPILocalRoot=(myLocalRank .EQ. 0)

IF (nProcessors.EQ.nLocalProcs) THEN
  SWRITE(UNIT_stdOUt,'(A,I0,A,I0,A)') ' | Starting gathered I/O communication with ',nLocalProcs,' procs in ',1,' group'
ELSE
  SWRITE(UNIT_stdOUt,'(A,I0,A,I0,A,I0,A)') ' | Starting gathered I/O communication with ',nLocalProcs,' procs each in ',&
                                                        nProcessors/nLocalProcs,' groups for a total number of ',&
                                                        nProcessors,' procs'
END IF

! now split global communicator into small group leaders and the others
MPI_COMM_LEADERS=MPI_COMM_NULL
MPI_COMM_WORKERS=MPI_COMM_NULL
myLeaderRank=-1
myWorkerRank=-1
IF(myLocalRank.EQ.0)THEN
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,0,0,MPI_COMM_LEADERS,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SPLIT',iError)
  CALL MPI_COMM_RANK( MPI_COMM_LEADERS,myLeaderRank,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_RANK',iError)
  CALL MPI_COMM_SIZE( MPI_COMM_LEADERS,nLeaderProcs,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SIZE',iError)
  nWorkerProcs=nProcessors-nLeaderProcs
ELSE
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,1,0,MPI_COMM_WORKERS,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SPLIT',iError)
  CALL MPI_COMM_RANK( MPI_COMM_WORKERS,myWorkerRank,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_RANK',iError)
  CALL MPI_COMM_SIZE( MPI_COMM_WORKERS,nWorkerProcs,iError)
  IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_SIZE',iError)
  nLeaderProcs=nProcessors-nWorkerProcs
END IF

END SUBROUTINE InitMPIvars


!==================================================================================================================================
!> Subroutine that controls the receive operations for the face data that has to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE StartReceiveMPIData(FaceData,DataSize,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: SendID                                   !< defines the send / receive direction -> 1=send MINE
                                                                          !< / receive YOUR, 2=send YOUR / receive MINE
INTEGER,INTENT(IN)            :: DataSize                                 !< size of one entry in array (e.g. one side:
                                                                          !< nVar*(N+1)**2
INTEGER,INTENT(IN)            :: LowerBound                               !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)            :: UpperBound                               !< upper side index for last dimension of FaceData
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs)                     !< communication handles
REAL,INTENT(OUT)              :: FaceData(DataSize,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iNBProc
!==================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =DataSize*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FaceData(:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_FLEXI,MPIRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPIData



!==================================================================================================================================
!> Subroutine that performs the send operations for the face data that has to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE StartSendMPIData(FaceData,DataSize,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: SendID                                   !< defines the send / receive direction -> 1=send MINE
                                                                          !< / receive YOUR, 2=send YOUR / receive MINE
INTEGER,INTENT(IN)            :: DataSize                                 !< size of one entry in array (e.g. one side:
                                                                          !< nVar*(N+1)*(N+1))
INTEGER,INTENT(IN)            :: LowerBound                               !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)            :: UpperBound                               !< upper side index for last dimension of FaceData
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs)                     !< communication handles
REAL,INTENT(IN)               :: FaceData(DataSize,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iNBProc
!==================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =DataSize*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FaceData(:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_FLEXI,MPIRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO ! iProc=1,nNBProcs

END SUBROUTINE StartSendMPIData


#if FV_ENABLED
!==================================================================================================================================
!> Subroutine that performs the send and receive operations for the FV_elems information at the face
!> that has to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE StartExchange_FV_Elems(FV_Elems,LowerBound,UpperBound,SendRequest,RecRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: SendID                          !< defines the send / receive direction -> 1=send MINE/receive YOUR,
                                                                 !< 2=send YOUR / receive MINE
INTEGER,INTENT(IN)            :: LowerBound                      !< lower side index for last dimension of FV_Elems
INTEGER,INTENT(IN)            :: UpperBound                      !< upper side index for last dimension of FV_Elems
TYPE(MPI_Request),INTENT(OUT) :: SendRequest(nNbProcs)           !< communicatio handles for send
TYPE(MPI_Request),INTENT(OUT) :: RecRequest(nNbProcs)            !< communicatio handles for receive
INTEGER,INTENT(INOUT)         :: FV_Elems(LowerBound:UpperBound) !< information about FV_Elems at faces to be communicated
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iNBProc
!==================================================================================================================================
DO iNbProc=1,nNbProcs
  ! Start send face data
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FV_Elems(SideID_start:SideID_end),nSendVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_FLEXI,SendRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    SendRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
  ! Start receive face data
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FV_Elems(SideID_start:SideID_end),nRecVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_FLEXI,RecRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    RecRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO ! iProc=1,nNBProcs

END SUBROUTINE StartExchange_FV_Elems
#endif


#if FV_ENABLED == 2
!==================================================================================================================================
!> Subroutine that performs the send and receive operations for the FV_elems information at the face
!> that has to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE StartExchange_FV_alpha(FV_alpha,LowerBound,UpperBound,SendRequest,RecRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: SendID                          !< defines the send / receive direction -> 1=send MINE/receive YOUR,
                                                                   !< 2=send YOUR / receive MINE
INTEGER,INTENT(IN)              :: LowerBound                      !< lower side index for last dimension of FV_Elems
INTEGER,INTENT(IN)              :: UpperBound                      !< upper side index for last dimension of FV_Elems
TYPE(MPI_Request),INTENT(OUT)   :: SendRequest(nNbProcs)           !< communicatio handles for send
TYPE(MPI_Request),INTENT(OUT)   :: RecRequest(nNbProcs)            !< communicatio handles for receive
REAL,INTENT(INOUT)              :: FV_alpha(LowerBound:UpperBound) !< information about FV_Elems at faces to be communicated
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iNBProc
!==================================================================================================================================
DO iNbProc=1,nNbProcs
  ! Start send face data
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FV_alpha(SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_FLEXI,SendRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    SendRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
  ! Start receive face data
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FV_alpha(SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_FLEXI,RecRequest(iNbProc),iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    RecRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO ! iProc=1,nNBProcs

END SUBROUTINE StartExchange_FV_alpha
#endif /*FV_ENABLED == 2*/


!==================================================================================================================================
!> We have to complete our non-blocking communication operations before we can (re)use the send / receive buffers
!==================================================================================================================================
SUBROUTINE FinishExchangeMPIData(nRequests,MPIRequest)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER          ,INTENT(IN)    :: nRequests             !< size of the handles
TYPE(MPI_Request),INTENT(INOUT) :: MPIRequest(nRequests) !< communication handles
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL MPI_WaitAll(nRequests,MPIRequest,MPI_STATUSES_IGNORE,iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_WaitAll',iError)

END SUBROUTINE FinishExchangeMPIData


!==================================================================================================================================
!> Deallocate MPI arrays
!==================================================================================================================================
SUBROUTINE FinalizeMPI()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(MPIRequest_U)
SDEALLOCATE(MPIRequest_Flux)
#if FV_ENABLED
SDEALLOCATE(MPIRequest_FV_Elems)
SDEALLOCATE(MPIRequest_FV_gradU)
#if ((FV_ENABLED == 2) && (PP_NodeType == 1))
SDEALLOCATE(MPIRequest_FV_U)
SDEALLOCATE(MPIRequest_FV_Flux)
#endif
#if FV_RECONSTRUCT
SDEALLOCATE(MPIRequest_Rec_MS)
SDEALLOCATE(MPIRequest_Rec_SM)
#endif
#endif
#if EDDYVISCOSITY
SDEALLOCATE(MPIRequest_SGS)
#endif
#if PARABOLIC
SDEALLOCATE(MPIRequest_gradU)
#endif /*PARABOLIC*/
SDEALLOCATE(NbProc)
SDEALLOCATE(nMPISides_Proc)
SDEALLOCATE(nMPISides_MINE_Proc)
SDEALLOCATE(nMPISides_YOUR_Proc)
SDEALLOCATE(offsetMPISides_MINE)
SDEALLOCATE(offsetMPISides_YOUR)
SDEALLOCATE(offsetElemMPI)
SDEALLOCATE(nMPISides_send)
SDEALLOCATE(nMPISides_rec)
SDEALLOCATE(OffsetMPISides_send)
SDEALLOCATE(OffsetMPISides_rec)

! Free MPI communicators
IF(MPI_COMM_WORKERS.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_WORKERS,iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_FREE',iError)
IF(MPI_COMM_LEADERS.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_LEADERS,iError)
IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_COMM_FREE',iError)

END SUBROUTINE FinalizeMPI
#endif /*USE_MPI*/

END MODULE MOD_MPI
