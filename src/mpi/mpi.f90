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
!> \brief Contains the routines that set up communicators and control non-blocking communication
!==================================================================================================================================
MODULE MOD_MPI
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitMPI
  MODULE PROCEDURE InitMPI
END INTERFACE

PUBLIC::InitMPI

#if USE_MPI
INTERFACE InitMPIvars
  MODULE PROCEDURE InitMPIvars
END INTERFACE

!INTERFACE StartReceiveMPIData
!  MODULE PROCEDURE StartReceiveMPIData
!END INTERFACE
!
!INTERFACE StartSendMPIData
!  MODULE PROCEDURE StartSendMPIData
!END INTERFACE

!INTERFACE FinishExchangeMPIData
!  MODULE PROCEDURE FinishExchangeMPIData
!END INTERFACE

#if FV_ENABLED
INTERFACE StartExchange_FV_Elems
  MODULE PROCEDURE StartExchange_FV_Elems
END INTERFACE
#endif

INTERFACE FinalizeMPI
  MODULE PROCEDURE FinalizeMPI
END INTERFACE

PUBLIC::InitMPIvars
PUBLIC::StartReceiveMPIData
PUBLIC::StartSendMPIData
#if FV_ENABLED
PUBLIC::StartExchange_FV_Elems
#endif
PUBLIC::FinishExchangeMPIData
PUBLIC::FinalizeMPI
#endif
!==================================================================================================================================

PUBLIC::DefineParametersMPI
CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersMPI()
! MODULES
USE MOD_ReadInTools,              ONLY: prms
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
!> Basic mpi initialization. Calls initialization routine of the mpi library and sets myRank, nProcessors and MPIRoot. If the code
!> is not compiled with mpi, InitMPI sets standard values for these variables.
!==================================================================================================================================
SUBROUTINE InitMPI(mpi_comm_IN)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL      :: mpi_comm_IN !< MPI communicator
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
LOGICAL :: initDone
!==================================================================================================================================
IF (PRESENT(mpi_comm_IN)) THEN
  MPI_COMM_FLEXI = mpi_comm_IN
ELSE
  CALL MPI_INIT(iError)
  CALL MPI_INITIALIZED(initDone,iError)
  IF(.NOT.initDone) CALL MPI_INIT(iError)
  IF(iError .NE. 0) &
    CALL Abort(__STAMP__,'Error in MPI_INIT',iError)
  MPI_COMM_FLEXI = MPI_COMM_WORLD
END IF

CALL MPI_COMM_RANK(MPI_COMM_FLEXI, myRank     , iError)
CALL MPI_COMM_SIZE(MPI_COMM_FLEXI, nProcessors, iError)
IF(iError .NE. 0) &
  CALL Abort(__STAMP__,'Could not get rank and number of processors',iError)
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
!> Initialize derived mpi variables used for communication
!==================================================================================================================================
SUBROUTINE InitMPIVars()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_ReadinTools,             ONLY: GETINT
USE MOD_Interpolation_Vars,      ONLY: InterpolationInitIsDone
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
DataSizeSide      =PP_nVar*(PP_N+1)*(PP_NZ+1)
DataSizeSidePrim  =PP_nVarPrim*(PP_N+1)*(PP_NZ+1)

! split communicator into smaller groups (e.g. for local nodes)
GroupSize=GETINT('GroupSize','0')
IF(GroupSize.LT.1)THEN ! group procs by node
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,myRank,myRank,MPI_COMM_NODE,iError)
ELSE ! use groupsize
  color=myRank/GroupSize
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,color,myRank,MPI_COMM_NODE,iError)
END IF
CALL MPI_COMM_RANK(MPI_COMM_NODE,myLocalRank,iError)
CALL MPI_COMM_SIZE(MPI_COMM_NODE,nLocalProcs,iError)
MPILocalRoot=(myLocalRank .EQ. 0)

! now split global communicator into small group leaders and the others
MPI_COMM_LEADERS=MPI_COMM_NULL
MPI_COMM_WORKERS=MPI_COMM_NULL
myLeaderRank=-1
myWorkerRank=-1
IF(myLocalRank.EQ.0)THEN
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,0,myRank,MPI_COMM_LEADERS,iError)
  CALL MPI_COMM_RANK( MPI_COMM_LEADERS,myLeaderRank,iError)
  CALL MPI_COMM_SIZE( MPI_COMM_LEADERS,nLeaderProcs,iError)
  nWorkerProcs=nProcessors-nLeaderProcs
ELSE
  CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI,1,myRank,MPI_COMM_WORKERS,iError)
  CALL MPI_COMM_RANK( MPI_COMM_WORKERS,myWorkerRank,iError)
  CALL MPI_COMM_SIZE( MPI_COMM_WORKERS,nWorkerProcs,iError)
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
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: SendID                                   !< defines the send / receive direction -> 1=send MINE
                                                                        !< / receive YOUR, 2=send YOUR / receive MINE
INTEGER,INTENT(IN)          :: DataSize                                 !< size of one entry in array (e.g. one side:
                                                                        !< nVar*(N+1)**2
INTEGER,INTENT(IN)          :: LowerBound                               !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)          :: UpperBound                               !< upper side index for last dimension of FaceData
INTEGER,INTENT(OUT)         :: MPIRequest(nNbProcs)                     !< communication handles
REAL,INTENT(OUT)            :: FaceData(DataSize,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
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
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: SendID                                   !< defines the send / receive direction -> 1=send MINE
                                                                        !< / receive YOUR, 2=send YOUR / receive MINE
INTEGER,INTENT(IN)          :: DataSize                                 !< size of one entry in array (e.g. one side:
                                                                        !< nVar*(N+1)*(N+1))
INTEGER,INTENT(IN)          :: LowerBound                               !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)          :: UpperBound                               !< upper side index for last dimension of FaceData
INTEGER,INTENT(OUT)         :: MPIRequest(nNbProcs)                     !< communication handles
REAL,INTENT(IN)             :: FaceData(DataSize,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
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
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
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
INTEGER,INTENT(IN)    :: SendID                          !< defines the send / receive direction -> 1=send MINE/receive YOUR,
                                                         !< 2=send YOUR / receive MINE
INTEGER,INTENT(IN)    :: LowerBound                      !< lower side index for last dimension of FV_Elems
INTEGER,INTENT(IN)    :: UpperBound                      !< upper side index for last dimension of FV_Elems
INTEGER,INTENT(OUT)   :: SendRequest(nNbProcs)           !< communicatio handles for send
INTEGER,INTENT(OUT)   :: RecRequest(nNbProcs)            !< communicatio handles for receive
INTEGER,INTENT(INOUT) :: FV_Elems(LowerBound:UpperBound) !< information about FV_Elems at faces to be communicated
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
  ELSE
    RecRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartExchange_FV_Elems
#endif



!==================================================================================================================================
!> We have to complete our non-blocking communication operations before we can (re)use the send / receive buffers
!==================================================================================================================================
SUBROUTINE FinishExchangeMPIData(nRequests,MPIRequest)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nRequests             !< size of the handles
INTEGER,INTENT(INOUT)       :: MPIRequest(nRequests) !< communication handles
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if LIBS_MPT
CALL MPI_WaitAll(nRequests,MPIRequest,MPI_STATUS_IGNORE,iError)
#else
CALL MPI_WaitAll(nRequests,MPIRequest,MPI_STATUSES_IGNORE,iError)
#endif
END SUBROUTINE FinishExchangeMPIData

!==================================================================================================================================
!> Deallocate MPI arrays
!==================================================================================================================================
SUBROUTINE FinalizeMPI()
! MODULES
USE MOD_MPI_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(MPIRequest_U)
SDEALLOCATE(MPIRequest_Flux)
#if FV_ENABLED
SDEALLOCATE(MPIRequest_FV_Elems)
SDEALLOCATE(MPIRequest_FV_gradU)
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
END SUBROUTINE FinalizeMPI

#endif /*USE_MPI*/

END MODULE MOD_MPI
