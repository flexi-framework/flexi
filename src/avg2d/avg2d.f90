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

!==================================================================================================================================
!>
!==================================================================================================================================
MODULE MOD_Avg2D
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE DefineParametersAvg2D
  MODULE PROCEDURE DefineParametersAvg2D
END INTERFACE

INTERFACE InitAvg2D
  MODULE PROCEDURE InitAvg2D
END INTERFACE

INTERFACE Avg2DSlave
  MODULE PROCEDURE Avg2DSlave
END INTERFACE

INTERFACE Avg2DMaster
  MODULE PROCEDURE Avg2DMaster
END INTERFACE

INTERFACE FinalizeAvg2D
  MODULE PROCEDURE FinalizeAvg2D
END INTERFACE

#if USE_MPI
INTERFACE SendAvg2DToMaster
  MODULE PROCEDURE SendAvg2DToMaster
END INTERFACE

INTERFACE RecvAvg2DOfMaster
  MODULE PROCEDURE RecvAvg2DOfMaster
END INTERFACE

INTERFACE SendAvg2DToSlave
  MODULE PROCEDURE SendAvg2DToSlave
END INTERFACE
#endif /* USE_MPI */

INTERFACE RecvAvg2DOfSlave
  MODULE PROCEDURE RecvAvg2DOfSlave
END INTERFACE

PUBLIC :: DefineParametersAvg2D
PUBLIC :: InitAvg2D
PUBLIC :: Avg2DSlave
PUBLIC :: Avg2DMaster
PUBLIC :: FinalizeAvg2D
#if USE_MPI
PUBLIC :: SendAvg2DToMaster
PUBLIC :: RecvAvg2DOfMaster
PUBLIC :: SendAvg2DToSlave
#endif /* USE_MPI */
PUBLIC :: RecvAvg2DOfSlave
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters needed for filtering
!==================================================================================================================================
SUBROUTINE DefineParametersAvg2D()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Avg2D")
CALL prms%CreateLogicalOption(   'doAvg2D',      "Turn-On to average the baseflow in spanwise direction")
CALL prms%CreateIntOption(       'Avg2DDir',     "Direction in which the averaging is performed","3")

END SUBROUTINE DefineParametersAvg2D

!==================================================================================================================================
!> Initialize all necessary information to perform averaging in z Direction
!==================================================================================================================================
SUBROUTINE InitAvg2D()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Avg2D_Vars
USE MOD_ReadInTools        ,ONLY:GETLOGICAL,GETINT
USE MOD_Mesh_Readin        ,ONLY:ReadIJKSorting
USE MOD_Mesh_Vars          ,ONLY:nElems,nElems_IJK,Elem_IJK,Elem_xGP
USE MOD_Interpolation_Vars ,ONLY:xGP,wGP
#if USE_MPI
USE MOD_Mesh_Vars          ,ONLY:nGlobalElems
USE MOD_MPI_Vars           ,ONLY:MPIRequest_Avg2DSend,MPIRequest_Avg2DRecv
USE MOD_MPI_Vars           ,ONLY:offsetElemMPI
#endif /* USE_MPI */
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iRank,i,j,k
#if USE_MPI
INTEGER             :: iElemToSend,nRecvVal,nSendVal
INTEGER             :: ElemRank
#endif /* USE_MPI */
INTEGER             :: iElemToRcv
INTEGER             :: recvProc,sendProc
INTEGER             :: iDiffRecvProc,maxIJ(2),nElemsloc_IJ(2)
INTEGER             :: iDiffSendProc,iElemK0,iProc,nMyElemsK0,nDiffElemsloc
INTEGER,ALLOCATABLE :: IJK_globalProc(:,:,:)
LOGICAL,ALLOCATABLE :: SendProcsList(:),RecvProcList(:)
INTEGER,ALLOCATABLE :: myElemsK0_temp(:),myElemsK0(:)
LOGICAL,ALLOCATABLE :: ijDone(:,:,:)
INTEGER,ALLOCATABLE :: posIProc_tmp(:),ElemsToSend_tmp(:,:,:,:),ElemsToRcv_tmp(:,:,:,:),posLocTmp(:,:)
REAL                :: dx
REAL,ALLOCATABLE    :: domainWidth(:)
REAL,ALLOCATABLE    :: SendBuffer(:,:),RcvBuffer(:,:)
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT AVG2D...'

doAvg2D = GETLOGICAL('doAvg2D','.FALSE.')
IF (.NOT. doAvg2D) THEN
  SWRITE(UNIT_stdOut,'(A)') ' | Baseflow is not averaged in spatial direction ...'
  SWRITE(UNIT_stdOut,'(A)') ' INIT AVG2D DONE...'
  SWRITE(UNIT_stdOut,'(132("-"))')
  RETURN
END IF

! Select averaging direction
Avg2DDir = GETINT('Avg2DDir')
IJK_Mask = PACK((/1,2,3/),(/1,2,3/).NE.Avg2DDir)

#if USE_MPI
! avoid MPI_ALLREDUCE due to excessive memory usage
CALL ReadIJKSorting(doGlobal=.TRUE.)
ALLOCATE(IJK_globalProc(nElems_IJK(1),nElems_IJK(2),nElems_IJK(3)))
IJK_globalProc=0
DO iElem=1,nGlobalElems
  DO iRank=0,nProcessors-1
    IF ((iElem.GT.offsetElemMPI(iRank)).AND.(iElem.LE.offsetElemMPI(iRank+1))) THEN
      ElemRank = iRank
      EXIT
    END IF
  END DO
  i=Elem_IJK(1,iElem)
  j=Elem_IJK(2,iElem)
  k=Elem_IJK(3,iElem)
  IJK_globalProc(i,j,k)=ElemRank
END DO
DEALLOCATE(Elem_IJK)
CALL ReadIJKSorting(doGlobal=.FALSE.)
#else
CALL ReadIJKSorting(doGlobal=.FALSE.)
ALLOCATE(IJK_globalProc(nElems_IJK(1),nElems_IJK(2),nElems_IJK(3)))
IJK_globalProc=0
#endif

!GET min/max i and j of elems of this proc
minIJ(1) = MINVAL(Elem_IJK(IJK_Mask(1),:))-1
minIJ(2) = MINVAL(Elem_IJK(IJK_Mask(2),:))-1
maxIJ(1) = MAXVAL(Elem_IJK(IJK_Mask(1),:))
maxIJ(2) = MAXVAL(Elem_IJK(IJK_Mask(2),:))
!Square of possible different Elements
nElemsloc_IJ(1) = maxIJ(1)-minIJ(1)
nElemsloc_IJ(2) = maxIJ(2)-minIJ(2)

! Count how many elements need to be received from different procs
ALLOCATE(myElemsK0_temp(nElems))
nMyElemsK0 = 0
DO iElem=1,nElems
  IF(Elem_IJK(Avg2DDir,iElem).EQ.1) THEN
    nMyElemsK0 = nMyElemsK0 + 1
    myElemsK0_temp(nMyElemsK0) = iElem
  END iF
END DO
IF(nMyElemsK0.GT.0) THEN
  ALLOCATE(myElemsK0(nMyElemsK0))
  myElemsK0=myElemsK0_temp(1:nMyElemsK0)
  SDEALLOCATE(myElemsK0_temp)
END IF

ALLOCATE(ElemsToRcv_tmp(2,nElemsloc_IJ(1),nElemsloc_IJ(2),0:nProcessors-1))
ElemsToRcv_tmp = -999
ALLOCATE(SendProcsList(0:nProcessors-1))
ALLOCATE(nElemsToRecv(0:nProcessors-1))
ALLOCATE(ijDone(nElems_IJK(IJK_Mask(1)),nElems_IJK(IJK_Mask(2)),0:nProcessors-1))
ijDone=.FALSE.
nElemsToRecv=0
SendProcsList=.FALSE.
DO iElemK0 = 1,nMyElemsK0
  iElem = myElemsK0(iElemK0)
  i=Elem_IJK(IJK_Mask(1),iElem)
  j=Elem_IJK(IJK_Mask(2),iElem)
  DO k=2,nElems_IJK(Avg2DDir)
    sendProc=IJK_globalProc(i,j,k)
    IF (sendProc .NE. myRank) THEN
      SendProcsList(sendProc)=.TRUE.
      IF (.NOT. ijDone(i,j,sendProc)) THEN
        ijDone(i,j,sendProc)=.TRUE.
        nElemsToRecv(sendProc)=nElemsToRecv(sendProc)+1
        ElemsToRcv_tmp(1,Elem_IJK(IJK_Mask(1),iElem)-minIJ(1),Elem_IJK(IJK_Mask(2),iElem)-minIJ(2),sendProc)=iElem
        ElemsToRcv_tmp(2,Elem_IJK(IJK_Mask(1),iElem)-minIJ(1),Elem_IJK(IJK_Mask(2),iElem)-minIJ(2),sendProc)=sendProc
      END IF
    END IF
  END DO
END DO
nSendProcs=COUNT(SendProcsList)
SDEALLOCATE(myElemsK0)
ALLOCATE(SendProcs(nSendProcs))
iDiffSendProc=1
DO iProc=0,nProcessors-1
  IF(SendProcsList(iProc)) THEN
    SendProcs(iDiffSendProc) = iProc
    iDiffSendProc=iDiffSendProc+1
  END IF
END DO

ALLOCATE(ElemsToRecvIJSorted(nSendProcs,MAXVAL(nElemsToRecv)))
ElemsToRecvIJSorted=0
ALLOCATE(posIProc_tmp(1:nSendProcs))
posIProc_tmp=0.
DO j = 1,nElemsloc_IJ(2)
  DO i = 1,nElemsloc_IJ(1)
    DO iDiffSendProc = 1,nSendProcs
      iProc=SendProcs(iDiffSendProc)
      IF (ElemsToRcv_tmp(2,i,j,iProc) .NE. -999) THEN
        IF (SendProcs(iDiffSendProc) .EQ. ElemsToRcv_tmp(2,i,j,iProc)) THEN
          posIProc_tmp(iDiffSendProc) = posIProc_tmp(iDiffSendProc)+1
          ElemsToRecvIJSorted(iDiffSendProc,posIProc_tmp(iDiffSendProc)) = ElemsToRcv_tmp(1,i,j,iProc)
          !EXIT
        END IF
      END IF
    END DO
  END DO
END DO
SDEALLOCATE(posIProc_tmp)
SDEALLOCATE(SendProcsList)
SDEALLOCATE(ElemsToRcv_tmp)

! Count how many elements need to be send to different procs
ALLOCATE(ElemsToSend_tmp(2,nElemsloc_IJ(1),nElemsloc_IJ(2),0:nProcessors-1))
ElemsToSend_tmp = -999
ALLOCATE(nElemsToSend(0:nProcessors-1))
nElemsToSend=0
ALLOCATE(RecvProcList(0:nProcessors-1))
RecvProcList=.FALSE.
ijDone=.FALSE.
DO iElem=1,nElems
  recvProc = IJK_globalProc(Elem_IJK(IJK_Mask(1),iElem),Elem_IJK(IJK_Mask(2),iElem),1)
  IF(Elem_IJK(Avg2DDir,iElem).NE.1 .AND. recvProc.NE.myRank) THEN
    RecvProcList(recvProc)=.TRUE.
    IF (.NOT. ijDone(Elem_IJK(IJK_Mask(1),iElem),Elem_IJK(IJK_Mask(2),iElem),recvProc)) THEN
      ijDone(Elem_IJK(IJK_Mask(1),iElem),Elem_IJK(IJK_Mask(2),iElem),recvProc)=.TRUE.
      nElemsToSend(recvProc) = nElemsToSend(recvProc)+1
      ElemsToSend_tmp(1,Elem_IJK(IJK_Mask(1),iElem)-minIJ(1),Elem_IJK(IJK_Mask(2),iElem)-minIJ(2),recvProc)=iElem
      ElemsToSend_tmp(2,Elem_IJK(IJK_Mask(1),iElem)-minIJ(1),Elem_IJK(IJK_Mask(2),iElem)-minIJ(2),recvProc) =recvProc
    END IF
  END iF
END DO
SDEALLOCATE(ijDone)
nRecvProcs=COUNT(RecvProcList)
ALLOCATE(RecvProcs(nRecvProcs))
iDiffRecvProc=1
DO iProc=0,nProcessors-1
  IF(RecvProcList(iProc)) THEN
    RecvProcs(iDiffRecvProc) = iProc
    iDiffRecvProc=iDiffRecvProc+1
  END IF
END DO

ALLOCATE(ElemsToSendIJSorted(nRecvProcs,MAXVAL(nElemsToSend)))
ElemsToSendIJSorted=0
ALLOCATE(posIProc_tmp(1:nRecvProcs))
posIProc_tmp=0.
DO j = 1,nElemsloc_IJ(2)
  DO i = 1,nElemsloc_IJ(1)
    DO iDiffRecvProc = 1,nRecvProcs
      iProc=RecvProcs(iDiffRecvProc)
      IF (ElemsToSend_tmp(2,i,j,iProc) .NE. -999) THEN
        IF (RecvProcs(iDiffRecvProc) .EQ. ElemsToSend_tmp(2,i,j,iProc)) THEN
          posIProc_tmp(iDiffRecvProc) = posIProc_tmp(iDiffRecvProc)+1
          ElemsToSendIJSorted(iDiffRecvProc,posIProc_tmp(iDiffRecvProc)) = ElemsToSend_tmp(1,i,j,iProc)
        END IF
      END IF
    END DO
  END DO
END DO
DEALLOCATE(IJK_globalProc)
DEALLOCATE(posIProc_tmp)
DEALLOCATE(RecvProcList)
DEALLOCATE(ElemsToSend_tmp)

!Build mapping to sort the globally averaged data which has been received by the master back into all i,j elements
ALLOCATE(MyElemsIJ(nElemsloc_IJ(1),nElemsloc_IJ(2),nElems))
ALLOCATE(nElemsIJ(nElemsloc_IJ(1),nElemsloc_IJ(2)))
MyElemsIJ=0
nElemsIJ=0
DO iElem = 1,nElems
  i=Elem_IJK(IJK_Mask(1),iElem)-minIJ(1)
  j=Elem_IJK(IJK_Mask(2),iElem)-minIJ(2)
  nElemsIJ(i,j)=nElemsIJ(i,j)+1
  MyElemsIJ(i,j,nElemsIJ(i,j))=iElem
END DO

nDiffElemsloc = 0
ALLOCATE(posLocTmp(nElemsloc_IJ(1),nElemsloc_IJ(2)))
posLocTmp=0
DO j = 1,nElemsloc_IJ(2)
  DO i = 1,nElemsloc_IJ(1)
    IF(nElemsIJ(i,j) .GT. 0 ) THEN
      nDiffElemsloc = nDiffElemsloc +1
      posLocTmp(i,j)=nDiffElemsloc
    END IF
  END DO
END DO
ALLOCATE(iDiffElem(nElems))
iDiffElem = -1
DO iElem=1,nElems
  i=Elem_IJK(IJK_Mask(1),iElem)-minIJ(1)
  j=Elem_IJK(IJK_Mask(2),iElem)-minIJ(2)
  IF(nElemsIJ(i,j) .GT. 0 ) THEN
    iDiffElem(iElem)=posLocTmp(i,j)
  END IF
END DO
SDEALLOCATE(posLocTmp)

! All Mappings are set up, now we can ALLOCATE all arrays to send and reive 2d averaged data
! Use only the conservative vector to average
nVarsAvg2D = PP_nVar

#if USE_MPI
ALLOCATE(SendBufferAvg2D(nVarsAvg2D,0:PP_N,0:PP_N,MAXVAL(nElemsToSend),nRecvProcs ))
ALLOCATE(RecvBufferAvg2D(nVarsAvg2D,0:PP_N,0:PP_N,MAXVAL(nElemsToRecv) ,nSendProcs))
#endif /* USE_MPI */
ALLOCATE(UAvg2D         (nVarsAvg2D,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(UAvg2DLocal    (nVarsAvg2D,0:PP_N,0:PP_N,nDiffElemsloc))
ALLOCATE(UAvg2DGlobal   (nVarsAvg2D,0:PP_N,0:PP_N,nDiffElemsloc))
#if USE_MPI
SendBufferAvg2D = 0.
RecvBufferAvg2D = 0.
#endif /* USE_MPI */
UAvg2D          = 0.
UAvg2DLocal     = 0.
UAvg2DGlobal    = 0.

!-------------------------------------------------------------------
! Compute the width of the domain

ALLOCATE(domainWidth(nDiffElemsloc))
domainWidth = 0.
! Every processor averages all its elements at all of its i,j position
DO iElem=1,nElems
  i=iDiffElem(iElem)
  dx = (Elem_xGP(3,0,0,PP_N,iElem) - Elem_xGP(3,0,0,0,iElem))*2.0/(xGP(PP_N)-xGP(0))
  DO k=0,PP_N
    domainWidth(i) = domainWidth(i) + wGP(k)/2.*dx * 1.0
  END DO
END DO

#if USE_MPI
ALLOCATE(MPIRequest_Avg2DSend(nRecvProcs))
ALLOCATE(MPIRequest_Avg2DRecv(nSendProcs))
ALLOCATE(SendBuffer(MAXVAL(nElemsToSend),nRecvProcs))
SendBuffer = 0.

! Start sending the domain width data
DO iDiffRecvProc = 1,nRecvProcs
  iProc    = RecvProcs(iDiffRecvProc)
  nSendVal = nElemsToSend(RecvProcs(iDiffRecvProc))
  DO iElemToSend = 1,nElemsToSend(RecvProcs(iDiffRecvProc))
    iElem = ElemsToSendIJSorted(iDiffRecvProc,iElemToSend)
    i=iDiffElem(iElem)
    SendBuffer(iElemToSend,iDiffRecvProc) = domainWidth(i)
  END DO
  CALL MPI_ISEND(SendBuffer(1:nElemsToSend(RecvProcs(iDiffRecvProc)),iDiffRecvProc),nSendVal,MPI_DOUBLE_PRECISION,  &
                 iProc,0,MPI_COMM_FLEXI,MPIRequest_Avg2DSend(iDiffRecvProc),iError)
END DO

! Start receiving the domain width data
ALLOCATE(RcvBuffer(MAXVAL(nElemsToRecv),nSendProcs))
RcvBuffer = 0.
DO iDiffSendProc = 1,nSendProcs
  iProc    = SendProcs(iDiffSendProc)
  nRecvVal = nElemsToRecv(SendProcs(iDiffSendProc))
  CALL MPI_IRecv(RcvBuffer(1:nElemsToRecv(SendProcs(iDiffSendProc)),iDiffSendProc),nRecvVal,MPI_DOUBLE_PRECISION,  &
                 iProc,0,MPI_COMM_FLEXI,MPIRequest_Avg2DRecv(iDiffSendProc),iError)
END DO
! Finish echange of domain width data
CALL MPI_WAITALL(nRecvProcs,MPIRequest_Avg2DSend(1:nRecvProcs),MPI_STATUSES_IGNORE,iError)
CALL MPI_WAITALL(nSendProcs,MPIRequest_Avg2DRecv(1:nSendProcs),MPI_STATUSES_IGNORE,iError)
#endif /* USE_MPI */

SpanWidth = -1
! Compute the global domain width
DO iDiffSendProc = 1,nSendProcs
  DO iElemToRcv = 1,nElemsToRecv(SendProcs(iDiffSendProc))
    iElem = ElemsToRecvIJSorted(iDiffSendProc,iElemToRcv)
    i=iDiffElem(iElem)
    domainWidth(i)=domainWidth(i)+RcvBuffer(iElemToRcv,iDiffSendProc)
  END DO
END DO
SpanWidth=MAXVAL(domainWidth(:))
SpanWidth=1./SpanWidth

SDEALLOCATE(domainWidth)
SDEALLOCATE(SendBuffer)
SDEALLOCATE(RcvBuffer)

SWRITE(UNIT_stdOut,'(A)') ' INIT AVG2D DONE...'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitAvg2D

!==================================================================================================================================
!> Each processor averages its own elements in z-Direction
!==================================================================================================================================
SUBROUTINE Avg2DSlave(UIn)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Avg2D_Vars
USE MOD_Mesh_Readin        ,ONLY:ReadIJKSorting
USE MOD_Mesh_Vars          ,ONLY:nElems,Elem_xGP
USE MOD_Interpolation_Vars ,ONLY:xGP,wGP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: UIn(nVarsAvg2D,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,k
REAL                :: dx
!==================================================================================================================================
! Every processor averages all its elements at all of its i,j position and multiplies it by its deltaX
UAvg2DLocal = 0.
DO iElem=1,nElems
  i  = iDiffElem(iElem)
  dx = (Elem_xGP(3,0,0,PP_N,iElem) - Elem_xGP(3,0,0,0,iElem))*2.0/(xGP(PP_N)-xGP(0))
  DO k=0,PP_N
    UAvg2DLocal(:,:,:,i) = UAvg2DLocal(:,:,:,i) + wGP(k)*0.5*dx * UIn(:,:,:,k,iElem)
  END DO
END DO
END SUBROUTINE Avg2DSlave

#if USE_MPI
!==================================================================================================================================
!> Each slave processor sends it pre-averaged elements to the master proc and opens its own receiving ports
!==================================================================================================================================
SUBROUTINE SendAvg2DToMaster()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Avg2D_Vars
USE MOD_MPI_Vars           ,ONLY:MPIRequest_Avg2DSend,MPIRequest_Avg2DRecv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,iDiffRecvProc,iDiffSendProc,iProc,nSendVal,nRecvVal,iElemToSend
!==================================================================================================================================
!Fill the send buffers and start sending
SendBufferAvg2D = 0.
DO iDiffRecvProc = 1,nRecvProcs
  iProc    = RecvProcs(iDiffRecvProc)
  nSendVal = nVarsAvg2D*(PP_N+1)*(PP_N+1)*nElemsToSend(RecvProcs(iDiffRecvProc))
  DO iElemToSend = 1,nElemsToSend(RecvProcs(iDiffRecvProc))
    iElem = ElemsToSendIJSorted(iDiffRecvProc,iElemToSend)
    i=iDiffElem(iElem)
    SendBufferAvg2D(:,:,:,iElemToSend,iDiffRecvProc) = UAvg2DLocal(:,:,:,i)
  END DO
  CALL MPI_ISEND(SendBufferAvg2D(:,:,:,1:nElemsToSend(RecvProcs(iDiffRecvProc)),iDiffRecvProc),nSendVal,MPI_DOUBLE_PRECISION,  &
                 iProc,0,MPI_COMM_FLEXI,MPIRequest_Avg2DSend(iDiffRecvProc),iError)
END DO

!Open the receive buffers
RecvBufferAvg2D = 0.
DO iDiffSendProc = 1,nSendProcs
  iProc    = SendProcs(iDiffSendProc)
  nRecvVal = nVarsAvg2D*(PP_N+1)*(PP_N+1)*nElemsToRecv(SendProcs(iDiffSendProc))
  CALL MPI_IRecv(RecvBufferAvg2D(:,:,:,1:nElemsToRecv(SendProcs(iDiffSendProc)),iDiffSendProc),nRecvVal,MPI_DOUBLE_PRECISION,  &
                 iProc,0,MPI_COMM_FLEXI,MPIRequest_Avg2DRecv(iDiffSendProc),iError)
END DO
END SUBROUTINE SendAvg2DToMaster

!==================================================================================================================================
!> Each master processor receives its i,j data to average 2d global
!==================================================================================================================================
SUBROUTINE RecvAvg2DOfMaster()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Avg2D_Vars
USE MOD_MPI_Vars           ,ONLY:MPIRequest_Avg2DSend,MPIRequest_Avg2DRecv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL MPI_WAITALL(nRecvProcs,MPIRequest_Avg2DSend(1:nRecvProcs),MPI_STATUSES_IGNORE,iError)
CALL MPI_WAITALL(nSendProcs,MPIRequest_Avg2DRecv(1:nSendProcs),MPI_STATUSES_IGNORE,iError)
END SUBROUTINE RecvAvg2DOfMaster
#endif /* USE_MPI */

!==================================================================================================================================
!> Each master processor sorts in the received data and performs averaging
!==================================================================================================================================
SUBROUTINE Avg2DMaster()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Avg2D_Vars
USE MOD_Mesh_Readin        ,ONLY:ReadIJKSorting
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,iElemToRecv,iDiffSendProc
!==================================================================================================================================
! Every processor averages all its elements at all of its i,j position and multiplies it by its deltaX
!UAvg2DGlobal = 0.
! Add its own locally averaged data
UAvg2DGlobal = UAvg2DLocal
! Add the received data
DO iDiffSendProc = 1,nSendProcs
  DO iElemToRecv = 1,nElemsToRecv(SendProcs(iDiffSendProc))
    iElem = ElemsToRecvIJSorted(iDiffSendProc,iElemToRecv)
    i=iDiffElem(iElem)
    UAvg2DGlobal(:,:,:,i)=UAvg2DGlobal(:,:,:,i)+RecvBufferAvg2D(:,:,:,iElemToRecv,iDiffSendProc)
  END DO
END DO
UAvg2DGlobal(:,:,:,:)=UAvg2DGlobal(:,:,:,:)*spanWidth
END SUBROUTINE Avg2DMaster

#if USE_MPI
!==================================================================================================================================
!> Each master processor sends it globally averaged elements to the slaves procs and opens its own receiving ports
!==================================================================================================================================
SUBROUTINE SendAvg2DToSlave()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Avg2D_Vars
USE MOD_MPI_Vars           ,ONLY:MPIRequest_Avg2DSend,MPIRequest_Avg2DRecv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,iDiffRecvProc,iDiffSendProc,iProc,nSendVal,nRecvVal,iElemToRecv
!==================================================================================================================================
!Everything has to be read in reverse to the inital sending and receive operation
!Open the receive buffers
SendBufferAvg2D = 0.
MPIRequest_Avg2DSend=0
DO iDiffRecvProc = 1,nRecvProcs
  iProc    = RecvProcs(iDiffRecvProc)
  nSendVal = nVarsAvg2D*(PP_N+1)*(PP_N+1)*nElemsToSend(RecvProcs(iDiffRecvProc))
  CALL MPI_IRecv(SendBufferAvg2D(:,:,:,1:nElemsToSend(RecvProcs(iDiffRecvProc)),iDiffRecvProc),nSendVal,MPI_DOUBLE_PRECISION,  &
                 iProc,0,MPI_COMM_FLEXI,MPIRequest_Avg2DSend(iDiffRecvProc),iError)
END DO

!Fill the send buffers and start sending
MPIRequest_Avg2DRecv=0
RecvBufferAvg2D = 0.
DO iDiffSendProc = 1,nSendProcs
  iProc    = SendProcs(iDiffSendProc)
  nRecvVal = nVarsAvg2D*(PP_N+1)*(PP_N+1)*nElemsToRecv(SendProcs(iDiffSendProc))
  DO iElemToRecv = 1,nElemsToRecv(SendProcs(iDiffSendProc))
    iElem = ElemsToRecvIJSorted(iDiffSendProc,iElemtoRecv)
    i=iDiffElem(iElem)
    RecvBufferAvg2D(:,:,:,iElemToRecv,iDiffSendProc) = UAvg2DGlobal(:,:,:,i)
  END DO
  CALL MPI_ISend(RecvBufferAvg2D(:,:,:,1:nElemsToRecv(SendProcs(iDiffSendProc)),iDiffSendProc),nRecvVal,MPI_DOUBLE_PRECISION,  &
                 iProc,0,MPI_COMM_FLEXI,MPIRequest_Avg2DRecv(iDiffSendProc),iError)
END DO
END SUBROUTINE SendAvg2DToSlave
#endif /* USE_MPI */

!==================================================================================================================================
!> Each slave processor receives its i,j data which is globally averaged and sorts it into its elements
!==================================================================================================================================
SUBROUTINE RecvAvg2DOfSlave(UOut)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Avg2D_Vars
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY:MPIRequest_Avg2DSend,MPIRequest_Avg2DRecv
#endif /* USE_MPI */
USE MOD_Mesh_Vars          ,ONLY:Elem_IJK,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)     :: UOut(nVarsAvg2D,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k,iDiffRecvProc,iElemToSend,iMyElemIJ,iElemIJ
!==================================================================================================================================
!Everything has to be read in reverse to the inital sending and receive operation
#if USE_MPI
CALL MPI_WAITALL(nRecvProcs,MPIRequest_Avg2DSend(1:nRecvProcs),MPI_STATUSES_IGNORE,iError)
CALL MPI_WAITALL(nSendProcs,MPIRequest_Avg2DRecv(1:nSendProcs),MPI_STATUSES_IGNORE,iError)
#endif /* USE_MPI */

! Sort in my own values, maybe so are actually local values but they are overriten by the global ones
UAvg2D = 0.
DO iElem=1,nElems
  i=iDiffElem(iElem)
  DO k=0,PP_N
    UAvg2D(:,:,:,k,iElem) = UAvg2DGlobal(:,:,:,i)
  END DO
END DO

!Everything has to be read in reverse to the inital sending
!Sort all averaged data back into the different elements
DO iDiffRecvProc = 1,nRecvProcs
  DO iElemToSend = 1,nElemsToSend(RecvProcs(iDiffRecvProc))
    iElem = ElemsToSendIJSorted(iDiffRecvProc,iElemToSend)
    i=Elem_IJK(IJK_Mask(1),iElem)-minIJ(1)
    j=Elem_IJK(IJK_Mask(2),iElem)-minIJ(2)
    !Loop over all my elements on this proc that are on the same i,j position
    DO iElemIJ = 1,nElemsIJ(i,j)
      iMyElemIJ = MyElemsIJ(i,j,iElemIJ)
      DO k=0,PP_N
        UAvg2D(:,:,:,k,iMyElemIJ) = SendBufferAvg2D(:,:,:,iElemToSend,iDiffRecvProc)
      END DO
    END DO
  END DO
END DO
UOut = UAvg2D
END SUBROUTINE RecvAvg2DOfSlave

!==================================================================================================================================
!> Deallocate Avg2D arrays
!==================================================================================================================================
SUBROUTINE FinalizeAvg2D()
! MODULES
USE MOD_Avg2D_Vars
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY:MPIRequest_Avg2DSend,MPIRequest_Avg2DRecv
#endif /* USE_MPI */
IMPLICIT NONE
!==================================================================================================================================
#if USE_MPI
SDEALLOCATE(MPIRequest_Avg2DRecv)
SDEALLOCATE(MPIRequest_Avg2DSend)
SDEALLOCATE(RecvBufferAvg2D)
SDEALLOCATE(SendBufferAvg2D)
SDEALLOCATE(RecvProcs)
SDEALLOCATE(SendProcs)
SDEALLOCATE(nElemsToRecv)
SDEALLOCATE(nElemsToSend)
#endif /* USE_MPI */
SDEALLOCATE(MyElemsIJ)
SDEALLOCATE(nElemsIJ)
SDEALLOCATE(UAvg2D)
SDEALLOCATE(UAvg2DLocal)
SDEALLOCATE(UAvg2DGlobal)
SDEALLOCATE(ElemsToSendIJSorted)
SDEALLOCATE(ElemsToRecvIJSorted)
SDEALLOCATE(iDiffElem)
END SUBROUTINE FinalizeAvg2D

END MODULE MOD_Avg2D
