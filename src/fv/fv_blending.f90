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
!> Module for the Finite Volume sub-cells shock capturing.
!>
!> DG elements, that are detected to contain a shock/high gradients/oscillations/..., can be switched to a Finite Volume scheme.
!> A DG element of polynomial degree N is subdivided into (N+1)^dim sub-cells (to each Gauss Point/DOF one FV sub-cell).
!> The FV sub-cells of such an element are updated using FV method with 2nd order TVD reconstruction (slope limiters).
!==================================================================================================================================
MODULE MOD_FV_Blending
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

#if FV_ENABLED == 2
PUBLIC:: FV_ExtendAlpha
#endif /*FV_ENABLED == 2*/
#if ((FV_ENABLED >= 2) && (PP_NodeType == 1))
PUBLIC:: FV_CommAlpha
#endif /*FV_ENABLED >= 2 && PP_NodeType == 1*/
#if FV_ENABLED == 2 || FV_ENABLED == 3
PUBLIC:: FV_Info
#endif /*FV_ENABLED == 2 || FV_ENABLED == 3*/
!==================================================================================================================================

CONTAINS

#if FV_ENABLED == 2
!==================================================================================================================================
!> Extend the blending coefficient FV_alpha
!==================================================================================================================================
SUBROUTINE FV_ExtendAlpha(FV_alpha)
! MODULES
USE MOD_PreProc
USE MOD_FV_Mortar        ,ONLY: FV_alpha_Mortar
USE MOD_FV_Vars          ,ONLY: FV_doExtendAlpha,FV_nExtendAlpha
USE MOD_FV_Vars          ,ONLY: FV_alpha_master,FV_alpha_slave,FV_alpha_extScale
USE MOD_Mesh_Vars        ,ONLY: nElems
#if USE_MPI
USE MOD_Mesh_Vars        ,ONLY: nSides
USE MOD_MPI              ,ONLY: StartExchange_FV_alpha,FinishExchangeMPIData
USE MOD_MPI_Vars         ,ONLY: MPIRequest_FV_Elems,nNbProcs
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: FV_alpha(nElems)    !< elementwise blending coefficient for DG/FV blending
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iElem
REAL               :: FV_alpha_current(nElems)
!==================================================================================================================================
IF (FV_doExtendAlpha) THEN
  ! Nullify arrays
  FV_alpha_master = 0.
  FV_alpha_slave  = 0.

  DO i=1,FV_nExtendAlpha
    ! Prolong blending factor to faces
    CALL FV_ProlongFValphaToFace(FV_alpha)

    ! TODO: You get here two times the network latency. Could be optimized
#if USE_MPI
    CALL FV_alpha_Mortar(FV_alpha_master,FV_alpha_slave,doMPISides=.TRUE.)
    CALL StartExchange_FV_alpha(FV_alpha_slave,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=2)
#endif
    CALL FV_alpha_Mortar(FV_alpha_master,FV_alpha_slave,doMPISides=.FALSE.)
#if USE_MPI
    CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems)

    CALL StartExchange_FV_alpha(FV_alpha_master,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=1)
    CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems)
#endif

    ! Compute maximum of neighbors across shared faces and scale with scaling factor
    FV_alpha_current = 0.
    CALL FV_ComputeExtendedAlpha(FV_alpha_master,FV_alpha_slave,FV_alpha_current,.FALSE.)
    DO iElem = 1,nElems
      FV_alpha(iElem) = MAX(FV_alpha(iElem),FV_alpha_ExtScale*FV_alpha_current(iElem))
    END DO
  END DO
END IF
END SUBROUTINE FV_ExtendAlpha


!==================================================================================================================================
!> Set FV_Alpha_slave and FV_Alpha_master information
!==================================================================================================================================
SUBROUTINE FV_ProlongFValphaToFace(FV_alpha)
! MODULES
USE MOD_FV_Vars         ,ONLY: FV_alpha_master,FV_alpha_slave
USE MOD_Mesh_Vars       ,ONLY: SideToElem,nSides,nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: FV_alpha(nElems)    !< elementwise blending coefficient for DG/FV blending
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSide,ElemID,nbElemID
!==================================================================================================================================
! array not allocated in postiMode
IF (.NOT.ALLOCATED(SideToElem)) RETURN

! set information whether elements adjacent to a side are DG or FV elements
DO iSide = 1,nSides
  ElemID    = SideToElem(S2E_ELEM_ID   ,iSide)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID,iSide)
  !master sides
  IF(ElemID  .GT.0) FV_alpha_master(iSide) = FV_alpha(ElemID)
  !slave side (ElemID,locSide and flip =-1 if not existing)
  IF(nbElemID.GT.0) FV_alpha_slave( iSide) = FV_alpha(nbElemID)
END DO
END SUBROUTINE FV_ProlongFValphaToFace


!==================================================================================================================================
!> Check for indicator value in neighboring elements
!==================================================================================================================================
SUBROUTINE FV_ComputeExtendedAlpha(FV_alpha_master,FV_alpha_slave,FV_alpha_current,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars ,ONLY: SideToElem,nSides
USE MOD_Mesh_Vars ,ONLY: firstMPISide_YOUR,lastMPISide_MINE
USE MOD_Mesh_Vars ,ONLY: nSides,nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN)   :: doMPISides  !<= .TRUE. only MPISides_YOUR+MPIMortar are filled

REAL,INTENT(IN)      :: FV_alpha_master( 1:nSides) !< (IN)  FV_alpha on master side
REAL,INTENT(IN)      :: FV_alpha_slave ( 1:nSides) !< (IN)  FV_alpha on slave side
REAL,INTENT(INOUT)   :: FV_alpha_current(1:nElems) !< (OUT) FV_alpha in elem
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,nbElemID,iSide
INTEGER            :: firstSideID,lastSideID
!==================================================================================================================================
! array not allocated in postiMode
IF (.NOT.ALLOCATED(SideToElem)) RETURN
! TODO: is this the correct behavior?

! set information whether elements adjacent to a side are DG or FV elements
DO iSide = 1,nSides
  ElemID    = SideToElem(S2E_ELEM_ID   ,iSide)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID,iSide)

  ! master sides
  IF (ElemID .GT. 0) &
    FV_alpha_current(ElemID)   = MAX(FV_alpha_current(ElemID)  ,MAX(FV_alpha_master(iSide),FV_alpha_slave(iSide)))

  ! slave sides
  IF (nbElemID .GT. 0) &
    FV_alpha_current(nbElemID) = MAX(FV_alpha_current(nbElemID),MAX(FV_alpha_master(iSide),FV_alpha_slave(iSide)))
END DO
END SUBROUTINE FV_ComputeExtendedAlpha
#endif


#if ((FV_ENABLED >= 2) && (PP_NodeType == 1))
!==================================================================================================================================
!> Extend the blending coefficient FV_alpha
!==================================================================================================================================
SUBROUTINE FV_CommAlpha(FV_alpha)
! MODULES
USE MOD_PreProc
USE MOD_FV_Mortar        ,ONLY: FV_alpha_Mortar
USE MOD_FV_Vars          ,ONLY: FV_alpha_master,FV_alpha_slave
#if USE_MPI
USE MOD_Mesh_Vars        ,ONLY: nSides,nElems
USE MOD_MPI              ,ONLY: StartExchange_FV_alpha,FinishExchangeMPIData
USE MOD_MPI_Vars         ,ONLY: MPIRequest_FV_Elems,nNbProcs
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: FV_alpha(nElems)    !< elementwise blending coefficient for DG/FV blending
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! TODO: You get here two times the network latency. Could be optimized
CALL FV_ProlongFValphaToFace(FV_alpha)
#if USE_MPI
! Prolong blending factor to faces
! CALL FV_ProlongFValphaToFace(doMPISides=.TRUE.)
CALL FV_alpha_Mortar(FV_alpha_master,FV_alpha_slave,doMPISides=.TRUE.)
CALL StartExchange_FV_alpha(FV_alpha_slave,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=2)
#endif
! CALL FV_ProlongFValphaToFace(doMPISides=.FALSE.)
CALL FV_alpha_Mortar(FV_alpha_master,FV_alpha_slave,doMPISides=.FALSE.)
#if USE_MPI
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems)

CALL StartExchange_FV_alpha(FV_alpha_master,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems)
#endif
END SUBROUTINE FV_CommAlpha
#endif /* FV_ENABLED */


#if FV_ENABLED == 2 || FV_ENABLED == 3
!==================================================================================================================================
!> Print information on the amount of FV blending
!==================================================================================================================================
SUBROUTINE FV_Info(iter)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: nGlobalElems
USE MOD_FV_Vars      ,ONLY: FV_alpha
USE MOD_Analyze_Vars ,ONLY: FV_totalAlpha
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER(KIND=DP),INTENT(IN) :: iter !< number of iterations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: FV_alpha_range(3)
!==================================================================================================================================
FV_alpha_range(1) = MINVAL(FV_alpha)
FV_alpha_range(2) = MAXVAL(FV_alpha)

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE    ,FV_alpha_range(1),1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE    ,FV_alpha_range(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE    ,FV_totalAlpha    ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(FV_alpha_range(1),0               ,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(FV_alpha_range(2),0               ,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
  CALL MPI_REDUCE(FV_totalAlpha    ,0               ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif /*USE_MPI*/

#if FV_ENABLED == 3
! DOF-wise FV_alpha
FV_totalAlpha = FV_totalAlpha / ((PP_N +1)*(PP_N +1)*(PP_NZ+1))
#endif /*FV_ENABLED*/

SWRITE(UNIT_stdOut,'(A,F8.3,A,F5.3,A,ES18.9)') ' FV_alpha    : ',FV_alpha_range(1),' - ',FV_alpha_range(2),&
                                              ', avg: '         ,FV_totalAlpha / REAL(nGlobalElems) / iter
FV_totalAlpha   = 0.
END SUBROUTINE FV_Info

#endif /* FV_ENABLED */
END MODULE MOD_FV_Blending
