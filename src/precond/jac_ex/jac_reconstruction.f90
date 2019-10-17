!!=================================================================================================================================
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
!> Contains the routines for the derivative of the fv_reconstruction procedure
!==================================================================================================================================
MODULE MOD_Jac_Reconstruction
#if FV_ENABLED && FV_RECONSTRUCT
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE Fill_ExtendedState
  MODULE PROCEDURE Fill_ExtendedState
END INTERFACE

INTERFACE FV_Reconstruction_Derivative
  MODULE PROCEDURE FV_Reconstruction_Derivative
END INTERFACE

INTERFACE FV_Reconstruction_Derivative_Surf
  MODULE PROCEDURE FV_Reconstruction_Derivative_Surf
END INTERFACE

#if PARABOLIC
INTERFACE JacFVGradients_Vol
  MODULE PROCEDURE JacFVGradients_Vol
END INTERFACE

INTERFACE JacFVGradients_nb
  MODULE PROCEDURE JacFVGradients_nb
END INTERFACE
#endif

PUBLIC::Fill_ExtendedState,FV_Reconstruction_Derivative,FV_Reconstruction_Derivative_Surf
#if PARABOLIC
PUBLIC::JacFVGradients_Vol,JacFVGradients_nb
#endif
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Fills extended arrays containing
!> - volume value per element
!> - values at the nodes of the first inner layer next to the DG element interface 
!> Those arrays are: UPrim_extended (containing primitive solution vector) and
!> FV_sdx_XI/ETA/ZETA_extended containing the inverse of the distance between subcell centers.
!==================================================================================================================================
SUBROUTINE Fill_ExtendedState(t,URec,URec_extended,FV_sdx_XI_extended,FV_sdx_ETA_extended &
#if PP_dim ==3 
                              ,FV_sdx_ZETA_extended &
#endif
                              )
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ProlongToFacePrim  ,ONLY: ProlongToFacePrim
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus
USE MOD_Mesh_Vars          ,ONLY: nElems,nSides,firstInnerSide,lastBCSide,lastMPISide_MINE
USE MOD_Mesh_Vars          ,ONLY: S2V,SideToElem
USE MOD_GetBoundaryFlux    ,ONLY: GetBoundaryState
USE MOD_Mesh_Vars          ,ONLY: NormVec,TangVec1,TangVec2,Face_xGP
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf
#if USE_MPI
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_MPI_Vars           ,ONLY: MPIRequest_Rec_SM,MPIRequest_Rec_MS,nNbProcs,DataSizeSidePrim
#endif
USE MOD_FV_Vars            ,ONLY: FV_Elems_Sum,FV_sdx_Face,FV_Vdm,FV_Elems
USE MOD_FV_Vars            ,ONLY: FV_sdx_XI,FV_sdx_ETA
#if PP_dim == 3
USE MOD_FV_Vars            ,ONLY: FV_sdx_ZETA
#endif
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)     :: t                                                                 !< physical time
REAL,INTENT(IN)     :: URec(         PP_nVarPrim, 0:PP_N,   0:PP_N,   0:PP_NZ, 1:nElems) !< rec volume solution
REAL,INTENT(OUT)    :: FV_sdx_XI_extended(        0:PP_N,   0:PP_NZ,  0:PP_N+1,1:nElems) !< FV_sdx_XI:storage order   (j,k,i,iElem)
REAL,INTENT(OUT)    :: FV_sdx_ETA_extended(       0:PP_N,   0:PP_NZ,  0:PP_N+1,1:nElems) !< FV_sdx_ETA:storage order  (i,k,j,iElem)
#if PP_dim == 3
REAL,INTENT(OUT)    :: FV_sdx_ZETA_extended(      0:PP_N,   0:PP_N,   0:PP_N+1,1:nElems) !< FV_sdx_ZETA:storage order (i,j,k,iElem)
REAL,INTENT(OUT)    :: URec_extended(PP_nVarPrim,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,1:nElems) !< extended rec volume solution
#else
REAL,INTENT(OUT)    :: URec_extended(PP_nVarPrim,-1:PP_N+1,-1:PP_N+1, 0:PP_NZ, 1:nElems) !< extended rec volume solution
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:ZDIM(PP_N),1:nSides)  :: URec_master,URec_slave
INTEGER                                                   :: p,q,l,i,j,k,ijk(3),locSideID,ElemID,SideID,flip,iElem
REAL,DIMENSION(0:PP_N,0:PP_NZ,1:nSides)                   :: FV_sdx_Face_loc
REAL,DIMENSION(1:PP_nVarPrim,0:PP_N,0:PP_NZ)              :: UPrim_Boundary
!==================================================================================================================================
! 1. Fill inner date with element volume data
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    URec_extended(:,i,j,k,iElem) = URec(:,i,j,k,iElem)
  END DO; END DO; END DO! i,j,k=0,PP_N
END DO ! iElem

! 2.   Add volume date of first layer of neighbouring cell
! 2.1. Copy first layer to surface and send data from slave to master and vice versa
#if USE_MPI
CALL StartReceiveMPIData(URec_slave ,DataSizeSidePrim,1,nSides,MPIRequest_Rec_SM(:,SEND),SendID=2) ! slave -> master
CALL StartReceiveMPIData(URec_master,DataSizeSidePrim,1,nSides,MPIRequest_Rec_MS(:,SEND),SendID=1) ! master -> slave
CALL ProlongToFacePrim(PP_N,URec,URec_master,URec_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL StartSendMPIData(   URec_slave ,DataSizeSidePrim,1,nSides,MPIRequest_Rec_SM(:,RECV),SendID=2) ! slave -> master
#endif
CALL ProlongToFacePrim(PP_N,URec,URec_master,URec_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
#if USE_MPI
! prolongation has to be finished before communication of URec_master as doMPISides=.FALSE. includes MPISidesMINE
CALL StartSendMPIData(   URec_master,DataSizeSidePrim,1,nSides,MPIRequest_Rec_MS(:,RECV),SendID=1) ! master -> slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Rec_SM)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Rec_MS)
#endif

! 3.   Fill information about distance between solution points at faces (FV_sdx_Face)
! 3.2. Send FV_sdx_Face from master to slave
DO SideID=firstInnerSide,lastMPISide_MINE
  SELECT CASE(FV_Elems_Sum(SideID))
  CASE(0) ! both DG
    CYCLE
  CASE(1) ! master=FV, slave=DG
    FV_sdx_Face_loc(:,:,SideID) = FV_sdx_Face(:,:,1,SideID)
  CASE(2) ! master=DG, slave=FV
    FV_sdx_Face_loc(:,:,SideID) = FV_sdx_Face(:,:,2,SideID)
  CASE(3) ! both FV
    FV_sdx_Face_loc(:,:,SideID) = FV_sdx_Face(:,:,3,SideID)
  CASE DEFAULT
    CALL Abort(__STAMP__, "FV_Elems_master/slave is not 0 or 1 somewhere!")
  END SELECT
END DO
#if USE_MPI
CALL StartReceiveMPIData(FV_sdx_Face_loc,(PP_N+1)*(PP_NZ+1),1,nSides,MPIRequest_Rec_MS(:,SEND),SendID=1) ! master -> slave
CALL StartSendMPIData(   FV_sdx_Face_loc,(PP_N+1)*(PP_NZ+1),1,nSides,MPIRequest_Rec_MS(:,RECV),SendID=1) ! master -> slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Rec_MS)
#endif
! 3.3. Fill FV_sdx_extended array with inner data
DO iElem=1,nElems
  DO p=0,PP_N; DO q=0,PP_NZ; DO l=1,PP_N
    FV_sdx_XI_extended(  p,q,l,iElem) = FV_sdx_XI(  p,q,l,iElem)
    FV_sdx_ETA_extended( p,q,l,iElem) = FV_sdx_ETA( p,q,l,iElem)
#if PP_dim == 3
    FV_sdx_ZETA_extended(p,q,l,iElem) = FV_sdx_ZETA(p,q,l,iElem)
#endif
  END DO; END DO; END DO
END DO

! 2.2. Switch Solution representation of URec_slave/master to FV if corresponding element is a DG element
! 2.3. Add first layer of URec of neighbouring cell to extended array
! 3.4. Add face data to extended FV_sdx array (take care of storage order of FV_sdx!)
! First, process the slave sides
DO SideID=firstInnerSide,nSides
  ! neighbor side !ElemID,locSideID and flip =-1 if not existing
  ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip      = SideToElem(S2E_FLIP,SideID)
  IF(ElemID.EQ.-1) CYCLE ! skip for firstMPISide_MINE:lastMPISide_MINE
  IF(FV_Elems(ElemID).EQ.0) CYCLE ! skip for DG elements
  IF(FV_Elems_Sum(SideID).EQ.2) CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_Vdm,URec_master(:,:,:,SideID)) ! switch master to FV
  DO q=0,PP_NZ; DO p=0,PP_N
    ijk=S2V(:,0,p,q,flip,locSideID)
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      URec_extended(:,     ijk(1)-1,ijk(2)  ,ijk(3)  ,ElemID) = URec_master(   :,p,q,SideID)
      FV_sdx_XI_extended(  ijk(2)  ,ijk(3)  ,ijk(1)  ,ElemID) = FV_sdx_face_loc( p,q,SideID)
    CASE(XI_PLUS)
      URec_extended(:,     ijk(1)+1,ijk(2)  ,ijk(3)  ,ElemID) = URec_master(   :,p,q,SideID)
      FV_sdx_XI_extended(  ijk(2)  ,ijk(3)  ,ijk(1)+1,ElemID) = FV_sdx_face_loc( p,q,SideID)
    CASE(ETA_MINUS)
      URec_extended(:,     ijk(1)  ,ijk(2)-1,ijk(3)  ,ElemID) = URec_master(   :,p,q,SideID)
      FV_sdx_ETA_extended( ijk(1)  ,ijk(3)  ,ijk(2)  ,ElemID) = FV_sdx_face_loc( p,q,SideID)
    CASE(ETA_PLUS)
      URec_extended(:,     ijk(1)  ,ijk(2)+1,ijk(3)  ,ElemID) = URec_master(   :,p,q,SideID)
      FV_sdx_ETA_extended( ijk(1)  ,ijk(3)  ,ijk(2)+1,ElemID) = FV_sdx_face_loc( p,q,SideID)
#if PP_dim == 3
    CASE(ZETA_MINUS)
      URec_extended(:,     ijk(1)  ,ijk(2)  ,ijk(3)-1,ElemID) = URec_master(   :,p,q,SideID)
      FV_sdx_ZETA_extended(ijk(1)  ,ijk(2)  ,ijk(3)  ,ElemID) = FV_sdx_face_loc( p,q,SideID)
    CASE(ZETA_PLUS)
      URec_extended(:,     ijk(1)  ,ijk(2)  ,ijk(3)+1,ElemID) = URec_master(   :,p,q,SideID)
      FV_sdx_ZETA_extended(ijk(1)  ,ijk(2)  ,ijk(3)+1,ElemID) = FV_sdx_face_loc( p,q,SideID)
#endif
    END SELECT
  END DO; END DO
END DO
! Second, process the master sides
DO SideID=firstInnerSide,lastMPISide_MINE
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  IF(FV_Elems(ElemID).EQ.0) CYCLE ! skip for DG elements
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  IF(FV_Elems_Sum(SideID).EQ.1) CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_Vdm,URec_slave(:,:,:,SideID)) ! switch slave to FV
  DO q=0,PP_NZ; DO p=0,PP_N
    ijk=S2V(:,0,p,q,0,locSideID)
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      URec_extended(:,     ijk(1)-1,ijk(2)  ,ijk(3)  ,ElemID) = URec_slave(   :,p,q,SideID)
      FV_sdx_XI_extended(  ijk(2)  ,ijk(3)  ,ijk(1)  ,ElemID) = FV_sdx_face_loc(p,q,SideID)
    CASE(XI_PLUS)
      URec_extended(:,     ijk(1)+1,ijk(2)  ,ijk(3)  ,ElemID) = URec_slave(   :,p,q,SideID)
      FV_sdx_XI_extended(  ijk(2)  ,ijk(3)  ,ijk(1)+1,ElemID) = FV_sdx_face_loc(p,q,SideID)
    CASE(ETA_MINUS)
      URec_extended(:,     ijk(1)  ,ijk(2)-1,ijk(3)  ,ElemID) = URec_slave(   :,p,q,SideID)
      FV_sdx_ETA_extended( ijk(1)  ,ijk(3)  ,ijk(2)  ,ElemID) = FV_sdx_face_loc(p,q,SideID)
    CASE(ETA_PLUS)
      URec_extended(:,     ijk(1)  ,ijk(2)+1,ijk(3)  ,ElemID) = URec_slave(   :,p,q,SideID)
      FV_sdx_ETA_extended( ijk(1)  ,ijk(3)  ,ijk(2)+1,ElemID) = FV_sdx_face_loc(p,q,SideID)
#if PP_dim == 3
    CASE(ZETA_MINUS)
      URec_extended(:,     ijk(1)  ,ijk(2)  ,ijk(3)-1,ElemID) = URec_slave(   :,p,q,SideID)
      FV_sdx_ZETA_extended(ijk(1)  ,ijk(2)  ,ijk(3)  ,ElemID) = FV_sdx_face_loc(p,q,SideID)
    CASE(ZETA_PLUS)
      URec_extended(:,     ijk(1)  ,ijk(2)  ,ijk(3)+1,ElemID) = URec_slave(   :,p,q,SideID)
      FV_sdx_ZETA_extended(ijk(1)  ,ijk(2)  ,ijk(3)+1,ElemID) = FV_sdx_face_loc(p,q,SideID)
#endif
    END SELECT
  END DO; END DO
END DO
! Third, do boundary sides
DO SideID=1,lastBCSide
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  IF(FV_Elems(ElemID).EQ.0) CYCLE ! skip for DG elements
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  ! set state with boundary condition
  CALL GetBoundaryState(SideID,t,PP_N,UPrim_boundary,URec_master(:,:,:,SideID),                    &
                        NormVec(:,:,:,1,SideID),TangVec1(:,:,:,1,SideID),TangVec2(:,:,:,1,SideID), &
                        Face_xGP(:,:,:,1,SideID))
  DO q=0,PP_NZ; DO p=0,PP_N
    ijk=S2V(:,0,p,q,0,locSideID)
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      URec_extended(:,     ijk(1)-1,ijk(2)  ,ijk(3)  ,ElemID) = UPrim_boundary(:,p,q         )
      FV_sdx_XI_extended(  ijk(2)  ,ijk(3)  ,ijk(1)  ,ElemID) = FV_sdx_face(     p,q,3,SideID)
    CASE(XI_PLUS)
      URec_extended(:,     ijk(1)+1,ijk(2)  ,ijk(3)  ,ElemID) = UPrim_boundary(:,p,q         )
      FV_sdx_XI_extended(  ijk(2)  ,ijk(3)  ,ijk(1)+1,ElemID) = FV_sdx_face(     p,q,3,SideID)
    CASE(ETA_MINUS)
      URec_extended(:,     ijk(1)  ,ijk(2)-1,ijk(3)  ,ElemID) = UPrim_boundary(:,p,q         )
      FV_sdx_ETA_extended( ijk(1)  ,ijk(3)  ,ijk(2)  ,ElemID) = FV_sdx_face(     p,q,3,SideID)
    CASE(ETA_PLUS)
      URec_extended(:,     ijk(1)  ,ijk(2)+1,ijk(3)  ,ElemID) = UPrim_boundary(:,p,q         )
      FV_sdx_ETA_extended( ijk(1)  ,ijk(3)  ,ijk(2)+1,ElemID) = FV_sdx_face(     p,q,3,SideID)
#if PP_dim == 3
    CASE(ZETA_MINUS)
      URec_extended(:,     ijk(1)  ,ijk(2)  ,ijk(3)-1,ElemID) = UPrim_boundary(:,p,q         )
      FV_sdx_ZETA_extended(ijk(1)  ,ijk(2)  ,ijk(3)  ,ElemID) = FV_sdx_face(     p,q,3,SideID)
    CASE(ZETA_PLUS)
      URec_extended(:,     ijk(1)  ,ijk(2)  ,ijk(3)+1,ElemID) = UPrim_boundary(:,p,q         )
      FV_sdx_ZETA_extended(ijk(1)  ,ijk(2)  ,ijk(3)+1,ElemID) = FV_sdx_face(     p,q,3,SideID)
#endif
    END SELECT
  END DO; END DO
END DO

END SUBROUTINE Fill_ExtendedState

!===================================================================================================================================
!> Calculates DOF-wise Jacobian of reconstruction procedure
!> dF/dUvol = dF/dU_LR * dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim * dUvol_prim/dUvol
!>               |            |                       |                   |
!>  FD in DGVolIntJac_FV    dCons/dPrim   derivative of reconstruction  dPrim/dCons
!===================================================================================================================================
SUBROUTINE FV_Reconstruction_Derivative(FV_sdx,FV_dx_L,FV_dx_R,UPrim_plus,UPrim_minus, &
                                        URec_extended,dUdUvol_plus,dUdUvol_minus)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_FV_Limiter     ,ONLY: FV_Limiter
USE MOD_FV_Vars        ,ONLY: LimiterType
USE MOD_Jacobian       ,ONLY: dPrimdCons,dConsdPrim
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: FV_sdx(0:PP_N+1)                          !< extended inverse of distance between centers of FV subcells
REAL,INTENT(IN)    :: FV_dx_L(0:PP_N)                           !< left distance between face and center of FV subcells
REAL,INTENT(IN)    :: FV_dx_R(0:PP_N)                           !< right distance between face and center of FV subcells
REAL,INTENT(IN)    :: UPrim_plus(PP_nVarPrim,0:PP_N)            !< primitive reconstructed value at plus side
REAL,INTENT(IN)    :: UPrim_minus(PP_nVarPrim,0:PP_N)           !< primitive reconstructed value at minus side
REAL,INTENT(IN)    :: URec_extended(PP_nVarPrim,-1:PP_N+1)      !< extended rec volume solution
REAL,INTENT(OUT)   :: dUdUvol_plus( PP_nVar,PP_nVar,0:PP_N,1:3) !< Jacobian of reconstruction procedure of left  state of interface
REAL,INTENT(OUT)   :: dUdUvol_minus(PP_nVar,PP_nVar,0:PP_N,1:3) !< Jacobian of reconstruction procedure of right state of interface
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                            :: iVar,i,ind
REAL,DIMENSION(PP_nVarPrim,0:PP_N)                 :: s_L,s_R,s_lim
REAL,DIMENSION(PP_nVar,PP_nVar, 0:PP_N  )          :: Jac_ConsPrim_plus,Jac_ConsPrim_minus
REAL,DIMENSION(PP_nVar,PP_nVar,-1:PP_N+1)          :: Jac_PrimCons
REAL,DIMENSION(PP_nVar,PP_nVar)                    :: matrix !< has to be used due to intel compiler error with matmul
REAL,DIMENSION(PP_nVar,PP_nVar)                    :: vector !< has to be used due to intel compiler error with matmul
!==================================================================================================================================
! storage order in dUdUvol(1:3):
! di/d(i-1), di/d(i), di/d(i+1)
dUdUvol_plus  = 0.
dUdUvol_minus = 0.

! 1. Do derivative of i_th surface data from cons to prim
DO i=0,PP_N
  CALL dConsdPrim(UPrim_plus( :,i),Jac_ConsPrim_plus( :,:,i))
  CALL dConsdPrim(UPrim_minus(:,i),Jac_ConsPrim_minus(:,:,i))
END DO
! 2. Do derivative of i-1:i+1 volume data from prim to cons 
DO i=-1,PP_N+1
  CALL dPrimdCons(URec_extended(:,i),Jac_PrimCons(:,:,i))
END DO

! 3. Calculate derivative of reconstruction
DO i=0,PP_N
  ! Calculate left and right gradients
  s_L(:,i) = (URec_extended(:,i  ) - URec_extended(:,i-1)) * FV_sdx(i  ) 
  s_R(:,i) = (URec_extended(:,i+1) - URec_extended(:,i  )) * FV_sdx(i+1)
  ! Limit gradients
  CALL FV_Limiter(s_L(:,i),s_R(:,i),s_lim(:,i)) ! only null- or minmod-limiter
  SELECT CASE(LimiterType)
  CASE(0) ! NullLimiter
    DO iVar=1,PP_nVar
      dUdUvol_minus(iVar,iVar,i,2) = 1.
      dUdUvol_plus( iVar,iVar,i,2) = 1.
    END DO !iVar
  CASE(1) ! MinMod
    DO iVar=1,PP_nVar
      IF(s_lim(iVar,i).EQ.0.)THEN ! first order
        dUdUvol_plus( iVar,iVar,i,2) = 1.
        dUdUvol_minus(iVar,iVar,i,2) = 1.
      ELSEIF(s_lim(iVar,i).EQ.s_L(iVar,i))THEN ! use left gradient
        dUdUvol_plus( iVar,iVar,i,1) = 0. - FV_sdx(i  ) * FV_dx_R(i)
        dUdUvol_plus( iVar,iVar,i,2) = 1. + FV_sdx(i  ) * FV_dx_R(i)

        dUdUvol_minus(iVar,iVar,i,1) = 0. + FV_sdx(i  ) * FV_dx_L(i)
        dUdUvol_minus(iVar,iVar,i,2) = 1. - FV_sdx(i  ) * FV_dx_L(i)
      ELSEIF(s_lim(iVar,i).EQ.s_R(iVar,i))THEN ! use right gradient
        dUdUvol_plus( iVar,iVar,i,2) = 1. - FV_sdx(i+1) * FV_dx_R(i)
        dUdUvol_plus( iVar,iVar,i,3) = 0. + FV_sdx(i+1) * FV_dx_R(i)

        dUdUvol_minus(iVar,iVar,i,2) = 1. + FV_sdx(i+1) * FV_dx_L(i)
        dUdUvol_minus(iVar,iVar,i,3) = 0. - FV_sdx(i+1) * FV_dx_L(i)
      ELSE
        CALL Abort(__STAMP__,'Slopes do not match with minmod in preconditioner!')
      END IF
    END DO !iVar
  CASE(9) ! Central
    DO iVar=1,PP_nVar
      dUdUvol_plus( iVar,iVar,i,1) = 0. + FV_dx_R(i) * 0.5*(-FV_sdx(i))
      dUdUvol_plus( iVar,iVar,i,2) = 1. + FV_dx_R(i) * 0.5*( FV_sdx(i)-FV_sdx(i+1))
      dUdUvol_plus( iVar,iVar,i,3) = 0. + FV_dx_R(i) * 0.5*( FV_sdx(i+1))

      dUdUvol_minus(iVar,iVar,i,1) = 0. - FV_dx_L(i) * 0.5*(-FV_sdx(i))
      dUdUvol_minus(iVar,iVar,i,2) = 1. - FV_dx_L(i) * 0.5*( FV_sdx(i) - FV_sdx(i+1))
      dUdUvol_minus(iVar,iVar,i,3) = 0. - FV_dx_L(i) * 0.5*( FV_sdx(i+1))
    END DO !iVar
  CASE DEFAULT 
    CALL Abort(__STAMP__,'No preconditioner for chosen limiter implemented!')
  END SELECT

  ! multiply: dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim
  DO ind=1,3
    dUdUvol_plus( :,:,i,ind) = MATMUL(Jac_ConsPrim_plus( :,:,i),dUdUvol_plus( :,:,i,ind))
    dUdUvol_minus(:,:,i,ind) = MATMUL(Jac_ConsPrim_minus(:,:,i),dUdUvol_minus(:,:,i,ind))
  END DO
  !multiply: (dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim) * dUvol_prim/dUvol
  matrix=dUdUvol_plus(:,:,i,1); vector=Jac_PrimCons(:,:,i-1)
  dUdUvol_plus( :,:,i,1) = MATMUL(matrix,vector)
  matrix=dUdUvol_plus(:,:,i,2); vector=Jac_PrimCons(:,:,i)
  dUdUvol_plus( :,:,i,2) = MATMUL(matrix,vector)
  matrix=dUdUvol_plus(:,:,i,3); vector=Jac_PrimCons(:,:,i+1)
  dUdUvol_plus( :,:,i,3) = MATMUL(matrix,vector)
  matrix=dUdUvol_minus(:,:,i,1); vector=Jac_PrimCons(:,:,i-1)
  dUdUvol_minus(:,:,i,1) = MATMUL(matrix,vector)
  matrix=dUdUvol_minus(:,:,i,2); vector=Jac_PrimCons(:,:,i)
  dUdUvol_minus(:,:,i,2) = MATMUL(matrix,vector)
  matrix=dUdUvol_minus(:,:,i,3); vector=Jac_PrimCons(:,:,i+1)
  dUdUvol_minus(:,:,i,3) = MATMUL(matrix,vector)
 
END DO !i

END SUBROUTINE FV_Reconstruction_Derivative

!===================================================================================================================================
!> Calculates Jacobian of reconstruction procedure of surface DOFs
!> dF/dUvol = dF/dU_LR * dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim * dUvol_prim/dUvol
!>               |            |                       |                   |
!>  FD in DGVolIntJac_FV    dCons/dPrim   derivative of reconstruction  dPrim/dCons
!===================================================================================================================================
SUBROUTINE FV_Reconstruction_Derivative_Surf(FV_sdx,FV_dx_L,FV_dx_R,FV_dx_L_nb,FV_dx_R_nb,UPrim_plus,UPrim_minus, &
                                             UPrim_plus_nb,UPrim_minus_nb,                                        &
                                             URec_extended,dUdUvol_plus,dUdUvol_minus)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_FV_Limiter     ,ONLY: FV_Limiter
USE MOD_FV_Vars        ,ONLY: LimiterType
USE MOD_Jacobian       ,ONLY: dPrimdCons,dConsdPrim
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: FV_sdx(0:PP_N+1)                               !< extended inverse of distance between centers of FV subcells
REAL,INTENT(IN)    :: FV_dx_L                                        !< left distance between face and center of FV subcells
REAL,INTENT(IN)    :: FV_dx_R                                        !< right distance between face and center of FV subcells
REAL,INTENT(IN)    :: FV_dx_L_nb                                     !< left dist. between face and center of FV subcells (nbElem)
REAL,INTENT(IN)    :: FV_dx_R_nb                                     !< right dist. between face and center of FV subcells (nbElem)
REAL,INTENT(IN)    :: UPrim_plus(PP_nVarPrim)                        !< primitive reconstructed value at plus side
REAL,INTENT(IN)    :: UPrim_minus(PP_nVarPrim)                       !< primitive reconstructed value at minus side
REAL,INTENT(IN)    :: UPrim_plus_nb(PP_nVarPrim)                     !< primitive reconstructed value at plus side (nbElem)
REAL,INTENT(IN)    :: UPrim_minus_nb(PP_nVarPrim)                    !< primitive reconstructed value at minus side (nbElem)
REAL,INTENT(IN)    :: URec_extended(PP_nVarPrim,-1:PP_N+1)           !< extended rec volume solution
REAL,INTENT(OUT)   :: dUdUvol_plus( PP_nVar,PP_nVar,PP_N:PP_N+1,1:2) !< Jacobian of reconstr. procedure at plus side of element
REAL,INTENT(OUT)   :: dUdUvol_minus(PP_nVar,PP_nVar,-1:0,2:3)        !< Jacobian of reconstr. procedure at minus side of element
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                            :: iVar,ind,i
REAL,DIMENSION(PP_nVarPrim)                        :: s_L_minus,s_R_minus,s_L_plus,s_R_plus,s_lim_minus,s_lim_plus
REAL,DIMENSION(PP_nVar,PP_nVar,PP_N:PP_N+1)        :: Jac_ConsPrim_plus
REAL,DIMENSION(PP_nVar,PP_nVar,-1:0)               :: Jac_ConsPrim_minus
REAL,DIMENSION(PP_nVar,PP_nVar, 0:1)               :: Jac_PrimCons_minus
REAL,DIMENSION(PP_nVar,PP_nVar,PP_N-1:PP_N)        :: Jac_PrimCons_plus
REAL,DIMENSION(PP_nVar,PP_nVar)                    :: matrix !< has to be used due to intel compiler error with matmul
REAL,DIMENSION(PP_nVar,PP_nVar)                    :: vector !< has to be used due to intel compiler error with matmul
!==================================================================================================================================
! storage order in dUdUvol(1:3):
! di/d(i-1), di/d(i), di/d(i+1)
dUdUvol_plus  = 0.
dUdUvol_minus = 0.

! 1. Do derivative of surface data from cons to prim
CALL dConsdPrim(UPrim_minus_nb,Jac_ConsPrim_minus(:,:,-1    ))
CALL dConsdPrim(UPrim_minus   ,Jac_ConsPrim_minus(:,:, 0    ))
CALL dConsdPrim(UPrim_plus    ,Jac_ConsPrim_plus( :,:,PP_N  ))
CALL dConsdPrim(UPrim_plus_nb ,Jac_ConsPrim_plus( :,:,PP_N+1))

! 2. Do derivative of interface and neighbouring volume data from prim to cons 
DO i=0,1
  CALL dPrimdCons(URec_extended(:,i),Jac_PrimCons_minus(:,:,i))
END DO
DO i=PP_N-1,PP_N
  CALL dPrimdCons(URec_extended(:,i),Jac_PrimCons_plus( :,:,i))
END DO

! 3. Calculate derivative of reconstruction
! Calculate left and right gradients
s_L_minus(:) = (URec_extended(:,0     ) - URec_extended(:,-1    )) * FV_sdx(0     ) 
s_R_minus(:) = (URec_extended(:,1     ) - URec_extended(:, 0    )) * FV_sdx(1     )
s_L_plus( :) = (URec_extended(:,PP_N  ) - URec_extended(:,PP_N-1)) * FV_sdx(PP_N  ) 
s_R_plus( :) = (URec_extended(:,PP_N+1) - URec_extended(:,PP_N  )) * FV_sdx(PP_N+1)
! Limit gradients
CALL FV_Limiter(s_L_minus(:),s_R_minus(:),s_lim_minus(:)) ! only null-, minmod- or central-limiter
CALL FV_Limiter(s_L_plus( :),s_R_plus( :),s_lim_plus( :)) ! only null-, minmod- or central-limiter
SELECT CASE(LimiterType)
CASE(0) ! NullLimiter
  DO iVar=1,PP_nVar
    dUdUvol_minus(iVar,iVar,0   ,2) = 1.
    dUdUvol_plus( iVar,iVar,PP_N,2) = 1.
  END DO !iVar
CASE(1) ! MinMod
  DO iVar=1,PP_nVar
    ! minus side of element
    ! derivatives of minus value at minus interface
    IF(s_lim_minus(iVar).EQ.0.)THEN ! first order
      dUdUvol_minus(iVar,iVar,0,2) = 1.
    ELSEIF(s_lim_minus(iVar).EQ.s_L_minus(iVar))THEN ! use left gradient
      dUdUvol_minus(iVar,iVar,0,2) = 1. - FV_sdx(0) * FV_dx_L
    ELSEIF(s_lim_minus(iVar).EQ.s_R_minus(iVar))THEN ! use right gradient
      dUdUvol_minus(iVar,iVar,0,2) = 1. + FV_sdx(1) * FV_dx_L
      dUdUvol_minus(iVar,iVar,0,3) = 0. - FV_sdx(1) * FV_dx_L
    ELSE
      CALL Abort(__STAMP__,'Slopes do not match with minmod in preconditioner!')
    END IF
    ! derivatives of plus value at minus interface, if reconstruction has been done with s_L_minus
    IF(ABS(UPrim_minus_nb(iVar)-(URec_extended(iVar,-1)+s_L_minus(iVar)*FV_dx_R_nb)).LE.1E-12)THEN 
      dUdUvol_minus(iVar,iVar,-1,3) = 0. + FV_dx_R_nb * ( FV_sdx(0))
    END IF
    ! plus side of element
    ! derivatives of plus values at plus interface
    IF(s_lim_plus(iVar).EQ.0.)THEN ! first order
      dUdUvol_plus(iVar,iVar,PP_N,2) = 1.
    ELSEIF(s_lim_plus(iVar).EQ.s_L_plus(iVar))THEN ! use left gradient
      dUdUvol_plus(iVar,iVar,PP_N,1) = 0. - FV_sdx(PP_N  ) * FV_dx_R
      dUdUvol_plus(iVar,iVar,PP_N,2) = 1. + FV_sdx(PP_N  ) * FV_dx_R
    ELSEIF(s_lim_plus(iVar).EQ.s_R_plus(iVar))THEN ! use right gradient
      dUdUvol_plus(iVar,iVar,PP_N,2) = 1. - FV_sdx(PP_N+1) * FV_dx_R
    ELSE
      CALL Abort(__STAMP__,'Slopes do not match with minmod in preconditioner!')
    END IF
    ! derivatives of minus values at plus interface, if reconstruction has been done with s_R_plus
    IF(ABS(UPrim_plus_nb(iVar)-(URec_extended(iVar,PP_N+1)-s_R_plus(iVar)*FV_dx_L_nb)).LE.1E-12)THEN 
      dUdUvol_plus( iVar,iVar,PP_N+1,1) = 0. - FV_dx_L_nb * (-FV_sdx(PP_N+1))
    END IF
  END DO !iVar
CASE(9) ! Central
  DO iVar=1,PP_nVar
    ! minus side of interface
    ! derivatives of minus value at minus interface
    dUdUvol_minus(iVar,iVar,0     ,2) = 1. - FV_dx_L    * 0.5*( FV_sdx(0) - FV_sdx(1))
    dUdUvol_minus(iVar,iVar,0     ,3) = 0. - FV_dx_L    * 0.5*( FV_sdx(1))
    ! derivatives of plus value at minus interface
    dUdUvol_minus(iVar,iVar,-1    ,3) = 0. + FV_dx_R_nb * 0.5*( FV_sdx(0))

    ! plus side of interface
    ! derivatives of plus values at plus interface
    dUdUvol_plus( iVar,iVar,PP_N  ,1) = 0. + FV_dx_R    * 0.5*(-FV_sdx(PP_N))
    dUdUvol_plus( iVar,iVar,PP_N  ,2) = 1. + FV_dx_R    * 0.5*( FV_sdx(PP_N)-FV_sdx(PP_N+1))
    ! derivatives of minus values at plus interface
    dUdUvol_plus( iVar,iVar,PP_N+1,1) = 0. - FV_dx_L_nb * 0.5*(-FV_sdx(PP_N+1))
  END DO !iVar
CASE DEFAULT 
  CALL Abort(__STAMP__,'No preconditioner for chosen limiter implemented!')
END SELECT

! multiply: dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim
DO ind=2,3
  dUdUvol_minus(:,:,0,ind) = MATMUL(Jac_ConsPrim_minus(:,:,0),dUdUvol_minus(:,:,0,ind))
END DO
dUdUvol_minus(:,:,-1,3) = MATMUL(Jac_ConsPrim_minus(:,:,-1),dUdUvol_minus(:,:,-1,3))
DO ind=1,2
  dUdUvol_plus( :,:,PP_N,ind) = MATMUL(Jac_ConsPrim_plus( :,:,PP_N),dUdUvol_plus( :,:,PP_N,ind))
END DO
dUdUvol_plus( :,:,PP_N+1,1) = MATMUL(Jac_ConsPrim_plus( :,:,PP_N+1),dUdUvol_plus( :,:,PP_N+1,1))

! multiply: (dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim) * dUvol_prim/dUvol
matrix=dUdUvol_plus(:,:,PP_N,1);   vector=Jac_PrimCons_plus(:,:,PP_N-1)
dUdUvol_plus( :,:,PP_N,1)   = MATMUL(matrix,vector)
matrix=dUdUvol_plus(:,:,PP_N,2);   vector=Jac_PrimCons_plus(:,:,PP_N) 
dUdUvol_plus( :,:,PP_N,2)   = MATMUL(matrix,vector)
matrix=dUdUvol_plus(:,:,PP_N+1,1); vector=Jac_PrimCons_plus(:,:,PP_N)
dUdUvol_plus( :,:,PP_N+1,1) = MATMUL(matrix,vector)

matrix=dUdUvol_minus(:,:,-1,3); vector=Jac_PrimCons_minus(:,:,0)
dUdUvol_minus(:,:,-1,3) = MATMUL(matrix,vector)
matrix=dUdUvol_minus(:,:,0,2);  vector=Jac_PrimCons_minus(:,:,0)
dUdUvol_minus(:,:,0,2)  = MATMUL(matrix,vector)
matrix=dUdUvol_minus(:,:,0,3);  vector=Jac_PrimCons_minus(:,:,1)
dUdUvol_minus(:,:,0,3)  = MATMUL(matrix,vector)

END SUBROUTINE FV_Reconstruction_Derivative_Surf

#if PARABOLIC
!===================================================================================================================================
!> Computes the Volume gradient Jacobian of the reconstruction procedure dQprim/dUprim (Q= Grad U) only w.r.t. the volume DOFs of
!> the current element.
!> As the reconstruction procedure is the same for all primitive variables and is done independently for each variable the Jacobian
!> has diagonal shape with the same coefficient on the diagonal. Therefore, we calculate only this scalar value and give it back.
!> The reconstruction procedure for the gradients takes the central limiter. Therefore the gradient at a DOF (i,j,k) is depending on
!> the left and right neighbour and the DOF itself.
!===================================================================================================================================
SUBROUTINE JacFVGradients_Vol(dir,iElem,Jac_reconstruct)
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars            ,ONLY: FV_sdx_XI,FV_sdx_ETA,FV_Metrics_fTilde_sJ,FV_Metrics_gTilde_sJ
USE MOD_Jac_Ex_Vars        ,ONLY: FV_sdx_XI_extended,FV_sdx_ETA_extended
#if PP_dim == 3     
USE MOD_FV_Vars            ,ONLY: FV_sdx_ZETA,FV_Metrics_hTilde_sJ
USE MOD_Jac_Ex_Vars        ,ONLY: FV_sdx_ZETA_extended
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                           :: iElem   !< current element index
INTEGER,INTENT(IN)                           :: dir     !< current physical direction
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                             :: Jac_reconstruct(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim)
                                                !< Jacobian of volume gradients in direction dir w.r.t. primitive volume solution
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: i,j,k,ll
!===================================================================================================================================
Jac_reconstruct = 0.
DO k=0,PP_NZ
  DO j=0,PP_N
    DO i=0,PP_N
      ! Contribution by gradient reconstruction procedure with central gradient (i,j,k) -> gradient, (ll) -> state 
      IF(i.GT.0)THEN
        Jac_reconstruct(i,j,k,i-1,1) = 0.5*(-FV_sdx_XI(j,k,i,iElem))
      END IF                     
      Jac_reconstruct(  i,j,k,i  ,1) = 0.5*FV_sdx_XI_extended(j,k,i,iElem) - 0.5*FV_sdx_XI_extended(j,k,i+1,iElem)
      IF(i.LT.PP_N)THEN          
        Jac_reconstruct(i,j,k,i+1,1) = 0.5*( FV_sdx_XI(j,k,i+1,iElem))
      END IF
      IF(j.GT.0)THEN
        Jac_reconstruct(i,j,k,j-1,2) = 0.5*(-FV_sdx_ETA(i,k,j,iElem))
      END IF                     
      Jac_reconstruct(  i,j,k,j  ,2) = 0.5*FV_sdx_ETA_extended(i,k,j,iElem) - 0.5*FV_sdx_ETA_extended(i,k,j+1,iElem)
      IF(j.LT.PP_N)THEN          
        Jac_reconstruct(i,j,k,j+1,2) = 0.5*( FV_sdx_ETA(i,k,j+1,iElem))
      END IF
#if PP_dim==3
      IF(k.GT.0)THEN
        Jac_reconstruct(i,j,k,k-1,3) = 0.5*(-FV_sdx_ZETA(i,j,k,iElem))
      END IF                     
      Jac_reconstruct(  i,j,k,k  ,3) = 0.5*FV_sdx_ZETA_extended(i,j,k,iElem) - 0.5*FV_sdx_ZETA_extended(i,j,k+1,iElem)
      IF(k.LT.PP_N)THEN          
        Jac_reconstruct(i,j,k,k+1,3) = 0.5*( FV_sdx_ZETA(i,j,k+1,iElem))
      END IF
#endif
      ! Contribution by lifting volume integral
      DO ll=0,PP_N
        Jac_reconstruct(i,j,k,ll,1) = Jac_reconstruct(i,j,k,ll,1)*FV_Metrics_fTilde_sJ(dir,ll,j,k,iElem)
        Jac_reconstruct(i,j,k,ll,2) = Jac_reconstruct(i,j,k,ll,2)*FV_Metrics_gTilde_sJ(dir,i,ll,k,iElem)
#if PP_dim==3
        Jac_reconstruct(i,j,k,ll,3) = Jac_reconstruct(i,j,k,ll,3)*FV_Metrics_hTilde_sJ(dir,i,j,ll,iElem)
#endif
      END DO ! ll
    END DO !i
  END DO !j
END DO !k
END SUBROUTINE JacFVGradients_Vol

!===================================================================================================================================
!> Computes the dependency of FV gradient of neighbouring element at the interface on the volume DOFs of the considered element.
!> Again, as the reconstruction of the gradients is independent of the variable (iVar) itself, the output of this routine is a
!> scalar.
!> For the second order reconstruction with central limiter the gradient of the neighbouring element only depends on the outer layer
!> of the current element.
!> |...:...:...: g | x :...:...: x | g :...:...:...| => gradients g only depent on the outer layer DOF g of the current element
!> |...:...:...:...|...:...:...:...|...:...:...:...|
!>    nbElem_minus       iElem        nbElem_plus
!===================================================================================================================================
SUBROUTINE JacFVGradients_nb(dir,iElem,dQ_dUVolOuter)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars                 ,ONLY: ElemToSide,nBCSides
#if FV_ENABLED
USE MOD_Jac_Ex_Vars               ,ONLY: FV_sdx_XI_extended,FV_sdx_ETA_extended
USE MOD_FV_Vars                   ,ONLY: FV_Metrics_fTilde_sJ,FV_Metrics_gTilde_sJ
#if PP_dim == 3
USE MOD_Jac_Ex_Vars               ,ONLY: FV_sdx_ZETA_extended
USE MOD_FV_Vars                   ,ONLY: FV_Metrics_hTilde_sJ
#endif
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: dir   !< considered direction (1,2,3)
INTEGER,INTENT(IN) :: iElem !< considered element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#if PP_dim == 3
REAL,INTENT(OUT) :: dQ_dUVolOuter(0:PP_N,0:PP_NZ,6,0:PP_N)!< Jacobian of surface gradients of the neighbour element w.r.t. primitive
                                                          !< solution of current element
#else
REAL,INTENT(OUT) :: dQ_dUVolOuter(0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iLocSide,SideID,i,j,k
!===================================================================================================================================
dQ_dUVolOuter=0.

#if PP_dim == 3
DO iLocSide=1,6
#else    
DO iLocSide=2,5
#endif    
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
  IF(SideID.LE.nBCSides) CYCLE  !for boundary conditions, dQ_dUVolOuter=0.

  SELECT CASE(iLocSide)
  CASE(XI_MINUS,XI_PLUS)
    DO k=0,PP_NZ
      DO j=0,PP_N
        IF(iLocSide.EQ.XI_MINUS)THEN
          dQ_dUVolOuter(j,k,XI_MINUS,   0) = 0.5*FV_sdx_XI_extended(j,k,     0,iElem) * FV_Metrics_fTilde_sJ(dir,   0,j,k,iElem)
        ELSE
          dQ_dUVolOuter(j,k,XI_PLUS ,PP_N) = 0.5*FV_sdx_XI_extended(j,k,PP_N+1,iElem) * FV_Metrics_fTilde_sJ(dir,PP_N,j,k,iElem)
        END IF
      END DO !j
    END DO !k
  CASE(ETA_MINUS,ETA_PLUS)
    DO k=0,PP_NZ
      DO i=0,PP_N
        IF(iLocSide.EQ.ETA_MINUS)THEN
          dQ_dUVolOuter(i,k,ETA_MINUS,   0) = 0.5*FV_sdx_ETA_extended(i,k,     0,iElem) * FV_Metrics_gTilde_sJ(dir,i,   0,k,iElem)
        ELSE
          dQ_dUVolOuter(i,k,ETA_PLUS ,PP_N) = 0.5*FV_sdx_ETA_extended(i,k,PP_N+1,iElem) * FV_Metrics_gTilde_sJ(dir,i,PP_N,k,iElem)
        END IF
      END DO !i
    END DO !k
#if PP_dim==3
  CASE(ZETA_MINUS,ZETA_PLUS)
    DO j=0,PP_N
      DO i=0,PP_N
        IF(iLocSide.EQ.ZETA_MINUS)THEN
          dQ_dUVolOuter(i,j,ZETA_MINUS,   0) = 0.5*FV_sdx_ZETA_extended(i,j,     0,iElem) * FV_Metrics_hTilde_sJ(dir,i,j,   0,iElem)
        ELSE
          dQ_dUVolOuter(i,j,ZETA_PLUS ,PP_N) = 0.5*FV_sdx_ZETA_extended(i,j,PP_N+1,iElem) * FV_Metrics_hTilde_sJ(dir,i,j,PP_N,iElem)
        END IF
      END DO !i
    END DO !j
#endif
  END SELECT
END DO !iLocSide
END SUBROUTINE JacFVGradients_nb
#endif /*PARABOLIC*/
#endif /*FV_ENABLED && FV_RECONSTRUCT*/

END MODULE MOD_Jac_Reconstruction
