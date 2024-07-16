!!=================================================================================================================================
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
#if FV_ENABLED
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Contains the main parts of the second order reconstruction of the FV subcells.
!==================================================================================================================================
MODULE MOD_FV_Reconstruction
! MODULES
IMPLICIT NONE
PRIVATE

#if FV_RECONSTRUCT
INTERFACE FV_PrepareSurfGradient
  MODULE PROCEDURE FV_PrepareSurfGradient
END INTERFACE

INTERFACE FV_CalcGradients
  MODULE PROCEDURE FV_CalcGradients
END INTERFACE

INTERFACE FV_SurfCalcGradients
  MODULE PROCEDURE FV_SurfCalcGradients
END INTERFACE

INTERFACE FV_SurfCalcGradients_BC
  MODULE PROCEDURE FV_SurfCalcGradients_BC
END INTERFACE

#if VOLINT_VISC
INTERFACE FV_SurfCalcGradients_Parabolic
  MODULE PROCEDURE FV_SurfCalcGradients_Parabolic
END INTERFACE
#endif

PUBLIC::FV_PrepareSurfGradient
PUBLIC::FV_SurfCalcGradients,FV_SurfCalcGradients_BC,FV_CalcGradients
#if VOLINT_VISC
PUBLIC::FV_SurfCalcGradients_Parabolic
#endif
#endif /* FV_RECONSTRUCT */
!==================================================================================================================================

CONTAINS

#if FV_RECONSTRUCT

!==================================================================================================================================
!> Fills FV_multi_master and FV_multi_slave arrays.
!> These arrays contain:
!>  - for DG: solution at the nodes of the first inner layer next to the DG element interface
!>  - for FV: slopes in normal direction to the DG element interfaces between first and second layer
!==================================================================================================================================
PPURE SUBROUTINE FV_PrepareSurfGradient(UPrim,FV_multi_master,FV_multi_slave,doMPIsides)
! MODULES                                                                                                                          !
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars ,ONLY: firstMPISide_YOUR,lastMPISide_YOUR,firstInnerSide,lastMPISide_MINE
USE MOD_Mesh_Vars ,ONLY: firstMortarMPISide,lastMortarMPISide,firstBCSide
USE MOD_Mesh_Vars ,ONLY: S2V,S2V2,SideToElem
USE MOD_Mesh_Vars ,ONLY: nElems,nSides
USE MOD_FV_Vars   ,ONLY: FV_Elems,FV_sdx_XI,FV_sdx_ETA
#if PP_dim == 3
USE MOD_FV_Vars   ,ONLY: FV_sdx_ZETA
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
! real_in, real_out, real_out, log_in
REAL,INTENT(IN)    :: UPrim          (PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< primitive volume solution
REAL,INTENT(OUT)   :: FV_multi_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)        !< DG: solution at first inner layer,
                                                                                  !< FV: slope between first and second layer
REAL,INTENT(OUT)   :: FV_multi_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)        !< DG: solution at first inner layer,
                                                                                  !< FV: slope between first and second layer
LOGICAL,INTENT(IN) :: doMPIsides                                       !< =.TRUE. only MPI sides are filled, =.FALSE. inner sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,i,j,k,locSideID,ElemID,SideID,flip,firstSideID,lastSideID,ijk0(3),ijk(2)
REAL               :: tmp(PP_nVarPrim,0:PP_N,0:PP_NZ)
!==================================================================================================================================

! First process the Slave sides
IF(doMPISides)THEN
  ! YOUR are filled
  firstSideID = firstMPISide_YOUR
  lastSideID  = lastMPISide_YOUR
ELSE
  ! (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstInnerSide
  lastSideID  = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ! neighbor side !ElemID,locSideID and flip =-1 if not existing
  ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip      = SideToElem(S2E_FLIP,SideID)
  IF (ElemID.EQ.-1) CYCLE

  IF (FV_Elems(ElemID).EQ.0) THEN ! DG Element
    ! for DG elements copy solution of first inner layer next to the interface
    DO q=0,PP_NZ; DO p=0,PP_N
      ijk0=S2V(:,0,p,q,flip,locSideID)
      FV_multi_slave(:,p,q,SideID) = UPrim(:,ijk0(1),ijk0(2),ijk0(3),ElemID)
    END DO; END DO
  ELSE ! FV Element
    ! for FV elements calculate slopes between first and second layer
    !    TODO: USE S2V directly in each of the following 6 cases and not S2V2 below!!!
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      DO k=0,PP_NZ; DO j=0,PP_N
        tmp(:,j,k) =  (UPrim(:,1,j,k,ElemID) - UPrim(:,0,j,k,ElemID)) * FV_sdx_XI(j,k,1,ElemID)
      END DO; END DO
    CASE(ETA_MINUS)
      DO k=0,PP_NZ; DO i=0,PP_N
        tmp(:,i,k) =  (UPrim(:,i,1,k,ElemID) - UPrim(:,i,0,k,ElemID)) * FV_sdx_ETA(i,k,1,ElemID)
      END DO; END DO
#if PP_dim == 3
    CASE(ZETA_MINUS)
      DO j=0,PP_NZ; DO i=0,PP_N
        tmp(:,i,j) =  (UPrim(:,i,j,1,ElemID) - UPrim(:,i,j,0,ElemID)) * FV_sdx_ZETA(i,j,1,ElemID)
      END DO; END DO
#endif
    CASE(XI_PLUS)
      DO k=0,PP_NZ; DO j=0,PP_N
        tmp(:,j,k) =  (UPrim(:,PP_N-1,j,k,ElemID) - UPrim(:,PP_N,j,k,ElemID)) * FV_sdx_XI(j,k,PP_N,ElemID)
      END DO; END DO
    CASE(ETA_PLUS)
      DO k=0,PP_NZ; DO i=0,PP_N
        tmp(:,i,k) =  (UPrim(:,i,PP_N-1,k,ElemID) - UPrim(:,i,PP_N,k,ElemID)) * FV_sdx_ETA(i,k,PP_N,ElemID)
      END DO; END DO
#if PP_dim == 3
    CASE(ZETA_PLUS)
      DO j=0,PP_NZ; DO i=0,PP_N
        tmp(:,i,j) =  (UPrim(:,i,j,PP_N-1,ElemID) - UPrim(:,i,j,PP_N,ElemID)) * FV_sdx_ZETA(i,j,PP_N,ElemID)
      END DO; END DO
#endif
    END SELECT

    DO q=0,PP_NZ; DO p=0,PP_N
      ijk=S2V2(:,p,q,flip,locSideID)
      FV_multi_slave(:,p,q,SideID) = tmp(:,ijk(1),ijk(2))
    END DO; END DO
  END IF
END DO

! Second process master sides, U_master is always MINE
! master side, flip=0
IF(doMPISides)THEN
  ! only MPI mortars
  firstSideID = firstMortarMPISide
   lastSideID =  lastMortarMPISide
ELSE
  ! BCSides, (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstBCSide
   lastSideID = lastMPISide_MINE
END IF
DO SideID=firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  IF (ElemID.EQ.-1) CYCLE

  IF (FV_Elems(ElemID).EQ.0) THEN ! DG Element
    ! for DG elements copy solution of first inner layer next to the interface
    DO q=0,PP_NZ; DO p=0,PP_N
      ijk0=S2V(:,0,p,q,0,locSideID)
      FV_multi_master(:,p,q,SideID) = UPrim(:,ijk0(1),ijk0(2),ijk0(3),ElemID)
    END DO; END DO
  ELSE ! FV Element
    ! for FV elements calculate slopes between first and second layer
    !    TODO: USE S2V directly in each of the following 6 cases and not S2V2 below!!!
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      DO k=0,PP_NZ; DO j=0,PP_N
        tmp(:,j,k) =  (UPrim(:,1,j,k,ElemID) - UPrim(:,0,j,k,ElemID)) * FV_sdx_XI(j,k,1,ElemID)
      END DO; END DO
    CASE(ETA_MINUS)
      DO k=0,PP_NZ; DO i=0,PP_N
        tmp(:,i,k) =  (UPrim(:,i,1,k,ElemID) - UPrim(:,i,0,k,ElemID)) * FV_sdx_ETA(i,k,1,ElemID)
      END DO; END DO
#if PP_dim == 3
    CASE(ZETA_MINUS)
      DO j=0,PP_NZ; DO i=0,PP_N
        tmp(:,i,j) =  (UPrim(:,i,j,1,ElemID) - UPrim(:,i,j,0,ElemID)) * FV_sdx_ZETA(i,j,1,ElemID)
      END DO; END DO
#endif
    CASE(XI_PLUS)
      DO k=0,PP_NZ; DO j=0,PP_N
        tmp(:,j,k) =  (UPrim(:,PP_N-1,j,k,ElemID) - UPrim(:,PP_N,j,k,ElemID)) * FV_sdx_XI(j,k,PP_N,ElemID)
      END DO; END DO
    CASE(ETA_PLUS)
      DO k=0,PP_NZ; DO i=0,PP_N
        tmp(:,i,k) =  (UPrim(:,i,PP_N-1,k,ElemID) - UPrim(:,i,PP_N,k,ElemID)) * FV_sdx_ETA(i,k,PP_N,ElemID)
      END DO; END DO
#if PP_dim == 3
    CASE(ZETA_PLUS)
      DO j=0,PP_NZ; DO i=0,PP_N
        tmp(:,i,j) =  (UPrim(:,i,j,PP_N-1,ElemID) - UPrim(:,i,j,PP_N,ElemID)) * FV_sdx_ZETA(i,j,PP_N,ElemID)
      END DO; END DO
#endif
    END SELECT

    DO q=0,PP_NZ; DO p=0,PP_N
      ijk=S2V2(:,p,q,0,locSideID)
      FV_multi_master(:,p,q,SideID) = tmp(:,ijk(1),ijk(2))
    END DO; END DO
  END IF
END DO !SideID
END SUBROUTINE FV_PrepareSurfGradient


!==================================================================================================================================
!> Calculate slopes across DG element interfaces (unlimited). Uses UPrim_master/UPrim_slave as well as
!> FV_multi_master/FV_multi_slave. They are limited in the FV_ProlongToDGFace routine.
!==================================================================================================================================
SUBROUTINE FV_SurfCalcGradients(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,&
        FV_surf_gradU,doMPIsides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars   ,ONLY: FV_sdx_Face,FV_Elems_Sum
USE MOD_Mesh_Vars ,ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars ,ONLY: nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)    :: UPrim_master   (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive master solution (without reconstruction)
REAL,INTENT(IN)    :: UPrim_slave    (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive slave solution (without reconstruction)
REAL,INTENT(IN)    :: FV_multi_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< DG solution at the first inner node next to the
                                                                           !< interface of the master element
REAL,INTENT(IN)    :: FV_multi_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< DG solution at the first inner node next to the
                                                                           !< interface of the slave element
LOGICAL,INTENT(IN) :: doMPIsides                                           !< =.TRUE. only MPI sides are filled, =.FALSE. inner sides
REAL,INTENT(OUT)   :: FV_surf_gradU  (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< slope over the interface
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: firstSideID,lastSideID,SideID,p,q
!==================================================================================================================================

IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
  lastSideID  = lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstInnerSide
  lastSideID  = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
  SELECT CASE(FV_Elems_Sum(SideID))
  CASE(0) ! both DG
    CYCLE
  CASE(1) ! master=FV, slave=DG
    DO q=0,PP_NZ; DO p=0,PP_N
      ! use FV_multi_slave (see FV_PrepareSurfGradient)
      FV_surf_gradU(:,p,q,SideID) = (UPrim_master(:,p,q,SideID) - FV_multi_slave(:,p,q,SideID)) * FV_sdx_Face(p,q,1,SideID)
    END DO; END DO ! p,q=0,PP_N
  CASE(2) ! master=DG, slave=FV
    DO q=0,PP_NZ; DO p=0,PP_N
      ! use FV_multi_master (see FV_PrepareSurfGradient)
      FV_surf_gradU(:,p,q,SideID) = (FV_multi_master(:,p,q,SideID) - UPrim_slave(:,p,q,SideID)) * FV_sdx_Face(p,q,2,SideID)
    END DO; END DO ! p,q=0,PP_N
  CASE(3) ! both FV
    DO q=0,PP_NZ; DO p=0,PP_N
      FV_surf_gradU(:,p,q,SideID) = (UPrim_master(:,p,q,SideID) - UPrim_slave(:,p,q,SideID)) * FV_sdx_Face(p,q,3,SideID)
    END DO; END DO ! p,q=0,PP_N
  CASE DEFAULT
    CALL Abort(__STAMP__, "FV_Elems is not 0 or 1 somewhere!")
  END SELECT

END DO
END SUBROUTINE FV_SurfCalcGradients


!==================================================================================================================================
!> Get slopes at the boundary conditions.
!==================================================================================================================================
SUBROUTINE FV_SurfCalcGradients_BC(UPrim_master,FV_surf_gradU,t)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars       ,ONLY: firstBCSide,lastBCSide,nSides
USE MOD_Mesh_Vars       ,ONLY: NormVec,TangVec1,TangVec2,Face_xGP
USE MOD_FV_Vars         ,ONLY: FV_Elems_master,FV_sdx_Face
USE MOD_GetBoundaryFlux ,ONLY: GetBoundaryFVgradient
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)    :: t                                                  !< physical time
REAL,INTENT(IN)    :: UPrim_master( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution
REAL,INTENT(INOUT) :: FV_surf_gradU(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< slope over the BC interface
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: SideID
!==================================================================================================================================
DO SideID=firstBCSide,lastBCSide
  IF (FV_Elems_master(SideID).EQ.0) CYCLE
  CALL GetBoundaryFVgradient(SideID,t,FV_surf_gradU(:,:,:,SideID),&
                                       UPrim_master(:,:,:,SideID),&
                                          NormVec(:,:,:,1,SideID),&
                                         TangVec1(:,:,:,1,SideID),&
                                         TangVec2(:,:,:,1,SideID),&
                                         Face_xGP(:,:,:,1,SideID),&
                                        FV_sdx_Face(:,:,:,SideID))
END DO
END SUBROUTINE FV_SurfCalcGradients_BC

!==================================================================================================================================
!> Calculate slopes in the inside and limit them using TVD limiters.
!> It is important to use physical distances (and not reference distances) to calculate the slope, since otherwise the
!> slope limiter can not be applied. (Scenario: inner-cell-stretching)
!> Additionally build central limited slopes for the computation of gradients used for the viscous fluxes.
!==================================================================================================================================
SUBROUTINE FV_CalcGradients(UPrim,FV_surf_gradU,gradUxi,gradUeta,gradUzeta &
#if VOLINT_VISC
        ,gradUxi_central,gradUeta_central,gradUzeta_central &
#endif
    )
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars        ,ONLY: FV_sdx_XI,FV_sdx_ETA,FV_Elems
USE MOD_FV_Limiter     ,ONLY: FV_Limiter
USE MOD_Mesh_Vars      ,ONLY: nElems,nSides
#if PP_dim == 3
USE MOD_FV_Vars        ,ONLY: FV_sdx_ZETA
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: UPrim            (PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< primitive volume solution
REAL,INTENT(IN)  :: FV_surf_gradU    (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides)      !< slopes over element interfaces
REAL,INTENT(OUT) :: gradUxi          (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems) !< physical slope in   xi-direction (TVD limited)
REAL,INTENT(OUT) :: gradUeta         (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems) !< physical slope in  eta-direction (TVD limited)
REAL,INTENT(OUT) :: gradUzeta        (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N,nElems) !< physical slope in zeta-direction (TVD limited)
#if VOLINT_VISC
REAL,INTENT(OUT) :: gradUxi_central  (PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< physical slope in   xi-direction (mean value)
REAL,INTENT(OUT) :: gradUeta_central (PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< physical slope in  eta-direction (mean value)
REAL,INTENT(OUT) :: gradUzeta_central(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< physical slope in zeta-direction (mean value)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUxi_tmp
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUeta_tmp
#if PP_dim == 3
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1) :: gradUzeta_tmp
#endif
INTEGER                                             :: iElem,l,iVar,p,q
#if VOLINT_VISC
INTEGER                                             :: i,j,k
#endif
!==================================================================================================================================
DO iElem=1,nElems
  IF (FV_Elems(iElem).EQ.0) CYCLE ! DG Element
  ! Strategy:
  ! 1. compute slopes between subcells in all 3 directions and store in temporary arrays
  !    Attention concerning the storage order: last index corresponds to the respective direction (xi/eta/zeta)!
  ! 2. copy the slopes over DG interfaces to the 0.th / PP_N+1.th index of the temporary arrays
  ! 3. limit the slopes from the temporary array and store in gradUxi/eta/zeta

  ! 1. gradients of inner subcells
  DO l=1,PP_N
    DO iVar=1,PP_nVarPrim
      gradUxi_tmp  (iVar,:,:,l) = (UPrim(iVar,l,:,:,iElem) - UPrim(iVar,l-1,:,:,iElem)) * FV_sdx_XI  (:,:,l,iElem)
      gradUeta_tmp (iVar,:,:,l) = (UPrim(iVar,:,l,:,iElem) - UPrim(iVar,:,l-1,:,iElem)) * FV_sdx_ETA (:,:,l,iElem)
#if PP_dim == 3
      gradUzeta_tmp(iVar,:,:,l) = (UPrim(iVar,:,:,l,iElem) - UPrim(iVar,:,:,l-1,iElem)) * FV_sdx_ZETA(:,:,l,iElem)
#endif
    END DO
  END DO
#if FV_ENABLED != 3
  ! 2. gradients of subcells at DG interface
  ! xi direction
  CALL CopySurfaceToVolume(FV_surf_gradU,gradUxi_tmp,iElem,XI_MINUS,0)
  CALL CopySurfaceToVolume(FV_surf_gradU,gradUxi_tmp,iElem,XI_PLUS ,PP_N+1)
  ! eta direction
  CALL CopySurfaceToVolume(FV_surf_gradU,gradUeta_tmp,iElem,ETA_MINUS,0)
  CALL CopySurfaceToVolume(FV_surf_gradU,gradUeta_tmp,iElem,ETA_PLUS ,PP_N+1)
  ! zeta direction
#if PP_dim == 3
  CALL CopySurfaceToVolume(FV_surf_gradU,gradUzeta_tmp,iElem,ZETA_MINUS,0)
  CALL CopySurfaceToVolume(FV_surf_gradU,gradUzeta_tmp,iElem,ZETA_PLUS ,PP_N+1)
#endif
#endif /*FV_ENABLED != 3*/

  ! 3. limit
#if FV_ENABLED != 3
  DO l=0,PP_N
  ! 3. Outer subcells are not limited but retain zero slope due to conservativity
#else
  DO l=1,PP_N-1
#endif /*FV_ENABLED != 3*/
   DO q=0,PP_NZ; DO p=0,PP_N
       CALL FV_Limiter(gradUxi_tmp  (:,p,q,l),gradUxi_tmp  (:,p,q,l+1),gradUxi  (:,p,q,l,iElem))
       CALL FV_Limiter(gradUeta_tmp (:,p,q,l),gradUeta_tmp (:,p,q,l+1),gradUeta (:,p,q,l,iElem))
#if PP_dim == 3
       CALL FV_Limiter(gradUzeta_tmp(:,p,q,l),gradUzeta_tmp(:,p,q,l+1),gradUzeta(:,p,q,l,iElem))
#endif
   END DO; END DO ! q, p
  END DO ! l
#if VOLINT_VISC
  CALL FV_CalcGradients_Parabolic(iElem,gradUxi_tmp,gradUeta_tmp &
#if PP_dim==3
                                       ,gradUzeta_tmp            &
#endif
                                 )
  CALL FV_PrepareSurfGradient_Parabolic(iElem,gradUxi_tmp,gradUeta_tmp &
#if PP_dim==3
                                             ,gradUzeta_tmp            &
#endif
                                       )
  ! limit with central limiter for viscous fluxes
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    gradUxi_central  (:,i,j,k,iElem) = 0.5*(gradUxi_tmp  (PRIM_LIFT,j,k,i)+gradUxi_tmp  (PRIM_LIFT,j,k,i+1))
    gradUeta_central (:,i,j,k,iElem) = 0.5*(gradUeta_tmp (PRIM_LIFT,i,k,j)+gradUeta_tmp (PRIM_LIFT,i,k,j+1))
#if PP_dim == 3
    gradUzeta_central(:,i,j,k,iElem) = 0.5*(gradUzeta_tmp(PRIM_LIFT,i,j,k)+gradUzeta_tmp(PRIM_LIFT,i,j,k+1))
#endif
  END DO; END DO; END DO! i,j,k=0,PP_N
#endif /*VOLINT_VISC*/
END DO

END SUBROUTINE FV_CalcGradients

!==================================================================================================================================
!> Copy surface data at the face specified by dir to the volume.
!==================================================================================================================================
PPURE SUBROUTINE CopySurfaceToVolume(surface,volume,iElem,dir,l)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars      ,ONLY: S2V2,nSides,ElemToSide
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)    :: surface(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< master surface data
INTEGER,INTENT(IN) :: iElem                                        !< element index
INTEGER,INTENT(IN) :: dir                                          !< face number of the element (1..6)
INTEGER,INTENT(IN) :: l                                            !< volume index where to store the volume data
REAL,INTENT(INOUT) :: volume(PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1)  !< output array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,ijk(2),SideID,flip
!==================================================================================================================================
SideID=ElemToSide(E2S_SIDE_ID,dir,iElem)
flip  =ElemToSide(E2S_FLIP,   dir,iElem)

SELECT CASE(dir)
CASE(XI_MINUS,ETA_MINUS,ZETA_MINUS)
  DO q=0,PP_NZ; DO p=0,PP_N
    ijk = S2V2(:,p,q,flip,dir)
    volume(:,ijk(1),ijk(2),l) = surface(:,p,q,SideID)
  END DO; END DO ! p,q=0,PP_N
CASE(XI_PLUS,ETA_PLUS,ZETA_PLUS)
  DO q=0,PP_NZ; DO p=0,PP_N
    ijk = S2V2(:,p,q,flip,dir)
    volume(:,ijk(1),ijk(2),l) = -surface(:,p,q,SideID)
  END DO; END DO ! p,q=0,PP_N
END SELECT

IF (flip.GT.0) THEN
  volume(:,:,:,l) = -volume(:,:,:,l)
END IF
END SUBROUTINE CopySurfaceToVolume

#if VOLINT_VISC
!==================================================================================================================================
!> Calculate slopes across DG element interfaces (unlimited). Uses UPrim_master/UPrim_slave as well as
!> FV_multi_master/FV_multi_slave. They are limited in the FV_ProlongToDGFace routine.
!==================================================================================================================================
SUBROUTINE FV_SurfCalcGradients_Parabolic()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars ,ONLY: gradUx_master,gradUx_slave,gradUy_master,gradUy_slave
USE MOD_FV_Vars      ,ONLY: FV_Elems_Sum,FV_surf_gradU_master,FV_surf_gradU_slave
!USE MOD_Mesh_Vars    ,ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars    ,ONLY: nSides,NormVec,TangVec1,firstMortarInnerSide
#if PP_dim==3
USE MOD_Mesh_Vars    ,ONLY: TangVec2
USE MOD_Lifting_Vars ,ONLY: gradUz_master,gradUz_slave
#endif /* PP_dim==3 */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: firstSideID,lastSideID,SideID,p,q
!==================================================================================================================================

!IF(doMPISides)THEN
!  ! fill only flux for MINE MPISides
!  firstSideID = firstMPISide_MINE
!  lastSideID  = lastMPISide_MINE
!ELSE
!  ! fill only InnerSides
!  firstSideID = firstInnerSide
!  lastSideID  = lastInnerSide
!END IF

DO SideID = firstMortarInnerSide, nSides
  SELECT CASE(FV_Elems_Sum(SideID))
  CASE(0) ! both DG
    CYCLE
  CASE(1) ! master=FV, slave=DG
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      gradUx_master(:,p,q,SideID) =0.5*(NormVec (1,p,q,1,SideID)*FV_surf_gradU_master(:,1,p,q,SideID)&
                                       +TangVec1(1,p,q,1,SideID)*FV_surf_gradU_master(:,2,p,q,SideID)&
                                       +TangVec2(1,p,q,1,SideID)*FV_surf_gradU_master(:,3,p,q,SideID)&
                                       +gradUx_slave(:,p,q,SideID))
      gradUy_master(:,p,q,SideID) =0.5*(NormVec (2,p,q,1,SideID)*FV_surf_gradU_master(:,1,p,q,SideID)&
                                       +TangVec1(2,p,q,1,SideID)*FV_surf_gradU_master(:,2,p,q,SideID)&
                                       +TangVec2(2,p,q,1,SideID)*FV_surf_gradU_master(:,3,p,q,SideID)&
                                       +gradUy_slave(:,p,q,SideID))
      gradUz_master(:,p,q,SideID) =0.5*(NormVec (3,p,q,1,SideID)*FV_surf_gradU_master(:,1,p,q,SideID)&
                                       +TangVec1(3,p,q,1,SideID)*FV_surf_gradU_master(:,2,p,q,SideID)&
                                       +TangVec2(3,p,q,1,SideID)*FV_surf_gradU_master(:,3,p,q,SideID)&
                                       +gradUz_slave(:,p,q,SideID))
#else
      gradUx_master(:,p,q,SideID) =0.5*(NormVec (1,p,q,1,SideID)*FV_surf_gradU_master(:,1,p,q,SideID)&
                                       +TangVec1(1,p,q,1,SideID)*FV_surf_gradU_master(:,2,p,q,SideID)&
                                       +gradUx_slave(:,p,q,SideID))
      gradUy_master(:,p,q,SideID) =0.5*(NormVec (2,p,q,1,SideID)*FV_surf_gradU_master(:,1,p,q,SideID)&
                                       +TangVec1(2,p,q,1,SideID)*FV_surf_gradU_master(:,2,p,q,SideID)&
                                       +gradUy_slave(:,p,q,SideID))
#endif /* PP_dim==3 */
    END DO; END DO ! p,q=0,PP_N
  CASE(2) ! master=DG, slave=FV
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      gradUx_master(:,p,q,SideID) =0.5*(NormVec (1,p,q,1,SideID)*FV_surf_gradU_slave(:,1,p,q,SideID)&
                                       +TangVec1(1,p,q,1,SideID)*FV_surf_gradU_slave(:,2,p,q,SideID)&
                                       +TangVec2(1,p,q,1,SideID)*FV_surf_gradU_slave(:,3,p,q,SideID)&
                                       +gradUx_master(:,p,q,SideID))
      gradUy_master(:,p,q,SideID) =0.5*(NormVec (2,p,q,1,SideID)*FV_surf_gradU_slave(:,1,p,q,SideID)&
                                       +TangVec1(2,p,q,1,SideID)*FV_surf_gradU_slave(:,2,p,q,SideID)&
                                       +TangVec2(2,p,q,1,SideID)*FV_surf_gradU_slave(:,3,p,q,SideID)&
                                       +gradUy_master(:,p,q,SideID))
      gradUz_master(:,p,q,SideID) =0.5*(NormVec (3,p,q,1,SideID)*FV_surf_gradU_slave(:,1,p,q,SideID)&
                                       +TangVec1(3,p,q,1,SideID)*FV_surf_gradU_slave(:,2,p,q,SideID)&
                                       +TangVec2(3,p,q,1,SideID)*FV_surf_gradU_slave(:,3,p,q,SideID)&
                                       +gradUz_master(:,p,q,SideID))
#else /* PP_dim==3 */
      gradUx_master(:,p,q,SideID) =0.5*(NormVec (1,p,q,1,SideID)*FV_surf_gradU_slave(:,1,p,q,SideID)&
                                       +TangVec1(1,p,q,1,SideID)*FV_surf_gradU_slave(:,2,p,q,SideID)&
                                       +gradUx_master(:,p,q,SideID))
      gradUy_master(:,p,q,SideID) =0.5*(NormVec (2,p,q,1,SideID)*FV_surf_gradU_slave(:,1,p,q,SideID)&
                                       +TangVec1(2,p,q,1,SideID)*FV_surf_gradU_slave(:,2,p,q,SideID)&
                                       +gradUy_master(:,p,q,SideID))
#endif /* PP_dim==3 */
    END DO; END DO ! p,q=0,PP_N
  CASE(3) ! both FV
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_dim==3
      gradUx_master(:,p,q,SideID) =NormVec (1,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,1,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,1,p,q,SideID))&
                                  +TangVec1(1,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,2,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,2,p,q,SideID))&
                                  +TangVec2(1,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,3,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,3,p,q,SideID))
      gradUy_master(:,p,q,SideID) =NormVec (2,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,1,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,1,p,q,SideID))&
                                  +TangVec1(2,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,2,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,2,p,q,SideID))&
                                  +TangVec2(2,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,3,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,3,p,q,SideID))
      gradUz_master(:,p,q,SideID) =NormVec (3,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,1,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,1,p,q,SideID))&
                                  +TangVec1(3,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,2,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,2,p,q,SideID))&
                                  +TangVec2(3,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,3,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,3,p,q,SideID))
#else /* PP_dim==3 */
      gradUx_master(:,p,q,SideID) = NormVec(1,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,1,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,1,p,q,SideID))&
                                 + TangVec1(1,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,2,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,2,p,q,SideID))
      gradUy_master(:,p,q,SideID) = NormVec(2,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,1,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,1,p,q,SideID))&
                                 + TangVec1(2,p,q,1,SideID)*0.5*(FV_surf_gradU_master(:,2,p,q,SideID) &
                                                                +FV_surf_gradU_slave (:,2,p,q,SideID))
#endif /* PP_dim==3 */

      END DO; END DO ! p,q=0,PP_N
  CASE DEFAULT
    CALL Abort(__STAMP__, "FV_Elems is not 0, 1, 2 or 3 somewhere!")
  END SELECT
  gradUx_slave(:,:,:,SideID) = gradUx_master(:,:,:,SideID)
  gradUy_slave(:,:,:,SideID) = gradUy_master(:,:,:,SideID)
#if PP_dim==3
  gradUz_slave(:,:,:,SideID) = gradUz_master(:,:,:,SideID)
#endif /* PP_dim==3 */
END DO ! SideID = 1, nSides
END SUBROUTINE FV_SurfCalcGradients_Parabolic

!==================================================================================================================================
!> Calculate gradients used for the viscous fluxes on faces
!> The calculation of the gradient which are tangential to their respective side is done with the following 6 Point Stencil:
!>           __________________________________
!>          |          |           |           |
!>          |          |           |           |
!>          |(i-1,j+1)-|> (i,j+1) -|> (i+1,j+1)|
!>          |          |           |           |
!>          |__________|____=>_____|___________|
!>          |          |    =>     |           |
!>    j^    |          |           |           |
!>     |    | (i-1,j) -|>  (i,j)  -|> (i+1,j)  |
!>     |    |          |           |           |
!>     |    |__________|___________|___________|
!>     |----------> i
!> => : Averaged Gradient to calulate, where the averaging is a simple arithmetic mean between all 4 FD stencils
!> -> : Simple FD Gradient between to cells
!==================================================================================================================================
SUBROUTINE FV_CalcGradients_Parabolic(iElem,gradUxi_tmp,gradUeta_tmp  &
#if PP_dim==3
                                           ,gradUzeta_tmp             &
#endif
                                     )
  ! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars  ,ONLY: FV_Metrics_fTilde_sJ_xi,FV_Metrics_gTilde_sJ_xi
USE MOD_FV_Vars  ,ONLY: FV_Metrics_fTilde_sJ_eta,FV_Metrics_gTilde_sJ_eta
USE MOD_FV_Vars  ,ONLY: gradUx_xi  ,gradUy_xi
USE MOD_FV_Vars  ,ONLY: gradUx_eta ,gradUy_eta
#if (PP_dim==3)
USE MOD_FV_Vars  ,ONLY: FV_Metrics_fTilde_sJ_zeta,FV_Metrics_gTilde_sJ_zeta
USE MOD_FV_Vars  ,ONLY: FV_Metrics_hTilde_sJ_xi
USE MOD_FV_Vars  ,ONLY: FV_Metrics_hTilde_sJ_eta
USE MOD_FV_Vars  ,ONLY: FV_Metrics_hTilde_sJ_zeta
USE MOD_FV_Vars  ,ONLY: gradUx_zeta ,gradUy_zeta ,gradUz_zeta,gradUz_xi,gradUz_eta
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: iElem
REAL,INTENT(IN)     :: gradUxi_tmp  (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1)
REAL,INTENT(IN)     :: gradUeta_tmp (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1)
#if PP_dim==3
REAL,INTENT(IN)     :: gradUzeta_tmp(PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
REAL                :: gradU_central_face_xi    (PP_dim,PP_nVarLifting,0:PP_N-1,0:PP_N  ,0:PP_NZ  )
REAL                :: gradU_central_face_eta   (PP_dim,PP_nVarLifting,0:PP_N  ,0:PP_N-1,0:PP_NZ  )
REAL                :: gradUxi_central_face     (       PP_nVarLifting,0:PP_N+1,0:PP_N  ,0:PP_NZ  )
REAL                :: gradUeta_central_face    (       PP_nVarLifting,0:PP_N  ,0:PP_N+1,0:PP_NZ  )
#if (PP_dim==3)
REAL                :: gradU_central_face_zeta  (PP_dim,PP_nVarLifting,0:PP_N  ,0:PP_N  ,0:PP_NZ-1)
REAL                :: gradUzeta_central_face   (       PP_nVarLifting,0:PP_N  ,0:PP_N  ,0:PP_NZ+1)
#endif
!==================================================================================================================================
! NOTE: Main steps:
! 1. map gradients form GRAD to PRIM to LIFT on faces for viscous fluxes
! 2. Calculate the gradients in xi-, eta-, zeta-direction on all FV inner sides via aritmetic mean of four neighbour cells
!    (use the "unlimited" central gradients and not the "limited" gradients)
! 3. Calculate the gradients in x-,  y- ,  z-   direction on all FV inner sides with corresponding FV-metrics
!----------------------------------------------------------------------------------------------------------------------------------

! 1. map gradients on faces for viscous fluxes
DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N+1
  gradUxi_central_face(:,i,j,k)   = gradUxi_tmp(PRIM_LIFT,j,k,i)
END DO; END DO; END DO! i,j,k=0,PP_N
DO k=0,PP_NZ; DO j=0,PP_N+1; DO i=0,PP_N
  gradUeta_central_face(:,i,j,k)  = gradUeta_tmp(PRIM_LIFT,i,k,j)
END DO; END DO; END DO! i,j,k=0,PP_N
#if PP_dim == 3
DO k=0,PP_NZ+1; DO j=0,PP_N; DO i=0,PP_N
  gradUzeta_central_face(:,i,j,k) = gradUzeta_tmp(PRIM_LIFT,i,j,k)
END DO; END DO; END DO! i,j,k=0,PP_N
#endif

! 2. Calculate the gradients in xi-, eta-, zeta-direction on all FV inner sides via aritmetic mean of four neighbour cells
! First the gradients normal to the gradient normal to the sides is copied
! From Green's theorem: du(xi)_dxi = 1/dxi (u_{i+1,j,k} - u_{i,j,k})
DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N-1
  gradU_central_face_xi(1,:,i,j,k)   = gradUxi_central_face(:,i+1,j,k)
END DO; END DO; END DO

DO k=0,PP_NZ; DO j=0,PP_N-1; DO i=0,PP_N
  gradU_central_face_eta(2,:,i,j,k)  = gradUeta_central_face(:,i,j+1,k)
END DO; END DO; END DO

#if (PP_dim==3)
DO k=0,PP_NZ-1; DO j=0,PP_N; DO i=0,PP_N
  gradU_central_face_zeta(3,:,i,j,k) = gradUzeta_central_face(:,i,j,k+1)
END DO; END DO; END DO
#endif

! Now the gradients along the sides are calculated
! From Green's theorem (FD gradient): du(eta)_dxi = 1/deta (u_{i+1/2,j+1/2,k} - u_{i+1/2,j+1/2,k})
! To prevent checkerboarding: u_{i+1/2,j+1/2,k} = 1/4 (u_{i,j,k}+u_{i,j+1,k}+u_{i+1,j,k}+u_{i+1,j+1,k})
DO k=0,PP_NZ; DO j=0,PP_N-1; DO i=0,PP_N
  gradU_central_face_eta(1,:,i,j,k)   = 0.25*(gradUxi_central_face(:,i  ,j  ,k  ) + &
                                              gradUxi_central_face(:,i+1,j  ,k  ) + &
                                              gradUxi_central_face(:,i  ,j+1,k  ) + &
                                              gradUxi_central_face(:,i+1,j+1,k  ))
END DO; END DO; END DO

DO k=0,PP_NZ; DO j=0,PP_N  ; DO i=0,PP_N-1
  gradU_central_face_xi(2,:,i,j,k)   = 0.25*(gradUeta_central_face(:,i  ,j  ,k  ) + &
                                             gradUeta_central_face(:,i+1,j  ,k  ) + &
                                             gradUeta_central_face(:,i  ,j+1,k  ) + &
                                             gradUeta_central_face(:,i+1,j+1,k  ))
END DO; END DO; END DO

#if (PP_dim==3)
DO k=0,PP_NZ-1; DO j=0,PP_N; DO i=0,PP_N
  gradU_central_face_zeta(1,:,i,j,k)  = 0.25*(gradUxi_central_face(:,i  ,j  ,k  ) + &
                                              gradUxi_central_face(:,i+1,j  ,k  ) + &
                                              gradUxi_central_face(:,i  ,j  ,k+1) + &
                                              gradUxi_central_face(:,i+1,j  ,k+1))
END DO; END DO; END DO

DO k=0,PP_NZ-1  ; DO j=0,PP_N; DO i=0,PP_N
  gradU_central_face_zeta(2,:,i,j,k) = 0.25*(gradUeta_central_face(:,i  ,j  ,k  ) + &
                                             gradUeta_central_face(:,i  ,j  ,k+1) + &
                                             gradUeta_central_face(:,i  ,j+1,k  ) + &
                                             gradUeta_central_face(:,i  ,j+1,k+1))
END DO; END DO; END DO

DO k=0,PP_NZ; DO j=0,PP_N  ; DO i=0,PP_N-1
  gradU_central_face_xi(3,:,i,j,k)  = 0.25*(gradUzeta_central_face(:,i  ,j  ,k  ) + &
                                            gradUzeta_central_face(:,i+1,j  ,k  ) + &
                                            gradUzeta_central_face(:,i  ,j  ,k+1) + &
                                            gradUzeta_central_face(:,i+1,j  ,k+1))
END DO; END DO; END DO

DO k=0,PP_NZ  ; DO j=0,PP_N-1; DO i=0,PP_N
  gradU_central_face_eta(3,:,i,j,k) = 0.25*(gradUzeta_central_face(:,i  ,j  ,k  ) + &
                                            gradUzeta_central_face(:,i  ,j  ,k+1) + &
                                            gradUzeta_central_face(:,i  ,j+1,k  ) + &
                                            gradUzeta_central_face(:,i  ,j+1,k+1))
END DO; END DO; END DO
#endif

! 3. Calculate the gradients in x-, y-, z-direction on all FV inner sides with corresponding FV-metrics
DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N-1
#if (PP_dim==3)
  gradUx_xi(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_xi(1,i,j,k,iElem)*gradU_central_face_xi(1,:,i,j,k) &
                           + FV_Metrics_gTilde_sJ_xi(1,i,j,k,iElem)*gradU_central_face_xi(2,:,i,j,k) &
                           + FV_Metrics_hTilde_sJ_xi(1,i,j,k,iElem)*gradU_central_face_xi(3,:,i,j,k)
  gradUy_xi(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_xi(2,i,j,k,iElem)*gradU_central_face_xi(1,:,i,j,k) &
                           + FV_Metrics_gTilde_sJ_xi(2,i,j,k,iElem)*gradU_central_face_xi(2,:,i,j,k) &
                           + FV_Metrics_hTilde_sJ_xi(2,i,j,k,iElem)*gradU_central_face_xi(3,:,i,j,k)
  gradUz_xi(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_xi(3,i,j,k,iElem)*gradU_central_face_xi(1,:,i,j,k) &
                           + FV_Metrics_gTilde_sJ_xi(3,i,j,k,iElem)*gradU_central_face_xi(2,:,i,j,k) &
                           + FV_Metrics_hTilde_sJ_xi(3,i,j,k,iElem)*gradU_central_face_xi(3,:,i,j,k)
#else
  gradUx_xi(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_xi(1,i,j,k,iElem)*gradU_central_face_xi(1,:,i,j,k) &
                           + FV_Metrics_gTilde_sJ_xi(1,i,j,k,iElem)*gradU_central_face_xi(2,:,i,j,k)
  gradUy_xi(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_xi(2,i,j,k,iElem)*gradU_central_face_xi(1,:,i,j,k) &
                           + FV_Metrics_gTilde_sJ_xi(2,i,j,k,iElem)*gradU_central_face_xi(2,:,i,j,k)
#endif
END DO; END DO; END DO! i,j,k=0,PP_N

DO k=0,PP_NZ; DO j=0,PP_N-1; DO i=0,PP_N
#if (PP_dim==3)
  gradUx_eta(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_eta(1,i,j,k,iElem)*gradU_central_face_eta(1,:,i,j,k) &
                            + FV_Metrics_gTilde_sJ_eta(1,i,j,k,iElem)*gradU_central_face_eta(2,:,i,j,k) &
                            + FV_Metrics_hTilde_sJ_eta(1,i,j,k,iElem)*gradU_central_face_eta(3,:,i,j,k)
  gradUy_eta(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_eta(2,i,j,k,iElem)*gradU_central_face_eta(1,:,i,j,k) &
                            + FV_Metrics_gTilde_sJ_eta(2,i,j,k,iElem)*gradU_central_face_eta(2,:,i,j,k) &
                            + FV_Metrics_hTilde_sJ_eta(2,i,j,k,iElem)*gradU_central_face_eta(3,:,i,j,k)
  gradUz_eta(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_eta(3,i,j,k,iElem)*gradU_central_face_eta(1,:,i,j,k) &
                            + FV_Metrics_gTilde_sJ_eta(3,i,j,k,iElem)*gradU_central_face_eta(2,:,i,j,k) &
                            + FV_Metrics_hTilde_sJ_eta(3,i,j,k,iElem)*gradU_central_face_eta(3,:,i,j,k)
#else
  gradUx_eta(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_eta(1,i,j,k,iElem)*gradU_central_face_eta(1,:,i,j,k) &
                            + FV_Metrics_gTilde_sJ_eta(1,i,j,k,iElem)*gradU_central_face_eta(2,:,i,j,k)
  gradUy_eta(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_eta(2,i,j,k,iElem)*gradU_central_face_eta(1,:,i,j,k) &
                            + FV_Metrics_gTilde_sJ_eta(2,i,j,k,iElem)*gradU_central_face_eta(2,:,i,j,k)
#endif
END DO; END DO; END DO! i,j,k=0,PP_N

#if (PP_dim==3)
DO k=0,PP_NZ-1; DO j=0,PP_N; DO i=0,PP_N
  gradUx_zeta(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_zeta(1,i,j,k,iElem)*gradU_central_face_zeta(1,:,i,j,k) &
                             + FV_Metrics_gTilde_sJ_zeta(1,i,j,k,iElem)*gradU_central_face_zeta(2,:,i,j,k) &
                             + FV_Metrics_hTilde_sJ_zeta(1,i,j,k,iElem)*gradU_central_face_zeta(3,:,i,j,k)
  gradUy_zeta(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_zeta(2,i,j,k,iElem)*gradU_central_face_zeta(1,:,i,j,k) &
                             + FV_Metrics_gTilde_sJ_zeta(2,i,j,k,iElem)*gradU_central_face_zeta(2,:,i,j,k) &
                             + FV_Metrics_hTilde_sJ_zeta(2,i,j,k,iElem)*gradU_central_face_zeta(3,:,i,j,k)
  gradUz_zeta(:,i,j,k,iElem) = FV_Metrics_fTilde_sJ_zeta(3,i,j,k,iElem)*gradU_central_face_zeta(1,:,i,j,k) &
                             + FV_Metrics_gTilde_sJ_zeta(3,i,j,k,iElem)*gradU_central_face_zeta(2,:,i,j,k) &
                             + FV_Metrics_hTilde_sJ_zeta(3,i,j,k,iElem)*gradU_central_face_zeta(3,:,i,j,k)
END DO; END DO; END DO! i,j,k=0,PP_N
#endif

END SUBROUTINE FV_CalcGradients_Parabolic

!==================================================================================================================================
!> Calculate gradients used for the viscous fluxes on faces
!==================================================================================================================================
SUBROUTINE FV_PrepareSurfGradient_Parabolic(iElem,gradUxi_tmp,gradUeta_tmp  &
#if PP_dim==3
                                           ,gradUzeta_tmp             &
#endif
                                     )
  ! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars  ,ONLY: FV_Metrics_fTilde_sJ_xi,FV_Metrics_gTilde_sJ_xi
USE MOD_FV_Vars  ,ONLY: FV_Metrics_fTilde_sJ_eta,FV_Metrics_gTilde_sJ_eta
USE MOD_FV_Vars  ,ONLY: gradUx_xi  ,gradUy_xi
USE MOD_FV_Vars  ,ONLY: gradUx_eta ,gradUy_eta
#if (PP_dim==3)
USE MOD_FV_Vars  ,ONLY: FV_Metrics_fTilde_sJ_zeta,FV_Metrics_gTilde_sJ_zeta
USE MOD_FV_Vars  ,ONLY: FV_Metrics_hTilde_sJ_xi
USE MOD_FV_Vars  ,ONLY: FV_Metrics_hTilde_sJ_eta
USE MOD_FV_Vars  ,ONLY: FV_Metrics_hTilde_sJ_zeta
USE MOD_FV_Vars  ,ONLY: gradUx_zeta ,gradUy_zeta ,gradUz_zeta,gradUz_xi,gradUz_eta
USE MOD_FV_Vars   ,ONLY: FV_Metrics_TangVec2_slave,FV_Metrics_TangVec2_master
#endif
USE MOD_FV_Vars   ,ONLY: FV_Metrics_NormVec_slave,FV_Metrics_TangVec1_slave
USE MOD_FV_Vars   ,ONLY: FV_Metrics_NormVec_master,FV_Metrics_TangVec1_master
USE MOD_FV_Vars   ,ONLY: FV_surf_gradU_master,FV_surf_gradU_slave
#if VOLINT_VISC
USE MOD_Mesh_Vars          ,ONLY: ElemToSide,S2V
#endif /* VOLINT_VISC */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: iElem
REAL,INTENT(IN)     :: gradUxi_tmp  (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1)
REAL,INTENT(IN)     :: gradUeta_tmp (PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1)
#if PP_dim==3
REAL,INTENT(IN)     :: gradUzeta_tmp(PP_nVarPrim,0:PP_N,0:PP_NZ,0:PP_N+1)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: locSideID,SideID,p,q,flip,ijk(3)
REAL                :: gradMapPrim_XI(PP_nVarLifting),gradMapPrim_ETA(PP_nVarLifting)
#if PP_dim==3
REAL                :: gradMapPrim_Zeta(PP_nVarLifting)
#endif /* PP_dim==3 */
!==================================================================================================================================
#if PP_dim==3
DO locSideID = 1, 6
#else
DO locSideID = 2, 5
#endif /* PP_dim==3 */
  SideID = ElemToSide(E2S_SIDE_ID,locSideID,iElem)
  flip   = ElemToSide(E2S_FLIP   ,locSideID,iElem)
  DO q=0,PP_NZ; DO p=0,PP_N
    ijk = S2V(:,0,p,q,flip,locSideID)
    ! calculate the local contribution depending on normal and tangential directiont
    SELECT CASE (locSideID)
      CASE(XI_MINUS, XI_PLUS)
        ! Across DG Side
        IF (locSideID .EQ. XI_MINUS) THEN
          gradMapPrim_XI = gradUxi_tmp(PRIM_LIFT,ijk(2),ijk(3),ijk(1))
        ELSE
          gradMapPrim_XI = gradUxi_tmp(PRIM_LIFT,ijk(2),ijk(3),ijk(1)+1)
        END IF
        ! ETA parallel to DG Side: standard arithmetic mean for central gradient
        gradMapPrim_Eta  = 0.5*(gradUeta_tmp(PRIM_LIFT,ijk(1),ijk(3),ijk(2)) &
                               +gradUeta_tmp(PRIM_LIFT,ijk(1),ijk(3),ijk(2)+1))
#if PP_dim==3
        ! ZETA parallel to DG Side: standard arithmetic mean for central gradient
        gradMapPrim_Zeta = 0.5*(gradUzeta_tmp(PRIM_LIFT,ijk(1),ijk(2),ijk(3)) &
                               +gradUzeta_tmp(PRIM_LIFT,ijk(1),ijk(2),ijk(3)+1))
#endif /* PP_dim==3 */
      CASE(ETA_MINUS, ETA_PLUS)
        ! XI parallel to DG Side: standard arithmetic mean for central gradient
        gradMapPrim_XI   = 0.5*(gradUxi_tmp(PRIM_LIFT,ijk(2),ijk(3),ijk(1)) &
                               +gradUxi_tmp(PRIM_LIFT,ijk(2),ijk(3),ijk(1)+1))
        ! Across DG Side
        IF (locSideID .EQ. ETA_MINUS) THEN
          gradMapPrim_Eta = gradUeta_tmp(PRIM_LIFT,ijk(1),ijk(3),ijk(2))
        ELSE
          gradMapPrim_Eta = gradUeta_tmp(PRIM_LIFT,ijk(1),ijk(3),ijk(2)+1)
        END IF
#if PP_dim==3
        ! ZETA parallel to DG Side: standard arithmetic mean for central gradient
        gradMapPrim_Zeta = 0.5*(gradUzeta_tmp(PRIM_LIFT,ijk(1),ijk(2),ijk(3)) &
                               +gradUzeta_tmp(PRIM_LIFT,ijk(1),ijk(2),ijk(3)+1))
#endif /* PP_dim==3 */
#if PP_dim==3
      CASE(ZETA_MINUS, ZETA_PLUS)
        ! XI parallel to DG Side: standard arithmetic mean for central gradient
        gradMapPrim_XI   = 0.5*(gradUxi_tmp(PRIM_LIFT,ijk(2),ijk(3),ijk(1)) &
                               +gradUxi_tmp(PRIM_LIFT,ijk(2),ijk(3),ijk(1)+1))
        ! ETA parallel to DG Side: standard arithmetic mean for central gradient
        gradMapPrim_Eta  = 0.5*(gradUeta_tmp(PRIM_LIFT,ijk(1),ijk(3),ijk(2)) &
                               +gradUeta_tmp(PRIM_LIFT,ijk(1),ijk(3),ijk(2)+1))
        ! Across DG Side
        IF (locSideID .EQ. ZETA_MINUS) THEN
          gradMapPrim_Zeta = gradUeta_tmp(PRIM_LIFT,ijk(1),ijk(2),ijk(3))
        ELSE
          gradMapPrim_Zeta = gradUeta_tmp(PRIM_LIFT,ijk(1),ijk(2),ijk(3)+1)
        END IF
#endif /* PP_dim==3 */
    END SELECT
    ! Rotate from element local xi/eta/zeta sytem to normal system of master side
    IF (flip.EQ.0) THEN
      FV_surf_gradU_master (:,1,p,q,SideID)= FV_Metrics_NormVec_master (1,p,q,SideID)*gradMapPrim_XI (:)&
                                           + FV_Metrics_NormVec_master (2,p,q,SideID)*gradMapPrim_ETA(:)
      FV_surf_gradU_master (:,2,p,q,SideID)= FV_Metrics_TangVec1_master(1,p,q,SideID)*gradMapPrim_XI (:)&
                                           + FV_Metrics_TangVec1_master(2,p,q,SideID)*gradMapPrim_ETA(:)
#if PP_dim==3
      FV_surf_gradU_master (:,1,p,q,SideID)= FV_surf_gradU_master (:,1,p,q,SideID) &
                                           + FV_Metrics_NormVec_master(3,p,q,SideID)*gradMapPrim_ZETA(:)
      FV_surf_gradU_master (:,2,p,q,SideID)= FV_surf_gradU_master (:,2,p,q,SideID) &
                                           + FV_Metrics_TangVec1_master(3,p,q,SideID)*gradMapPrim_ZETA(:)
      FV_surf_gradU_master (:,3,p,q,SideID)= FV_Metrics_TangVec2_master(1,p,q,SideID)*gradMapPrim_XI  (:)&
                                           + FV_Metrics_TangVec2_master(2,p,q,SideID)*gradMapPrim_ETA (:)&
                                           + FV_Metrics_TangVec2_master(3,p,q,SideID)*gradMapPrim_ZETA(:)
#endif /* PP_dim==3 */
    ELSE
      FV_surf_gradU_slave (:,1,p,q,SideID)= FV_Metrics_NormVec_slave (1,p,q,SideID)*gradMapPrim_XI  (:)&
                                          + FV_Metrics_NormVec_slave (2,p,q,SideID)*gradMapPrim_ETA (:)
      FV_surf_gradU_slave (:,2,p,q,SideID)= FV_Metrics_TangVec1_slave(1,p,q,SideID)*gradMapPrim_XI  (:)&
                                          + FV_Metrics_TangVec1_slave(2,p,q,SideID)*gradMapPrim_ETA (:)
#if PP_dim==3
      FV_surf_gradU_slave (:,1,p,q,SideID)= FV_surf_gradU_slave (:,1,p,q,SideID) &
                                          + FV_Metrics_NormVec_slave(3,p,q,SideID)*gradMapPrim_ZETA (:)
      FV_surf_gradU_slave (:,2,p,q,SideID)= FV_surf_gradU_slave (:,2,p,q,SideID) &
                                          + FV_Metrics_TangVec1_slave(3,p,q,SideID)*gradMapPrim_ZETA(:)
      FV_surf_gradU_slave (:,3,p,q,SideID)= FV_Metrics_TangVec2_slave(1,p,q,SideID)*gradMapPrim_XI  (:)&
                                          + FV_Metrics_TangVec2_slave(2,p,q,SideID)*gradMapPrim_ETA (:)&
                                          + FV_Metrics_TangVec2_slave(3,p,q,SideID)*gradMapPrim_ZETA(:)
#endif /* PP_dim==3 */
    END IF ! flip
  END DO; END DO ! q,p = 0, PP_N
END DO ! locSide
END SUBROUTINE FV_PrepareSurfGradient_Parabolic
#endif /* VOLINT_VISC */
#endif /* FV_RECONSTRUCT */


END MODULE MOD_FV_Reconstruction
#endif /* FV_ENABLED */
