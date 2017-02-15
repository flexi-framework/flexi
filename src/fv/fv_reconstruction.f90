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
#if FV_ENABLED
#include "flexi.h"

!==================================================================================================================================
!> Contains the main parts of the second order reconstruction of the FV subcells.
!==================================================================================================================================
MODULE MOD_FV_Reconstruction
! MODULES
! IMPLICIT VARIABLE HANDLING
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

PUBLIC::FV_PrepareSurfGradient
PUBLIC::FV_SurfCalcGradients,FV_SurfCalcGradients_BC,FV_CalcGradients
#endif /* FV_RECONSTRUCT */
!==================================================================================================================================

CONTAINS

#if FV_RECONSTRUCT  


!==================================================================================================================================
!> Fills FV_multi_master and FV_multi_slave arrays. 
!> These arrays contain:
!>  - for DG: solution of first inner layer next to the DG element interface 
!>  - for FV: slopes in normal direction to the DG element interfaces between first and second layer 
!==================================================================================================================================
SUBROUTINE FV_PrepareSurfGradient(UPrim,FV_multi_master,FV_multi_slave,doMPIsides) 
! MODULES                                                                                                                          !
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars ,ONLY: firstMPISide_YOUR,lastMPISide_YOUR,firstInnerSide,lastMPISide_MINE
USE MOD_Mesh_Vars ,ONLY: firstMortarMPISide,lastMortarMPISide,firstBCSide
USE MOD_Mesh_Vars ,ONLY: S2V,S2V2,SideToElem
USE MOD_Mesh_Vars ,ONLY: nElems,nSides
USE MOD_FV_Vars   ,ONLY: FV_Elems,FV_sdx_XI,FV_sdx_ETA,FV_sdx_ZETA
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
! real_in, real_out, real_out, log_in
REAL,INTENT(IN)    :: UPrim(1:PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL,INTENT(OUT)   :: FV_multi_master(1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(OUT)   :: FV_multi_slave (1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
LOGICAL,INTENT(IN) :: doMPIsides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: p,q,locSideID,ElemID,SideID,flip,firstSideID,lastSideID,ijk0(3),ijk(2)
REAL    :: tmp(1:PP_nVarPrim,0:PP_N,0:PP_N)
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
    DO q=0,PP_N; DO p=0,PP_N
      ijk0=S2V(:,0,p,q,flip,locSideID)
      FV_multi_slave(:,p,q,SideID) = UPrim(:,ijk0(1),ijk0(2),ijk0(3),ElemID)
    END DO; END DO
  ELSE ! FV Element
    ! for FV elements calculate slopes between first and second layer 
    !    TODO: USE S2V directly in each of the following 6 cases and not S2V2 below!!!
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,1,p,q,ElemID) - UPrim(:,0,p,q,ElemID)) * FV_sdx_XI(1,p,q,ElemID)
      END DO; END DO
    CASE(ETA_MINUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,p,1,q,ElemID) - UPrim(:,p,0,q,ElemID)) * FV_sdx_ETA(p,1,q,ElemID)
      END DO; END DO
    CASE(ZETA_MINUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,p,q,1,ElemID) - UPrim(:,p,q,0,ElemID)) * FV_sdx_ZETA(p,q,1,ElemID)
      END DO; END DO
    CASE(XI_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,PP_N-1,p,q,ElemID) - UPrim(:,PP_N,p,q,ElemID)) * FV_sdx_XI(PP_N,p,q,ElemID)
      END DO; END DO
    CASE(ETA_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,p,PP_N-1,q,ElemID) - UPrim(:,p,PP_N,q,ElemID)) * FV_sdx_ETA(p,PP_N,q,ElemID)
      END DO; END DO
    CASE(ZETA_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,p,q,PP_N-1,ElemID) - UPrim(:,p,q,PP_N,ElemID)) * FV_sdx_ZETA(p,q,PP_N,ElemID)
      END DO; END DO
    END SELECT

    DO q=0,PP_N; DO p=0,PP_N
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
    DO q=0,PP_N; DO p=0,PP_N
      ijk0=S2V(:,0,p,q,0,locSideID)
      FV_multi_master(:,p,q,SideID) = UPrim(:,ijk0(1),ijk0(2),ijk0(3),ElemID)
    END DO; END DO
  ELSE ! FV Element
    ! for FV elements calculate slopes between first and second layer 
    !    TODO: USE S2V directly in each of the following 6 cases and not S2V2 below!!!
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,1,p,q,ElemID) - UPrim(:,0,p,q,ElemID)) * FV_sdx_XI(1,p,q,ElemID)
      END DO; END DO
    CASE(ETA_MINUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,p,1,q,ElemID) - UPrim(:,p,0,q,ElemID)) * FV_sdx_ETA(p,1,q,ElemID)
      END DO; END DO
    CASE(ZETA_MINUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,p,q,1,ElemID) - UPrim(:,p,q,0,ElemID)) * FV_sdx_ZETA(p,q,1,ElemID)
      END DO; END DO
    CASE(XI_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,PP_N-1,p,q,ElemID) - UPrim(:,PP_N,p,q,ElemID)) * FV_sdx_XI(PP_N,p,q,ElemID)
      END DO; END DO
    CASE(ETA_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,p,PP_N-1,q,ElemID) - UPrim(:,p,PP_N,q,ElemID)) * FV_sdx_ETA(p,PP_N,q,ElemID)
      END DO; END DO
    CASE(ZETA_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        tmp(:,p,q) =  (UPrim(:,p,q,PP_N-1,ElemID) - UPrim(:,p,q,PP_N,ElemID)) * FV_sdx_ZETA(p,q,PP_N,ElemID)
      END DO; END DO
    END SELECT

    DO q=0,PP_N; DO p=0,PP_N
      ijk=S2V2(:,p,q,0,locSideID)
      FV_multi_master(:,p,q,SideID) = tmp(:,ijk(1),ijk(2)) 
    END DO; END DO
  END IF
END DO !SideID
END SUBROUTINE FV_PrepareSurfGradient


!==================================================================================================================================
!> Calculate slopes across DG element interfaces (unlimited). Uses UPrim_master/UPrim_slave as well as FV_multi_master/FV_multi_slave.
!> They are limited in the FV_ProlongToDGFace routine.
!==================================================================================================================================
SUBROUTINE FV_SurfCalcGradients(UPrim_master,UPrim_slave,FV_multi_master,FV_multi_slave,&
        FV_surf_gradU_master,FV_surf_gradU_slave,doMPIsides)
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
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(IN)    :: UPrim_slave (PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(IN)    :: FV_multi_master(1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(IN)    :: FV_multi_slave (1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
LOGICAL,INTENT(IN) :: doMPIsides
REAL,INTENT(OUT)   :: FV_surf_gradU_master(1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(OUT)   :: FV_surf_gradU_slave (1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
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
    DO q=0,PP_N; DO p=0,PP_N
      ! use FV_multi_slave (see FV_PrepareSurfGradient)
      FV_surf_gradU_slave (:,p,q,SideID) = (FV_multi_slave(:,p,q,SideID) - UPrim_master(:,p,q,SideID)) * FV_sdx_Face(p,q,1,SideID)
      FV_surf_gradU_master(:,p,q,SideID) = -FV_surf_gradU_slave(:,p,q,SideID)
    END DO; END DO ! p,q=0,PP_N
  CASE(2) ! master=DG, slave=FV
    DO q=0,PP_N; DO p=0,PP_N
      ! use FV_multi_master (see FV_PrepareSurfGradient)
      FV_surf_gradU_slave (:,p,q,SideID) = (UPrim_slave(:,p,q,SideID) - FV_multi_master(:,p,q,SideID)) * FV_sdx_Face(p,q,2,SideID)
      FV_surf_gradU_master(:,p,q,SideID) = -FV_surf_gradU_slave(:,p,q,SideID)
    END DO; END DO ! p,q=0,PP_N
  CASE(3) ! both FV
    DO q=0,PP_N; DO p=0,PP_N
      FV_surf_gradU_slave (:,p,q,SideID) = (UPrim_slave(:,p,q,SideID) - UPrim_master(:,p,q,SideID)) * FV_sdx_Face(p,q,3,SideID) 
      FV_surf_gradU_master(:,p,q,SideID) = -FV_surf_gradU_slave(:,p,q,SideID) 
    END DO; END DO ! p,q=0,PP_N
  CASE DEFAULT
    CALL Abort(__STAMP__, "FV_Elems is not 0 or 1 somewhere!")
  END SELECT

END DO
END SUBROUTINE FV_SurfCalcGradients


!==================================================================================================================================
!> Get slopes at the boundary conditions.
!==================================================================================================================================
SUBROUTINE FV_SurfCalcGradients_BC(UPrim_master,FV_surf_gradU_master,t)
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
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: UPrim_master(        PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(INOUT) :: FV_surf_gradU_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: SideID
!==================================================================================================================================
DO SideID=firstBCSide,lastBCSide
  IF (FV_Elems_master(SideID).EQ.0) CYCLE
  CALL GetBoundaryFVgradient(SideID,t,FV_surf_gradU_master(:,:,:,SideID),&
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
SUBROUTINE FV_CalcGradients(UPrim,FV_surf_gradU_master,FV_surf_gradU_slave,gradUxi,gradUeta,gradUzeta,&
        gradUxi_central,gradUeta_central,gradUzeta_central)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars        ,ONLY: FV_sdx_XI,FV_sdx_ETA,FV_sdx_ZETA,FV_Elems
USE MOD_FV_Limiter     ,ONLY: FV_Limiter
USE MOD_Mesh_Vars      ,ONLY: nElems,nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)        :: UPrim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems)              !< (IN) volume solution
REAL,INTENT(IN)        :: FV_surf_gradU_master(1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)   !< (IN) slopes over DG element interfaces
REAL,INTENT(IN)        :: FV_surf_gradU_slave (1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)   !< 
REAL,INTENT(OUT)       :: gradUxi  (1:PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL,INTENT(OUT)       :: gradUeta (1:PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL,INTENT(OUT)       :: gradUzeta(1:PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL,INTENT(OUT)       :: gradUxi_central  (1:PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL,INTENT(OUT)       :: gradUeta_central (1:PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL,INTENT(OUT)       :: gradUzeta_central(1:PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N+1) :: gradUxi_tmp
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N+1) :: gradUeta_tmp
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N+1) :: gradUzeta_tmp
INTEGER                                            :: iElem,l,i,j,k,iVar,p,q
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
      gradUxi_tmp  (iVar,:,:,l) = (UPrim(iVar,l,:,:,iElem) - UPrim(iVar,l-1,:,:,iElem)) * FV_sdx_XI  (l,:,:,iElem)
      gradUeta_tmp (iVar,:,:,l) = (UPrim(iVar,:,l,:,iElem) - UPrim(iVar,:,l-1,:,iElem)) * FV_sdx_ETA (:,l,:,iElem)
      gradUzeta_tmp(iVar,:,:,l) = (UPrim(iVar,:,:,l,iElem) - UPrim(iVar,:,:,l-1,iElem)) * FV_sdx_ZETA(:,:,l,iElem)
    END DO
  END DO 
  ! 2. gradients of subcells at DG interface
  ! xi direction
  CALL CopySurfaceToVolume(FV_surf_gradU_master,FV_surf_gradU_slave,gradUxi_tmp,iElem,XI_MINUS,0)
  CALL CopySurfaceToVolume(FV_surf_gradU_master,FV_surf_gradU_slave,gradUxi_tmp,iElem,XI_PLUS ,PP_N+1)
  ! eta direction
  CALL CopySurfaceToVolume(FV_surf_gradU_master,FV_surf_gradU_slave,gradUeta_tmp,iElem,ETA_MINUS,0)
  CALL CopySurfaceToVolume(FV_surf_gradU_master,FV_surf_gradU_slave,gradUeta_tmp,iElem,ETA_PLUS ,PP_N+1)
  ! zeta direction
  CALL CopySurfaceToVolume(FV_surf_gradU_master,FV_surf_gradU_slave,gradUzeta_tmp,iElem,ZETA_MINUS,0)
  CALL CopySurfaceToVolume(FV_surf_gradU_master,FV_surf_gradU_slave,gradUzeta_tmp,iElem,ZETA_PLUS ,PP_N+1)

  ! 3. limit
  DO l=0,PP_N
   DO q=0,PP_N; DO p=0,PP_N
       CALL FV_Limiter(gradUxi_tmp  (:,p,q,l),gradUxi_tmp  (:,p,q,l+1),gradUxi  (:,p,q,l,iElem))
       CALL FV_Limiter(gradUeta_tmp (:,p,q,l),gradUeta_tmp (:,p,q,l+1),gradUeta (:,p,q,l,iElem))
       CALL FV_Limiter(gradUzeta_tmp(:,p,q,l),gradUzeta_tmp(:,p,q,l+1),gradUzeta(:,p,q,l,iElem))
   END DO; END DO ! q, p
  END DO ! l
  ! limit with central limiter for viscous fluxes
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    gradUxi_central  (:,i,j,k,iElem) = 0.5*(gradUxi_tmp  (:,j,k,i)+gradUxi_tmp  (:,j,k,i+1))
    gradUeta_central (:,i,j,k,iElem) = 0.5*(gradUeta_tmp (:,i,k,j)+gradUeta_tmp (:,i,k,j+1))
    gradUzeta_central(:,i,j,k,iElem) = 0.5*(gradUzeta_tmp(:,i,j,k)+gradUzeta_tmp(:,i,j,k+1))
  END DO; END DO; END DO! i,j,k=0,PP_N
END DO
END SUBROUTINE FV_CalcGradients

!==================================================================================================================================
!> Copy surface data at the face specified by dir to the volume. 
!==================================================================================================================================
PURE SUBROUTINE CopySurfaceToVolume(surface_master,surface_slave,volume,iElem,dir,l)
! MODULES
USE MOD_PreProc        ,ONLY: PP_N
USE MOD_Mesh_Vars      ,ONLY: V2S,nSides,ElemToSide
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)         :: surface_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< master surface data 
REAL,INTENT(IN)         :: surface_slave (PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< slave surface data
INTEGER,INTENT(IN)      :: iElem                                             !< element index 
INTEGER,INTENT(IN)      :: dir                                               !< face number of the element (1..6)
INTEGER,INTENT(IN)      :: l                                                 !< volume index where to store the volume data
REAL,INTENT(INOUT)      :: volume(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_N+1)        !< output array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,pq(2),SideID,flip
!================================================================================================================================== 
SideID=ElemToSide(E2S_SIDE_ID,dir,iElem)
flip  =ElemToSide(E2S_FLIP,   dir,iElem)
SELECT CASE(dir)
CASE(XI_MINUS)
  IF (flip.EQ.0) THEN
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,0,i,j,flip,dir)
      volume(:,i,j,l) = surface_master(:,pq(1),pq(2),SideID)
    END DO; END DO
  ELSE
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,0,i,j,flip,dir)
      volume(:,i,j,l) = surface_slave (:,pq(1),pq(2),SideID)
    END DO; END DO
  END IF
CASE(XI_PLUS)
  IF (flip.EQ.0) THEN
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,0,i,j,flip,dir)
      volume(:,i,j,l) = -surface_master(:,pq(1),pq(2),SideID)
    END DO; END DO
  ELSE
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,0,i,j,flip,dir)
      volume(:,i,j,l) = -surface_slave (:,pq(1),pq(2),SideID)
    END DO; END DO
  END IF

CASE(ETA_MINUS)
  IF (flip.EQ.0) THEN
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,i,0,j,flip,dir)
      volume(:,i,j,l) = surface_master(:,pq(1),pq(2),SideID)
    END DO; END DO
  ELSE
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,i,0,j,flip,dir)
      volume(:,i,j,l) = surface_slave (:,pq(1),pq(2),SideID)
    END DO; END DO
  END IF
CASE(ETA_PLUS)
  IF (flip.EQ.0) THEN
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,i,0,j,flip,dir)
      volume(:,i,j,l) = -surface_master(:,pq(1),pq(2),SideID)
    END DO; END DO
  ELSE
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,i,0,j,flip,dir)
      volume(:,i,j,l) = -surface_slave (:,pq(1),pq(2),SideID)
    END DO; END DO
  END IF

CASE(ZETA_MINUS)
  IF (flip.EQ.0) THEN
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,i,j,0,flip,dir)
      volume(:,i,j,l) = surface_master(:,pq(1),pq(2),SideID)
    END DO; END DO
  ELSE
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,i,j,0,flip,dir)
      volume(:,i,j,l) = surface_slave (:,pq(1),pq(2),SideID)
    END DO; END DO
  END IF
CASE(ZETA_PLUS)
  IF (flip.EQ.0) THEN
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,i,j,0,flip,dir)
      volume(:,i,j,l) = -surface_master(:,pq(1),pq(2),SideID)
    END DO; END DO
  ELSE
    DO i=0,PP_N; DO j=0,PP_N
      pq=V2S(1:2,i,j,0,flip,dir)
      volume(:,i,j,l) = -surface_slave (:,pq(1),pq(2),SideID)
    END DO; END DO
  END IF
END SELECT
END SUBROUTINE CopySurfaceToVolume

#endif /* FV_RECONSTRUCT */


END MODULE MOD_FV_Reconstruction
#endif /* FV_ENABLED */
