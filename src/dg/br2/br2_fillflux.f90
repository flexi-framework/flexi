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
#if PARABOLIC
#include "flexi.h"

!==================================================================================================================================
!> Fills the inner, periodic and bc fluxes for the DG gradients at the interfaces
!==================================================================================================================================
MODULE MOD_Lifting_FillFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_FillFlux
  MODULE PROCEDURE Lifting_FillFlux
END INTERFACE

INTERFACE Lifting_FillFlux_BC
  MODULE PROCEDURE Lifting_FillFlux_BC
END INTERFACE

PUBLIC::Lifting_FillFlux,Lifting_FillFlux_BC
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> \brief Computes the BR2 Surface Fluxes in direction "dir"
!>
!> The routine fills the flux arrays for the sides ranging from firstSideID to lastSideID using the BR2 approximation of surface
!> fluxes, Surfelem contribution is considered as well
!> In case one of the elements contains a FV solution, the FV integral means from the FV element have to be transformed onto the DG
!> nodes set.
!==================================================================================================================================
PPURE SUBROUTINE Lifting_FillFlux(dir,UPrimface_master,UPrimface_slave,Flux,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,   lastInnerSide
USE MOD_Mesh_Vars,       ONLY: firstMPISide_MINE,lastMPISide_MINE
#if FV_ENABLED
USE MOD_FV_Vars         ,ONLY: FV_Elems_Sum,FV_sVdm
USE MOD_ChangeBasisByDim,ONLY: ChangeBasisSurf
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: dir                        !< direction (x,y,z)
LOGICAL,INTENT(IN) :: doMPISides                 !< =.TRUE. only MINE MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(INOUT) :: UPrimface_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< solution on the master sides
REAL,INTENT(INOUT) :: UPrimface_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< solution on the slave sides
REAL,INTENT(OUT)   :: Flux(1:PP_nVarPrim,0:PP_N,0:PP_NZ,nSides)             !< surface flux contribution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            ::SideID,p,q,firstSideID,lastSideID
#if FV_ENABLED
REAL               :: UPrim_glob(1:PP_nVarPrim,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstInnerSide
   lastSideID =  lastInnerSide
END IF

DO SideID = firstSideID,lastSideID
#if FV_ENABLED
  SELECT CASE(FV_Elems_Sum(SideID))
  CASE(0) ! both DG
#endif
    ! BR2 uses strong form, i.e. subtract the solution from the inside
    !BR2: Flux = 1/2(UR+UL)-UL=1/2(UR-UL)
    Flux(:,:,:,SideID)=0.5*(UPrimface_slave(:,:,:,SideID)-UPrimface_master(:,:,:,SideID))
#if FV_ENABLED
  CASE(1) ! master=FV, slave=DG
    CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_sVdm,UPrimface_master(:,:,:,SideID),UPrim_glob)
    Flux(:,:,:,SideID)=0.5*(UPrimface_slave(:,:,:,SideID)-UPrim_glob(:,:,:))
  CASE(2) ! master=DG, slave=FV
    CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_sVdm,UPrimface_slave(:,:,:,SideID),UPrim_glob)
    Flux(:,:,:,SideID)=0.5*(UPrim_glob(:,:,:)-UPrimface_master(:,:,:,SideID))
  CASE(3) ! both FV
    CYCLE
  END SELECT
#endif

  ! multiply by surface area
  DO q=0,PP_NZ
    DO p=0,PP_N
      Flux(:,p,q,SideID)=Flux(:,p,q,SideID)*NormVec(dir,p,q,0,SideID)*SurfElem(p,q,0,SideID)
    END DO
  END DO
END DO ! SideID

END SUBROUTINE Lifting_FillFlux



!==================================================================================================================================
!> \brief Computes the BR2 Surface Fluxes for Boundary Conditions in all three spatial directions.
!> Surfelem contribution is considered as well
!>
!> The routine fills the flux arrays for the BC sides using the BR2 approximation of surface fluxes, Surfelem contribution is
!> considered as well. In case one of the elements contains a FV solution, the FV integral means from the FV element have to be
!> transformed onto the DG nodes set.
!==================================================================================================================================
SUBROUTINE Lifting_FillFlux_BC(t,UPrim_master,FluxX,FluxY,FluxZ)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,TangVec1,TangVec2,Face_xGP,SurfElem,nSides,nBCSides
USE MOD_GetBoundaryFlux, ONLY: Lifting_GetBoundaryFlux
#if FV_ENABLED
USE MOD_FV_Vars         ,ONLY: FV_Elems_master
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: t                                                !< Current time
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< Primitive solution from the inside
REAL,INTENT(OUT)   :: FluxX(       PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< gradient flux in x-dir
REAL,INTENT(OUT)   :: FluxY(       PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< gradient flux in y-dir
REAL,INTENT(OUT)   :: FluxZ(       PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< gradient flux in z-dir
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q
!==================================================================================================================================
! fill flux for boundary sides

! only compute boundary flux once use FluxZ as temp storage
DO SideID=1,nBCSides
  IF (FV_Elems_master(SideID).GT.0) CYCLE
  CALL Lifting_GetBoundaryFlux(SideID,t,UPrim_master(:,:,:,SideID),FluxZ(:,:,:,SideID),&
       NormVec(:,:,:,0,SideID),TangVec1(:,:,:,0,SideID),TangVec2(:,:,:,0,SideID),Face_xGP(:,:,:,0,SideID),SurfElem(:,:,0,SideID))
  DO q=0,PP_NZ; DO p=0,PP_N
    FluxX(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(1,p,q,0,SideID)
    FluxY(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(2,p,q,0,SideID)
#if PP_dim == 3
    FluxZ(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(3,p,q,0,SideID)
#else
    FluxZ(:,p,q,SideID)=0.
#endif
  END DO; END DO
END DO

END SUBROUTINE Lifting_FillFlux_BC

END MODULE MOD_Lifting_FillFlux
#endif /* PARABOLIC */
