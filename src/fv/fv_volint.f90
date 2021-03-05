!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
!> This module contains only the volume operator of the FV sub-cells method.
!==================================================================================================================================
MODULE MOD_FV_VolInt
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE FV_VolInt
  MODULE PROCEDURE FV_VolInt
END INTERFACE

PUBLIC::FV_VolInt
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Volume operator of the FV sub-cells method.
!> The following steps are performed direction-by-direction (XI/ETA/ZETA) for every inner slice in the respective direction:
!> - reconstruct solution at the sub-cell interfaces
!> - evaluate Riemann solver at the slices
!> - apply fluxes to the left and right sub-cell of the slice
!> are evaluated in the volume integral of the lifting procedure.
!==================================================================================================================================
SUBROUTINE FV_VolInt(UPrim,Ut_FV)
! MODULES
USE MOD_PreProc                         ! all PP_*** variables
USE MOD_Flux         ,ONLY: EvalFlux3D  ! 3D fluxes
USE MOD_FV_Vars
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: gradUx,gradUy,gradUz
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
#endif
USE MOD_Riemann      ,ONLY: Riemann
USE MOD_Mesh_Vars    ,ONLY: nElems
USE MOD_EOS          ,ONLY: PrimToCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)    :: UPrim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  !< solution vector of primitive variables
REAL,INTENT(INOUT) :: Ut_FV(PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  !< time derivative of conservative solution vector
                                                                        !< for FV elements
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                       :: i,j,k,p,q,iElem
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_NZ)    :: UPrim_L,UPrim_R
REAL,DIMENSION(PP_nVar    ,0:PP_N,0:PP_NZ)    :: UCons_L,UCons_R
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ,0:PP_N) :: diffFlux_x,diffFlux_y,diffFlux_z
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ)        :: Fvisc_FV
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: f,g,h      !< viscous volume fluxes at GP
#endif
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ)        :: F_FV
!==================================================================================================================================

! This routine works as follows:
! The tensor product stucture is used to evaluate the fluxes first for all interfaces/slices in xi-direction, then in eta- and
! at last in zeta-direction.
!   0. The viscous fluxes in all sub-cells are calculated.
! For each direction the following steps are performed.
!   1. copy viscous flux
!   2. reconstruct primitive solution at the cell interfaces
!   3. calculate conservative solution at the interfaces (from the reconstructed primitive solution)
!   4. evaluate Riemann solver to get the advective flux
!   5. calculate viscous flux in normal direction at the interface (mean value of the adjacent viscous fluxes)
!   6. add flux to the left and right sub-cell of the interface

DO iElem=1,nElems
  IF (FV_Elems(iElem).EQ.0) CYCLE ! DG Element

#if PARABOLIC
  ! 0. Eval viscous flux for all sub-cells
  CALL EvalDiffFlux3D(UPrim(:,:,:,:,iElem), gradUx(:,:,:,:,iElem), gradUy(:,:,:,:,iElem), gradUz(:,:,:,:,iElem),f,g,h,iElem)
#endif

  ! === Xi-Direction ============
#if PARABOLIC
  ! 1. copy viscous fluxes to temporary array (for performance)
  DO i=0,PP_N
    diffFlux_x(:,:,:,i) = f(:,i,:,:)
    diffFlux_y(:,:,:,i) = g(:,i,:,:)
    diffFlux_z(:,:,:,i) = h(:,i,:,:)
  END DO ! i=0,PP_N
#endif
  ! This is the only nullification we have to do, the rest will be overwritten accordingly. Nullify here to catch only the FV
  ! elements. Might be smarter to do this in a single, long operation?
  Ut_FV(:,0,:,:,iElem) = 0.
  ! iterate over all inner slices in xi-direction
  DO i=1,PP_N
    DO q=0,PP_NZ; DO p=0,PP_N
      ! 2. reconstruct solution at left and right side of the interface/slice
#if FV_RECONSTRUCT
      UPrim_L(:,p,q) = UPrim(:,i-1,p,q,iElem) + gradUxi(:,p,q,i-1,iElem) * FV_dx_XI_R(p,q,i-1,iElem)
      UPrim_R(:,p,q) = UPrim(:,i  ,p,q,iElem) - gradUxi(:,p,q,i  ,iElem) * FV_dx_XI_L(p,q,i  ,iElem)
#else
      UPrim_L(:,p,q) = UPrim(:,i-1,p,q,iElem)
      UPrim_R(:,p,q) = UPrim(:,i  ,p,q,iElem)
#endif
    END DO; END DO ! p,q=0,PP_N
    ! 3. convert primitve solution to conservative
    CALL PrimToCons(PP_N,UPrim_L, UCons_L)
    CALL PrimToCons(PP_N,UPrim_R, UCons_R)

    ! 4. calculate advective part of the flux
    CALL Riemann(PP_N,F_FV,UCons_L,UCons_R,UPrim_L,UPrim_R,          &
        FV_NormVecXi (:,:,:,i,iElem), FV_TangVec1Xi(:,:,:,i,iElem), FV_TangVec2Xi(:,:,:,i,iElem),.FALSE.)

#if PARABOLIC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_NZ; DO p=0,PP_N
      Fvisc_FV(:,p,q)=0.5*(FV_NormVecXi(1,p,q,i,iElem)*(diffFlux_x(:,p,q,i-1)+diffFlux_x(:,p,q,i)) &
                          +FV_NormVecXi(2,p,q,i,iElem)*(diffFlux_y(:,p,q,i-1)+diffFlux_y(:,p,q,i)) &
                          +FV_NormVecXi(3,p,q,i,iElem)*(diffFlux_z(:,p,q,i-1)+diffFlux_z(:,p,q,i)))
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*PARABOLIC*/

    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO k=0,PP_NZ; DO j=0,PP_N
      Ut_FV(:,i-1,j,k,iElem) = Ut_FV(:,i-1,j,k,iElem) + F_FV(:,j,k) * FV_SurfElemXi_sw(j,k,i,iElem)
      ! During our first sweep, the DOF here has never been touched and can thus be overwritten
      Ut_FV(:,i  ,j,k,iElem) =                     -1.* F_FV(:,j,k) * FV_SurfElemXi_sw(j,k,i,iElem)
    END DO; END DO
  END DO ! i

  ! === Eta-Direction ===========
#if PARABOLIC
  ! 1. copy fluxes to temporary array (for performance)
  DO j=0,PP_N
    diffFlux_x(:,:,:,j) = f(:,:,j,:)
    diffFlux_y(:,:,:,j) = g(:,:,j,:)
    diffFlux_z(:,:,:,j) = h(:,:,j,:)
  END DO ! j=0,PP_N
#endif

  ! iterate over all inner slices in eta-direction
  DO j=1,PP_N
    DO q=0,PP_NZ; DO p=0,PP_N
      ! 2. reconstruct solution at left and right side of the interface/slice
#if FV_RECONSTRUCT
      UPrim_L(:,p,q) = UPrim(:,p,j-1,q,iElem) + gradUeta(:,p,q,j-1,iElem) * FV_dx_ETA_R(p,q,j-1,iElem)
      UPrim_R(:,p,q) = UPrim(:,p,j  ,q,iElem) - gradUeta(:,p,q,j  ,iElem) * FV_dx_ETA_L(p,q,j  ,iElem)
#else
      UPrim_L(:,p,q) = UPrim(:,p,j-1,q,iElem)
      UPrim_R(:,p,q) = UPrim(:,p,j  ,q,iElem)
#endif
    END DO; END DO ! p,q=0,PP_N
    ! 3. convert primitve solution to conservative
    CALL PrimToCons(PP_N,UPrim_L, UCons_L)
    CALL PrimToCons(PP_N,UPrim_R, UCons_R)

    ! 4. calculate advective part of the flux
    CALL Riemann(PP_N,F_FV,UCons_L,UCons_R,UPrim_L,UPrim_R,          &
        FV_NormVecEta (:,:,:,j,iElem), FV_TangVec1Eta(:,:,:,j,iElem), FV_TangVec2Eta(:,:,:,j,iElem),.FALSE.)
#if PARABOLIC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_NZ; DO p=0,PP_N
      Fvisc_FV(:,p,q)=0.5*(FV_NormVecEta(1,p,q,j,iElem)*(diffFlux_x(:,p,q,j-1)+diffFlux_x(:,p,q,j)) &
                          +FV_NormVecEta(2,p,q,j,iElem)*(diffFlux_y(:,p,q,j-1)+diffFlux_y(:,p,q,j)) &
                          +FV_NormVecEta(3,p,q,j,iElem)*(diffFlux_z(:,p,q,j-1)+diffFlux_z(:,p,q,j)))
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*PARABOLIC*/

    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO k=0,PP_NZ; DO i=0,PP_N
      Ut_FV(:,i,j-1,k,iElem) = Ut_FV(:,i,j-1,k,iElem) + F_FV(:,i,k) * FV_SurfElemEta_sw(i,k,j,iElem)
      Ut_FV(:,i,j  ,k,iElem) = Ut_FV(:,i,j  ,k,iElem) - F_FV(:,i,k) * FV_SurfElemEta_sw(i,k,j,iElem)
    END DO; END DO
  END DO ! j

#if PP_dim == 3
  ! === Zeta-Direction ============
  ! 1. no copy of viscous fluxes to diffFlux_x/y/z required, since f,g,h already have the correct memory layout

  ! iterate over all inner slices in zeta-direction
  DO k=1,PP_N
    DO q=0,PP_N; DO p=0,PP_N
      ! 2. reconstruct solution at left and right side of the interface/slice
#if FV_RECONSTRUCT
      UPrim_L(:,p,q) = UPrim(:,p,q,k-1,iElem) + gradUzeta(:,p,q,k-1,iElem) * FV_dx_ZETA_R(p,q,k-1,iElem)
      UPrim_R(:,p,q) = UPrim(:,p,q,k  ,iElem) - gradUzeta(:,p,q,k  ,iElem) * FV_dx_ZETA_L(p,q,k  ,iElem)
#else
      UPrim_L(:,p,q) = UPrim(:,p,q,k-1,iElem)
      UPrim_R(:,p,q) = UPrim(:,p,q,k  ,iElem)
#endif
    END DO; END DO ! p,q=0,PP_N
    ! 3. convert primitve solution to conservative
    CALL PrimToCons(PP_N,UPrim_L, UCons_L)
    CALL PrimToCons(PP_N,UPrim_R, UCons_R)

    ! 4. calculate advective part of the flux
    CALL Riemann(PP_N,F_FV,UCons_L,UCons_R,UPrim_L,UPrim_R,          &
        FV_NormVecZeta (:,:,:,k,iElem), FV_TangVec1Zeta(:,:,:,k,iElem), FV_TangVec2Zeta(:,:,:,k,iElem),.FALSE.)
#if PARABOLIC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_N; DO p=0,PP_N
      Fvisc_FV(:,p,q)=0.5*(FV_NormVecZeta(1,p,q,k,iElem)*(f(:,p,q,k-1)+f(:,p,q,k)) &
                          +FV_NormVecZeta(2,p,q,k,iElem)*(g(:,p,q,k-1)+g(:,p,q,k)) &
                          +FV_NormVecZeta(3,p,q,k,iElem)*(h(:,p,q,k-1)+h(:,p,q,k)))
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*PARABOLIC*/

    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO j=0,PP_N; DO i=0,PP_N
      Ut_FV(:,i,j,k-1,iElem) = Ut_FV(:,i,j,k-1,iElem) + F_FV(:,i,j) * FV_SurfElemZeta_sw(i,j,k,iElem)
      Ut_FV(:,i,j,k  ,iElem) = Ut_FV(:,i,j,k  ,iElem) - F_FV(:,i,j) * FV_SurfElemZeta_sw(i,j,k,iElem)
    END DO; END DO
  END DO ! k
#endif /* PP_dim == 3 */

END DO ! iElem
END SUBROUTINE FV_VolInt

END MODULE MOD_FV_VolInt
#endif /* FV_ENABLED */
