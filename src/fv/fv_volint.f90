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
SUBROUTINE FV_VolInt(UPrim,Ut)
! MODULES
USE MOD_PreProc                         ! all PP_*** variables
USE MOD_Flux         ,ONLY: EvalFlux3D  ! 3D fluxes
USE MOD_FV_Vars
#if VOLINT_VISC
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
REAL,INTENT(INOUT) :: Ut(   PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  !< time derivative of conservative solution vector
                                                                         !< for FV elements
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                         :: i,j,k,p,q,iElem
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_NZ)      :: UPrim_L,UPrim_R
REAL,DIMENSION(PP_nVar    ,0:PP_N,0:PP_NZ)      :: UCons_L,UCons_R
#if VOLINT_VISC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ,0:PP_N-1) :: diffFlux_x,diffFlux_y
REAL,DIMENSION(PP_nVar,0:PP_N-1,0:PP_N,0:PP_NZ) :: f_xi,g_xi,h_xi
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N-1,0:PP_NZ) :: f_eta,g_eta,h_eta
REAL                                            :: UPrim_xi  (PP_nVarPrim,0:PP_N-1,0:PP_N  ,0:PP_NZ  )
REAL                                            :: UPrim_eta (PP_nVarPrim,0:PP_N  ,0:PP_N-1,0:PP_NZ  )
#if PP_dim == 3
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ,0:PP_N-1) :: diffFlux_z
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ-1) :: f_zeta,g_zeta,h_zeta
REAL                                            :: UPrim_zeta(PP_nVarPrim,0:PP_N  ,0:PP_N  ,0:PP_NZ-1)
#endif /*PP_dim == 3*/
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ)          :: Fvisc_FV
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)   :: f,g,h      !< viscous volume fluxes at GP
#endif
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ)        :: F_FV
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: Ut_FV
#if FV_ENABLED==3
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: Flux_xi
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: Flux_eta
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: Ut_FV_xi
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: Ut_FV_eta
#if PP_dim == 3
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: Flux_zeta
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: Ut_FV_zeta
#endif /*PP_dim == 3*/
#endif /*FV_ENABLED==3*/
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
#if FV_ENABLED == 1
  IF (FV_Elems(iElem) .EQ. 0) CYCLE ! DG Element
#elif FV_ENABLED == 2
  IF (FV_alpha(iElem) .LT. EPSILON(0.)) CYCLE
#endif

#if VOLINT_VISC
  ! 0. Eval viscous flux for all sub-cells
  ! 0.1 Calculate mean states on the FV element faces
  UPrim_xi   = 0.5*(UPrim(:,0:PP_N-1,:,:,iElem)+UPrim(:,1:PP_N,:,:,iElem))
  UPrim_eta  = 0.5*(UPrim(:,:,0:PP_N-1,:,iElem)+UPrim(:,:,1:PP_N,:,iElem))
#if (PP_dim==3)
  UPrim_zeta = 0.5*(UPrim(:,:,:,0:PP_NZ-1,iElem)+UPrim(:,:,:,1:PP_NZ,iElem))
#endif
  ! 0.2 Calculate the viscous fluxes on all inner faces, here three times: xi-faces, eta-faces, zeta-faces
  CALL EvalDiffFlux3D(UPrim_xi  ,gradUx_xi  (:,:,:,:,iElem),gradUy_xi  (:,:,:,:,iElem),gradUz_xi  (:,:,:,:,iElem),f_xi  ,g_xi  ,h_xi  ,iElem,PP_N-1,PP_N  ,PP_NZ  )
  CALL EvalDiffFlux3D(UPrim_eta ,gradUx_eta (:,:,:,:,iElem),gradUy_eta (:,:,:,:,iElem),gradUz_eta (:,:,:,:,iElem),f_eta ,g_eta ,h_eta ,iElem,PP_N  ,PP_N-1,PP_NZ  )
#if PP_dim==3
  CALL EvalDiffFlux3D(UPrim_zeta,gradUx_zeta(:,:,:,:,iElem),gradUy_zeta(:,:,:,:,iElem),gradUz_zeta(:,:,:,:,iElem),f_zeta,g_zeta,h_zeta,iElem,PP_N  ,PP_N  ,PP_NZ-1)
#endif
#endif /*VOLINT_VISC*/

  ! === Xi-Direction ============
#if VOLINT_VISC
  ! 1. copy viscous fluxes to temporary array (for performance)
  DO i=0,PP_N-1
    diffFlux_x(:,:,:,i) = f_xi(:,i,:,:)
    diffFlux_y(:,:,:,i) = g_xi(:,i,:,:)
#if PP_dim == 3
    diffFlux_z(:,:,:,i) = h_xi(:,i,:,:)
#endif
  END DO ! i=0,PP_N
#endif
  ! This is the only nullification we have to do, the rest will be overwritten accordingly. Nullify here to catch only the FV
  ! elements. Might be smarter to do this in a single, long operation?
  Ut_FV(:,0,:,:) = 0.
  ! For FV subblend only
#if FV_ENABLED == 3
  Ut_FV_xi  (:,:,:,:) = 0.
  Ut_FV_eta (:,:,:,:) = 0.
#if PP_dim == 3
  Ut_FV_zeta(:,:,:,:) = 0.
#endif /* PP_dim == 3 */
#endif /*FV_ENABLED*/
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

#if VOLINT_VISC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_NZ; DO p=0,PP_N
      Fvisc_FV(:,p,q)=FV_NormVecXi(1,p,q,i,iElem)*diffFlux_x(:,p,q,i-1) &
#if PP_dim == 3
                     +FV_NormVecXi(3,p,q,i,iElem)*diffFlux_z(:,p,q,i-1) &
#endif
                     +FV_NormVecXi(2,p,q,i,iElem)*diffFlux_y(:,p,q,i-1)
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*VOLINT_VISC*/

#if FV_ENABLED <= 2
    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO k=0,PP_NZ; DO j=0,PP_N
      Ut_FV(:,i-1,j,k) = Ut_FV(:,i-1,j,k) + F_FV(:,j,k) * FV_SurfElemXi_sw(j,k,i,iElem) * FV_w_inv(i-1)
      ! During our first sweep, the DOF here has never been touched and can thus be overwritten
      Ut_FV(:,i  ,j,k) =               -1.* F_FV(:,j,k) * FV_SurfElemXi_sw(j,k,i,iElem) * FV_w_inv(i)
    END DO; END DO
#elif FV_ENABLED == 3
    DO k=0,PP_NZ; DO j=0,PP_N
      Ut_FV_xi  (:,i-1,j,k) = Ut_FV_xi  (  :,i-1,j,k)  + F_FV(:,j,k) * FV_SurfElemXi_sw(j,k,i,iElem) * FV_w_inv(i-1)
      Ut_FV_xi  (:,i  ,j,k) =                          - F_FV(:,j,k) * FV_SurfElemXi_sw(j,k,i,iElem) * FV_w_inv(i)
    END DO; END DO
#endif /*FV_ENABLED*/
END DO ! i
#if FV_ENABLED == 3
  DO i=1,PP_N-1
    Ut_FV_xi  (:,i,:,:) = Ut_FV_xi  (:,i,:,:) + Ut_FV_xi  (:,i-1,:,:)
  END DO ! i
#endif /*FV_ENABLED*/

  ! === Eta-Direction ===========
#if VOLINT_VISC
  ! 1. copy fluxes to temporary array (for performance)
  DO j=0,PP_N-1
    diffFlux_x(:,:,:,j) = f_eta(:,:,j,:)
    diffFlux_y(:,:,:,j) = g_eta(:,:,j,:)
#if PP_dim == 3
    diffFlux_z(:,:,:,j) = h_eta(:,:,j,:)
#endif
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
#if VOLINT_VISC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_NZ; DO p=0,PP_N
      Fvisc_FV(:,p,q)=FV_NormVecEta(1,p,q,j,iElem)*diffFlux_x(:,p,q,j-1) &
#if PP_dim == 3
                     +FV_NormVecEta(3,p,q,j,iElem)*diffFlux_z(:,p,q,j-1) &
#endif
                     +FV_NormVecEta(2,p,q,j,iElem)*diffFlux_y(:,p,q,j-1)
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*VOLINT_VISC*/

#if FV_ENABLED <= 2
    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO k=0,PP_NZ; DO i=0,PP_N
      Ut_FV(:,i,j-1,k) = Ut_FV(:,i,j-1,k) + F_FV(:,i,k) * FV_SurfElemEta_sw(i,k,j,iElem) * FV_w_inv(j-1)
      Ut_FV(:,i,j  ,k) = Ut_FV(:,i,j  ,k) - F_FV(:,i,k) * FV_SurfElemEta_sw(i,k,j,iElem) * FV_w_inv(j)
    END DO; END DO
#elif FV_ENABLED == 3
    DO k=0,PP_NZ; DO i=0,PP_N
      Ut_FV_eta (:,i,j-1,k) =  Ut_FV_eta (:,i,j-1,k)   + F_FV(:,i,k) * FV_SurfElemEta_sw(i,k,j,iElem) * FV_w_inv(j-1)
      Ut_FV_eta (:,i,j  ,k) =                          - F_FV(:,i,k) * FV_SurfElemEta_sw(i,k,j,iElem) * FV_w_inv(j)
    END DO; END DO
#endif /*FV_ENABLED*/
END DO ! j
#if FV_ENABLED == 3
  DO j=1,PP_N-1
    Ut_FV_eta (:,:,j,:) = Ut_FV_eta (:,:,j,:) + Ut_FV_eta (:,:,j-1,:)
  END DO ! j
#endif /*FV_ENABLED*/

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
#if VOLINT_VISC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_N; DO p=0,PP_N
      Fvisc_FV(:,p,q)=FV_NormVecZeta(1,p,q,k,iElem)*f_zeta(:,p,q,k-1) &
                     +FV_NormVecZeta(2,p,q,k,iElem)*g_zeta(:,p,q,k-1) &
                     +FV_NormVecZeta(3,p,q,k,iElem)*h_zeta(:,p,q,k-1)
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*VOLINT_VISC*/

#if FV_ENABLED <= 2
    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO j=0,PP_N; DO i=0,PP_N
      Ut_FV(:,i,j,k-1) = Ut_FV(:,i,j,k-1) + F_FV(:,i,j) * FV_SurfElemZeta_sw(i,j,k,iElem) * FV_w_inv(k-1)
      Ut_FV(:,i,j,k  ) = Ut_FV(:,i,j,k  ) - F_FV(:,i,j) * FV_SurfElemZeta_sw(i,j,k,iElem) * FV_w_inv(k)
    END DO; END DO
#elif FV_ENABLED == 3
    DO j=0,PP_N; DO i=0,PP_N
      Ut_FV_zeta(:,i,j,k-1) = Ut_FV_zeta(:,i,j,k-1)    + F_FV(:,i,j) * FV_SurfElemZeta_sw(i,j,k,iElem) * FV_w_inv(k-1)
      Ut_FV_zeta(:,i,j,k  ) =                          - F_FV(:,i,j) * FV_SurfElemZeta_sw(i,j,k,iElem) * FV_w_inv(k)
    END DO; END DO
#endif /*FV_ENABLED*/
  END DO ! k
#if FV_ENABLED == 3
  DO k=1,PP_N-1
    Ut_FV_zeta(:,:,:,k) = Ut_FV_zeta(:,:,:,k) + Ut_FV_zeta(:,:,:,k-1)
  END DO ! k
#endif /*FV_ENABLED*/
#endif /* PP_dim == 3 */

#if FV_ENABLED == 3
  ! 2. Blend the DG and the FV volume integral                     (Rueda-Ramírez, 2022; formula: 18)
  ! 3. Calculate the volume integral from the blended inner fluxes (Rueda-Ramírez, 2022; formula: 17)

  ! u_t = - (f^hat_{i,i+1}-f^hat_{i-1,i})^DG - (f^hat_{i,i+1}-f^hat_{i-1,i})^FV
  ! f^hat_{a,b} = (1-\alpha_{a,b}) f^hat_{a,b}^DG + \alpha_{a,b} f^hat_{a,b}^FV
  ! f^hat_{i,i+1}^DG = f^hat_{i-1,i}^DG + f^hat_{i,i+1}^DG
  DO k=0,PP_NZ; DO j=0,PP_N
    ! Blend in xi   direction
    Flux_xi  (:,0,j,k)       = Ut_xi    (:,0,j,k,iElem)
    Ut_xi    (:,0,j,k,iElem) = Flux_xi  (:,0,j,k      ) * (1.-FV_alpha(FV_int(1),0,j,k,iElem))  + Ut_FV_xi  (:,0,j,k) * FV_alpha(FV_int(1),0,j,k,iElem)
    Ut       (:,0,j,k,iElem) =                            Ut_xi  (             :,0,j,k,iElem)
    DO i=1,PP_N-1
      Flux_xi  (:,i,j,k)       = Flux_xi  (:,i-1,j,k) + Ut_xi  (     :,i,j  ,k,iElem)
      Ut_xi    (:,i,j,k,iElem) = Flux_xi  (:,i  ,j,k) * (1.-FV_alpha(FV_int(1),i,j  ,k,iElem)) + Ut_FV_xi  (:,i,j,k) * FV_alpha(FV_int(1),i,j,k,iElem)
      Ut       (:,i,j,k,iElem) =                        (-Ut_xi  (           :,i-1,j,k,iElem)  + Ut_xi     (:,i,j,k,iElem))
    END DO ! i
    ! Calculate the volume integral from the blended inner fluxes (Rueda-Ramírez, 2022; formula: 17)
    Flux_xi  (:,PP_N,j,k)       = Ut_xi  (  :,PP_N,j,k,iElem)
    Ut_xi    (:,PP_N,j,k,iElem) = Flux_xi  (:,PP_N  ,j,k) * (1.-FV_alpha(FV_int(1),PP_N,j,k,iElem))  + Ut_FV_xi  (:,PP_N,j,k) * FV_alpha(FV_int(1),PP_N,j,k,iElem)
    Ut       (:,PP_N,j,k,iElem) =                               Ut_xi  (         :,PP_N,j,k,iElem)
  END DO; END DO ! j,k

  DO k=0,PP_NZ; DO i=0,PP_N
    ! Blend in eta  direction
    Flux_eta (:,i,0,k)       = Ut_eta   (:,i,0,k,iElem)
    Ut_eta   (:,i,0,k,iElem) = Flux_eta (:,i,0,k      ) * (1.-FV_alpha(FV_int(1),i,0,k,iElem))  + Ut_FV_eta (:,i,0,k) * FV_alpha(FV_int(1),i,0,k,iElem)
    Ut       (:,i,0,k,iElem) = Ut       (:,i,0,k,iElem) + Ut_eta (             :,i,0,k,iElem)
    DO j=1,PP_N-1
      Flux_eta (:,i,j,k)       = Flux_eta (:,i,j-1,k)     + Ut_eta (     :,i,j  ,k,iElem)
      Ut_eta   (:,i,j,k,iElem) = Flux_eta (:,i,j  ,k)     * (1.-FV_alpha(FV_int(1),i,j  ,k,iElem)) + Ut_FV_eta (:,i,j,k) * FV_alpha(FV_int(1),i,j,k,iElem)
      Ut       (:,i,j,k,iElem) = Ut       (:,i,j,k,iElem) + (-Ut_eta (           :,i,j-1,k,iElem)  + Ut_eta    (:,i,j,k,iElem))
    END DO ! i
    ! Calculate the volume integral from the blended inner fluxes (Rueda-Ramírez, 2022; formula: 17)
    Flux_eta (:,i,PP_N,k)       = Ut_eta (  :,i,PP_N,k,iElem)
    Ut_eta   (:,i,PP_N,k,iElem) = Flux_eta (:,i,PP_N  ,k) * (1.-FV_alpha(FV_int(1),i,PP_N,k,iElem)) + Ut_FV_eta (:,i,PP_N,k) * FV_alpha(FV_int(1),i,PP_N,k,iElem)
    Ut       (:,i,PP_N,k,iElem) = Ut(       :,i,PP_N,k,iElem) + Ut_eta (         :,i,PP_N,k,iElem)
  END DO; END DO ! j,k
#if PP_dim == 3
  DO j=0,PP_N; DO i=0,PP_N
    ! Blend in zeta direction
    Flux_zeta(:,i,j,0)       = Ut_zeta  (:,i,j,0,iElem)
    Ut_zeta  (:,i,j,0,iElem) = Flux_zeta(:,i,j,0      ) * (1.-FV_alpha(FV_int(1),i,j,0,iElem)) + Ut_FV_zeta(:,i,j,0) * FV_alpha(FV_int(1),i,j,0,iElem)
    Ut       (:,i,j,0,iElem) = Ut       (:,i,j,0,iElem) + Ut_zeta(             :,i,j,0,iElem)
    DO k=1,PP_NZ-1
      Flux_zeta(:,i,j,k)       = Flux_zeta(:,i,j,k-1)     + Ut_zeta(     :,i,j,k  ,iElem)
      Ut_zeta  (:,i,j,k,iElem) = Flux_zeta(:,i,j,k  )     * (1.-FV_alpha(FV_int(1),i,j,k  ,iElem)) + Ut_FV_zeta(:,i,j,k) * FV_alpha(FV_int(1),i,j,k,iElem)
      Ut       (:,i,j,k,iElem) = Ut       (:,i,j,k,iElem) + (-Ut_zeta(           :,i,j,k-1,iElem)  + Ut_zeta(   :,i,j,k,iElem))
    END DO ! i
    ! Calculate the volume integral from the blended inner fluxes (Rueda-Ramírez, 2022; formula: 17)
    Flux_zeta(:,i,j,PP_N)       = Ut_zeta(  :,i,j,PP_N,iElem)
    Ut_zeta  (:,i,j,PP_N,iElem) = Flux_zeta(:,i,j,PP_N  ) * (1.-FV_alpha(FV_int(1),i,j,PP_N,iElem)) + Ut_FV_zeta(:,i,j,PP_N) * FV_alpha(FV_int(1),i,j,PP_N,iElem)
    Ut       (:,i,j,PP_N,iElem) = Ut(       :,i,j,PP_N,iElem) + Ut_zeta(         :,i,j,PP_N,iElem)
  END DO; END DO ! j,k
#endif /* PP_dim == 3 */
#elif FV_ENABLED == 2
  ! Blend the solutions together
  Ut(:,:,:,:,iElem) = (1 - FV_alpha(iElem)) * Ut(:,:,:,:,iElem) + FV_alpha(iElem)*Ut_FV
#else
  Ut(:,:,:,:,iElem) = Ut_FV
#endif /*FV_BLENDING*/
END DO ! iElem
END SUBROUTINE FV_VolInt

END MODULE MOD_FV_VolInt
#endif /* FV_ENABLED */
