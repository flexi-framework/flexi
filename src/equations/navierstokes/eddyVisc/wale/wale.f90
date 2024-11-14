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

!===================================================================================================================================
!> Subroutines necessary for calculating the eddy-viscosity in the Wall-Adapting Local Eddy-Viscosity (WALE) model, derived in:
!>   - Nicoud, F., Ducros, F. Subgrid-Scale Stress Modelling Based on the Square of the Velocity Gradient Tensor. Flow, Turbulence
!>     and Combustion 62, 183–200 (1999). https://doi.org/10.1023/A:1009995426001
!==================================================================================================================================
MODULE MOD_WALE
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitWALE
  MODULE PROCEDURE InitWALE
END INTERFACE

INTERFACE WALE
  MODULE PROCEDURE WALE_Point
  MODULE PROCEDURE WALE_Volume
END INTERFACE

INTERFACE FinalizeWALE
  MODULE PROCEDURE FinalizeWALE
END INTERFACE

PUBLIC::InitWALE, WALE_Volume, FinalizeWALE
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get model parameters and initialize WALE model
!===================================================================================================================================
SUBROUTINE InitWALE()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EddyVisc_Vars
USE MOD_ReadInTools        ,ONLY: GETREAL,GETLOGICAL
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone,wGP
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone,nElems,sJ
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: CellVol
!===================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.WALEInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitWALE not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT WALE...'

! Read model coefficient
CS = GETREAL('CS')

! Allocate precomputed (model constant*filter width)**2
ALLOCATE(CSdeltaS2(nElems))

! WALE: (CS*deltaS)**2 * beta * dens
! Calculate the filter width deltaS := (Cell volume)^(1/3) / (PP_N+1)
DO iElem=1,nElems
  CellVol = 0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CellVol = CellVol + wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
  DeltaS(iElem)    = CellVol**(1./3.) / (REAL(PP_N)+1.)
  CsDeltaS2(iElem) = (DeltaS(iElem)*CS)**2.
END DO

WALEInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT WALE DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitWALE


!===================================================================================================================================
!> Compute WALE eddy-visosity
!===================================================================================================================================
PPURE SUBROUTINE WALE_Point(gradUx,gradUy,gradUz,dens,CsDeltaS2,muSGS)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx, gradUy, gradUz   !> Gradients in x,y,z directions
REAL                          ,INTENT(IN)  :: dens       !> pointwise density
REAL                          ,INTENT(IN)  :: CsDeltaS2  !> constant factor (CS*deltaS)**2
REAL                          ,INTENT(OUT) :: muSGS      !> pointwise eddyviscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: g(3,3)      !> Velocity gradient tensor (g_ij)
REAL                                 :: g_sq(3,3)   !> Squared velocity gradient tensor (g^2_ij = g_ik * g_kj)
REAL                                 :: tr_g_sq     !> Volumetric part of squared velocity gradient tensor (g^2_kk/3.)
REAL                                 :: norm_S      !> Norm of strain rate tensor (S_ij*S_ij)
                                                    !>   with: S_ij = 1/2*(g_ij + g_ji)
REAL                                 :: norm_S_d    !> Norm of symmetric, deviatoric part of g_sq (S^d_ij*S^d_ij)
                                                    !>   with: S^d_ij = 1/2*(g^2_ij + g^2_ji) - delta_ij*g^2_kk/3.   (Eq.10)
REAL,PARAMETER                       :: eps = 1.e-9 !> Small number to avoid division by zero
!===================================================================================================================================
! Build velocity gradient tensor (g_ij)
g(:,1)=gradUx(LIFT_VELV)
g(:,2)=gradUy(LIFT_VELV)
g(:,3)=gradUz(LIFT_VELV)

! Squared gradient tensor (g^2_ij = g_ik * g_kj)
g_sq(:,:) = MATMUL(g,g)

! Volumetric contribution of squared gradient tensor (g_kk/3.).
tr_g_sq = (g_sq(1,1) + g_sq(2,2) + g_sq(3,3))/3.

! Compute contractions of type (A_ij*A_ij) of S and S_d, respectively
! (both matrices are symmetric => just take off-diagonal entries from upper triangle twice!)
! 1.) S_ij*S_ij
norm_S   = 0.5*( (g(1,2)+g(2,1))**2 + (g(2,3)+g(3,2))**2 + (g(3,1)+g(1,3))**2 ) & ! Off-diagonal elements
         + g(1,1)**2 + g(2,2)**2 + g(3,3)**2                                      ! Diagonal elements
! 2.) S^d_ij*S^d_ij (Account for trace as in Eq. (10))
norm_S_d = 0.5*( (g_sq(1,2)+g_sq(2,1))**2 + (g_sq(2,3)+g_sq(3,2))**2 + (g_sq(3,1)+g_sq(1,3))**2 ) & ! Off-diagonal
         + (g_sq(1,1)-tr_g_sq)**2 + (g_sq(2,2)-tr_g_sq)**2 + (g_sq(3,3)-tr_g_sq)**2                 ! Diagonal

! WALE model: mu = rho * (DeltaS*C_w)**2 * beta
! Check if gradients are small, since then quotient A/B tends towards 0/0, which is undefined.
! Set mu_sgs=0 manually in this case instead. The limits here are chosen in an ad hoc manner.
IF (norm_S_d .LT. eps) THEN
  muSGS = 0.
ELSE
  muSGS = dens * CsDeltaS2 * norm_S_d**(3./2.) / (norm_S**(5./2.) + norm_S_d**(5./4.)) ! Eq. (13)
END IF
END SUBROUTINE WALE_Point


!===================================================================================================================================
!> Compute WALE Eddy-Visosity for the volume
!===================================================================================================================================
SUBROUTINE WALE_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems
USE MOD_EddyVisc_Vars,     ONLY: CsDeltaS2, muSGS
USE MOD_Lifting_Vars,      ONLY: gradUx, gradUy, gradUz
USE MOD_DG_Vars,           ONLY: U
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
!===================================================================================================================================
DO iElem = 1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CALL WALE_Point(gradUx(   :,i,j,k,iElem), gradUy(:,i,j,k,iElem), gradUz(:,i,j,k,iElem), &
                         U(DENS,i,j,k,iElem), CsDeltaS2(     iElem),  muSGS(1,i,j,k,iElem))
  END DO; END DO; END DO ! i,j,k
END DO
END SUBROUTINE WALE_Volume


!===============================================================================================================================
!> Deallocate arrays and finalize variables used by WALE SGS model
!===============================================================================================================================
SUBROUTINE FinalizeWALE()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===============================================================================================================================
SDEALLOCATE(CsDeltaS2)
WALEInitIsDone = .FALSE.

END SUBROUTINE FinalizeWALE

END MODULE MOD_WALE
