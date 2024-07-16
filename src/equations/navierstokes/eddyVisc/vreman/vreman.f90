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
!> Subroutines necessary for calculating Vreman's Eddy-Viscosity, originally derived in
!>   - Vreman, A. W. "An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory and applications."
!>     Physics of fluids 16.10 (2004): 3670-3681.
!===================================================================================================================================
MODULE MOD_Vreman
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitVreman
   MODULE PROCEDURE InitVreman
END INTERFACE

INTERFACE Vreman
   MODULE PROCEDURE Vreman_Point
   MODULE PROCEDURE Vreman_Volume
END INTERFACE

INTERFACE FinalizeVreman
   MODULE PROCEDURE FinalizeVreman
END INTERFACE

PUBLIC::InitVreman, Vreman_Volume, FinalizeVreman
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get model parameters and initialize Vreman's model
!===================================================================================================================================
SUBROUTINE InitVreman()
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
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.VremanInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitVreman not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT VREMAN...'

! Read model coefficient
! Vreman model, paper CS=Smagorinsky constant: 0.18
! Vreman model, paper CS=Smagorinsky constant for FLEXI: 0.11
CS        = GETREAL('CS')

! Allocate precomputed (model constant*filter width)**2
ALLOCATE(CSdeltaS2(nElems))

! Vreman: (CS*deltaS)**2 * SQRT(B/A) * dens
! Precompute first term and store in damp
! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )
DO iElem=1,nElems
  CellVol = 0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CellVol = CellVol + wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
  DeltaS(iElem)    = CellVol**(1./3.)  / (REAL(PP_N)+1.)
  CSdeltaS2(iElem) = 2.5*(CS * DeltaS(iElem))**2
END DO

VremanInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT VREMAN DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitVreman


!===================================================================================================================================
!> Compute Vreman Eddy-Visosity
!===================================================================================================================================
PPURE SUBROUTINE Vreman_Point(gradUx,gradUy,gradUz,dens,CSdeltaS2,muSGS)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx, gradUy, gradUz   !> Gradients in x,y,z directions
REAL                          ,INTENT(IN)  :: dens       !> pointwise density
REAL                          ,INTENT(IN)  :: CSdeltaS2  !> filter width
REAL                          ,INTENT(OUT) :: muSGS      !> pointwise eddyviscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: i,j,k
REAL         :: alpha(3,3),beta(3,3)
REAL         :: A,B
!===================================================================================================================================
alpha(1,:)=gradUx(LIFT_VELV)
alpha(2,:)=gradUy(LIFT_VELV)
alpha(3,:)=gradUz(LIFT_VELV)
beta=0.
DO j=1,3; DO i=1,3; DO k=1,3
  beta(i,j)=beta(i,j)+alpha(k,i)*alpha(k,j)
END DO; END DO; END DO! i,j,k=1,3
B=beta(1,1)*beta(2,2)-beta(1,2)**2+beta(1,1)*beta(3,3)-beta(1,3)**2+beta(2,2)*beta(3,3)-beta(2,3)**2

A=0
DO j=1,3; DO i=1,3
  A=A+alpha(i,j)*alpha(i,j)
END DO; END DO! i,j=1,3

! Vreman: 2.5*(CS * deltaS)**2 * SQRT(B/A) * dens
! Check if gradients are small, since then quotient A/B tends towards 0/0, which is undefined.
! Set mu_sgs=0 manually in this case instead. The limits here are chosen in an ad hoc manner.
IF (B .LT. 1d-12 .OR. A .LT. 1d-5) THEN
  muSGS = 0.
ELSE
  muSGS = CSdeltaS2 * SQRT(B/A) * dens
END IF
END SUBROUTINE Vreman_Point


!===================================================================================================================================
!> Compute Vreman Eddy-Visosity for the volume
!===================================================================================================================================
SUBROUTINE Vreman_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems
USE MOD_EddyVisc_Vars,     ONLY: CSdeltaS2, muSGS
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
    CALL Vreman_Point(gradUx(   :,i,j,k,iElem), gradUy(:,i,j,k,iElem), gradUz(:,i,j,k,iElem), &
                           U(DENS,i,j,k,iElem),      CSdeltaS2(iElem),  muSGS(1,i,j,k,iElem))
  END DO; END DO; END DO ! i,j,k
END DO
END SUBROUTINE Vreman_Volume


!===============================================================================================================================
!> Deallocate arrays and finalize variables used by Vreman SGS model
!===============================================================================================================================
SUBROUTINE FinalizeVreman()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===============================================================================================================================
SDEALLOCATE(CSdeltaS2)
VremanInitIsDone = .FALSE.

END SUBROUTINE FinalizeVreman

END MODULE MOD_Vreman
