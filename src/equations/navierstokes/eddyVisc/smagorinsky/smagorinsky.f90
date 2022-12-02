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
#include "flexi.h"
#include "eos.h"

!===================================================================================================================================
!> Subroutines necessary for calculating Smagorinsky Eddy-Viscosity, originally derived in
!>   - Smagorinsky, Joseph. "General circulation experiments with the primitive equations: I. The basic experiment."
!>     Monthly weather review 91.3 (1963): 99-164.
!>
!> The Van-Driest damping for the Smagorinsky model for channel flow is originally used in
!>   - Moin, Parviz, and John Kim. "Numerical investigation of turbulent channel flow."
!>     Journal of fluid mechanics 118 (1982): 341-377.
!===================================================================================================================================
MODULE MOD_Smagorinsky
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitSmagorinsky
   MODULE PROCEDURE InitSmagorinsky
END INTERFACE

INTERFACE Smagorinsky
   MODULE PROCEDURE Smagorinsky_Point
   MODULE PROCEDURE Smagorinsky_Volume
END INTERFACE

INTERFACE FinalizeSmagorinsky
   MODULE PROCEDURE FinalizeSmagorinsky
END INTERFACE

PUBLIC::InitSmagorinsky, Smagorinsky_Volume, FinalizeSmagorinsky
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get model parameters and initialize Smagorinsky model
!===================================================================================================================================
SUBROUTINE InitSmagorinsky()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EddyVisc_Vars
USE MOD_ReadInTools        ,ONLY: GETREAL,GETLOGICAL
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone,wGP
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone,nElems,sJ,Elem_xGP
USE MOD_EOS_Vars           ,ONLY: mu0
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: CellVol,yPlus
!===================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.SmagorinskyInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitSmagorinsky not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SMAGORINSKY...'

! Read model coefficient
CS = GETREAL('CS')
! Do Van Driest style damping or not (only for channel!)
VanDriest = GETLOGICAL('VanDriest')

ALLOCATE(damp(1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
damp = 1.

! Smago: (damp*CS*deltaS)**2 * S_eN * dens
! Precompute first term and store in damp
! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )
DO iElem=1,nElems
  CellVol = 0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CellVol = CellVol + wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
    IF (VanDriest)THEN
      yPlus = (1. - ABS(Elem_xGP(2,i,j,k,iElem))) / mu0   ! y-dir
      damp(1,i,j,k,iElem) = 1. - EXP(-yPlus/26.) ! Van Driest damping factor
    END IF
  END DO; END DO; END DO
  DeltaS(iElem) = CellVol**(1./3.)  / (REAL(PP_N)+1.)

  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    damp(1,i,j,k,iElem) = (damp(1,i,j,k,iElem) * CS * deltaS(iElem))**2
  END DO; END DO; END DO
END DO

SmagorinskyInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SMAGORINSKY DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitSmagorinsky

!===================================================================================================================================
!> Compute Smagorinsky Eddy-Visosity
!===================================================================================================================================
PPURE SUBROUTINE Smagorinsky_Point(gradUx,gradUy,gradUz,dens,damp,muSGS)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx, gradUy, gradUz   !> Gradients in x,y,z directions
REAL                          ,INTENT(IN)  :: dens    !> pointwise density
REAL                          ,INTENT(IN)  :: damp    !> constant factor (damp*CS*deltaS)**2
REAL                          ,INTENT(OUT) :: muSGS   !> pointwise eddyviscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: S_eN
!===================================================================================================================================
! Already take the square root of 2 into account here
S_eN = SQRT ( 2.*(gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2 + gradUz(LIFT_VEL3)**2) &
              + ( gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2) )**2                    &
              + ( gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3) )**2                    &
              + ( gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3) )**2 )
! Smagorinsky model: (damp * CS * deltaS)**2 * S_eN * rho
! we store the first constant term in damp
muSGS = damp * S_eN * dens
END SUBROUTINE Smagorinsky_Point

!===================================================================================================================================
!> Compute Smagorinsky Eddy-Visosity for the volume
!===================================================================================================================================
SUBROUTINE Smagorinsky_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems
USE MOD_EddyVisc_Vars,     ONLY: damp, muSGS
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
    CALL Smagorinsky_Point(gradUx(:   ,i,j,k,iElem), gradUy(:,i,j,k,iElem), gradUz(:,i,j,k,iElem), &
                                U(DENS,i,j,k,iElem),   damp(1,i,j,k,iElem),  muSGS(1,i,j,k,iElem))
  END DO; END DO; END DO ! i,j,k
END DO
END SUBROUTINE Smagorinsky_Volume

!===============================================================================================================================
!> Deallocate arrays and finalize variables used by Smagorinsky SGS model
!===============================================================================================================================
SUBROUTINE FinalizeSmagorinsky()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===============================================================================================================================
SmagorinskyInitIsDone = .FALSE.
END SUBROUTINE FinalizeSmagorinsky

END MODULE MOD_Smagorinsky
