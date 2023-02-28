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
!> Subroutines for the sigma-model proposed in
!>  - Nicoud, Franck, et al. "Using singular values to build a subgrid-scale model for large eddy simulations."
!>    Physics of fluids 23.8 (2011): 085106.
!===================================================================================================================================
MODULE MOD_SigmaModel
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitSigmaModel
   MODULE PROCEDURE InitSigmaModel
END INTERFACE

INTERFACE SigmaModel
   MODULE PROCEDURE SigmaModel_Point
   MODULE PROCEDURE SigmaModel_Volume
END INTERFACE

INTERFACE FinalizeSigmaModel
   MODULE PROCEDURE FinalizeSigmaModel
END INTERFACE

PUBLIC::InitSigmaModel,SigmaModel_Volume,FinalizeSigmaModel
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get model parameters and initialize sigma-model
!===================================================================================================================================
SUBROUTINE InitSigmaModel()
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
INTEGER :: i,iElem,j,k
REAL    :: CellVol
!===================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.SigmaModelInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitSigmaModel not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SIGMA-MODEL...'

! Read model coefficient
CS = GETREAL('CS')

! Allocate precomputed (model constant*filter width)**2
ALLOCATE(CSdeltaS2(nElems))

! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )
DO iElem=1,nElems
  CellVol = 0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CellVol = CellVol +wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
  DeltaS(iElem)    = CellVol**(1./3.)  / (REAL(PP_N)+1.)
  CSdeltaS2(iElem) = (CS * deltaS(iElem))**2
END DO

SigmaModelInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SIGMA-MODEL DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitSigmaModel

!===================================================================================================================================
!> Compute sigma-Model Eddy-Visosity
!===================================================================================================================================
SUBROUTINE SigmaModel_Point(gradUx,gradUy,gradUz,dens,CSdeltaS2,muSGS)
! MODULES
USE MOD_EddyVisc_Vars,     ONLY:CS
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx, gradUy, gradUz   !> Gradients in x,y,z directions
REAL                          ,INTENT(IN)  :: dens       !> pointwise density
REAL                          ,INTENT(IN)  :: CSdeltaS2  !> filter width
REAL                          ,INTENT(OUT) :: muSGS      !> pointwise eddyviscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DSYEV
! LOCAL VARIABLES
INTEGER            :: info
REAL               :: d_model
REAL               :: sigma(3), lambda(3)
REAL               :: G_mat(3,3)
REAL               :: work(9)  ! lapack work array
!===================================================================================================================================
G_Mat(1,1) = gradUx(LIFT_VEL1)*gradUx(LIFT_VEL1) + gradUx(LIFT_VEL2)*gradUx(LIFT_VEL2) + gradUx(LIFT_VEL3)*gradUx(LIFT_VEL3)
G_Mat(1,2) = gradUx(LIFT_VEL1)*gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2)*gradUy(LIFT_VEL2) + gradUx(LIFT_VEL3)*gradUy(LIFT_VEL3)
G_Mat(1,3) = gradUx(LIFT_VEL1)*gradUz(LIFT_VEL1) + gradUx(LIFT_VEL2)*gradUz(LIFT_VEL2) + gradUx(LIFT_VEL3)*gradUz(LIFT_VEL3)
G_Mat(2,1) = G_Mat(1,2)
G_Mat(2,2) = gradUy(LIFT_VEL1)*gradUy(LIFT_VEL1) + gradUy(LIFT_VEL2)*gradUy(LIFT_VEL2) + gradUy(LIFT_VEL3)*gradUy(LIFT_VEL3)
G_Mat(2,3) = gradUy(LIFT_VEL1)*gradUz(LIFT_VEL1) + gradUy(LIFT_VEL2)*gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3)*gradUz(LIFT_VEL3)
G_Mat(3,1) = G_Mat(1,3)
G_Mat(3,2) = G_Mat(2,3)
G_Mat(3,3) = gradUz(LIFT_VEL1)*gradUz(LIFT_VEL1) + gradUz(LIFT_VEL2)*gradUz(LIFT_VEL2) + gradUz(LIFT_VEL3)*gradUz(LIFT_VEL3)

! LAPACK
CALL DSYEV('N','U',3,G_Mat,3,lambda,work,9,info)
IF(info .NE. 0) THEN
  WRITE(*,*)'Eigenvalue Computation failed 3D',info
  d_model = 0.
ELSE
  sigma = SQRT(MAX(0.,lambda)) ! ensure were not negative
  d_model = (sigma(1)*(sigma(3)-sigma(2))*(sigma(2)-sigma(1)))/(sigma(3)**2)
END IF
! Sigma-Model
muSGS = CSdeltaS2 * d_model * dens
END SUBROUTINE SigmaModel_Point

!===================================================================================================================================
!> Compute SigmaModel Eddy-Visosity for the volume
!===================================================================================================================================
SUBROUTINE SigmaModel_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems
USE MOD_EddyVisc_Vars,     ONLY: muSGS, CSdeltaS2
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
    CALL SigmaModel_Point(gradUx(    :,i,j,k,iElem), gradUy(:,i,j,k,iElem), gradUz(:,i,j,k,iElem), &
                               U(DENS,i,j,k,iElem),       CSdeltaS2(iElem),  muSGS(1,i,j,k,iElem))
  END DO; END DO; END DO ! i,j,k
END DO
END SUBROUTINE SigmaModel_Volume

!===============================================================================================================================
!> Deallocate arrays and finalize variables used by SigmaModel SGS model
!===============================================================================================================================
SUBROUTINE FinalizeSigmaModel()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===============================================================================================================================
SDEALLOCATE(CSdeltaS2)
SigmaModelInitIsDone = .FALSE.

END SUBROUTINE FinalizeSigmaModel

END MODULE MOD_SigmaModel
