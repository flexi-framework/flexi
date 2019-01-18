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

!===================================================================================================================================
!> Subroutines necessary for calculating SigmaModel Eddy-Viscosity
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
!> Get some parameters needed by SigmaModel modules and initialize SigmaModel
!===================================================================================================================================
SUBROUTINE InitSigmaModel()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EddyVisc_Vars
USE MOD_ReadInTools        ,ONLY: GETREAL,GETLOGICAL
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone,wGP
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone,nElems,sJ
USE MOD_Testcase_Vars      ,ONLY: testcase
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
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SigmaModel...'

! Read the variables used for LES model
! SigmaModel model
CS     = GETREAL('CS')
! Do Van Driest style damping or not
VanDriest = GETLOGICAL('VanDriest','.FALSE.')

! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )
DO iElem=1,nElems
  CellVol = 0.
  DO i=0,PP_N
    DO j=0,PP_N
      DO k=0,PP_N
        CellVol = CellVol +wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
      END DO
    END DO
  END DO
  DeltaS(iElem) = ( CellVol)**(1./3.)  / (REAL(PP_N)+1.)
END DO

SigmaModelInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SigmaModel DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitSigmaModel

!===================================================================================================================================
!> Compute SigmaModel Eddy-Visosity
!===================================================================================================================================
SUBROUTINE SigmaModel_Point(gradUx,gradUy,gradUz,dens,deltaS,muSGS)
! MODULES
USE MOD_EddyVisc_Vars,     ONLY:CS
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: gradUx, gradUy, gradUz
REAL                       ,INTENT(IN)  :: dens, deltaS
REAL                       ,INTENT(OUT) :: muSGS
!-----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DSYEV
! LOCAL VARIABLES
REAL               :: sigma(3), lambda(3)
REAL               :: G_mat(3,3)
REAL               :: work(9)  !lapack work array
INTEGER            :: info
REAL               :: d_model
!===================================================================================================================================
G_Mat(1,1) = gradUx(2)*gradUx(2) + gradUx(3)*gradUx(3) + gradUx(4)*gradUx(4)
G_Mat(1,2) = gradUx(2)*gradUy(2) + gradUx(3)*gradUy(3) + gradUx(4)*gradUy(4)
G_Mat(1,3) = gradUx(2)*gradUz(2) + gradUx(3)*gradUz(3) + gradUx(4)*gradUz(4)
G_Mat(2,1) = G_Mat(1,2)
G_Mat(2,2) = gradUy(2)*gradUy(2) + gradUy(3)*gradUy(3) + gradUy(4)*gradUy(4)
G_Mat(2,3) = gradUy(2)*gradUz(2) + gradUy(3)*gradUz(3) + gradUy(4)*gradUz(4)
G_Mat(3,1) = G_Mat(1,3)
G_Mat(3,2) = G_Mat(2,3)
G_Mat(3,3) = gradUz(2)*gradUz(2) + gradUz(3)*gradUz(3) + gradUz(4)*gradUz(4)

! LAPACK
CALL DSYEV('N','U',3,G_Mat,3,lambda,work,9,info)
IF(info .NE. 0) THEN
  WRITE(*,*)'Eigenvalue Computation failed 3D',info
  d_model = 0.
ELSE
  sigma = SQRT(MAX(0.,lambda)) ! ensure were not negative
  d_model = (sigma(3)*(sigma(1)-sigma(2))*(sigma(2)-sigma(3)))/(sigma(1)**2)
END IF
! SigmaModel model
muSGS = (CS*deltaS)**2. * d_model * dens
END SUBROUTINE SigmaModel_Point

!===================================================================================================================================
!> Compute SigmaModel Eddy-Visosity for the volume
!===================================================================================================================================
SUBROUTINE SigmaModel_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems
USE MOD_EddyVisc_Vars,     ONLY: muSGS, deltaS
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
  DO k = 0,PP_NZ; DO j = 0,PP_N; DO i = 0,PP_N
    CALL SigmaModel_Point(gradUx(:,i,j,k,iElem), gradUy(:,i,j,k,iElem), gradUz(:,i,j,k,iElem), &
                                 U(1,i,j,k,iElem),        deltaS(iElem),  muSGS(1,i,j,k,iElem))
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
!===============================================================================================================================
SigmaModelInitIsDone = .FALSE.
END SUBROUTINE FinalizeSigmaModel

END MODULE MOD_SigmaModel
