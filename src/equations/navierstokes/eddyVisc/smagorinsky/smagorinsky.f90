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
!> Subroutines necessary for calculating Smagorinsky Eddy-Viscosity
!===================================================================================================================================
MODULE MOD_Smagorinsky
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitSmagorinsky
   MODULE PROCEDURE InitSmagorinsky
END INTERFACE

INTERFACE FinalizeSmagorinsky
   MODULE PROCEDURE FinalizeSmagorinsky
END INTERFACE

PUBLIC::InitSmagorinsky,Smagorinsky,FinalizeSmagorinsky
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get some parameters needed by Smagorinsky modules and initialize Smagorinsky
!===================================================================================================================================
SUBROUTINE InitSmagorinsky()
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
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.SmagorinskyInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitSmagorinsky not ready to be called or already called.")
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SMAGORINSKY...'

! Read the variables used for LES model
! Smagorinsky model
CS     = GETREAL('CS')
IF(testcase.EQ."channel") THEN
  ! Do Van Driest style damping or not
  VanDriest = GETLOGICAL('VanDriest','.FALSE.')
END IF

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

SmagorinskyInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SMAGORINSKY DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitSmagorinsky

!===================================================================================================================================
!> Compute Smagorinsky Eddy-Visosity at a given point in the volume
!===================================================================================================================================
SUBROUTINE Smagorinsky(iElem,i,j,k,muSGS)
! MODULES
USE MOD_PreProc
USE MOD_EddyVisc_Vars,     ONLY:deltaS,CS,VanDriest,muSGSmax
USE MOD_EOS_Vars,          ONLY:mu0
USE MOD_Mesh_Vars,         ONLY:Elem_xGP
USE MOD_Lifting_Vars,      ONLY:gradUx,gradUy,gradUz
USE MOD_DG_Vars,           ONLY:U
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem             !< index of current element
!> indices of the current volume point
INTEGER,INTENT(IN)                        :: i,j,k
!> gradients of the velocities w.r.t. all directions
REAL,INTENT(INOUT)                        :: muSGS             !< local SGS viscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: S_eN
REAL                :: yPlus,damp
!===================================================================================================================================
! Already take the square root of 2 into account here
S_eN = 2*(gradUx(2,i,j,k,iElem)**2. + gradUy(3,i,j,k,iElem)**2. + gradUz(4,i,j,k,iElem)**2.)
S_eN = S_eN + ( gradUy(2,i,j,k,iElem) + gradUx(3,i,j,k,iElem) )**2.
S_eN = S_eN + ( gradUz(2,i,j,k,iElem) + gradUx(4,i,j,k,iElem) )**2.
S_eN = S_eN + ( gradUz(3,i,j,k,iElem) + gradUy(4,i,j,k,iElem) )**2.
S_eN = SQRT(S_eN)
! Smagorinsky model
IF(.NOT.VanDriest)THEN
  damp=1.
ELSE
  yPlus = (1. - ABS(Elem_xGP(2,i,j,k,iElem)))/mu0
  damp = 1. - EXP(-yPlus/26.) ! Van Driest damping factor
END IF
muSGS= (damp**2*CS*deltaS(iElem))**2. * S_eN*U(1,i,j,k,iElem)
muSGSmax(iElem) = MAX(muSGS,muSGSmax(iElem))
END SUBROUTINE Smagorinsky

!===============================================================================================================================
!> Deallocate arrays and finalize variables used by Smagorinsky SGS model
!===============================================================================================================================
SUBROUTINE FinalizeSmagorinsky()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!===============================================================================================================================
SmagorinskyInitIsDone = .FALSE.
END SUBROUTINE FinalizeSmagorinsky

END MODULE MOD_Smagorinsky
