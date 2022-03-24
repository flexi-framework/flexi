!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
!> Module containing the main procedures needed to build the connection between the WM interface and the boundary points.
!===================================================================================================================================
MODULE MOD_Newton
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE Newton
  MODULE PROCEDURE Newton_2D
  MODULE PROCEDURE Newton_3D
END INTERFACE

PUBLIC :: Newton

CONTAINS

!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE Newton_2D(NIn,xIn,dxIn,Xi_CLN,wBary_CLN,xZeroIn,Xi,epsOut,LagOut,FOut)
! MODULES                                                                                                                          !
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys
USE MOD_Mathtools,         ONLY: INVERSE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: NIn
REAL,INTENT(IN)           :: xIn(        :)
REAL,INTENT(IN)           :: dxIn(  1:2,1:2,0:NIn,0:Nin)
REAL,INTENT(IN)           :: Xi_CLN(     :)
REAL,INTENT(IN)           :: wBary_CLN(  :)
REAL,INTENT(IN)           :: xZeroIn(   1:2,0:NIn,0:NIn)
REAL,INTENT(INOUT)        :: Xi(    1:2)
REAL,INTENT(OUT),OPTIONAL :: EpsOut
REAL,INTENT(OUT),OPTIONAL :: LagOut(1:2,0:NIn)
REAL,INTENT(OUT),OPTIONAL :: FOut(  1:2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER         :: N_dim  = 2
REAL,PARAMETER            :: epsTol = 1.E-8
INTEGER                   :: i,j,iter
REAL                      :: eps_F
REAL                      :: F(   1:N_dim)
REAL                      :: Jac( 1:N_dim,1:N_dim)
REAL                      :: sJac(1:N_dim,1:N_dim)
REAL                      :: LagVol( 0:NIn,0:NIn)
REAL                      :: Lag(1:2,0:NIn)
!===================================================================================================================================

CALL LagrangeInterpolationPolys(Xi(1),NIn,Xi_CLN,wBary_CLN,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NIn,Xi_CLN,wBary_CLN,Lag(2,:))

F = -xIn

DO j=0,NIn; DO i=0,NIn
  LagVol(i,j) = Lag(1,i)*Lag(2,j)
  F = F + xZeroIn(1:2,i,j)*LagVol(i,j)
END DO; END DO

eps_F = epsTol*SUM(F*F) ! relative error to initial guess
iter  = 0
IF (PRESENT(epsOut)) epsOut = eps_F

DO WHILE ((SUM(F*F).GT.eps_F).AND.(iter.LT.100))
  iter = iter+1
  ! Compute F Jacobian dx/dXi
  Jac  = 0.

  DO j=0,NIn; DO i=0,NIn
    Jac = Jac + dxIn(1:2,1:2,i,j)*LagVol(i,j)
  END DO; END DO

  ! Compute inverse of Jacobian
  sJac = INVERSE(Jac)

  ! Iterate Xi using Newton step
  Xi = Xi - MATMUL(sJac,F)
  ! if Newton gets outside reference space range [-1,1], exit.
  ! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
  IF((iter.GT.3).AND.(ANY(ABS(Xi).GT.1.2))) EXIT

  ! Compute function value
  CALL LagrangeInterpolationPolys(Xi(1),NIn,Xi_CLN,wBary_CLN,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),NIn,Xi_CLN,wBary_CLN,Lag(2,:))

  ! F(xi) = x(xi) - xIn
  F = -xIn

  DO j=0,NIn; DO i=0,NIn
    LagVol(i,j) = Lag(1,i)*Lag(2,j)
    F = F + xZeroIn(1:2,i,j)*LagVol(i,j)
  END DO; END DO
END DO

IF (PRESENT(LagOut)) LagOut = Lag
IF (PRESENT(FOut))   FOut   = F

END SUBROUTINE Newton_2D


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE Newton_3D(NIn,xIn,dxIn,Xi_CLN,wBary_CLN,xZeroIn,Xi,epsOut,LagOut,FOut)
! MODULES                                                                                                                          !
use mod_globals
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys
USE MOD_Mathtools,         ONLY: INVERSE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: NIn
REAL,INTENT(IN)           :: xIn(        :)
REAL,INTENT(IN)           :: dxIn(  1:3,1:3,0:NIn,0:Nin,0:NIn)
REAL,INTENT(IN)           :: Xi_CLN(     :)
REAL,INTENT(IN)           :: wBary_CLN(  :)
REAL,INTENT(IN)           :: xZeroIn(   1:3,0:NIn,0:NIn,0:NIn)
REAL,INTENT(INOUT)        :: Xi(    1:3)
REAL,INTENT(OUT),OPTIONAL :: EpsOut
REAL,INTENT(OUT),OPTIONAL :: LagOut(1:3,0:NIn)
REAL,INTENT(OUT),OPTIONAL :: FOut(  1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER         :: N_dim  = 3
REAL,PARAMETER            :: epsTol = 1.E-8
INTEGER                   :: i,j,k,iter
REAL                      :: eps_F
REAL                      :: F(   1:N_dim)
REAL                      :: Jac( 1:N_dim,1:N_dim)
REAL                      :: sJac(1:N_dim,1:N_dim)
REAL                      :: LagVol( 0:NIn,0:NIn,0:NIn)
REAL                      :: Lag(1:3,0:NIn)
!===================================================================================================================================

CALL LagrangeInterpolationPolys(Xi(1),NIn,Xi_CLN,wBary_CLN,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NIn,Xi_CLN,wBary_CLN,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),NIn,Xi_CLN,wBary_CLN,Lag(3,:))

F = -xIn

DO k=0,NIn; DO j=0,NIn; DO i=0,NIn
  LagVol(i,j,k) = Lag(1,i)*Lag(2,j)*Lag(3,k)
  F = F + xZeroIn(1:3,i,j,k)*LagVol(i,j,k)
END DO; END DO; END DO

eps_F = epsTol*SUM(F*F) ! relative error to initial guess
iter  = 0
IF (PRESENT(epsOut)) epsOut = eps_F

DO WHILE ((SUM(F*F).GT.eps_F).AND.(iter.LT.100))
  iter = iter+1
  ! Compute F Jacobian dx/dXi
  Jac  = 0.

  DO k=0,NIn; DO j=0,NIn; DO i=0,NIn
    Jac = Jac + dxIn(1:3,1:3,i,j,k)*LagVol(i,j,k)
  END DO; END DO; END DO

  ! Compute inverse of Jacobian
  sJac = INVERSE(Jac)

  ! Iterate Xi using Newton step
  Xi = Xi - MATMUL(sJac,F)
  ! if Newton gets outside reference space range [-1,1], exit.
  ! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
  IF((iter.GT.3).AND.(ANY(ABS(Xi).GT.1.2))) EXIT

  ! Compute function value
  CALL LagrangeInterpolationPolys(Xi(1),NIn,Xi_CLN,wBary_CLN,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),NIn,Xi_CLN,wBary_CLN,Lag(2,:))
  CALL LagrangeInterpolationPolys(Xi(3),NIn,Xi_CLN,wBary_CLN,Lag(3,:))

  ! F(xi) = x(xi) - xIn
  F = -xIn

  DO k=0,NIn; DO j=0,NIn; DO i=0,NIn
    LagVol(i,j,k) = Lag(1,i)*Lag(2,j)*Lag(3,k)
    F = F + xZeroIn(1:3,i,j,k)*LagVol(i,j,k)
  END DO; END DO; END DO
END DO

IF (PRESENT(LagOut)) LagOut = Lag
IF (PRESENT(FOut))   FOut   = F

END SUBROUTINE Newton_3D

END MODULE MOD_Newton
