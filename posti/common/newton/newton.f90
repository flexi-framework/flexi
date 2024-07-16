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

!===================================================================================================================================
!> Module containing a Newton algorithm to find the coordinates in reference space for a point in physical coordinates
!===================================================================================================================================
MODULE MOD_Newton
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Newton
  MODULE PROCEDURE Newton_2D
  MODULE PROCEDURE Newton_3D
END INTERFACE

PUBLIC :: Newton
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Employ a Newton algorithm to find the parametric coordinates xi of an interpolation point in an element defined by CL points on
!> NIn. We try to solve F(xi) = x(xi) - x_InterpolationPoint = 0 for the parametric coordinates in 2D.
!> Newton iteration: xi_(n+1) = xi_n - (J(xi_n))^(-1)*F(xi_n), the Jacobian is the derivative of the mesh coordinates w.r.t. the
!> parametric coordinates.
!===================================================================================================================================
SUBROUTINE Newton_2D(NIn,xIn,dxIn,Xi_CLN,wBary_CLN,xZeroIn,Xi,epsOut,LagOut,FOut)
! MODULES                                                                                                                          !
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys
USE MOD_Mathtools,         ONLY: INVERSE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: NIn                         !> Polynomial degree of coordinate mapping
REAL,INTENT(IN)           :: xIn(        :)              !> Physical coordinates of point that is searched
REAL,INTENT(IN)           :: dxIn(  1:2,1:2,0:NIn,0:Nin) !> Jacobian of coordinate mapping: dx_i/dxi_j
REAL,INTENT(IN)           :: Xi_CLN(     :)              !> Interpolation points of Lagrange polynomials of degree NIn
REAL,INTENT(IN)           :: wBary_CLN(  :)              !> Barycentric  weights of Lagrange polynomials of degree NIn
REAL,INTENT(IN)           :: xZeroIn(1:2,0:NIn,0:NIn)    !> Physical coordinates at interpolation points Xi_CLN
REAL,INTENT(INOUT)        :: Xi(    1:2)                 !> Input:  Intial guess of reference coordiantes for Newton
                                                         !> Output: Found reference coordinates of point XIn by Newton
REAL,INTENT(OUT),OPTIONAL :: EpsOut                      !> relative error bound for convergence
REAL,INTENT(OUT),OPTIONAL :: LagOut(1:2,0:NIn)           !> Lagrange polynomials evaluated at Xi
REAL,INTENT(OUT),OPTIONAL :: FOut(  1:2)                 !> L2-error of final Xi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER         :: nDim    = 2
INTEGER,PARAMETER         :: maxIter = 100
REAL,PARAMETER            :: epsTol  = 1.E-8
INTEGER                   :: i,j,iter
REAL                      :: eps_F
REAL                      :: F(   1:nDim)
REAL                      :: Jac( 1:nDim,1:nDim)
REAL                      :: LagVol( 0:NIn,0:NIn)
REAL                      :: Lag(1:2,0:NIn)
!===================================================================================================================================

! Evaluate Lagrange polynomials at xi
CALL LagrangeInterpolationPolys(Xi(1),NIn,Xi_CLN,wBary_CLN,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NIn,Xi_CLN,wBary_CLN,Lag(2,:))

! F(xi) = x(xi) - xIn
F = -xIn
DO j=0,NIn; DO i=0,NIn
  LagVol(i,j) = Lag(1,i)*Lag(2,j)
  F = F + xZeroIn(1:2,i,j)*LagVol(i,j)
END DO; END DO

eps_F = epsTol*SUM(F*F) ! relative error to initial guess
iter  = 0
IF (PRESENT(epsOut)) epsOut = eps_F

! While L2 error exceeds eps_F and max. iteration number not reached
DO WHILE ((SUM(F*F).GT.eps_F).AND.(iter.LT.maxIter))
  iter = iter+1

  ! Compute F Jacobian dx/dXi
  Jac  = 0.
  DO j=0,NIn; DO i=0,NIn
    Jac = Jac + dxIn(1:2,1:2,i,j)*LagVol(i,j)
  END DO; END DO

  ! Perform Newton step
  Xi = Xi - MATMUL(INVERSE(Jac),F)

  ! if Newton gets outside reference space range [-1,1], exit.
  ! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
  IF((iter.GT.3).AND.(ANY(ABS(Xi).GT.1.2))) EXIT

  ! Evaluate Lagrange polynomials at new xi
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
!> Employ a Newton algorithm to find the parametric coordinates xi of an interpolation point in an element defined by CL points on
!> NIn. We try to solve F(xi) = x(xi) - x_InterpolationPoint = 0 for the parametric coordinates in 3D.
!> Newton iteration: xi_(n+1) = xi_n - (J(xi_n))^(-1)*F(xi_n), the Jacobian is the derivative of the mesh coordinates w.r.t. the
!> parametric coordinates.
!===================================================================================================================================
SUBROUTINE Newton_3D(NIn,xIn,dxIn,Xi_CLN,wBary_CLN,xZeroIn,Xi,epsOut,LagOut,FOut)
! MODULES                                                                                                                          !
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys
USE MOD_Mathtools,         ONLY: INVERSE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: NIn                               !> Polynomial degree of coordinate mapping
REAL,INTENT(IN)           :: xIn(        :)                    !> Physical coordinates of point that is searched
REAL,INTENT(IN)           :: dxIn(  1:3,1:3,0:NIn,0:Nin,0:NIn) !> Jacobian of coordinate mapping: dx_i/dxi_j
REAL,INTENT(IN)           :: Xi_CLN(     :)                    !> Interpolation points of Lagrange polynomials of degree NIn
REAL,INTENT(IN)           :: wBary_CLN(  :)                    !> Barycentric  weights of Lagrange polynomials of degree NIn
REAL,INTENT(IN)           :: xZeroIn(   1:3,0:NIn,0:NIn,0:NIn) !> Physical coordinates at interpolation points Xi_CLN
REAL,INTENT(INOUT)        :: Xi(    1:3)                       !> Input:  Intial guess of reference coordiantes for Newton
                                                               !> Output: Found reference coordinates of point XIn by Newton
REAL,INTENT(OUT),OPTIONAL :: EpsOut                            !> relative error bound for convergence
REAL,INTENT(OUT),OPTIONAL :: LagOut(1:3,0:NIn)                 !> Lagrange polynomials evaluated at Xi
REAL,INTENT(OUT),OPTIONAL :: FOut(  1:3)                       !> L2-error of final Xi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER         :: nDim    = 3
INTEGER,PARAMETER         :: maxIter = 100
REAL,PARAMETER            :: epsTol  = 1.E-8
INTEGER                   :: i,j,k,iter
REAL                      :: eps_F
REAL                      :: F(   1:nDim)
REAL                      :: Jac( 1:nDim,1:nDim)
REAL                      :: LagVol( 0:NIn,0:NIn,0:NIn)
REAL                      :: Lag(1:3,0:NIn)
!===================================================================================================================================

! Evaluate Lagrange polynomials at xi
CALL LagrangeInterpolationPolys(Xi(1),NIn,Xi_CLN,wBary_CLN,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NIn,Xi_CLN,wBary_CLN,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),NIn,Xi_CLN,wBary_CLN,Lag(3,:))

! F(xi) = x(xi) - xIn
F = -xIn
DO k=0,NIn; DO j=0,NIn; DO i=0,NIn
  LagVol(i,j,k) = Lag(1,i)*Lag(2,j)*Lag(3,k)
  F = F + xZeroIn(1:3,i,j,k)*LagVol(i,j,k)
END DO; END DO; END DO

eps_F = epsTol*SUM(F*F) ! relative error to initial guess
iter  = 0
IF (PRESENT(epsOut)) epsOut = eps_F

! While L2 error exceeds eps_F and max. iteration number not reached
DO WHILE ((SUM(F*F).GT.eps_F).AND.(iter.LT.maxIter))
  iter = iter+1

  ! Compute F Jacobian dx/dXi
  Jac  = 0.
  DO k=0,NIn; DO j=0,NIn; DO i=0,NIn
    Jac = Jac + dxIn(1:3,1:3,i,j,k)*LagVol(i,j,k)
  END DO; END DO; END DO

  ! Perform Newton step
  Xi = Xi - MATMUL(INVERSE(Jac),F)

  ! if Newton gets outside reference space range [-1,1], exit.
  ! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
  IF((iter.GT.3).AND.(ANY(ABS(Xi).GT.1.2))) EXIT

  ! Evaluate Lagrange polynomials at new xi
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
