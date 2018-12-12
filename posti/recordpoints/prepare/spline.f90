#include "flexi.h"

!===================================================================================================================================
!> Module to handle the Recordpoints
!===================================================================================================================================
MODULE MOD_Spline
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE GetSpline
  MODULE PROCEDURE GetSpline
END INTERFACE

INTERFACE GetEquiPoints
  MODULE PROCEDURE GetEquiPoints
END INTERFACE

INTERFACE EvalEquiError
  MODULE PROCEDURE EvalEquiError
END INTERFACE

INTERFACE EvalSpline
  MODULE PROCEDURE EvalSpline
END INTERFACE

INTERFACE EvalSplineDeriv
  MODULE PROCEDURE EvalSplineDeriv
END INTERFACE

PUBLIC :: GetSpline,GetEquiPoints,EvalSpline,EvalSplineDeriv,EvalEquiError
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Define a spline in ndim dimensions using nCP controlpoints stored in xCP. Return the coefficients for all nCP-1 segments and
!> the parameterization variable s_out.
!> If s_in is specified, s_out=s_in. If not, the sum of the linear distances between the xCP is used as s.
!===================================================================================================================================
SUBROUTINE GetSpline(ndim,nCP,xCP,coeff,s_out,s_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: ndim,nCP
REAL,INTENT(IN)                 :: xCP(ndim,nCP)
REAL,INTENT(IN),OPTIONAL        :: s_in(nCP)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: s_out(nCP)
REAL,INTENT(OUT)                :: coeff(ndim,4,nCp-1)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nSp,dim,i
REAL,DIMENSION(nCP)             :: a,c
REAL,DIMENSION(nCP-1)           :: b,d,h
REAL,DIMENSION(nCP-2)           :: r,a_,b_
REAL,DIMENSION(nCP-3)           :: c_
!===================================================================================================================================

nSp=nCP-1 ! number of piecewise polynomials

!if a parametrization coordinate is present, use that one
IF(PRESENT(s_in)) THEN
  s_out=s_in
ELSE
  ! initial parametrization with independent variable s
  s_out(1)=0.
  DO i=2,nSp+1
    s_out(i)=s_out(i-1)+SQRT(SUM((xCP(1:ndim,i)-xCP(1:ndim,i-1))**2))
  END DO! i=2,nSp+1
END IF!(PRESENT(s_in))
DO i=2,nSp+1
  h(i-1)=s_out(i)-s_out(i-1)
  IF(ABS(h(i-1)).LT.1e-10) &
    STOP 'points on spline appear to be the same!!'
END DO! i=2,nSp+1

!get coeffients dimension-wise
DO dim=1,ndim
  c(1)=0.
  c(nSp+1)=0.
  DO i=1,nSp+1
    a(i)=xCP(dim,i)
  END DO
  ! LGS matrix

  !diagonal
  DO i=1,nSp-1
    b_(i)=2*(h(i)+h(i+1))
  END DO
  !right hand side
  DO i=1,nSp-1
    r(i)=3/h(i+1)*(a(i+2)-a(i+1)) -3/h(i)*(a(i+1)-a(i))
  END DO
  IF(nSp .GT. 2) THEN
    !upper diagonal
    DO i=1,nSp-2
        c_(i)=h(i)
    END DO
    !lower diagonal
    a_(1)=0.
    DO i=2,nSp-1
        a_(i)=h(i-1)
    END DO
    !solve
    !thomas alg.
    CALL thomas( nSp-1,a_,b_,c_,r,c(2:nSP) )
  ELSEIF(nSp.EQ.2)THEN
    c(2)=r(1)/b_(1)
  END IF

  DO i=1,nSp
    b(i)=1/h(i)*(a(i+1)-a(i)) - h(i)/3*(c(i+1) + 2*c(i))
    d(i)=1/(3*h(i))*(c(i+1)-c(i))
  END DO
  coeff(dim,1,:)=a(1:nSp)
  coeff(dim,2,:)=b(1:nSp)
  coeff(dim,3,:)=c(1:nSp)
  coeff(dim,4,:)=d(1:nSp)
END DO! dim=1,ndim

END SUBROUTINE GetSpline



!===================================================================================================================================
!> Given a set of nP_in points xP_in in ndim dimension, return  the desired nP_out equidistant points over the arc length of
!> the spline that connects the input points xP_in. t_equi contains the corresponding coordinates in the parameter space of the
!> spline.
!===================================================================================================================================
SUBROUTINE GetEquiPoints(ndim,nP_in,nP_out,xP_in,xP_out,t_equi)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: ndim,nP_in,nP_out
REAL,INTENT(IN)                 :: xP_in(ndim,nP_in)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: xP_out(ndim,nP_out),t_equi(nP_out)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,iP,iSp,nSuper,nSp,iter
REAL                            :: t_ini(nP_in),t_tmp(1,nP_in),s(nP_in),s_tmp(nP_in)
REAL                            :: t_loc,s_loc,x_loc(ndim),x_loc_old(ndim)
REAL                            :: coeff(ndim,4,nP_in-1)
REAL                            :: t_coeff(1,4,nP_in-1)
REAL                            :: ds_equi,ds(nP_in-1),EquiErr!,ds_measured,s_measured
!===================================================================================================================================
nSp=nP_in-1
!iterative procedure acc. to Hindenlang Diss p.83p
xP_out = xP_in
DO iter=1,10
  ! 1. get spline through the given points
  CALL GetSpline(3,nP_in,xP_out,coeff,t_ini)

  ! 2. calculate arclength on supersampled points s(t)
  nSuper=10
  s=0.
  ds=0.
  x_loc=xP_out(:,1)
  DO iSp=1,nSp
    DO i=2,nSuper
      t_loc=t_ini(iSP)+REAL(i-1)/REAL(nSuper-1)*(t_ini(iSp+1)-t_ini(iSp))
      x_loc_old=x_loc
      x_loc(:)=coeff(:,1,iSp)+coeff(:,2,iSp)*(t_loc-t_ini(iSp)) &
              +coeff(:,3,iSp)*(t_loc-t_ini(iSp))**2 + coeff(:,4,iSp)*(t_loc-t_ini(iSp))**3
      ds(iSp)=ds(iSp)+SQRT(SUM((x_loc(1:ndim)-x_loc_old(1:ndim))**2))
    END DO
    s(iSp+1)=s(iSp)+ds(iSp)
  END DO
  ds_equi=1./REAL(nP_out-1)*s(nP_in)
  ! evaluate the error (departure from equidistant spacing in s)
  EquiErr=0.
  DO iSp=1,nSp
    EquiErr=MAX(EquiErr,ABS(ds(iSp)-ds_equi))
  END DO
  !WRITE(*,*)'iter, EquiErr= ',iter, EquiErr
  IF(EquiErr.LT.2e-16) RETURN

  ! 3. create inverse mapping t(s) as spline
  t_tmp(1,:)=t_ini
  s_tmp = s
  CALL GetSpline(1,nP_in,t_tmp,t_coeff,s,s_in=s_tmp)
  !get new t at equidistant s arclength points
  DO iP=1,nP_out
    s_loc=REAL(iP-1)*ds_equi
    CALL EvalSpline(1,nP_in,s_loc,s,t_coeff,t_tmp(:,iP)) ! t(s) at the equidistant point
  END DO! i=1,nP_out
  t_equi(:)=t_tmp(1,:)

  ! 4. evaluate the mapping at equidistant points x(t(s))
  DO iP=1,nP_out
    t_loc=t_equi(iP)
    CALL EvalSpline(3,nP_in,t_loc,t_ini,coeff,xP_out(:,iP)) ! x(t(s)) at the equidistant point
  END DO! i=1,nP_out
END DO !iter
WRITE(UNIT_StdOut,'(A,F12.5)')' Warning: equidistant routine did not converge! , error = ',EquiErr
END SUBROUTINE GetEquiPoints


!===================================================================================================================================
!> Given a set of nP_in points xP_in in ndim dimension, return the desired nP_out equidistant points over the arc length of
!> the spline that connects the input points xP_in. t_equi contains the corresponding coordinates in the parameter space of the
!> spline.
!===================================================================================================================================
SUBROUTINE EvalEquiError(ndim,nP_in,xP_in,EquiErr)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: ndim,nP_in
REAL,INTENT(IN)                 :: xP_in(ndim,nP_in)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: EquiErr
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,iSp,nSuper,nSp
REAL                            :: t_ini(nP_in),s(nP_in)
REAL                            :: t_loc,x_loc(ndim),x_loc_old(ndim)
REAL                            :: coeff(ndim,4,nP_in-1)
REAL                            :: ds(nP_in-1),ds_equi
!===================================================================================================================================
nSp=nP_in-1
! 1. get spline through the given points
CALL GetSpline(3,nP_in,xP_in,coeff,t_ini)

! 2. calculate arclength on supersampled points s(t)
nSuper=10
s=0.
ds=0.
x_loc=xP_in(:,1)
DO iSp=1,nSp
  DO i=2,nSuper
    t_loc=t_ini(iSP)+REAL(i-1)/REAL(nSuper-1)*(t_ini(iSp+1)-t_ini(iSp))
    x_loc_old=x_loc
    x_loc(:)=coeff(:,1,iSp)+coeff(:,2,iSp)*(t_loc-t_ini(iSp)) &
            +coeff(:,3,iSp)*(t_loc-t_ini(iSp))**2 + coeff(:,4,iSp)*(t_loc-t_ini(iSp))**3
    ds(iSp)=ds(iSp)+SQRT(SUM((x_loc(1:ndim)-x_loc_old(1:ndim))**2))
  END DO
  s(iSp+1)=s(iSp)+ds(iSp)
END DO
ds_equi=1./REAL(nP_in-1)*s(nP_in)
! evaluate the error (departure from equidistant spacing in s)
EquiErr=0.
DO iSp=1,nSp
  EquiErr=MAX(EquiErr,ABS(ds(iSp)-ds_equi))
END DO
END SUBROUTINE EvalEquiError


!===================================================================================================================================
!> Evaluate a spline at the parametrization coordinate position s_in and return the solution in x_loc.
!> The spline is defined using nCP,s and coeff (can be calculated using GetSpline)
!===================================================================================================================================
SUBROUTINE EvalSpline(ndim,nCP,s_in,s,coeff,x_loc)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: ndim,nCP
REAL,INTENT(IN)                 :: s_in,s(nCP),coeff(ndim,4,nCp-1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: x_loc(ndim)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSp
REAL                            :: s_loc
!===================================================================================================================================
!find the interval s_i<s_loc<s_i+1
s_loc=MAX(s_in,s(1))
s_loc=MIN(s_loc,s(nCP))
iSp=1
DO WHILE(s_loc .GT. s(iSp+1))
  iSp=iSp+1
END DO
!iSp=MAX(1,iSp-1)

x_loc(:)=coeff(:,1,iSp)+coeff(:,2,iSp)*(s_loc-s(iSp)) &
          +coeff(:,3,iSp)*(s_loc-s(iSp))**2 + coeff(:,4,iSp)*(s_loc-s(iSp))**3

END SUBROUTINE EvalSpline


!===================================================================================================================================
!> Evaluate the derivative of a given spline (ndim,nCP,s,coeff) at the parameter space location s_in
!===================================================================================================================================
SUBROUTINE EvalSplineDeriv(ndim,nCP,s_in,s,coeff,dx_loc)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: ndim,nCP
REAL,INTENT(IN)                 :: s_in,s(nCP),coeff(ndim,4,nCp-1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: dx_loc(ndim)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSp
REAL                            :: s_loc
!===================================================================================================================================
!find the interval s_i<s_loc<s_i+1
s_loc=MAX(s_in,s(1))
s_loc=MIN(s_loc,s(nCP))
iSp=1
DO WHILE(s_loc .GT. s(iSp+1))
  iSp=iSp+1
END DO

dx_loc(:)=coeff(:,2,iSp) &
          +2.*coeff(:,3,iSp)*(s_loc-s(iSp)) + 3.*coeff(:,4,iSp)*(s_loc-s(iSp))**2
END SUBROUTINE EvalSplineDeriv


!===================================================================================================================================
!> Solve the tridiagonal equation system of the form
!> | b_1 c_1                              .  |
!> | a_2 b_2 c_2                          .  |
!> |  .                                   .  |
!> |  .                                      |
!> |  .                                      |
!> |                a_(m-1) b_(m-1) c_(m-1)  |
!> |                        a_m     b_m      |
!===================================================================================================================================
SUBROUTINE thomas( m,a_,b_,c_,r,x )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: m
REAL,INTENT(INOUT),DIMENSION(m) :: a_,b_,r
REAL,INTENT(IN),DIMENSION(m-1)  :: c_
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: x(m)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i
REAL                            :: q
!===================================================================================================================================
DO i=2,m
  q=a_(i)/b_(i-1)
  b_(i)=b_(i)-q*c_(i-1)
  r(i)=r(i)-q*r(i-1)
END DO
x(m)=r(m)/b_(m)
DO i=m-1,1,-1
  x(i)=(r(i)-c_(i)*x(i+1))/b_(i)
END DO
END SUBROUTINE thomas

END MODULE MOD_Spline
