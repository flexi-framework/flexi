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
!> Contains the routines for the calculation of the analytical flux jacobians of the Navier-Stokes equations
!===================================================================================================================================
MODULE MOD_Jacobian
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE dConsdPrimTemp
  MODULE PROCEDURE dConsdPrimTemp
END INTERFACE

INTERFACE dPrimTempdCons
  MODULE PROCEDURE dPrimTempdCons
END INTERFACE

INTERFACE EvalAdvFluxJacobianPoint
  MODULE PROCEDURE EvalAdvFluxJacobianPoint
END INTERFACE

PUBLIC::EvalAdvFluxJacobian,EvalAdvFluxJacobianPoint
#if PARABOLIC
PUBLIC::EvalDiffFluxJacobian
PUBLIC::EvalFluxGradJacobian
#endif
PUBLIC::dConsdPrimTemp,dPrimTempdCons
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Navier-Stokes-Problem:
!> The Jacobian of the advective Flux with respect to the conservative variables U
!===================================================================================================================================
SUBROUTINE EvalAdvFluxJacobian(U,UPrim,fJac,gJac,hJac)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars           ,ONLY:nDOFElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------  
REAL,DIMENSION(PP_nVar,nDOFElem),INTENT(IN)              :: U               !< Conservative solution
REAL,DIMENSION(PP_nVarPrim,nDOFElem),INTENT(IN)          :: UPrim           !< Primitive solution
REAL,DIMENSION(PP_nVar,PP_nVar,nDOFElem),INTENT(OUT)     :: fJac,gJac,hJac  !< Derivative of the physical fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!===================================================================================================================================
DO i=1,nDOFElem
  CALL EvalAdvFluxJacobianPoint(U(:,i),UPrim(:,i),fJac(:,:,i),gJac(:,:,i),hJac(:,:,i))
END DO !i
END SUBROUTINE EvalAdvFluxJacobian

!===================================================================================================================================
!> Navier-Stokes-Problem:
!> The Jacobian of the advective Flux with respect to the conservative variables U
!===================================================================================================================================
PPURE SUBROUTINE EvalAdvFluxJacobianPoint(U,UPrim,fJac,gJac,hJac)
! MODULES
USE MOD_EOS_Vars          ,ONLY:Kappa,KappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------  
REAL,DIMENSION(PP_nVar),INTENT(IN)              :: U               !< Conservative solution
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)          :: UPrim           !< Primitive solution
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: fJac,gJac,hJac  !< Derivative of the physical fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: KappaM2
REAL    :: uv,uu,vv,absu,v1,v2,srho
REAL    :: a1,phi
#if PP_dim==3
REAL    :: uw,vw,ww,v3
#endif
!===================================================================================================================================
KappaM2   = Kappa-2.

srho = 1./UPrim(1)
v1=UPrim(2)
v2=UPrim(3)
uv=UPrim(2)*UPrim(3)
uu=UPrim(2)*UPrim(2)
vv=UPrim(3)*UPrim(3)
#if PP_dim==3
v3=UPrim(4)
uw=UPrim(2)*UPrim(4)
vw=UPrim(3)*UPrim(4)
ww=UPrim(4)*UPrim(4)
absu=uu+vv+ww
phi  = kappaM1*0.5*absu
a1   = kappa * U(5)*sRho - phi


fJac(1,1:5)= (/          0.,             1.,          0.,           0.,        0. /)
fJac(2,1:5)= (/      phi-uu, (1-kappaM2)*v1, -kappaM1*v2,  -kappaM1*v3,   kappaM1 /)
fJac(3,1:5)= (/         -uv,             v2,          v1,           0.,        0. /)
fJac(4,1:5)= (/         -uw,             v3,          0.,           v1,        0. /)
fJac(5,1:5)= (/ v1*(phi-a1),  a1-kappaM1*uu, -kappaM1*uv,  -kappaM1*uw,  kappa*v1 /)


gJac(1,1:5)= (/          0.,           0.,              1.,          0.,       0. /)
gJac(2,1:5)= (/         -uv,           v2,              v1,          0.,       0. /)
gJac(3,1:5)= (/      phi-vv,  -kappaM1*v1, (1.-kappaM2)*v2, -kappaM1*v3,  kappaM1 /)
gJac(4,1:5)= (/         -vw,           0.,              v3,          v2,       0. /)
gJac(5,1:5)= (/ v2*(phi-a1),  -kappaM1*uv,   a1-kappaM1*vv, -kappaM1*vw, kappa*v2 /)


hJac(1,1:5)= (/          0.,          0.,           0.,              1.,       0. /)
hJac(2,1:5)= (/         -uw,          v3,           0.,              v1,       0. /)
hJac(3,1:5)= (/         -vw,          0.,           v3,              v2,       0. /)
hJac(4,1:5)= (/      phi-ww, -kappaM1*v1,  -kappaM1*v2, (1.-kappaM2)*v3,  kappaM1 /)
hJac(5,1:5)= (/ v3*(phi-a1), -kappaM1*uw,  -kappaM1*vw,   a1-kappaM1*ww, kappa*v3 /)
#else
absu=uu+vv
phi  = kappaM1*0.5*absu
a1   = kappa * U(5)*sRho - phi


fJac(1,1:5)= (/          0.,             1.,          0.,    0.,        0. /)
fJac(2,1:5)= (/      phi-uu, (1-kappaM2)*v1, -kappaM1*v2,    0.,   kappaM1 /)
fJac(3,1:5)= (/         -uv,             v2,          v1,    0.,        0. /)
fJac(4,1:5)= (/          0.,             0.,          0.,    0.,        0. /)
fJac(5,1:5)= (/ v1*(phi-a1),  a1-kappaM1*uu, -kappaM1*uv,    0.,  kappa*v1 /)


gJac(1,1:5)= (/          0.,          0.,              1.,   0.,        0. /)
gJac(2,1:5)= (/         -uv,          v2,              v1,   0.,        0. /)
gJac(3,1:5)= (/      phi-vv, -kappaM1*v1, (1.-kappaM2)*v2,   0.,   kappaM1 /)
gJac(4,1:5)= (/          0.,          0.,              0.,   0.,        0. /)
gJac(5,1:5)= (/ v2*(phi-a1), -kappaM1*uv,   a1-kappaM1*vv,   0.,  kappa*v2 /)

hJac(:,:)=0.
#endif
END SUBROUTINE EvalAdvFluxJacobianPoint

#if PARABOLIC
!===================================================================================================================================
!> The Jacobian of the diffusion flux with respect to the conservative variables U
!===================================================================================================================================
SUBROUTINE EvalDiffFluxJacobian(nDOF_loc,U,UPrim,gradUx,gradUy,gradUz,fJac,gJac,hJac &
#if EDDYVISCOSITY
                                ,muSGS &
#endif
                                )
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars     ,ONLY:s23,s43
USE MOD_Viscosity
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------  
!----------------------------------------------------------------------------------------------------------------------------------  
INTEGER,INTENT(IN)                                   :: nDOF_loc             !< number of degrees of freedom
REAL,DIMENSION(PP_nVar        ,nDOF_loc),INTENT(IN)  :: U                    !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim    ,nDOF_loc),INTENT(IN)  :: UPrim                !< solution in primitive variables
REAL,DIMENSION(PP_nVarPrim    ,nDOF_loc),INTENT(IN)  :: gradUx,gradUy,gradUz !< primitive gradients
REAL,DIMENSION(PP_nVar,PP_nVar,nDOF_loc),INTENT(OUT) :: fJac,gJac,hJac       !< Derivative of the physical diffusive fluxes
#if EDDYVISCOSITY
REAL,DIMENSION(1              ,nDOF_loc),INTENT(IN)  :: muSGS                !< eddy viscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: dir,i
REAL                :: muS
#if PP_dim==3
REAL                :: tau(3,3)
#else
REAL                :: tau(2,2)
#endif
!===================================================================================================================================
fJac = 0.
gJac = 0.
hJac = 0.

DO i=1,nDOF_loc
  muS = VISCOSITY_PRIM(UPrim(:,i))
#if EDDYVISCOSITY
  muS = muS    + muSGS(1,i)
#endif

#if PP_dim==3
  tau(1,1) = muS * ( s43 * gradUx(2,i) - s23 * gradUy(3,i) - s23 * gradUz(4,i)) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
  tau(2,2) = muS * (-s23 * gradUx(2,i) + s43 * gradUy(3,i) - s23 * gradUz(4,i)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
  tau(3,3) = muS * (-s23 * gradUx(2,i) - s23 * gradUy(3,i) + s43 * gradUz(4,i)) !-2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
  tau(1,2) = muS * (gradUy(2,i) + gradUx(3,i))               !mu*(u_y+v_x)
  tau(2,1) = tau(1,2)
  tau(1,3) = muS * (gradUz(2,i) + gradUx(4,i))               !mu*(u_z+w_x)
  tau(3,1) = tau(1,3)
  tau(2,3) = muS * (gradUz(3,i) + gradUy(4,i))               !mu*(y_z+w_y)
  tau(3,2) = tau(2,3)
#else
  tau(1,1) = muS * ( s43 * gradUx(2,i) - s23 * gradUy(3,i))  ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
  tau(2,2) = muS * (-s23 * gradUx(2,i) + s43 * gradUy(3,i))  !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
  tau(1,2) = muS * (gradUy(2,i) + gradUx(3,i))               !mu*(u_y+v_x)
  tau(2,1) = tau(1,2)
#endif

  ! dF^d(1:4)/dU(:) = 0, only the energy equation is directly depending on the conservative solution!
  ! dF^d(5)/dU(1) = (u*tau_(1,1) + v*tau_(1,2) + w*tau_(1,3))/rho
  ! dF^d(5)/dU(2) = - tau_(1,1)/rho
  ! dF^d(5)/dU(3) = - tau_(1,2)/rho
  ! dF^d(5)/dU(4) = - tau_(1,3)/rho
  ! dF^d(5)/dU(5) = 0.
  DO dir=1,PP_dim
    fJac(5,1,i) = fJac(5,1,i) + tau(1,dir)*UPrim(1+dir,i)
  END DO
  fJac(5,1,i) = fJac(5,1,i)/ U(1,i)
  fJac(5,2,i) = -tau(1,1)  / U(1,i)
  fJac(5,3,i) = -tau(1,2)  / U(1,i)
#if PP_dim==3
  fJac(5,4,i) = -tau(1,3)  / U(1,i)
#endif

  ! dG^d(1:4)/dU(:) = 0, only the energy equation is directly depending on the conservative solution!
  ! dG^d(5)/dU(1) = (u*tau_(2,1) + v*tau_(2,2) + w*tau_(2,3))/rho
  ! dG^d(5)/dU(2) = - tau_(2,1)/rho
  ! dG^d(5)/dU(3) = - tau_(2,2)/rho
  ! dG^d(5)/dU(4) = - tau_(2,3)/rho
  ! dG^d(5)/dU(5) = 0.
  DO dir=1,PP_dim
    gJac(5,1,i) = gJac(5,1,i) + tau(2,dir)*UPrim(1+dir,i)
  END DO
  gJac(5,1,i) = gJac(5,1,i)/ U(1,i)
  gJac(5,2,i) = -tau(2,1)  / U(1,i)
  gJac(5,3,i) = -tau(2,2)  / U(1,i)
#if PP_dim==3
  gJac(5,4,i) = -tau(2,3)  / U(1,i)
#endif

#if PP_dim==3
  ! dH^d(1:4)/dU(:) = 0, only the energy equation is directly depending on the conservative solution!
  ! dH^d(5)/dU(1) = (u*tau_(3,1) + v*tau_(3,2) + w*tau_(3,3))/rho
  ! dH^d(5)/dU(2) = - tau_(3,1)/rho
  ! dH^d(5)/dU(3) = - tau_(3,2)/rho
  ! dH^d(5)/dU(4) = - tau_(3,3)/rho
  ! dH^d(5)/dU(5) = 0.
  DO dir=1,PP_dim
    hJac(5,1,i) = hJac(5,1,i) + tau(3,dir)*UPrim(1+dir,i)
  END DO
  hJac(5,1,i) = hJac(5,1,i)/ U(1,i)
  hJac(5,2,i) = -tau(3,1)  / U(1,i)
  hJac(5,3,i) = -tau(3,2)  / U(1,i)
  hJac(5,4,i) = -tau(3,3)  / U(1,i)
#endif
END DO

END SUBROUTINE EvalDiffFluxJacobian

!===================================================================================================================================
!> Computes the volume derivative of the analytical diffusive flux with respect to the gradient of U: d(F^v)/dQ, Q=grad U
!===================================================================================================================================
SUBROUTINE EvalFluxGradJacobian(nDOF_loc,U,UPrim,fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz &
#if EDDYVISCOSITY
                               ,muSGS &
#endif
                               )
! MODULES
USE MOD_PreProc
USE MOD_Viscosity
USE MOD_Equation_Vars,ONLY:s43,s23
USE MOD_EOS_Vars,     ONLY:cp,Pr
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: PrSGS
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                       :: nDOF_loc !< number of degrees of freedom
REAL,DIMENSION(PP_nVar    ,nDOF_loc),INTENT(IN)          :: U        !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim,nDOF_loc),INTENT(IN)          :: UPrim    !< solution in primitive variables
#if EDDYVISCOSITY
REAL,DIMENSION(1          ,nDOF_loc),INTENT(IN)          :: muSGS    !< eddy viscosity
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVarPrim,nDOF_loc),INTENT(OUT) :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz !<
                                                            !> Gradient of the diffusive Cartesian fluxes (iVar,i,j,k) w.r.t. grad
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
REAL                :: muS,lambda
!===================================================================================================================================
DO i=1,nDOF_loc
  muS    = VISCOSITY_PRIM(UPrim(:,i))
  lambda = THERMAL_CONDUCTIVITY_H(muS)
  !Add turbulent sub grid scale viscosity to mu
#if EDDYVISCOSITY
  muS    = muS    + muSGS(1,i)
  lambda = lambda + muSGS*cp/PrSGS
#endif
#if PP_dim==3
  ! derivatives of diffusive flux in x-direction
  fJacQx(1,1:6,i) = 0.
  fJacQx(2,1:6,i) = (/ 0.,           -muS*s43,                 0.,                 0., 0.,      0./)
  fJacQx(3,1:6,i) = (/ 0.,                 0.,               -muS,                 0., 0.,      0./)
  fJacQx(4,1:6,i) = (/ 0.,                 0.,                 0.,               -muS, 0.,      0./)
  fJacQx(5,1:6,i) = (/ 0., -muS*s43*UPrim(2,i),   -muS*UPrim(3,i),    -muS*UPrim(4,i), 0., -lambda/)

  fJacQy(1,1:6,i) = 0.
  fJacQy(2,1:6,i) = (/ 0.,                 0.,            muS*s23,                 0., 0.,      0./)
  fJacQy(3,1:6,i) = (/ 0.,               -muS,                 0.,                 0., 0.,      0./)
  fJacQy(4,1:6,i) = 0.
  fJacQy(5,1:6,i) = (/ 0.,    -muS*UPrim(3,i), muS*s23*UPrim(2,i),                 0., 0.,      0./)

  fJacQz(1,1:6,i) = 0.
  fJacQz(2,1:6,i) = (/ 0.,                 0.,                 0.,            muS*s23, 0.,      0./)
  fJacQz(3,1:6,i) = 0.
  fJacQz(4,1:6,i) = (/ 0.,               -muS,                 0.,                 0., 0.,      0./)
  fJacQz(5,1:6,i) = (/ 0.,    -muS*UPrim(4,i),                 0., muS*s23*UPrim(2,i), 0.,      0./)


  ! derivatives of diffusive flux in y-direction
  gJacQx(1,1:6,i) = 0.
  gJacQx(2,1:6,i) = (/ 0.,                 0.,               -muS,                 0., 0.,      0./)
  gJacQx(3,1:6,i) = (/ 0.,            muS*s23,                 0.,                 0., 0.,      0./)
  gJacQx(4,1:6,i) = 0.
  gJacQx(5,1:6,i) = (/ 0., muS*s23*UPrim(3,i),    -muS*UPrim(2,i),                 0., 0.,      0./)

  gJacQy(1,1:6,i) = 0.
  gJacQy(2,1:6,i) = (/ 0.,               -muS,                 0.,                 0., 0.,      0./)
  gJacQy(3,1:6,i) = (/ 0.,                 0.,           -muS*s43,                 0., 0.,      0./)
  gJacQy(4,1:6,i) = (/ 0.,                 0.,                 0.,               -muS, 0.,      0./)
  gJacQy(5,1:6,i) = (/ 0.,    -muS*UPrim(2,i),-muS*s43*UPrim(3,i),    -muS*UPrim(4,i), 0., -lambda/)

  gJacQz(1,1:6,i) = 0.
  gJacQz(2,1:6,i) = 0.
  gJacQz(3,1:6,i) = (/ 0.,                 0.,                 0.,            muS*s23, 0.,      0./)
  gJacQz(4,1:6,i) = (/ 0.,                 0.,               -muS,                 0., 0.,      0./)
  gJacQz(5,1:6,i) = (/ 0.,                 0.,    -muS*UPrim(4,i), muS*s23*UPrim(3,i), 0.,      0./)

  ! derivatives of diffusive flux in z-direction
  hJacQx(1,1:6,i) = 0.
  hJacQx(2,1:6,i) = (/ 0.,                 0.,                 0.,               -muS, 0.,      0./)
  hJacQx(3,1:6,i) = 0.
  hJacQx(4,1:6,i) = (/ 0.,            muS*s23,                 0.,                 0., 0.,      0./)
  hJacQx(5,1:6,i) = (/ 0., muS*s23*UPrim(4,i),                 0.,    -muS*UPrim(2,i), 0.,      0./)

  hJacQy(1,1:6,i) = 0.
  hJacQy(2,1:6,i) = 0.
  hJacQy(3,1:6,i) = (/ 0.,                 0.,                 0.,               -muS, 0.,      0./)
  hJacQy(4,1:6,i) = (/ 0.,                 0.,            mu0*s23,                 0., 0.,      0./)
  hJacQy(5,1:6,i) = (/ 0.,                 0., muS*s23*UPrim(4,i),    -muS*UPrim(3,i), 0.,      0./)

  hJacQz(1,1:6,i) = 0.
  hJacQz(2,1:6,i) = (/ 0.,               -muS,                 0.,                 0., 0.,      0./)
  hJacQz(3,1:6,i) = (/ 0.,                 0.,               -muS,                 0., 0.,      0./)
  hJacQz(4,1:6,i) = (/ 0.,                 0.,                 0.,           -muS*s43, 0.,      0./)
  hJacQz(5,1:6,i) = (/ 0.,    -muS*UPrim(2,i),    -muS*UPrim(3,i),-muS*s43*UPrim(4,i), 0., -lambda/)
#else
  ! derivatives of diffusive flux in x-direction
  fJacQx(1,1:6,i) = 0.
  fJacQx(2,1:6,i) = (/ 0.,           -muS*s43,                 0.,                 0., 0.,      0./)
  fJacQx(3,1:6,i) = (/ 0.,                 0.,               -muS,                 0., 0.,      0./)
  fJacQx(4,1:6,i) = 0.
  fJacQx(5,1:6,i) = (/ 0., -muS*s43*UPrim(2,i),   -muS*UPrim(3,i),                 0., 0., -lambda/)

  fJacQy(1,1:6,i) = 0.
  fJacQy(2,1:6,i) = (/ 0.,                 0.,            muS*s23,                 0., 0.,      0./)
  fJacQy(3,1:6,i) = (/ 0.,               -muS,                 0.,                 0., 0.,      0./)
  fJacQy(4,1:6,i) = 0.
  fJacQy(5,1:6,i) = (/ 0.,    -muS*UPrim(3,i), muS*s23*UPrim(2,i),                 0., 0.,      0./)

  fJacQz(:,:,i) = 0.

  ! derivatives of diffusive flux in y-direction
  gJacQx(1,1:6,i) = 0.
  gJacQx(2,1:6,i) = (/ 0.,                 0.,               -muS,                 0., 0.,      0./)
  gJacQx(3,1:6,i) = (/ 0.,            muS*s23,                 0.,                 0., 0.,      0./)
  gJacQx(4,1:6,i) = 0.
  gJacQx(5,1:6,i) = (/ 0., muS*s23*UPrim(3,i),    -muS*UPrim(2,i),                 0., 0.,      0./)

  gJacQy(1,1:6,i) = 0.
  gJacQy(2,1:6,i) = (/ 0.,               -muS,                 0.,                 0., 0.,      0./)
  gJacQy(3,1:6,i) = (/ 0.,                 0.,           -muS*s43,                 0., 0.,      0./)
  gJacQy(4,1:6,i) = 0.
  gJacQy(5,1:6,i) = (/ 0.,    -muS*UPrim(2,i),-muS*s43*UPrim(3,i),                 0., 0., -lambda/)

  gJacQz(:,:,i) = 0.

  ! derivatives of diffusive flux in z-direction
  hJacQx(:,:,i) = 0.
  hJacQy(:,:,i) = 0.
  hJacQz(:,:,i) = 0.
#endif
END DO
END SUBROUTINE EvalFluxGradJacobian
#endif /*PARABOLIC*/

!===================================================================================================================================
!> The Jacobian of the transformation from primitive to conservative variables
!===================================================================================================================================
SUBROUTINE dConsdPrim(UPrim,Jac)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars                 ,ONLY: KappaM1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: Jac      !< cons to prim Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                            :: UE(PP_2Var),dpdrho,dedrho
!===================================================================================================================================
UE(PRIM) = UPrim
UE(DENS) = UPrim(1)
UE(SRHO) = 1./UE(DENS)


dpdrho = KappaM1*0.5*SUM(UE(VELV)*UE(VELV))
#if PP_dim == 3
dedrho = (UE(VEL1)**2+UE(VEL2)**2+UE(VEL3)**2) - dpdrho / KappaM1
#else
dedrho = (UE(VEL1)**2+UE(VEL2)**2            ) - dpdrho / KappaM1
#endif

Jac(1,1:5)= (/      1.,                0.,                0.,                0.,         0. /)
Jac(2,1:5)= (/UE(VEL1),          UE(DENS),                0.,                0.,         0. /)
Jac(3,1:5)= (/UE(VEL2),                0.,          UE(DENS),                0.,         0. /)
#if PP_dim == 3
Jac(4,1:5)= (/UE(VEL3),                0.,                0.,          UE(DENS),         0. /)
Jac(5,1:5)= (/  dedrho, UE(DENS)*UE(VEL1), UE(DENS)*UE(VEL2), UE(DENS)*UE(VEL3), 1./KappaM1 /)
#else
Jac(4,1:5)= 0.
Jac(5,1:5)= (/  dedrho, UE(DENS)*UE(VEL1), UE(DENS)*UE(VEL2),                0., 1./KappaM1 /)
#endif

END SUBROUTINE dConsdPrim

!===================================================================================================================================
!> The Jacobian of the transformation from primitive (including temperature) to conservative variables
!===================================================================================================================================
SUBROUTINE dConsdPrimTemp(UPrim,Jac)
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVarPrim),INTENT(OUT) :: Jac      !< cons to prim Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! fill jacobian without temperature
CALL dConsdPrim(UPrim,Jac(1:5,1:5))
Jac(1:5,6) = 0.
END SUBROUTINE dConsdPrimTemp

!===================================================================================================================================
!> The Jacobian of the transformation from conservative to primitive variables
!===================================================================================================================================
SUBROUTINE dPrimdCons(UPrim,Jac)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars                 ,ONLY:KappaM1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: Jac      !< prim to cons Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                            :: UE(PP_2Var)
!===================================================================================================================================
UE(PRIM) = UPrim
UE(SRHO) = 1./UPrim(1)

Jac(1,1:5)= (/                                1.,                0.,                0.,                0.,      0. /)
Jac(2,1:5)= (/                -UE(VEL1)*UE(SRHO),          UE(SRHO),                0.,                0.,      0. /)
Jac(3,1:5)= (/                -UE(VEL2)*UE(SRHO),                0.,          UE(SRHO),                0.,      0. /)
#if PP_dim == 3
Jac(4,1:5)= (/                -UE(VEL3)*UE(SRHO),                0.,                0.,          UE(SRHO),      0. /)
Jac(5,1:5)= (/KappaM1*0.5*SUM(UE(VELV)*UE(VELV)), -UE(VEL1)*KappaM1, -UE(VEL2)*KappaM1, -UE(VEL3)*KappaM1, KappaM1 /)
#else
Jac(4,1:5)= 0.
Jac(5,1:5)= (/KappaM1*0.5*SUM(UE(VELV)*UE(VELV)), -UE(VEL1)*KappaM1, -UE(VEL2)*KappaM1,                0., KappaM1 /)
#endif
END SUBROUTINE dPrimdCons

!===================================================================================================================================
!> The Jacobian of the transformation from conservative to primitive (including temperature) variables
!===================================================================================================================================
SUBROUTINE dPrimTempdCons(UPrim,Jac)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars                 ,ONLY:KappaM1,R
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim,PP_nVar),INTENT(OUT) :: Jac      !< prim to cons Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                            :: sRhoR,dpdU(5)
!===================================================================================================================================
! fill jacobian without temperature
CALL dPrimdCons(UPrim,Jac(1:5,1:5))

! fill jacobian of transformation to temperature
#if PP_dim==3
dpdU(1)   =  KappaM1*0.5*(UPrim(2)**2+UPrim(3)**2+UPrim(4)**2)
#else
dpdU(1)   =  KappaM1*0.5*(UPrim(2)**2+UPrim(3)**2)
#endif
dpdU(2:4) = -KappaM1*UPrim(2:4)
dpdU(5)   =  KappaM1
sRhoR     =  1./(R*UPrim(1))

Jac(6,1)   = dpdU(1  )*sRhoR-UPrim(5)*sRhoR/UPrim(1)
Jac(6,2:4) = dpdU(2:4)*sRhoR
Jac(6,5)   = dpdU(5  )*sRhoR

END SUBROUTINE dPrimTempdCons

END MODULE MOD_Jacobian
