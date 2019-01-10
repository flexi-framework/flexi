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

!==================================================================================================================================
!> Contains the routine EvalFlux3D which computes the complete NSE flux f,g,h for all DOFs in one Element: used in volume integral
!> Contains the routine EvalFlux1D_Adv which computes the Euler flux f for all DOFs of one Element side: used in Riemann_Adv
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE EvalFlux3D
  MODULE PROCEDURE EvalFlux3D_Point
  MODULE PROCEDURE EvalFlux3D_Volume
END INTERFACE

INTERFACE EvalEulerFlux1D
  MODULE PROCEDURE EvalEulerFlux1D
END INTERFACE

INTERFACE EvalEulerFlux1D_fast
  MODULE PROCEDURE EvalEulerFlux1D_fast
END INTERFACE

#if PARABOLIC
INTERFACE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D_Point
  MODULE PROCEDURE EvalDiffFlux3D_Surface
  MODULE PROCEDURE EvalDiffFlux3D_Volume
END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::EvalFlux3D, EvalEulerFlux1D, EvalEulerFlux1D_fast
#if PARABOLIC
PUBLIC::EvalDiffFlux3D
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute Navier-Stokes fluxes using the conservative variables and derivatives
!==================================================================================================================================
PPURE SUBROUTINE EvalFlux3D_Point(U,UPrim,f,g,h)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U        !< Conservative solution
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim    !< Primitive solution
!> Physical fluxes in x/y/z direction
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: f,g,h
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Ep
!==================================================================================================================================
! auxiliary variables
Ep   = U(5) + UPrim(5)
#if PP_dim==3
! Euler part
! Euler fluxes x-direction
f(1) = U(2)                       ! rho*u
f(2) = U(2) * UPrim(2) + UPrim(5) ! rho*u²+p
f(3) = U(2) * UPrim(3)            ! rho*u*v
f(4) = U(2) * UPrim(4)            ! rho*u*w
f(5) = Ep * UPrim(2)              ! (rho*e+p)*u
! Euler fluxes y-direction
g(1) = U(3)                       ! rho*v
g(2) = f(3)                       ! rho*u*v
g(3) = U(3) * UPrim(3) + UPrim(5) ! rho*v²+p
g(4) = U(3) * UPrim(4)            ! rho*v*w
g(5) = Ep * UPrim(3)              ! (rho*e+p)*v
! Euler fluxes z-direction
h(1) = U(4)                       ! rho*v
h(2) = f(4)                       ! rho*u*w
h(3) = g(4)                       ! rho*v*w
h(4) = U(4) * UPrim(4) + UPrim(5) ! rho*v²+p
h(5) = Ep * UPrim(4)              ! (rho*e+p)*w

#else

! Euler part
! Euler fluxes x-direction
f(1) = U(2)                       ! rho*u
f(2) = U(2)*UPrim(2)+UPrim(5)     ! rho*u²+p
f(3) = U(2)*UPrim(3)              ! rho*u*v
f(4) = 0.
f(5) = Ep*UPrim(2)                ! (rho*e+p)*u
! Euler fluxes y-direction
g(1)= U(3)                        ! rho*v
g(2)= f(3)                        ! rho*u*v
g(3)= U(3)*UPrim(3)+UPrim(5)      ! rho*v²+p
g(4)= 0.
g(5)= Ep*UPrim(3)                 ! (rho*e+p)*v
! Euler fluxes z-direction
h   = 0.
#endif
END SUBROUTINE EvalFlux3D_Point

PPURE SUBROUTINE EvalFlux3D_Volume(Nloc,U,UPrim,f,g,h)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER                                               ,INTENT(IN)  :: Nloc     !< Polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U        !< Conservative solution
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim    !< Primitive solution
!> Physical fluxes in x,y,z directions
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f,g,h
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,ZDIM(Nloc);  DO j=0,Nloc; DO i=0,Nloc
  CALL EvalFlux3D_Point(U(:,i,j,k),UPrim(:,i,j,k),f(:,i,j,k),g(:,i,j,k),h(:,i,j,k))
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalFlux3D_Volume


#if PARABOLIC
!==================================================================================================================================
!> Compute Navier-Stokes diffusive flux using the conservative variables and derivatives.
!==================================================================================================================================
PPURE SUBROUTINE EvalDiffFlux3D_Point(UPrim,gradUx,gradUy,gradUz,f,g,h &
#if EDDYVISCOSITY
                                      ,muSGS &
#endif
)
! MODULES
USE MOD_Equation_Vars,ONLY: s23,s43
USE MOD_EOS_Vars,     ONLY: cp,Pr
USE MOD_Viscosity
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: PrSGS
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim                !< Solution vector
!> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: gradUx,gradUy,gradUz
!> Physical fluxes in x,y,z directions
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: f,g,h
#if EDDYVISCOSITY
REAL                       ,INTENT(IN)  :: muSGS                !< SGS viscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: muS,lambda
REAL                :: tau_xx,tau_yy,tau_xy
#if PP_dim==3
REAL                :: tau_zz,tau_xz,tau_yz
#endif
!==================================================================================================================================
! Viscous part
! ideal gas law
muS    = VISCOSITY_PRIM(UPrim)
lambda = THERMAL_CONDUCTIVITY_H(muS)
!Add turbulent sub grid scale viscosity to mu
#if EDDYVISCOSITY
muS    = muS    + muSGS
lambda = lambda + muSGS*cp/PrSGS
#endif

ASSOCIATE( v1     => UPrim(2),  v2     => UPrim(3),  v3     => UPrim(4), &
           gradT1 => GradUx(6), gradT2 => GradUy(6), gradT3 => GradUz(6) )
#if PP_dim==3
! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dT /)

! viscous fluxes in x-direction
tau_xx = muS * ( s43 * gradUx(2) - s23 * gradUy(3) - s23 * gradUz(4)) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = muS * (-s23 * gradUx(2) + s43 * gradUy(3) - s23 * gradUz(4)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_zz = muS * (-s23 * gradUx(2) - s23 * gradUy(3) + s43 * gradUz(4)) !-2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
tau_xy = muS * (gradUy(2) + gradUx(3))               !mu*(u_y+v_x)
tau_xz = muS * (gradUz(2) + gradUx(4))               !mu*(u_z+w_x)
tau_yz = muS * (gradUz(3) + gradUy(4))               !mu*(y_z+w_y)

f(1) = 0.
f(2) = -tau_xx                                       ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(3) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
f(4) = -tau_xz                                       ! F_euler-mu*(u_z+w_x)
f(5) = -tau_xx*v1-tau_xy*v2-tau_xz*v3-lambda*gradT1  ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
! viscous fluxes in y-direction
g(1) = 0.
g(2) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
g(3) = -tau_yy                                       ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(4) = -tau_yz                                       ! F_euler-mu*(y_z+w_y)
g(5) = -tau_xy*v1-tau_yy*v2-tau_yz*v3-lambda*gradT2  ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
! viscous fluxes in z-direction
h(1) = 0.
h(2) = -tau_xz                                       ! F_euler-mu*(u_z+w_x)
h(3) = -tau_yz                                       ! F_euler-mu*(y_z+w_y)
h(4) = -tau_zz                                       ! F_euler-4/3*mu*w_z+2/3*mu*(u_x+v_y)
h(5) = -tau_xz*v1-tau_yz*v2-tau_zz*v3-lambda*gradT3  ! F_euler-(tau_zx*u+tau_zy*v+tau_zz*w-q_z) q_z=-lambda*T_z
#else
! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dT /)

! viscous fluxes in x-direction
tau_xx = muS * ( s43 * gradUx(2) - s23 * gradUy(3))  ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = muS * (-s23 * gradUx(2) + s43 * gradUy(3))  !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_xy = muS * (gradUy(2) + gradUx(3))               !mu*(u_y+v_x)

f(1) = 0.
f(2) = -tau_xx                                       ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(3) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
f(4) = 0.
f(5) = -tau_xx*v1-tau_xy*v2-lambda*gradT1            ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
! viscous fluxes in y-direction
g(1) = 0.
g(2) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
g(3) = -tau_yy                                       ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(4) = 0.
g(5) = -tau_xy*v1-tau_yy*v2-lambda*gradT2  ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
! viscous fluxes in z-direction
h    = 0.
#endif
END ASSOCIATE
END SUBROUTINE EvalDiffFlux3D_Point

!==================================================================================================================================
!> Compute Navier-Stokes diffusive flux for the volume
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Volume(UPrim,gradUx,gradUy,gradUz,f,g,h,iElem)
! MODULES
USE MOD_PreProc
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: UPrim                !< Solution vector
!> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: gradUx,gradUy,gradUz
!> Physical fluxes in x,y,z directions
REAL,DIMENSION(PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT) :: f,g,h
INTEGER, INTENT(IN)                                           :: iElem                !< element index in global array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,PP_NZ;  DO j=0,PP_N; DO i=0,PP_N
  CALL EvalDiffFlux3D_Point(Uprim(:,i,j,k),gradUx(:,i,j,k),gradUy(:,i,j,k),gradUz(:,i,j,k), &
                                                 f(:,i,j,k),     g(:,i,j,k),     h(:,i,j,k)  &
#if EDDYVISCOSITY
                            ,muSGS(1,i,j,k,iElem)&
#endif
                            )
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalDiffFlux3D_Volume

!==================================================================================================================================
!> Compute Navier-Stokes diffusive flux for the surface
!==================================================================================================================================
PPURE SUBROUTINE EvalDiffFlux3D_Surface(Nloc,UPrim,gradUx,gradUy,gradUz,f,g,h &
#if EDDYVISCOSITY
                                 ,muSGS &
#endif
                                 )
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER                                        ,INTENT(IN)  :: Nloc                 !< Polynomial degree of input solution
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim                !< Solution vector
!> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx,gradUy,gradUz
!> Physical fluxes in x,y,z directions
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f,g,h
#if EDDYVISCOSITY
REAL,DIMENSION(1          ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: muSGS                !< SGS viscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j
!==================================================================================================================================
DO j=0,ZDIM(Nloc); DO i=0,PP_N
  CALL EvalDiffFlux3D_Point(Uprim(:,i,j),gradUx(:,i,j),gradUy(:,i,j),gradUz(:,i,j), &
                                               f(:,i,j),     g(:,i,j),     h(:,i,j)  &
#if EDDYVISCOSITY
                            ,muSGS(1,i,j) &
#endif
                            )
END DO; END DO ! i,j
END SUBROUTINE EvalDiffFlux3D_Surface
#endif /*PARABOLIC*/


!==================================================================================================================================
!> Computes 1D Euler flux using the conservative variables.
!==================================================================================================================================
PPURE SUBROUTINE EvalEulerFlux1D(U,F)
! MODULES
USE MOD_EOS_Vars ,ONLY:KappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_nVar) !< vector of conservative variables
REAL,INTENT(OUT)    :: F(PP_nVar) !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: UE(PP_2Var)      ! auxiliary variables
!==================================================================================================================================
! auxiliary variables
! TODO: ATTENTION: Temperature of UE not filled!!!
UE(CONS)=U
UE(SRHO)=1./UE(DENS)
UE(VELV)=VELOCITY_HE(UE)
UE(PRES)=PRESSURE_HE(UE)
! Euler fluxes x-direction
F(1)= U(MOM1)                     ! rho*u
F(2)= U(MOM1)*UE(VEL1)+UE(PRES)   ! rho*u²+p
F(3)= U(MOM1)*UE(VEL2)            ! rho*u*v
#if PP_dim==3
F(4)= U(MOM1)*UE(VEL3)            ! rho*u*w
#else
F(4)=0.
#endif
F(5)=(U(ENER)+UE(PRES))*UE(VEL1)  ! (rho*e+p)*u
END SUBROUTINE EvalEulerFlux1D


!==================================================================================================================================
!> Computes 1D Euler flux using the conservative and primitive variables (for better performance)
!==================================================================================================================================
PPURE SUBROUTINE EvalEulerFlux1D_fast(U,F)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_2Var) !< vector of conservative and primitive variables
REAL,INTENT(OUT)    :: F(PP_nVar) !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Euler fluxes x-direction
F(1)= U(MOM1)                 ! rho*u
F(2)= U(MOM1)*U(VEL1)+U(PRES) ! rho*u²+p
F(3)= U(MOM1)*U(VEL2)         ! rho*u*v
#if PP_dim==3
F(4)= U(MOM1)*U(VEL3)         ! rho*u*w
#else
F(4)= 0.
#endif
F(5)=(U(ENER)+U(PRES))*U(VEL1)! (rho*e+p)*u
END SUBROUTINE EvalEulerFlux1D_fast

END MODULE MOD_Flux
