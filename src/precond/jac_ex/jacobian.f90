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
#if EQNSYSNR==2
#include "eos.h"
#endif

!===================================================================================================================================
!> Contains the routines for the calculation of the analytical flux jacobians of the different equation systems
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
PUBLIC::EvalAdvFluxJacobian
#if EQNSYSNR==2
#if PARABOLIC
PUBLIC::EvalDiffFluxJacobian
#endif
#endif
!===================================================================================================================================

CONTAINS

#if EQNSYSNR==1
!===================================================================================================================================
!> Linear Scalar Advection Difussion:
!> The Jacobian of the analytical advective Flux with respect to the Variable U
!===================================================================================================================================
SUBROUTINE EvalAdvFluxJacobian(U,UPrim,fJac,gJac,hJac)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY:AdvVel
USE MOD_DG_Vars           ,ONLY:nDOFElem,imex
USE MOD_TimeDisc_Vars     ,ONLY:TimeDiscMode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar,nDOFElem),INTENT(IN)              :: U
REAL,DIMENSION(PP_nVarPrim,nDOFElem),INTENT(IN)          :: UPrim
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar,nDOFElem),INTENT(OUT) :: fJac,gJac,hJac             ! Cartesian fluxes (iVar,i,j,k)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
fJac = AdvVel(1)
gJac = AdvVel(2)
#if PP_dim==3
hJac = AdvVel(3)
#endif
END SUBROUTINE EvalAdvFluxJacobian
#endif

#if EQNSYSNR==2
SUBROUTINE EvalAdvFluxJacobian(U,UPrim,fJac,gJac,hJac)
!===================================================================================================================================
! Navier-Stokes-Problem:
! The Jacobian of the advective Flux with respect to the primitive Variable U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars          ,ONLY:Kappa,KappaM1
USE MOD_DG_Vars           ,ONLY:nDOFElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------  
REAL,DIMENSION(PP_nVar,nDOFElem),INTENT(IN)              :: U
REAL,DIMENSION(PP_nVarPrim,nDOFElem),INTENT(IN)          :: UPrim
REAL,DIMENSION(PP_nVar,PP_nVar,nDOFElem),INTENT(OUT)     :: fJac,gJac,hJac  ! Derivative of theCartesian fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: KappaM2
REAL    :: uv,uu,vv,absu,v1,v2,srho
INTEGER :: i
REAL    :: a1,phi
#if PP_dim==3
REAL    :: uw,vw,ww,v3
#endif
!===================================================================================================================================
KappaM2   = Kappa-2.

DO i=1,nDOFElem
  srho = 1./UPrim(1,i)
  v1=UPrim(2,i)
  v2=UPrim(3,i)
  uv=UPrim(2,i)*UPrim(3,i)
  uu=UPrim(2,i)*UPrim(2,i)
  vv=UPrim(3,i)*UPrim(3,i)
#if PP_dim==3
  v3=UPrim(4,i)
  uw=UPrim(2,i)*UPrim(4,i)
  vw=UPrim(3,i)*UPrim(4,i)
  ww=UPrim(4,i)*UPrim(4,i)
  absu=uu+vv+ww
  phi  = kappaM1*0.5*absu
  a1   = kappa * U(5,i)*sRho - phi


  fJac(1,1:5,i)= (/          0.,             1.,          0.,           0.,        0. /)
  fJac(2,1:5,i)= (/      phi-uu, (1-kappaM2)*v1, -kappaM1*v2,  -kappaM1*v3,   kappaM1 /)
  fJac(3,1:5,i)= (/         -uv,             v2,          v1,           0.,        0. /)
  fJac(4,1:5,i)= (/         -uw,             v3,          0.,           v1,        0. /)
  fJac(5,1:5,i)= (/ v1*(phi-a1),  a1-kappaM1*uu, -kappaM1*uv,  -kappaM1*uw,  kappa*v1 /)


  gJac(1,1:5,i)= (/          0.,           0.,              1.,          0.,       0. /)
  gJac(2,1:5,i)= (/         -uv,           v2,              v1,          0.,       0. /)
  gJac(3,1:5,i)= (/      phi-vv,  -kappaM1*v1, (1.-kappaM2)*v2, -kappaM1*v3,  kappaM1 /)
  gJac(4,1:5,i)= (/         -vw,           0.,              v3,          v2,       0. /)
  gJac(5,1:5,i)= (/ v2*(phi-a1),  -kappaM1*uv,   a1-kappaM1*vv, -kappaM1*vw, kappa*v2 /)


  hJac(1,1:5,i)= (/          0.,          0.,           0.,              1.,       0. /)
  hJac(2,1:5,i)= (/         -uw,          v3,           0.,              v1,       0. /)
  hJac(3,1:5,i)= (/         -vw,          0.,           v3,              v2,       0. /)
  hJac(4,1:5,i)= (/      phi-ww, -kappaM1*v1,  -kappaM1*v2, (1.-kappaM2)*v3,  kappaM1 /)
  hJac(5,1:5,i)= (/ v3*(phi-a1), -kappaM1*uw,  -kappaM1*vw,   a1-kappaM1*ww, kappa*v3 /)
#else
  absu=uu+vv
  phi  = kappaM1*0.5*absu
  a1   = kappa * U(5,i)*sRho - phi


  fJac(1,1:5,i)= (/          0.,             1.,          0.,    0.,        0. /)
  fJac(2,1:5,i)= (/      phi-uu, (1-kappaM2)*v1, -kappaM1*v2,    0.,   kappaM1 /)
  fJac(3,1:5,i)= (/         -uv,             v2,          v1,    0.,        0. /)
  fJac(4,1:5,i)= (/          0.,             0.,          0.,    0.,        0. /)
  fJac(5,1:5,i)= (/ v1*(phi-a1),  a1-kappaM1*uu, -kappaM1*uv,    0.,  kappa*v1 /)


  gJac(1,1:5,i)= (/          0.,          0.,              1.,   0.,        0. /)
  gJac(2,1:5,i)= (/         -uv,          v2,              v1,   0.,        0. /)
  gJac(3,1:5,i)= (/      phi-vv, -kappaM1*v1, (1.-kappaM2)*v2,   0.,   kappaM1 /)
  gJac(4,1:5,i)= (/          0.,          0.,              0.,   0.,        0. /)
  gJac(5,1:5,i)= (/ v2*(phi-a1), -kappaM1*uv,   a1-kappaM1*vv,   0.,  kappa*v2 /)

  hJac(:,:,i)=0.
#endif
END DO !i
END SUBROUTINE EvalAdvFluxJacobian

#if PARABOLIC
SUBROUTINE EvalDiffFluxJacobian(Nloc,U,UPrim,gradUx,gradUy,gradUz,fJac,gJac,hJac)
!===================================================================================================================================
! The Jacobian of the diffusion Flux with respect to the primitive Variable U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars          ,ONLY:sKappaM1,KappasPr,kappaM1,R
USE MOD_Equation_Vars     ,ONLY:s23,s43
USE MOD_EOS_Vars          ,ONLY:cp,Pr
USE MOD_Viscosity
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------  
!----------------------------------------------------------------------------------------------------------------------------------  
INTEGER,INTENT(IN)                                       :: Nloc
REAL,DIMENSION(PP_nVar,Nloc),INTENT(IN)              :: U
REAL,DIMENSION(PP_nVarPrim,Nloc),INTENT(IN)          :: UPrim
REAL,DIMENSION(PP_nVarPrim,Nloc),INTENT(IN)          :: gradUx,gradUy,gradUz
REAL,DIMENSION(PP_nVar,PP_nVar,Nloc),INTENT(OUT)     :: fJac,gJac,hJac  ! Derivative of theCartesian fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: v1v1,v2v2,v1v2
REAL                :: srho,e,v1,v2,srhomu
REAL                :: rhox,rhoy
REAL                :: gradv1x,gradv2x
REAL                :: gradv1y,gradv2y
REAL                :: gradex,gradey
REAL                :: muS
REAL                :: lambda
REAL                :: s43mKappasPr                            ! (4/3-kappa/Pr)
REAL                :: mKappasPr                               ! (1-kappa/Pr)
REAL                :: KappaM1sR
INTEGER             :: i
!REAL                :: tau_xx,tau_yy,tau_xy
#if PP_dim==3
REAL                :: gradv3x,gradv3y
REAL                :: v3,v3v3,v2v3,v1v3
REAL                :: rhoz
REAL                :: gradv1z,gradv2z,gradv3z
REAL                :: gradez
#endif
!===================================================================================================================================
s43mKappasPr = s43-KappasPr
mKappasPr    = 1.-KappasPr
KappaM1sR=KappaM1/R

fJac=0.
gJac=0.
hJac=0.

DO i=1,Nloc
  srho = 1. / U(1,i) ! 1/rho
  v1   = UPrim(2,i) 
  v2   = UPrim(3,i) 
  v1v1 = v1*v1
  v2v2 = v2*v2
  v1v2 = v1*v2

  e=U(5,i)*srho                              ! e...specific energy
  !ideal gas law
  muS=VISCOSITY_PRIM(UPrim(:,i))
  lambda=THERMAL_CONDUCTIVITY_H(muS)

  srhomu  = srho*muS !add viscosity
  rhox    = srho*srhomu*gradUx(1,i)    !ATTENTION: rhox=mu0*drhodx/rho^2 !!!! 
  rhoy    = srho*srhomu*gradUy(1,i) 

  gradv1x = srhomu*gradUx(2,i)    !gradu(1,1)= u_x, gradu(1,2)=u_y, gradu(1,3)=u_z
  gradv2x = srhomu*gradUx(3,i)    !gradu(2,1)= v_x, gradu(2,2)=v_y, gradu(2,3)=v_z
  gradex  = srhomu*(sKappaM1*srho*(gradUx(5,i)-UPrim(5,i)*gradUx(1,i)*srho) + SUM(UPrim(2:PP_dim+1,i)*gradUx(2:PP_dim+1,i)))
  gradv1y = srhomu*gradUy(2,i)
  gradv2y = srhomu*gradUy(3,i)
  gradey  = srhomu*(sKappaM1*srho*(gradUy(5,i)-UPrim(5,i)*gradUy(1,i)*srho) + SUM(UPrim(2:PP_dim+1,i)*gradUy(2:PP_dim+1,i)))

#if PP_dim==3
  v3   = UPrim(4,i) 
  v3v3 = v3*v3
  v1v3 = v1*v3
  v2v3 = v2*v3

  rhoz   = srho*srhomu*gradUz(1,i) 

  ! compute derivatives via product rule (a*b)'=a'*b+a*b'
  gradv3x = srhomu*gradUx(4,i)    !gradu(3,1)= w_x, gradu(3,2)=w_y, gradu(3,3)=w_z
  gradv3y = srhomu*gradUy(4,i)
  gradv1z = srhomu*gradUz(2,i)
  gradv2z = srhomu*gradUz(3,i)
  gradv3z = srhomu*gradUz(4,i)
  gradez  = srhomu*(sKappaM1*srho*(gradUz(5,i)-UPrim(5,i)*gradUz(1,i)*srho) + SUM(UPrim(2:4,i)*gradUz(2:4,i)))
  ! viscous fluxes in x-direction      
  ! fJac(i,j)=d f(U_i) / d U_j
  ! Attention: Use product rule, since gradv1x= srho*(gradUx2 - (rho*v1)/rho *gradUx1)
  fJac(2,1  ,i) = ( s43*(gradv1x-v1*rhox)-s23*(gradv2y-v2*rhoy + gradv3z-v3*rhoz) )
  fJac(3,1  ,i) = ( (gradv1y-v1*rhoy) +(gradv2x-v2*rhox) ) 
  fJac(4,1  ,i) = ( (gradv1z-v1*rhoz) +(gradv3x-v3*rhox) ) 
  fJac(5,1  ,i) = ( v1*(s43mKappasPr*(2*gradv1x-v1*rhox)-s23*(2*(gradv2y+gradv3z)-v2*rhoy -v3*rhoz) ) &
                              +mKappasPr*(2*(v2*gradv2x+v3*gradv3x)-(v2v2+v3v3)*rhox) -v1v2*rhoy-v1v3*rhoz &
                              +2*(v2*gradv1y +v3*gradv1z)  &
                              + KappasPr*(gradex-e*rhox) ) 
  fJac(2,2:4,i) = (/ s43*rhox, -s23*rhoy, -s23*rhoz /) 
  fJac(3,2:4,i) = (/     rhoy,      rhox,        0. /) 
  fJac(4,2:4,i) = (/     rhoz,        0.,      rhox /) 
  fJac(5,2  ,i) = ( -s43mKappasPr*(gradv1x-v1*rhox) +s23*(gradv2y+gradv3z) +v2*rhoy +v3*rhoz )
  fJac(5,3  ,i) = ( -s23*v1*rhoy -gradv1y -mKappasPr*(gradv2x-v2*rhox) )
  fJac(5,4  ,i) = ( -s23*v1*rhoz -gradv1z -mKappasPr*(gradv3x-v3*rhox) )
  fJac(5,5  ,i) = ( KappasPr*rhox  ) 

  ! viscous fluxes in y-direction      
  ! gJac(i,j)=d g(U_i) / d U_j
  gJac(2,1  ,i) = ( (gradv1y-v1*rhoy) +(gradv2x-v2*rhox) ) 
  gJac(3,1  ,i) = ( s43*(gradv2y-v2*rhoy)-s23*(gradv1x-v1*rhox + gradv3z-v3*rhoz) ) 
  gJac(4,1  ,i) = ( (gradv2z-v2*rhoz) +(gradv3y-v3*rhoy) ) 
  gJac(5,1  ,i) = ( v2*(s43mKappasPr*(2*gradv2y-v2*rhoy)-s23*(2*(gradv1x+gradv3z)-v1*rhox -v3*rhoz) ) &
                              +mKappasPr*(2*(v1*gradv1y+v3*gradv3y)-(v1v1+v3v3)*rhoy) -v1v2*rhox -v2v3*rhoz &
                              +2*(v1*gradv2x+v3*gradv2z)  &
                              + KappasPr*(gradey-e*rhoy) ) 
  gJac(2,2:4,i) = (/     rhoy,      rhox,        0. /) 
  gJac(3,2:4,i) = (/-s23*rhox,  s43*rhoy, -s23*rhoz /) 
  gJac(4,2:4,i) = (/       0.,      rhoz,      rhoy /) 
  gJac(5,2  ,i) = ( -s23*v2*rhox -gradv2x -mKappasPr*(gradv1y-v1*rhoy) )
  gJac(5,3  ,i) = ( -s43mKappasPr*(gradv2y-v2*rhoy) +s23*(gradv1x+gradv3z) +v1*rhox +v3*rhoz )
  gJac(5,4  ,i) = ( -s23*v2*rhoz -gradv2z -mKappasPr*(gradv3y-v3*rhoy) )
  gJac(5,5  ,i) = ( KappasPr*rhoy  ) 

  ! viscous fluxes in z-direction      
  ! hJac(i,j)=d h(U_i) / d U_j
  hJac(2,1  ,i) = ( (gradv1z-v1*rhoz) +(gradv3x-v3*rhox) )
  hJac(3,1  ,i) = ( (gradv2z-v2*rhoz) +(gradv3y-v3*rhoy) ) 
  hJac(4,1  ,i) = ( s43*(gradv3z-v3*rhoz)-s23*(gradv2y-v2*rhoy + gradv1x-v1*rhox) ) 
  hJac(5,1  ,i) = ( v3*(s43mKappasPr*(2*gradv3z-v3*rhoz)-s23*(2*(gradv2y+gradv1x)-v2*rhoy -v1*rhox) ) &
                              +mKappasPr*(2*(v1*gradv1z+v2*gradv2z)-(v1v1+v2v2)*rhoz) -v1v3*rhox - v2v3*rhoy   &
                              +2*(v1*gradv3x +v2*gradv3y) &
                              + KappasPr*(gradez-e*rhoz) ) 
  hJac(2,2:4,i) = (/     rhoz,        0.,      rhox /) 
  hJac(3,2:4,i) = (/       0.,      rhoz,      rhoy /) 
  hJac(4,2:4,i) = (/-s23*rhox, -s23*rhoy,  s43*rhoz /) 
  hJac(5,2  ,i) = ( -s23*v3*rhox -gradv3x -mKappasPr*(gradv1z-v1*rhoz) )
  hJac(5,3  ,i) = ( -s23*v3*rhoy -gradv3y -mKappasPr*(gradv2z-v2*rhoz) )
  hJac(5,4  ,i) = ( -s43mKappasPr*(gradv3z-v3*rhoz) +s23*(gradv2y+gradv1x) +v1*rhox +v2*rhoy )
  hJac(5,5  ,i) = ( KappasPr*rhoz  ) 
#else
  ! viscous fluxes in x-direction      
  ! fJac(i,j)=d f(U_i) / d U_j

  !tau_xx=muS*( s43*gradUx(2,i)-s23*gradUy(3,i))
  !tau_yy=muS*(-s23*gradUx(2,i)+s43*gradUy(3,i)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
  !tau_xy=muS*(gradUy(2,i)+gradUx(3,i))

  fJac(2,1  ,i) = ( s43*(gradv1x-v1*rhox)-s23*(gradv2y-v2*rhoy) )
  fJac(3,1  ,i) = ( (gradv1y-v1*rhoy) +(gradv2x-v2*rhox) ) 
  fJac(4,:  ,i) = 0. 
  fJac(5,1  ,i) = ( v1*(s43mKappasPr*(2*gradv1x-v1*rhox)-s23*(2*gradv2y-v2*rhoy) ) &
                              +mKappasPr*(2*(v2*gradv2x)-(v2v2)*rhox) -v1v2*rhoy &
                              +2*(v2*gradv1y)  &
                              + KappasPr*(gradex-e*rhox) ) 
  ! not working
  !fJac(5,1  ,i) = v1*fJac(2,1,i)+tau_xx*v1*srho + v2*fJac(3,1,i)+tau_xy*v2*srho &
                  !-KappasPr*(-gradex+e*rhox+v1*(gradv1x-v1*rhox)+v2*(gradv2x-v2*rhox))

  fJac(2,2:4,i) = (/ s43*rhox, -s23*rhoy, 0. /) 
  fJac(3,2:4,i) = (/     rhoy,      rhox,        0. /) 
  fJac(5,2  ,i) = ( -s43mKappasPr*(gradv1x-v1*rhox) +s23*(gradv2y) +v2*rhoy )
  !the same as below
  !fJac(5,2  ,i) = fJac(2,2,i)*v1 - tau_xx*srho + fJac(3,2,i)*v2 &
                  !-KappasPr*(-gradv1x+v1*rhox)
  fJac(5,3  ,i) = ( -s23*v1*rhoy -gradv1y -mKappasPr*(gradv2x-v2*rhox) )
  !the same as below
  !fJac(5,3  ,i) = fJac(2,3,i)*v1 +fJac(3,3,i)*v2 - tau_xy*srho &
                  !-KappasPr*(-gradv2x+v2*rhox)
  fJac(5,4  ,i) = 0.
  fJac(5,5  ,i) = ( KappasPr*rhox  ) 

  ! viscous fluxes in y-direction      
  ! gJac(i,j)=d g(U_i) / d U_j
  gJac(2,1  ,i) = fJac(3,1,i) 
  gJac(3,1  ,i) = ( s43*(gradv2y-v2*rhoy)-s23*(gradv1x-v1*rhox ) ) 
  gJac(4,:  ,i) = 0. 
  gJac(5,1  ,i) = ( v2*(s43mKappasPr*(2*gradv2y-v2*rhoy)-s23*(2*(gradv1x)-v1*rhox) ) &
                              +mKappasPr*(2*(v1*gradv1y)-(v1v1)*rhoy) -v1v2*rhox  &
                              +2*(v1*gradv2x)  &
                              + KappasPr*(gradey-e*rhoy) ) 
  ! not working
  !gJac(5,1  ,i) = v1*gJac(2,1,i)+tau_xy*v1*srho + v2*gJac(3,1,i)+tau_yy*v2*srho &
                  !-KappasPr*(-gradey+e*rhoy+v1*(gradv1y-v1*rhoy)+v2*(gradv2y-v2*rhoy))
  gJac(2,2:4,i) = (/     rhoy,      rhox,        0. /) 
  gJac(3,2:4,i) = (/-s23*rhox,  s43*rhoy,        0. /) 
  gJac(5,2  ,i) = ( -s23*v2*rhox -gradv2x -mKappasPr*(gradv1y-v1*rhoy) )
  !the same as below
  !gJac(5,2  ,i) =  gJac(2,2,i)*v1 - tau_xy*srho + gJac(3,2,i)*v2 &
                   !-KappasPr*(-gradv1y+v1*rhoy)
  gJac(5,3  ,i) = ( -s43mKappasPr*(gradv2y-v2*rhoy) +s23*(gradv1x) +v1*rhox )
  !the same as below
  !gJac(5,3  ,i) =  gJac(2,3,i)*v1 + gJac(3,3,i)*v2 - tau_yy*srho &
                   !-KappasPr*(-gradv2y+v2*rhoy)
  gJac(5,4  ,i) = 0.
  gJac(5,5  ,i) = ( KappasPr*rhoy  ) 

  ! viscous fluxes in z-direction      
  ! hJac(i,j)=d h(U_i) / d U_j
  hJac(:,:  ,i) = 0.
#endif
  

END DO !i
END SUBROUTINE EvalDiffFluxJacobian
#endif /*parabolic*/
#endif /*EQNSYS*/

END MODULE MOD_Jacobian
