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
INTERFACE dConsdPrimTemp
  MODULE PROCEDURE dConsdPrim
END INTERFACE

INTERFACE dPrimTempdCons
  MODULE PROCEDURE dPrimdCons
END INTERFACE

PUBLIC::EvalAdvFluxJacobian
#if PARABOLIC
PUBLIC::EvalDiffFluxJacobian
PUBLIC::EvalFluxGradJacobian
#endif
PUBLIC::dConsdPrimTemp,dPrimTempdCons
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> RANS-SA equations:
!> The Jacobian of the analytical advective Flux with respect to the Variable U
!===================================================================================================================================
SUBROUTINE EvalAdvFluxJacobian(U,UPrim,fJac,gJac,hJac)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars           ,ONLY:nDOFElem
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
INTEGER :: i
!===================================================================================================================================
DO i=1,nDOFElem
  CALL EvalAdvFluxJacobianPoint(U(:,i),UPrim(:,i),fJac(:,:,i),gJac(:,:,i),hJac(:,:,i))
END DO !i
END SUBROUTINE EvalAdvFluxJacobian

!===================================================================================================================================
!> RANS-SA equations:
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
REAL    :: a1,phi,muT
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
muT  = U(6)


fJac(1,1:6)= (/          0.,             1.,          0.,           0.,        0.,    0. /)
fJac(2,1:6)= (/      phi-uu, (1-kappaM2)*v1, -kappaM1*v2,  -kappaM1*v3,   kappaM1,    0. /)
fJac(3,1:6)= (/         -uv,             v2,          v1,           0.,        0.,    0. /)
fJac(4,1:6)= (/         -uw,             v3,          0.,           v1,        0.,    0. /)
fJac(5,1:6)= (/ v1*(phi-a1),  a1-kappaM1*uu, -kappaM1*uv,  -kappaM1*uw,  kappa*v1,    0. /)
fJac(6,1:6)= (/-muT*v1*srho,       muT*srho,          0.,           0.,        0.,    v1 /)


gJac(1,1:6)= (/          0.,           0.,              1.,          0.,       0.,    0. /)
gJac(2,1:6)= (/         -uv,           v2,              v1,          0.,       0.,    0. /)
gJac(3,1:6)= (/      phi-vv,  -kappaM1*v1, (1.-kappaM2)*v2, -kappaM1*v3,  kappaM1,    0. /)
gJac(4,1:6)= (/         -vw,           0.,              v3,          v2,       0.,    0. /)
gJac(5,1:6)= (/ v2*(phi-a1),  -kappaM1*uv,   a1-kappaM1*vv, -kappaM1*vw, kappa*v2,    0. /)
gJac(6,1:6)= (/-muT*v2*srho,           0.,        muT*srho,          0.,       0.,    v2 /)


hJac(1,1:6)= (/          0.,          0.,           0.,              1.,       0.,    0. /)
hJac(2,1:6)= (/         -uw,          v3,           0.,              v1,       0.,    0. /)
hJac(3,1:6)= (/         -vw,          0.,           v3,              v2,       0.,    0. /)
hJac(4,1:6)= (/      phi-ww, -kappaM1*v1,  -kappaM1*v2, (1.-kappaM2)*v3,  kappaM1,    0. /)
hJac(5,1:6)= (/ v3*(phi-a1), -kappaM1*uw,  -kappaM1*vw,   a1-kappaM1*ww, kappa*v3,    0. /)
gJac(6,1:6)= (/-muT*v3*srho,           0.,          0.,        muT*srho,       0.,    v3 /)
#else
absu=uu+vv
phi  = kappaM1*0.5*absu
a1   = kappa * U(5)*sRho - phi


fJac(1,1:6)= (/          0.,             1.,          0.,    0.,        0.,    0. /)
fJac(2,1:6)= (/      phi-uu, (1-kappaM2)*v1, -kappaM1*v2,    0.,   kappaM1,    0. /)
fJac(3,1:6)= (/         -uv,             v2,          v1,    0.,        0.,    0. /)
fJac(4,1:6)= (/          0.,             0.,          0.,    0.,        0.,    0. /)
fJac(5,1:6)= (/ v1*(phi-a1),  a1-kappaM1*uu, -kappaM1*uv,    0.,  kappa*v1,    0. /)
fJac(6,1:6)= (/-muT*v1*srho,       muT*srho,          0.,    0.,        0.,    v1 /)


gJac(1,1:6)= (/          0.,          0.,              1.,   0.,        0.,    0. /)
gJac(2,1:6)= (/         -uv,          v2,              v1,   0.,        0.,    0. /)
gJac(3,1:6)= (/      phi-vv, -kappaM1*v1, (1.-kappaM2)*v2,   0.,   kappaM1,    0. /)
gJac(4,1:6)= (/          0.,          0.,              0.,   0.,        0.,    0. /)
gJac(5,1:6)= (/ v2*(phi-a1), -kappaM1*uv,   a1-kappaM1*vv,   0.,  kappa*v2,    0. /)
gJac(6,1:6)= (/-muT*v2*srho,           0.,        muT*srho,  0.,        0.,    v2 /)

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
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------  
!----------------------------------------------------------------------------------------------------------------------------------  
INTEGER,INTENT(IN)                                   :: nDOF_loc             !< number of degrees of freedom
REAL,DIMENSION(PP_nVar        ,nDOF_loc),INTENT(IN)  :: U                    !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim    ,nDOF_loc),INTENT(IN)  :: UPrim                !< solution in primitive variables
REAL,DIMENSION(PP_nVarPrim    ,nDOF_loc),INTENT(IN)  :: gradUx,gradUy,gradUz !< primitive gradients
REAL,DIMENSION(PP_nVar,PP_nVar,nDOF_loc),INTENT(OUT) :: fJac,gJac,hJac       !< Derivative of the Cartesian fluxes (iVar,i,j,k)
#if EDDYVISCOSITY
REAL,DIMENSION(1              ,nDOF_loc),INTENT(IN)  :: muSGS                !< eddyviscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL abort(__STAMP__,'DiffFlux Jacobian for RANS_SA equations not yet implemented!')
fJac = 0.
gJac = 0.
hJac = 0.
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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                              :: nDOF_loc !< number of degrees of freedom
REAL,DIMENSION(PP_nVar    ,nDOF_loc),INTENT(IN) :: U        !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim,nDOF_loc),INTENT(IN) :: UPrim    !< solution in primitive variables
#if EDDYVISCOSITY
REAL,DIMENSION(1          ,nDOF_loc),INTENT(IN) :: muSGS    !< eddyviscosity
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar,nDOF_loc),INTENT(OUT) :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz !<
                                                        !> Jacobian of the diffusive Cartesian fluxes (iVar,i,j,k) w.r.t gradients
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
fJacQx = 0.
fJacQy = 0.
fJacQz = 0.

gJacQx = 0.
gJacQy = 0.
gJacQz = 0.

#if PP_dim==3
hJacQx = 0.
hJacQy = 0.
hJacQz = 0.
#endif
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
REAL,DIMENSION(PP_nVar,PP_nVarPrim),INTENT(OUT) :: Jac      !< cons to prim Jacobian
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
! dependency on temperature
Jac(1:5,6) = 0.
! dependency on SA viscosity
Jac(1,7)   = UE(NUSA)
Jac(2:5,7) = 0.
Jac(6,7)   = UE(DENS)
END SUBROUTINE dConsdPrim

!===================================================================================================================================
!> The Jacobian of the transformation from conservative to primitive variables
!===================================================================================================================================
SUBROUTINE dPrimdCons(UPrim,Jac)
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
REAL                                            :: UE(PP_2Var)
REAL                                            :: sRhoR,dpdU(5)
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

! kinematic SA viscosity
Jac(1:6,6) = 0.
Jac(7,1)   = -UE(NUSA) * UE(SRHO)
Jac(7,2:5) = 0.
Jac(7,6)   = UE(SRHO)
END SUBROUTINE dPrimdCons

END MODULE MOD_Jacobian
