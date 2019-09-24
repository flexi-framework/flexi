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
!> Contains the routines to calculate the analytical derivatives of the viscous flux with respect to the gradients
!===================================================================================================================================
MODULE MOD_GradJacobian
#if PARABOLIC   
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

PUBLIC::EvalFluxGradJacobian
!===================================================================================================================================

CONTAINS


#if EQNSYSNR==1
!===================================================================================================================================
!> Computes the volume derivative of the analytical diffusive flux with respect to the gradient of U: d(F^v)/dQ, Q=grad U
!===================================================================================================================================
SUBROUTINE EvalFluxGradJacobian(U,UPrim,fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz &
#if EDDYVISCOSITY
                               ,muSGS &
#endif
                               )
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,      ONLY:nDOFElem
USE MOD_Equation_Vars,ONLY:DiffC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar    ,nDOFElem),INTENT(IN) :: U
REAL,DIMENSION(PP_nVarPrim,nDOFElem),INTENT(IN) :: UPrim
#if EDDYVISCOSITY
REAL,DIMENSION(1          ,nDOFElem),INTENT(IN) :: muSGS
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
                                                        ! Gradient of the diffusive Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,PP_nVar,Nloc),INTENT(OUT) :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
fJacQx = -DiffC
fJacQy = 0.
fJacQz = 0.

gJacQx = 0.
gJacQy = -DiffC
gJacQz = 0.

#if PP_dim==3
hJacQx = 0.
hJacQy = 0.
hJacQz = -DiffC
#endif
END SUBROUTINE EvalFluxGradJacobian
#endif /*linearscalaradvection*/

#if EQNSYSNR==2
!===================================================================================================================================
!> Computes the volume derivative of the analytical diffusive flux with respect to the gradient of U: d(F^v)/dQ, Q=grad U
!===================================================================================================================================
SUBROUTINE EvalFluxGradJacobian(U,UPrim,fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz &
#if EDDYVISCOSITY
                               ,muSGS &
#endif
                               )
! MODULES
USE MOD_PreProc
USE MOD_Viscosity
USE MOD_DG_Vars,      ONLY:nDOFElem
USE MOD_Equation_Vars,ONLY:s43,s23
USE MOD_EOS_Vars,     ONLY:cp,Pr
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: PrSGS
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar    ,nDOFElem),INTENT(IN)          :: U
REAL,DIMENSION(PP_nVarPrim,nDOFElem),INTENT(IN)          :: UPrim
#if EDDYVISCOSITY
REAL,DIMENSION(1          ,nDOFElem),INTENT(IN)          :: muSGS
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
                                                        ! Gradient of the diffusive Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,PP_nVarPrim,nDOFElem),INTENT(OUT) :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
REAL                :: muS,lambda
!===================================================================================================================================
DO i=1,nDOFElem
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
#endif /*navierstokes*/

#endif /*PARABOLIC*/
END MODULE MOD_GradJacobian
