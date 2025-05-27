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
!> Contains the routine EvalFlux3D which computes the complete flux f,g,h for all DOFs in one Element: used in volume integral
!==================================================================================================================================
#include "flexi.h"
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC::EvalFlux3D
#if PARABOLIC
PUBLIC::EvalDiffFlux2D
PUBLIC::EvalDiffFlux3D
#endif
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for a single element.
!==================================================================================================================================
SUBROUTINE EvalFlux3D(Nloc,ULoc,dummy,f,g,h)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                                 :: Nloc     !< Polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: ULoc     !< Solution
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: dummy    !< primitive solution (useless here)
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f        !< Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: g        !< Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: h        !< Cartesian fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
f(1,:,:,:) = 0.5*Uloc(1,:,:,:)*Uloc(1,:,:,:)
f(2,:,:,:) =     Uloc(1,:,:,:)*Uloc(2,:,:,:)

g(1,:,:,:) =     Uloc(2,:,:,:)*Uloc(1,:,:,:)
g(2,:,:,:) = 0.5*Uloc(2,:,:,:)*Uloc(2,:,:,:)

#if PP_dim==3
f(3,:,:,:) =     Uloc(1,:,:,:)*Uloc(3,:,:,:)
g(3,:,:,:) =     Uloc(2,:,:,:)*Uloc(3,:,:,:)
h(1,:,:,:) =     Uloc(3,:,:,:)*Uloc(1,:,:,:)
h(2,:,:,:) =     Uloc(3,:,:,:)*Uloc(2,:,:,:)
h(3,:,:,:) = 0.5*Uloc(3,:,:,:)*Uloc(3,:,:,:)
#else
f(3,:,:,:) = 0.
g(3,:,:,:) = 0.
h          = 0.
#endif

END SUBROUTINE EvalFlux3D

#if PARABOLIC
!==================================================================================================================================
!> Compute diffusive fluxes with diffusion coefficient DiffC using the conservative variables for a single side with a
!> variable polynomial degree.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux2D(Nloc,f,g,h,U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                              :: Nloc          !< Polynomial degree
REAL,DIMENSION(PP_nVarPrim,   0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_Face        !< Solution
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx_Face   !< Gradient in x-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUy_Face   !< Gradient in y-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUz_Face   !< Gradient in z-direction
REAL,DIMENSION(PP_nVar,       0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f             !< Cartesian fluxes (iVar,i,j)
REAL,DIMENSION(PP_nVar,       0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: g             !< Cartesian fluxes (iVar,i,j)
REAL,DIMENSION(PP_nVar,       0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: h             !< Cartesian fluxes (iVar,i,j)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
f(1,:,:) = -DiffC*gradUx_Face(1,:,:)
f(2,:,:) = -DiffC*gradUx_Face(2,:,:)

g(1,:,:) = -DiffC*gradUy_Face(1,:,:)
g(2,:,:) = -DiffC*gradUy_Face(2,:,:)

#if PP_dim==3
f(3,:,:) = -DiffC*gradUx_Face(3,:,:)
g(3,:,:) = -DiffC*gradUy_Face(3,:,:)
h(1,:,:) = -DiffC*gradUz_Face(1,:,:)
h(2,:,:) = -DiffC*gradUz_Face(2,:,:)
h(3,:,:) = -DiffC*gradUz_Face(3,:,:)
#else
f(3,:,:) = 0.
g(3,:,:) = 0.
h        = 0.
#endif

END SUBROUTINE EvalDiffFlux2D

!==================================================================================================================================
!> Compute linear scalar diffusion fluxes with diffusion coefficient DiffC using the conservative
!> variables for a single volume element.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D(U,gradUx,gradUy,gradUz,f,g,h,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: U             !< Solution
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: gradUx        !< Gradient in x-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: gradUy        !< Gradient in y-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)  :: gradUz        !< Gradient in z-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT) :: f             !< Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT) :: g             !< Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT) :: h             !< Cartesian fluxes (iVar,i,j,k)
INTEGER, INTENT(IN)                                       :: iELem         !< Element index, not nedded in LinAdv but for NSE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
f(1,:,:,:) = -DiffC*gradUx(1,:,:,:)
f(2,:,:,:) = -DiffC*gradUx(2,:,:,:)

g(1,:,:,:) = -DiffC*gradUy(1,:,:,:)
g(2,:,:,:) = -DiffC*gradUy(2,:,:,:)

#if PP_dim==3
f(3,:,:,:) = -DiffC*gradUx(3,:,:,:)
g(3,:,:,:) = -DiffC*gradUy(3,:,:,:)
h(1,:,:,:) = -DiffC*gradUz(1,:,:,:)
h(2,:,:,:) = -DiffC*gradUz(2,:,:,:)
h(3,:,:,:) = -DiffC*gradUz(3,:,:,:)
#else
f(3,:,:,:) = 0.
g(3,:,:,:) = 0.
h          = 0.
#endif

END SUBROUTINE EvalDiffFlux3D

#endif /*PARABOLIC*/

END MODULE MOD_Flux
