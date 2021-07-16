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
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE EvalFlux3D
  MODULE PROCEDURE EvalFlux3D
END INTERFACE

#if PARABOLIC
INTERFACE EvalDiffFlux2D
  MODULE PROCEDURE EvalDiffFlux2D_overwrite
  MODULE PROCEDURE EvalDiffFlux2D
END INTERFACE
INTERFACE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D_overwrite
END INTERFACE
PUBLIC::EvalDiffFlux2D
PUBLIC::EvalDiffFlux3D
#endif

PUBLIC::EvalFlux3D
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for a single element.
!==================================================================================================================================
SUBROUTINE EvalFlux3D(Nloc,ULoc,dummy,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel
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
f = AdvVel(1)*Uloc(:,:,:,:)
g = AdvVel(2)*Uloc(:,:,:,:)
#if PP_dim==3
h = AdvVel(3)*Uloc(:,:,:,:)
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
f = -DiffC*gradUx_Face
g = -DiffC*gradUy_Face
#if PP_dim==3
h = -DiffC*gradUz_Face
#else
h = 0.
#endif
END SUBROUTINE EvalDiffFlux2D

!==================================================================================================================================
!> Compute diffusive fluxes with diffusion coefficient DiffC using the conservative variables for a single side, input will
!> be overwritten.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux2D_overwrite(U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim,   0:PP_N,0:PP_NZ),INTENT(IN)       :: U_Face         !< Solution
REAL,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_NZ),INTENT(INOUT) :: gradUx_Face    !< Gradient in x-direction
REAL,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_NZ),INTENT(INOUT) :: gradUy_Face    !< Gradient in y-direction
REAL,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_NZ),INTENT(INOUT) :: gradUz_Face    !< Gradient in z-direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
gradUx_Face = -DiffC*gradUx_Face
gradUy_Face = -DiffC*gradUy_Face
#if PP_dim==3
gradUz_Face = -DiffC*gradUz_Face
#else
gradUz_Face = 0.
#endif
END SUBROUTINE EvalDiffFlux2D_overwrite

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
f = -DiffC*gradUx(:,:,:,:)
g = -DiffC*gradUy(:,:,:,:)
#if PP_dim==3
h = -DiffC*gradUz(:,:,:,:)
#else
h = 0.
#endif
END SUBROUTINE EvalDiffFlux3D

!==================================================================================================================================
!> Compute linear scalar diffusion fluxes with diffusion coefficient DiffC using the conservative
!> variables for a single volume element, input will be overwritten.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_overwrite(ULoc,gradUx,gradUy,gradUz)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)    :: ULoc      !< Solution
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: gradUx    !< Gradient in x-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: gradUy    !< Gradient in y-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: gradUz    !< Gradient in z-direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
gradUx = -DiffC*gradUx(:,:,:,:)
gradUy = -DiffC*gradUy(:,:,:,:)
#if PP_dim==3
gradUz = -DiffC*gradUz(:,:,:,:)
#else
gradUz = 0.
#endif
END SUBROUTINE EvalDiffFlux3D_overwrite
#endif /*PARABOLIC*/

END MODULE MOD_Flux
