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
!==================================================================================================================================
!> Contains the routine EvalFlux3D which computes the complete flux f,g,h for all DOFs in one Element: used in volume integral
!==================================================================================================================================
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
!> Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFlux3D(NLoc,ULoc,dummy,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                 :: NLoc     !< Polynomial degree
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(IN)  :: ULoc     ! Solution
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(IN)  :: dummy    ! primitive solution (useless here)
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(OUT) :: f        ! Cartesian fluxe (iVar,i,j,k)
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(OUT) :: g        ! Cartesian fluxe (iVar,i,j,k)
REAL,DIMENSION(1,0:NLoc,0:NLoc,0:NLoc),INTENT(OUT) :: h        ! Cartesian fluxe (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
f = AdvVel(1)*Uloc(:,:,:,:)
g = AdvVel(2)*Uloc(:,:,:,:)
h = AdvVel(3)*Uloc(:,:,:,:)
END SUBROUTINE EvalFlux3D

#if PARABOLIC
SUBROUTINE EvalDiffFlux2D(Nloc,f,g,h,U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
!==================================================================================================================================
! Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for every volume Gauss point.
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: Nloc
REAL,DIMENSION(1,0:Nloc,0:Nloc),INTENT(IN)  :: U_Face
REAL,DIMENSION(1,0:Nloc,0:Nloc),INTENT(IN)  :: gradUx_Face,gradUy_Face,gradUz_Face
REAL,DIMENSION(1,0:Nloc,0:Nloc),INTENT(OUT) :: f,g,h
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
f = -DiffC*gradUx_Face
g = -DiffC*gradUy_Face
h = -DiffC*gradUz_Face
END SUBROUTINE EvalDiffFlux2D

SUBROUTINE EvalDiffFlux2D_overwrite(U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
!==================================================================================================================================
! Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for every volume Gauss point.
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(1,0:PP_N,0:PP_N),INTENT(IN)    :: U_Face
REAL,DIMENSION(1,0:PP_N,0:PP_N),INTENT(INOUT) :: gradUx_Face,gradUy_Face,gradUz_Face
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
gradUx_Face = -DiffC*gradUx_Face
gradUy_Face = -DiffC*gradUy_Face
gradUz_Face = -DiffC*gradUz_Face
END SUBROUTINE EvalDiffFlux2D_overwrite

!==================================================================================================================================
!> Compute linear scalar diffusion fluxes with diffusion coefficient DiffC using the conservative
!> variables for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D(U,gradUx,gradUy,gradUz,f,g,h,iElem)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)  :: U             !< Solution
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)  :: gradUx        !< Gradient in x-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)  :: gradUy        !< Gradient in y-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)  :: gradUz        !< Gradient in z-direction
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: f             !< Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: g             !< Cartesian fluxes (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: h             !< Cartesian fluxes (iVar,i,j,k)
INTEGER, INTENT(IN)                                      :: iELem         !< Element index, not nedded in LinAdv but for Navier-
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
f = -DiffC*gradUx(:,:,:,:)
g = -DiffC*gradUy(:,:,:,:)
h = -DiffC*gradUz(:,:,:,:)
END SUBROUTINE EvalDiffFlux3D

SUBROUTINE EvalDiffFlux3D_overwrite(ULoc,gradUx,gradUy,gradUz)
!==================================================================================================================================
! Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for every volume Gauss point.
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN)    :: ULoc
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: gradUx
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: gradUy
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: gradUz
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
gradUx = -DiffC*gradUx(:,:,:,:)
gradUy = -DiffC*gradUy(:,:,:,:)
gradUz = -DiffC*gradUz(:,:,:,:)
END SUBROUTINE EvalDiffFlux3D_overwrite
#endif /*PARABOLIC*/

END MODULE MOD_Flux
