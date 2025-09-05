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

#if PARABOLIC
INTERFACE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux2D_overwrite
  MODULE PROCEDURE EvalDiffFlux2D_Side
  MODULE PROCEDURE EvalDiffFlux2D_Point
  MODULE PROCEDURE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D_overwrite
#if FV_ENABLED
  MODULE PROCEDURE EvalDiffFlux3D_Volume_FV
#endif /*FV_ENABLED*/
END INTERFACE
#endif

PUBLIC:: EvalFlux3D
#if PARABOLIC
PUBLIC:: EvalDiffFlux3D
#endif
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute linear scalar advection fluxes with velocity AdvVel(3) using the conservative variables for a single element.
!==================================================================================================================================
SUBROUTINE EvalFlux3D(Nloc,ULoc,dummy,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel
! IMPLICIT VARIABLE HANDLING
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
SUBROUTINE EvalDiffFlux2D_Point(f,g,h,U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: U_Face        !< Solution
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx_Face   !< Gradient in x-direction
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUy_Face   !< Gradient in y-direction
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUz_Face   !< Gradient in z-direction
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: f             !< Cartesian fluxes (iVar,i,j)
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: g             !< Cartesian fluxes (iVar,i,j)
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: h             !< Cartesian fluxes (iVar,i,j)
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
END SUBROUTINE EvalDiffFlux2D_Point


!==================================================================================================================================
!> Compute diffusive fluxes with diffusion coefficient DiffC using the conservative variables for a single side with a
!> variable polynomial degree.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux2D_Side(Nloc,f,g,h,U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
! IMPLICIT VARIABLE HANDLING
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
END SUBROUTINE EvalDiffFlux2D_Side


!==================================================================================================================================
!> Compute diffusive fluxes with diffusion coefficient DiffC using the conservative variables for a single side, input will
!> be overwritten.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux2D_overwrite(U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
! IMPLICIT VARIABLE HANDLING
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
! IMPLICIT VARIABLE HANDLING
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
! IMPLICIT VARIABLE HANDLING
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


#if FV_ENABLED
!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Volume_FV(UPrim,gradUx,gradUy,gradUz,f,g,h,iElem,PP_N_xi,PP_N_eta,PP_N_zeta)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:DiffC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                                   :: iElem                !< element index in global array
INTEGER,INTENT(IN)                                                   :: PP_N_xi
INTEGER,INTENT(IN)                                                   :: PP_N_eta
INTEGER,INTENT(IN)                                                   :: PP_N_zeta
REAL,DIMENSION(PP_nVar,0:PP_N_xi,0:PP_N_eta,0:PP_N_zeta),INTENT(IN)  :: UPrim                !< Solution vector
!> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar,0:PP_N_xi,0:PP_N_eta,0:PP_N_zeta),INTENT(IN)  :: gradUx,gradUy,gradUz
!> Physical fluxes in x,y,z directions
REAL,DIMENSION(PP_nVar,0:PP_N_xi,0:PP_N_eta,0:PP_N_zeta),INTENT(OUT) :: f,g,h
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
END SUBROUTINE EvalDiffFlux3D_Volume_FV
#endif /*FV_ENABLED*/
#endif /*PARABOLIC*/

END MODULE MOD_Flux
