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

!===================================================================================================================================
!> Contains the computation of surface Jacobians on big mortar sides from the Jacobians on the small mortar sides.
!===================================================================================================================================
MODULE MOD_Jac_Ex_MortarU
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Jacobian_MortarU
  MODULE PROCEDURE Jacobian_Mortar
END INTERFACE

PUBLIC::Jacobian_MortarU
!===================================================================================================================================

INTEGER,PARAMETER :: TP_nVar = PP_nVar
CONTAINS

#include "jac_ex_mortar.t90"

END MODULE MOD_Jac_Ex_MortarU

MODULE MOD_Jac_Ex_MortarGrad
#if PARABOLIC
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Jacobian_MortarGrad
  MODULE PROCEDURE Jacobian_Mortar
END INTERFACE

PUBLIC::Jacobian_MortarGrad
!===================================================================================================================================

INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim
CONTAINS

#include "jac_ex_mortar.t90"

#endif
END MODULE MOD_Jac_Ex_MortarGrad

