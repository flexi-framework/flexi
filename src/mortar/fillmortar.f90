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

!==================================================================================================================================
!> \brief Routines that perform the projection operation between nonconforming interfaces using the operators set up in module
!> mortar
!> 
!> Contains the routines to
!> - interpolate the solution at the large sides to the small ones, which are used for flux computation
!> - project the flux from the small sides back to the large ones
!==================================================================================================================================
MODULE MOD_FillMortar
IMPLICIT NONE
PRIVATE

#define TP_NZ PP_NZ
#define WITHnVar 1

INTERFACE U_Mortar
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_Mortar
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

PUBLIC::U_Mortar,Flux_Mortar

CONTAINS
#include "fillmortar.t90"
END MODULE MOD_FillMortar

!==================================================================================================================================

MODULE MOD_FillMortarCons
IMPLICIT NONE
PRIVATE

#undef TP_NZ
#define TP_NZ PP_NZ
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE U_MortarCons
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_MortarCons
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

PUBLIC::U_MortarCons,Flux_MortarCons

CONTAINS
#include "fillmortar.t90"
END MODULE MOD_FillMortarCons

!==================================================================================================================================

MODULE MOD_FillMortarPrim
IMPLICIT NONE
PRIVATE

#undef TP_NZ
#define TP_NZ PP_NZ
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE U_MortarPrim
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_MortarPrim
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

PUBLIC::U_MortarPrim,Flux_MortarPrim

CONTAINS
#include "fillmortar.t90"
END MODULE MOD_FillMortarPrim

!==================================================================================================================================

MODULE MOD_FillMortar1
IMPLICIT NONE
PRIVATE

#undef TP_NZ
#define TP_NZ PP_NZ
INTEGER,PARAMETER :: TP_nVar = 1

INTERFACE U_Mortar1
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_Mortar1
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

PUBLIC::U_Mortar1,Flux_Mortar1

CONTAINS
#include "fillmortar.t90"
END MODULE MOD_FillMortar1

!==================================================================================================================================

MODULE MOD_FillMortar1_3D
IMPLICIT NONE
PRIVATE

#undef TP_NZ
#define TP_NZ PP_N
INTEGER,PARAMETER :: TP_nVar = 1

INTERFACE U_Mortar1
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_Mortar1
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

PUBLIC::U_Mortar1,Flux_Mortar1

CONTAINS
#include "fillmortar.t90"
END MODULE MOD_FillMortar1_3D

