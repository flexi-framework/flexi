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
!> Contains routines to interpolate the interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFace
IMPLICIT NONE
PRIVATE
#define WITHnVar 1

INTERFACE ProlongToFace
  MODULE PROCEDURE ProlongToFace
END INTERFACE

INTERFACE EvalElemFace
  MODULE PROCEDURE EvalElemFaceG
  MODULE PROCEDURE EvalElemFaceGL
END INTERFACE

PUBLIC::ProlongToFace,EvalElemFace

CONTAINS
#include "prolongtoface.t90"
END MODULE MOD_ProlongToFace

!==================================================================================================================================

MODULE MOD_ProlongToFaceCons
IMPLICIT NONE
PRIVATE
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE ProlongToFaceCons
  MODULE PROCEDURE ProlongToFace
END INTERFACE

PUBLIC::ProlongToFaceCons

CONTAINS
#include "prolongtoface.t90"
END MODULE MOD_ProlongToFaceCons

!==================================================================================================================================

MODULE MOD_ProlongToFacePrim
IMPLICIT NONE
PRIVATE
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE ProlongToFacePrim
  MODULE PROCEDURE ProlongToFace
END INTERFACE

PUBLIC::ProlongToFacePrim

CONTAINS
#include "prolongtoface.t90"
END MODULE MOD_ProlongToFacePrim

!==================================================================================================================================

MODULE MOD_ProlongToFace1
! MODULES
IMPLICIT NONE
PRIVATE
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = 1

INTERFACE ProlongToFace1
  MODULE PROCEDURE ProlongToFace
END INTERFACE

PUBLIC::ProlongToFace1

CONTAINS
#include "prolongtoface.t90"
END MODULE MOD_ProlongToFace1

!==================================================================================================================================

