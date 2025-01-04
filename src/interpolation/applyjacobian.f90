!==================================================================================================================================
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
!==================================================================================================================================
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> \brief Module containing routines that changes scalar solution representation between physical and reference coordinates
!> Depending on interface parameters it distinguishes between:
!> - Conversion with separate input and output variables
!> - Conversion where input will be overwritten by the transformed output
!> - Selective conversion only of DG elements where the input will be overwritten by the transformed output
!==================================================================================================================================
MODULE MOD_ApplyJacobian
IMPLICIT NONE
PRIVATE

#define WITHnVar 1

INTERFACE ApplyJacobian
   MODULE PROCEDURE ApplyJacobian
   MODULE PROCEDURE ApplyJacobian_local
   MODULE PROCEDURE ApplyJacobian_select
END INTERFACE

PUBLIC:: ApplyJacobian

CONTAINS
#include "applyjacobian.t90"
END MODULE MOD_ApplyJacobian

!==================================================================================================================================
!> \brief Module containing routines that changes conservative solution representation between physical and reference coordinates
!> Depending on interface parameters it distinguishes between:
!> - Conversion with separate input and output variables
!> - Conversion where input will be overwritten by the transformed output
!> - Selective conversion only of DG elements where the input will be overwritten by the transformed output
!==================================================================================================================================
MODULE MOD_ApplyJacobianCons
IMPLICIT NONE
PRIVATE

#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE ApplyJacobianCons
   MODULE PROCEDURE ApplyJacobian
   MODULE PROCEDURE ApplyJacobian_local
   MODULE PROCEDURE ApplyJacobian_select
END INTERFACE

PUBLIC:: ApplyJacobianCons

CONTAINS
#include "applyjacobian.t90"
END MODULE MOD_ApplyJacobianCons

!==================================================================================================================================
!> \brief Module containing routines that change primitive solution representation between physical and reference coordinates
!> Depending on interface parameters it distinguishes between:
!> - Conversion with separate input and output variables
!> - Conversion where input will be overwritten by the transformed output
!> - Selective conversion only of DG elements where the input will be overwritten by the transformed output
!==================================================================================================================================
MODULE MOD_ApplyJacobianPrim
IMPLICIT NONE
PRIVATE

#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE ApplyJacobianPrim
   MODULE PROCEDURE ApplyJacobian
   MODULE PROCEDURE ApplyJacobian_local
   MODULE PROCEDURE ApplyJacobian_select
END INTERFACE

PUBLIC:: ApplyJacobianPrim

CONTAINS
#include "applyjacobian.t90"
END MODULE MOD_ApplyJacobianPrim

!==================================================================================================================================
!> \brief Module containing routines that change primitive solution representation between physical and reference coordinates
!> Depending on interface parameters it distinguishes between:
!> - Conversion with separate input and output variables
!> - Conversion where input will be overwritten by the transformed output
!> - Selective conversion only of DG elements where the input will be overwritten by the transformed output
!==================================================================================================================================
MODULE MOD_ApplyJacobianLifting
IMPLICIT NONE
PRIVATE

#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVarLifting

INTERFACE ApplyJacobianLifting
   MODULE PROCEDURE ApplyJacobian
   MODULE PROCEDURE ApplyJacobian_local
   MODULE PROCEDURE ApplyJacobian_select
END INTERFACE

PUBLIC:: ApplyJacobianLifting

CONTAINS
#include "applyjacobian.t90"
END MODULE MOD_ApplyJacobianLifting
