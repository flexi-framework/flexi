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
!> Contains analyze routines specific to the linear scalar advection equation
!==================================================================================================================================
MODULE MOD_AnalyzeEquation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: AnalyzeEquation, InitAnalyzeEquation, FinalizeAnalyzeEquation
PUBLIC::DefineParametersAnalyzeEquation
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeEquation()
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE DefineParametersAnalyzeEquation


!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> Equation specific analyze routine
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

END SUBROUTINE AnalyzeEquation


!==================================================================================================================================
!> Finalizes variables necessary for analyze subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
