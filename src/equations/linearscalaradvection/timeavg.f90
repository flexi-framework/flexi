!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
!> Module to handle time averaging - not used in linear advection diffusion case!
!==================================================================================================================================
MODULE MOD_TimeAverage
! MODULES
USE MOD_Analyze_Vars
IMPLICIT NONE
PRIVATE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitCalcTimeAverage
  MODULE PROCEDURE InitCalcTimeAverage
END INTERFACE

INTERFACE FinalizeTimeAverage
  MODULE PROCEDURE FinalizeTimeAverage
END INTERFACE

INTERFACE CalcTimeAverage
  MODULE PROCEDURE CalcTimeAverage
END INTERFACE

PUBLIC::InitCalcTimeAverage, FinalizeTimeAverage, CalcTimeAverage
!==================================================================================================================================
CONTAINS



!==================================================================================================================================
!> Initializes the time averaging variable
!==================================================================================================================================
SUBROUTINE InitCalcTimeAverage()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE InitCalcTimeAverage


!==================================================================================================================================
!> Computes time averages
!==================================================================================================================================
SUBROUTINE CalcTimeAverage(Finalize,dt,t)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN)              :: Finalize
REAL,INTENT(IN)                 :: dt
REAL,INTENT(IN),OPTIONAL        :: t
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CalcTimeAverage



!==================================================================================================================================
!> Finalizes the time averaging routines
!==================================================================================================================================
SUBROUTINE FinalizeTimeAverage()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeTimeAverage

END MODULE MOD_TimeAverage
