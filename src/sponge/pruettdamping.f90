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
!> Subroutines needed for the Sponge base flow based on a moving average of the instationary flow field, also known as Pruett
!> damping. See "The temporally filtered Navier–Stokes equations: Properties of the residual stress" for details.
!==================================================================================================================================
MODULE MOD_PruettDamping
! MODULES
IMPLICIT NONE
PRIVATE
REAL :: tempFilterWidth !< Temporal filter width used for moving average
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitPruettDamping
  MODULE PROCEDURE InitPruettDamping
END INTERFACE

INTERFACE FinalizePruettDamping
  MODULE PROCEDURE FinalizePruettDamping
END INTERFACE

INTERFACE TempFilterTimeDeriv
  MODULE PROCEDURE TempFilterTimeDeriv
END INTERFACE

PUBLIC :: InitPruettDamping,  FinalizePruettDamping
PUBLIC :: TempFilterTimeDeriv
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Read in parameters needed by Pruett damping.
!==================================================================================================================================
SUBROUTINE InitPruettDamping()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,ONLY:GETREAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PruettDamping...'

tempFilterWidth=GETREAL('tempFilterWidth')

SWRITE(UNIT_stdOut,'(A)')' INIT PruettDamping DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitPruettDamping

!==================================================================================================================================
!>Integrate the Pruett baseflow (time-filtered solution) in time using a simple Euler forward approach:
!>\f$ \frac{d}{dt} \bar{u} \approx \frac{\bar{u}^{n+1} - \bar{u}^{n}}{\Delta t}= \frac{u^n-\bar{u}^n}{\Delta} \f$
!==================================================================================================================================
SUBROUTINE TempFilterTimeDeriv(UIn,dt)
! MODULES
USE MOD_PreProc
USE MOD_Sponge_Vars,ONLY:SpBaseFlow
USE MOD_Mesh_Vars,  ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: UIn(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< Global solution array
REAL,INTENT(IN) :: dt                                       !< Current timestep
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: fac
!==================================================================================================================================
fac=dt/tempFilterWidth
SpBaseFlow = SpBaseFlow+(UIn - SpBaseFlow)*fac
END SUBROUTINE TempFilterTimeDeriv


!==================================================================================================================================
!>Finalizes variables necessary for Pruett damping.
!==================================================================================================================================
SUBROUTINE FinalizePruettDamping()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizePruettDamping


END MODULE MOD_PruettDamping
