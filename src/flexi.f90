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
!> Control program of the Flexi code.
!==================================================================================================================================
PROGRAM Flexi
! MODULES
USE MOD_Globals
USE MOD_Flexi
USE MOD_TimeDisc,          ONLY:TimeDisc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Initialize
CALL InitFlexi(0)

! Run Simulation
CALL TimeDisc()

! Finalize
CALL FinalizeFlexi()
#if USE_MPI
! we also have to finalize MPI itself here
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif
END PROGRAM Flexi
