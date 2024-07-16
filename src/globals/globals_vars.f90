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

!==================================================================================================================================
!> Provides parameters, used globally (please use EXTREMLY carefully!)
!==================================================================================================================================
MODULE MOD_Globals_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=6),PARAMETER :: ProgramName  = 'FLEXI'               !> name of this program
REAL                       :: StartTime                            !< start time of the simulation
REAL                       :: WallTime                             !> Wall time needed by a simulation (is not reset by
                                                                   !> performing a load balance step, only by user restart)
REAL                       :: InitializationWallTime               !> Wall time needed to initialize a simulation (or
                                                                   !> re-initialize a simulation by performing a load balance
                                                                   !>  step)
REAL                       :: ReadMeshWallTime                     !> Wall time needed to read the mesh (SUBROUTINE ReadMesh)
REAL                       :: DomainDecompositionWallTime          !> Wall time needed for domain decomposition
REAL                       :: CommMeshReadinWallTime               !> Shared memory mesh communication
REAL                       :: SimulationEfficiency                 !> relates the simulated time to the used CPUh (SIMULATION TIME PER
                                                                   !> CALCULATION in [s]/[CPUh])
REAL                       :: StartT                               !> Timer start
REAL,PARAMETER             :: epsMach    = EPSILON(0.)
!===================================================================================================================================

END MODULE MOD_Globals_Vars
