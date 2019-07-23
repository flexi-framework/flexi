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
!==================================================================================================================================
!> Contains global variables used by the equation specific analyze modules.
!==================================================================================================================================
MODULE MOD_AnalyzeEquation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! precomputed variables
LOGICAL,ALLOCATABLE  :: isWall(:)                         !< marks wheter a boundary condition is a wall boundary
INTEGER              :: maxlen                            !< max length of BCName

! Variables for the specific analyze routines
LOGICAL              :: doCalcBodyForces    =.FALSE.      !< marks if body forces at walls shall be computed
LOGICAL              :: doCalcBulkState     =.FALSE.      !< marks if bulk state shall be computed
LOGICAL              :: doCalcMeanFlux      =.FALSE.      !< marks if mean flux shall be computed
LOGICAL              :: doCalcTotalStates   =.FALSE.      !< marks if total states (pt, Tt) shall be computed
LOGICAL              :: doCalcWallVelocity  =.FALSE.      !< marks if wall velocity shall be computed
LOGICAL              :: doCalcResiduals     =.FALSE.      !< marks if residuals shall be computed
LOGICAL              :: doWriteBodyForces   =.FALSE.      !< marks if body forces at walls shall be written to a file
LOGICAL              :: doWriteBulkState    =.FALSE.      !< marks if bulk state shall be written to a file
LOGICAL              :: doWriteMeanFlux     =.FALSE.      !< marks if mean flux shall be written to a file
LOGICAL              :: doWriteTotalStates  =.FALSE.      !< marks if total states (pt, Tt) shall be written to a file
LOGICAL              :: doWriteWallVelocity =.FALSE.      !< marks if wall velocity shall be written to a file
LOGICAL              :: doWriteResiduals    =.FALSE.      !< marks if residuals shall be written to a file
CHARACTER(LEN=255),ALLOCATABLE :: Filename_BodyForce(:)   !< output files for bodyforces per BC
CHARACTER(LEN=255)             :: Filename_Bulk           !< output file  for bulk velocity
CHARACTER(LEN=255),ALLOCATABLE :: Filename_MeanFlux(:)    !< output files for mean flux per BC
CHARACTER(LEN=255),ALLOCATABLE :: Filename_TotalStates(:) !< output files for total states per BC
CHARACTER(LEN=255),ALLOCATABLE :: Filename_WallVel(:)     !< output files for wall velocities per BC
CHARACTER(LEN=255)             :: Filename_Residuals      !< output file  for residuals

! Time averaging and fluctuation variables

LOGICAL              :: doCalcTimeAverage   =.FALSE.      !< marks if time averaging should be performed
LOGICAL              :: doCalcFluctuations  =.FALSE.      !< marks if time fluctuations should be computed
REAL   ,ALLOCATABLE  :: UAvg(:,:,:,:,:)                   !< time averaged solution U
REAL   ,ALLOCATABLE  :: UFluc(:,:,:,:,:)                  !< time averaged solution squared (U^2)
LOGICAL,ALLOCATABLE  :: CalcAvg(:)                        !< variables for which time averages should be computed (global indexing)
LOGICAL,ALLOCATABLE  :: CalcFluc(:)                       !< variables for which fluctuations should be computed (global indexing)
INTEGER,ALLOCATABLE  :: iAvg(:)                           !< map from (global) VariableList to index in UAvg array
INTEGER,ALLOCATABLE  :: iFluc(:)                          !< map from (global) VariableList to index in UFluc array
INTEGER,ALLOCATABLE  :: FlucAvgMap(:,:)                   !< map from index in UFluc array to index in UAvg array
                                                          !< (e.g. for mixed term uv: iFluc(1,1) -> u iFluc(2,1) -> v)
INTEGER              :: nVarAvg                           !< number of time averag variables
INTEGER              :: nVarFluc                          !< number of fluctuation variables
INTEGER              :: nVarFlucHasAvg                    !< number of fluctuations depending only on one time average
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAvgOut(:)       !< time averaged variable names
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesFlucOut(:)      !< fluctuation variable names
REAL                 :: dtAvg                             !< sum of timesteps
REAL                 :: dtOld                             !< dt from previous iteration
!==================================================================================================================================
END MODULE MOD_AnalyzeEquation_Vars
