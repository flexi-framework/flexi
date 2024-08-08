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
!==================================================================================================================================
!> Contains global variables used by the sponge modules.
!==================================================================================================================================
MODULE MOD_Sponge_Vars
! MODULES
USE MOD_Areas_Vars
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                        :: doSponge                  !< Turn on to employ sponge regions for reducing reflections at boundaries
INTEGER                        :: nSpongeRamps              !< number of sponge ramps
TYPE(tArea),ALLOCATABLE,TARGET :: Sponges(:)                !< array containing all sponge ramps
LOGICAL                        :: SpongeViz                 !< Turn on to write a visualization file of the sponge region and strength
LOGICAL                        :: WriteSponge               !< Turn on to write the sponge region and strength to the state file
INTEGER                        :: nSpongeElems              !< number of elements for which sponge is applied
INTEGER,ALLOCATABLE            :: spongeMap(:)              !< mapping from Elem -> spongElem
REAL,ALLOCATABLE,TARGET        :: SpBaseFlow(:,:,:,:,:)     !< precompute global reference state for whole field
REAL,POINTER                   :: SpBaseFlow_p(:,:,:,:,:)   !< Ponter to SpBaseFlow
INTEGER                        :: SpBaseFlowType            !< Specifies the type of baseflow
REAL,ALLOCATABLE               :: damping(:)                !< Strength of damping per face
REAL,ALLOCATABLE               :: SpongeMat(:,:,:,:)        !< precomputed sponge functions per DOF and sponge elem
REAL,ALLOCATABLE               :: SpongeMat_Out(:,:,:,:,:)  !< precomputed sponge functions per DOF and sponge elem
REAL,ALLOCATABLE               :: tempFilterWidthSp(:)      !< Filter width of each sponge region
REAL,ALLOCATABLE               :: SpongeDistance(:)         !< Array containing the distance of the ramping of the sponge
CHARACTER(LEN=255)             :: SpBaseFlowFile            !< File contiaining the sponge baseflow
!==================================================================================================================================
END MODULE MOD_Sponge_Vars
