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
!> Contains global variables used by the DG modules.
!==================================================================================================================================
MODULE MOD_Sponge_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: doSponge              !< Turn on to employ sponge regions for reducing reflections at boundaries
LOGICAL          :: SpongeViz             !< Turn on to write a visualization file of the sponge region and strength
LOGICAL          :: CalcPruettDamping=.FALSE. !< true if temporally varying, solution adaptive Pruett baseflow is used
INTEGER          :: nSpongeELems          !< number of elements for which sponge is applied
INTEGER,ALLOCATABLE :: spongeMap(:)       !< mapping from Elem -> spongElem
REAL             :: damping               !< Strenght of damping per face
REAL             :: xStart(3)             !< Coordinates of start postion of sponge ramp (shape=1) or center (shape=2)
REAL             :: SpAxis(3)             !< Axis vector of cylindrical sponge (SpongeShape=2)
REAL             :: SpVec(3)              !< Direction vector of the spong ramp (SpongeShape=1)
REAL             :: SpDistance            !< Length of sponge ramp.
                                          !< The sponge will have maximum strength at the end of the ramp and after that point
REAL             :: SpRadius              !< Radius of the sponge zone (SpongeShape=2)
REAL,ALLOCATABLE :: SpongeMat(:,:,:,:)    !< precomputed sponge functions per DOF and sponge elem
REAL,ALLOCATABLE :: SpBaseFlow(:,:,:,:,:) !< precompute global reference state for whole field
!==================================================================================================================================
END MODULE MOD_Sponge_Vars
