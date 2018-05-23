!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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

!===================================================================================================================================
!> Contains global variables provided by the walldistance routines
!===================================================================================================================================
MODULE MOD_Walldistance_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER             :: NSuper                    !< Polynomial degree used for supersampling for coarse search
INTEGER             :: NVisu                     !< Polynomial degree used for visualization
REAL,ALLOCATABLE    :: xSuper_Face(:,:,:,:)      !< Face coordinates (of BCs) on super sampling points
REAL,ALLOCATABLE    :: distance(:,:,:,:)         !< wall distance
INTEGER,ALLOCATABLE :: nearestFace(:,:,:,:,:)    !< Index of nearest supersampling point and side
REAL,ALLOCATABLE    :: D(:,:)                    !< Derivative matrix

! Vandermonde matrices
REAL,ALLOCATABLE    :: Vdm_GaussN_EquiNSuper(:,:)   !< Vandermonde from Gauss mesh points to equdistant mesh points in supersampling
!===================================================================================================================================
END MODULE MOD_Walldistance_Vars
