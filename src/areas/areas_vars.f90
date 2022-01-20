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
!> Contains global variables used by the sponge modules.
!==================================================================================================================================
MODULE MOD_Areas_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

!-----------------------------------------------------------------------------------------------------------------------------------
! PARAMETER DEFINITONS
!-----------------------------------------------------------------------------------------------------------------------------------

INTEGER,PARAMETER :: SHAPE_RAMP         = 1
INTEGER,PARAMETER :: SHAPE_CYLINDRICAL  = 2
INTEGER,PARAMETER :: SHAPE_CUBOID       = 3
INTEGER,PARAMETER :: SHAPE_POLYGON      = 4

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPE DEFINITONS
!-----------------------------------------------------------------------------------------------------------------------------------

! Type used for general shape definition
TYPE tShape
  SEQUENCE
  REAL                                    :: xStart(3)           = 0.    !< Starting Point for Area Ramp
  REAL                                    :: xEnd(3)             = 0.    !< End Point of Cuboit
  REAL                                    :: ArVec(3)            = 0.    !< Vector defining the ramp direction
  REAL                                    :: cubeoidCenter(3)    = 0.    !< Center of catresian aligned cube
  REAL                                    :: ArDistance          = 0.    !< Distance of the area
  REAL                                    :: ArRadius            = 0.    !< Radius of the cylindrical (3D) / radial (2D) Area
#if(PP_dim==3)
  REAL                                    :: ArAxis(3)           = 0.    !< Axis of the cylindrical Area (only 3D)
#endif
  INTEGER                                 :: nAreaVertices       = 0.    !< Number of Vertices for each polygon Area
  REAL,ALLOCATABLE                        :: AreaVertex(:,:)             !< Vertices defining polygon Area region
END TYPE tShape

! Type used to store the relevant information about an area
TYPE tArea
  SEQUENCE
  CHARACTER(LEN=255)                  :: AreaStr
  INTEGER                             :: AreaShape
  INTEGER                             :: nAreaElems         !< number of elements for which sponge is applied
  INTEGER,ALLOCATABLE                 :: AreaMap(:)         !< mapping from Elem -> spongElem
  TYPE(tShape)                        :: Shape
END TYPE tArea

! This type is used to create an array of pointers, pointing to each of the created areas (to loop over them)
TYPE tAreaList
  TYPE(tArea),POINTER                 :: pArea              !< Pointer to the area
END TYPE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                        :: doArea                !< Turn on to employ sponge regions for reducing reflections at boundaries
LOGICAL                        :: AreaViz               !< Turn on to write a visualization file of the sponge region and strength

INTEGER                        :: nAreas=0              !< Number of interfaces
TYPE(tAreaList),ALLOCATABLE    :: Areas(:)              !< List of all active interfaces
!==================================================================================================================================
END MODULE MOD_Areas_Vars
