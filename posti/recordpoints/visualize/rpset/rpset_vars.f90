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
!===================================================================================================================================
!> Contains global variables used for/by the RPSet
!===================================================================================================================================
MODULE MOD_RPSetVisuVisu_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: nRP_global              !< Global number of record points
INTEGER                         :: nRP_HDF5                !< Number of record points in the HDF5 file
INTEGER                         :: nGroups                 !< Number of groups from RP set definition
CHARACTER(LEN=255),ALLOCATABLE  :: GroupNames(:)           !< Names of all groups
INTEGER                         :: nLines                  !< Number of lines for output
INTEGER                         :: nPoints                 !< Number of points for output
INTEGER                         :: nPlanes                 !< Number of planes for output
LOGICAL,ALLOCATABLE             :: OutputGroup(:)          !< Indicating groups for output
INTEGER,ALLOCATABLE             :: Points_IDlist(:)        !< IDs of all points for output
INTEGER,ALLOCATABLE             :: Points_GroupIDlist(:)   !< List of group IDs of all points for output
REAL,ALLOCATABLE                :: x_RP(:,:)               !< Coodinates of record points (may be local)
REAL,ALLOCATABLE                :: xF_RP(:,:)              !< Physical coordinates of the record points
INTEGER,ALLOCATABLE             :: RPOutMap(:)             !< Mapping from output RPs to global RPs

TYPE tLine !< Type used to organize a line set
  CHARACTER(LEN=255)            :: Name                    !< Name of the line
  INTEGER                       :: GroupID                 !< ID of the group the line belongs to
  INTEGER                       :: nRP                     !< Number of RPs in the line
  INTEGER,ALLOCATABLE           :: IDlist(:)               !< List of IDs of the RPs that make up the line
  REAL,ALLOCATABLE              :: LocalCoord(:)           !< Local coordinates of the RPs
  REAL,ALLOCATABLE              :: Tmat(:,:)               !< Transformation matrix
END TYPE tLine

TYPE tPlane !< Type used to organize a plane set
  CHARACTER(LEN=255)            :: Name                    !< Name of the line
  INTEGER                       :: GroupID                 !< ID of the group the line belongs to
  INTEGER                       :: Type=0                  !< 0 - standard, 1 - sphere, 2 - BLPlane
  INTEGER                       :: nRP(2)                  !< Number of RPs in the two plane directions 
  INTEGER,ALLOCATABLE           :: IDlist(:,:)             !< List of IDs of the RPs that make up the line
  REAL,ALLOCATABLE              :: NormVec(:,:)            !< Normal vector
  REAL,ALLOCATABLE              :: TangVec(:,:)            !< Tangential vector
  REAL,ALLOCATABLE              :: LocalCoord(:,:,:)       !< Local coordinates of the RPs
  REAL,ALLOCATABLE              :: BLProps(:,:)            !< Boundary layer properties for BLPlanes
END TYPE tPlane

TYPE(tLine),POINTER             :: Lines(:)                !< Pointer to all output lines
TYPE(tPlane),POINTER            :: Planes(:)               !< Pointer to all output planes

LOGICAL                         :: RPSetInitIsDone = .FALSE.
!===================================================================================================================================
END MODULE MOD_RPSetVisuVisu_Vars
