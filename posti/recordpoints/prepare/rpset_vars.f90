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
MODULE MOD_RPSet_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                         :: RPSetInitIsDone=.FALSE.
INTEGER                         :: nRP_global                     !< Total number of record points
INTEGER                         :: nGroups                        !< Number of defined groups
INTEGER                         :: nPoints                        !< Number of single points
INTEGER                         :: nLines                         !< Number of line type structures
INTEGER                         :: nPlanes                        !< Number of plane type structures
INTEGER                         :: nBoxes                         !< Number of box type structures
INTEGER,ALLOCATABLE             :: OffsetRP(:,:)                  !< Offest array used for parallel input/output

TYPE tGroup                                                       !< Data type representing a single group
  CHARACTER(LEN=255)               :: Name         !< Name of the group
  INTEGER                          :: ID           !< ID of the group
  INTEGER                          :: nRP          !< Number of record points in the group
  TYPE(tRP_Ptr),POINTER            :: RP_ptr(:)    !< Array of pointers to the record points of the group
END TYPE tGroup

TYPE tPoint                                                       !< Data type representing a single point
  INTEGER                          :: GroupID       !< ID of the group the point belongs to
  TYPE(tRP),POINTER                :: RP            !< Pointer to the record points that makes up the single point
END TYPE tPoint

TYPE tLine                                                         !< Data type representing a single line
  CHARACTER(LEN=255)               :: Name         !< Name of the line
  INTEGER                          :: GroupID      !< ID of the group the line belongs to
  REAL                             :: xStart(3)    !< Start coordinate of the line
  REAL                             :: xEnd(3)      !< End coordinate of the line
  INTEGER                          :: nRP          !< Number of record points in the group
  TYPE(tRP_Ptr),POINTER            :: RP_ptr(:)    !< Array of pointers to the record points of the line
END TYPE tLine

TYPE tPlane                                                         !< Data type representing a single plane
  CHARACTER(LEN=255)               :: Name         !< Name of the plane
  INTEGER                          :: GroupID      !< ID of the group the plane belongs to
  REAL                             :: x(3,4)       !< 4 corner points of the plane
  INTEGER                          :: nRP(2)       !< Number of recordpoints in the two plane dimensions
  TYPE(tRP_Ptr),POINTER            :: RP_ptr(:,:)  !< Array of pointers to the record points of the group
  REAL,ALLOCATABLE                 :: NormVec(:,:) !< Normal vector at each of the points
  REAL,ALLOCATABLE                 :: TangVec(:,:) !< Tangential vector at each of the points
END TYPE tPlane

TYPE tBox                                                           !< Data type representing a single box
  CHARACTER(LEN=255)               :: Name            !< Name of the plane
  INTEGER                          :: GroupID         !< ID of the group the plane belongs to
  REAL                             :: x(3,8)          !< 4 corner points of the plane
  INTEGER                          :: nRP(3)          !< Number of recordpoints in the two plane dimensions
  TYPE(tRP_Ptr),POINTER            :: RP_ptr(:,:,:)   !< Array of pointers to the record points of the group
  REAL,ALLOCATABLE                 :: NormVec(:,:,:)  !< Normal vector at each of the points
  REAL,ALLOCATABLE                 :: TangVec1(:,:,:) !< Tangential vector at each of the points
  REAL,ALLOCATABLE                 :: TangVec2(:,:,:) !< Tangential vector at each of the points
END TYPE tBox

TYPE tRP                                                            !< Data type representing a single record point
  INTEGER                          :: ID           !< Identifier for the record point
  INTEGER                          :: ElemID       !< ID of the element in the mesh that the record point belongs to
  INTEGER                          :: GroupID      !< ID of the group the line belongs to
  REAL                             :: xi(3)        !< Reference coordinates in the mesh of the record point
  REAL                             :: x(3)         !< Physical coordinates of the record point
  REAL                             :: xF(3)        !< Temporary value used in Newton
END TYPE tRP

TYPE tRPlist                                                        !< Data structure to store list of all the record points
  TYPE(tRP),POINTER                :: RP           !< Pointer to the record point data structure
END TYPE tRPlist

TYPE tRP_Ptr                                                        !< Data structure for pointers to record points
  TYPE(tRP),POINTER                :: RP           !< node pointer
END TYPE tRP_Ptr



TYPE(tGroup),POINTER            :: Groups(:)                        !< Array of pointers to all groups
TYPE(tPoint),POINTER            :: Points(:)                        !< Array of pointers to all single points (not all recordpoints)
TYPE(tLine),POINTER             :: Lines(:)                         !< Array of pointers to all lines
TYPE(tPlane),POINTER            :: Planes(:)                        !< Array of pointers to all planes
TYPE(tBox),POINTER              :: Boxes(:)                         !< Array of pointers to all planes
TYPE(tRPlist),POINTER           :: RPlist(:)                        !< Array of pointers to all record points

INTERFACE GetNewRP
  MODULE PROCEDURE GetNewRP
END INTERFACE

PUBLIC:: GetNewRP
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Will create a single new record point data structure, defined by the group id and the coordinates
!===================================================================================================================================
SUBROUTINE GetNewRP(RP,GroupID,x)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tRP),POINTER ::  RP
INTEGER           ::  GroupID
REAL              ::  x(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ALLOCATE(RP)
nRP_global=nRP_global+1
RP%ID=nRP_global
RP%ElemID=-1
RP%GroupID=GroupID
Groups(GroupID)%nRP= Groups(GroupID)%nRP+1
RP%xi=HUGE(9.)
RP%xF=HUGE(9.)
RP%x(1:3)=x(1:3)
!WRITE(*,*) 'id,groupid',RP%ID,RP%GroupID

END SUBROUTINE GetNewRP

END MODULE MOD_RPSet_Vars
