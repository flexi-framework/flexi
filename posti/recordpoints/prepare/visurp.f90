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

!===================================================================================================================================
!> Module containing Record point visualization routines
!===================================================================================================================================
MODULE MOD_VisuRP
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE VisuRP
  MODULE PROCEDURE VisuRP
END INTERFACE

PUBLIC:: VisuRP
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Visualize the recordpoints as VTK structured grids. This routines prepares types that are used to store the points, planes
!> and lines and then calls the common routine to write those to the VTK file format.
!===================================================================================================================================
SUBROUTINE VisuRP()
! MODULES
USE MOD_Globals
USE MOD_Parameters
USE MOD_Output_Vars     ,ONLY: ProjectName
USE MOD_RPSet_Vars
USE MOD_VTKStructuredOutput
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(RPPoint)                :: RPPoints
TYPE(RPLine),ALLOCATABLE     :: RPLines(:)
TYPE(RPPlane),ALLOCATABLE    :: RPPlanes(:)
TYPE(RPBox),ALLOCATABLE      :: RPBoxes(:)
TYPE(tLine),POINTER          :: Line
TYPE(tPlane),POINTER         :: Plane
TYPE(tBox),POINTER           :: Box
CHARACTER(LEN=255)           :: GroupName
INTEGER                      :: iVar,i,j,k,iPlane,iBox,iLine
!===================================================================================================================================
IF(doVisuRP) THEN
  ! Prepare points structure
  RPPoints%nRPs = nPoints
  ALLOCATE(RPPoints%Coords(3,nPoints))
  DO iVar=1,3
    DO i=1,nPoints
      RPPoints%Coords(iVar,i)=Points(i)%RP%xF(iVar)
    END DO ! i
  END DO !iVar
  ! Prepare lines structure
  ALLOCATE(RPLines(nLines))
  DO iLine=1,nLines
    Line=>Lines(iLine)
    RPLines(iLine)%nRPs=Line%nRP
    GroupName=Groups(Line%GroupID)%Name
    RPLines(iLine)%name = TRIM(GroupName)//'_'//TRIM(Line%Name)
    ALLOCATE(RPLines(iLine)%Coords(3,Line%nRP))
    DO iVar=1,3
      DO i=1,Line%nRP
        RPLines(iLine)%Coords(iVar,i)=Line%RP_ptr(i)%RP%xF(iVar)
      END DO ! i
    END DO !iVar
  END DO
  ! Prepare plane structure
  ALLOCATE(RPPlanes(nPlanes))
  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    RPPlanes(iPlane)%nRPs=Plane%nRP
    GroupName=Groups(Plane%GroupID)%Name
    RPPlanes(iPlane)%name = TRIM(GroupName)//'_'//TRIM(Plane%Name)
    ALLOCATE(RPPlanes(iPlane)%Coords(3,Plane%nRP(1),Plane%nRP(2)))
    DO iVar=1,3
      DO j=1,Plane%nRP(2)
        DO i=1,Plane%nRP(1)
          RPPlanes(iPlane)%Coords(iVar,i,j)=Plane%RP_ptr(i,j)%RP%xF(iVar)
        END DO ! i
      END DO ! j
    END DO !iVar
  END DO
  ! Prepare box structure
  ALLOCATE(RPBoxes(nBoxes))
  DO iBox=1,nBoxes
    Box=>Boxes(iBox)
    RPBoxes(iBox)%nRPs=Box%nRP
    GroupName=Groups(Box%GroupID)%Name
    RPBoxes(iBox)%name = TRIM(GroupName)//'_'//TRIM(Box%Name)
    ALLOCATE(RPBoxes(iBox)%Coords(3,Box%nRP(1),Box%nRP(2),Box%nRP(3)))
    DO iVar=1,3
      DO k=1,Box%nRP(3)
        DO j=1,Box%nRP(2)
          DO i=1,Box%nRP(1)
            RPBoxes(iBox)%Coords(iVar,i,j,k)=Box%RP_ptr(i,j,k)%RP%xF(iVar)
          END DO ! i
        END DO ! j
      END DO ! k
    END DO !iVar
  END DO
  CALL WriteStructuredDataToVTK(ProjectName,nLines,nPlanes,nBoxes,RPPoints,RPLines,RPPlanes,RPBoxes,.FALSE.,0)
  WRITE(UNIT_StdOut,'(132("-"))')
END IF
END SUBROUTINE VisuRP

END MODULE MOD_VisuRP
