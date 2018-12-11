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
TYPE(tPlane),POINTER         :: Plane
TYPE(tLine),POINTER          :: Line
CHARACTER(LEN=255)           :: GroupName
INTEGER                      :: iVar,i,j,iPlane,iLine,nSets,iSet
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
  CALL WriteStructuredDataToVTK(ProjectName,nLines,nPlanes,RPPoints,RPLines,RPPlanes,withData=.FALSE.)
  WRITE(UNIT_StdOut,'(132("-"))')
END IF
END SUBROUTINE VisuRP

END MODULE MOD_VisuRP
