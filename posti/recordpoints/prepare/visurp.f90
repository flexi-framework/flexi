#include "flexi.h"

!===================================================================================================================================
!> Module conaining Record point visualization routines
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
!> Call of software specific output routines
!===================================================================================================================================
SUBROUTINE VisuRP()
! MODULES
USE MOD_Globals
USE MOD_Parameters
USE MOD_Output_Vars     ,ONLY: ProjectName
USE MOD_RPSet_Vars      ,ONLY: nRP_global
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: FileName
CHARACTER(LEN=255)            :: strOutputFile
!===================================================================================================================================
IF(doVisuRP) THEN
  nCoords=4
  ALLOCATE(CoordNames(4)) ! Time,x,y,z 
  CoordNames(1)='Time'
  CoordNames(2)='CoordinateX'
  CoordNames(3)='CoordinateY'
  CoordNames(4)='CoordinateZ'
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RPvisu'
  strOutputFile=TRIM(FileName)//'.plt'

  WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING RP POSITIONS TO ',strOutputFile
#ifdef HASTECPLOT
  CALL WriteDataToTecplotBinary(nRP_global,strOutputFile)
#else
  STOP 'TECPLOT visualization is not available.'
#endif
  WRITE(UNIT_StdOut,'(132("-"))')
END IF
END SUBROUTINE VisuRP


#ifdef WITHTECPLOT
!===================================================================================================================================
!> Subroutine to write record point data to Tecplot format
!===================================================================================================================================
SUBROUTINE WriteRPDataToTecplotBinary(nRP,FileString)
! MODULES
USE MOD_Globals
USE MOD_Parameters
USE MOD_RPSet_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nRP        !< Number of RP to be visualized 
CHARACTER(LEN=*),INTENT(IN)   :: FileString !< Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar
INTEGER             :: TECINI112,TECZNE112,TECDAT112,TECEND112
INTEGER             :: iStat
INTEGER             :: offset
INTEGER,ALLOCATABLE :: array_0(:),array_1(:),passiveVars(:)
CHARACTER(LEN=255)  :: Format_Title
CHARACTER(LEN=35)   :: VarString
INTEGER             :: iPoint,iLine,i,j,iPlane,GroupID
CHARACTER(LEN=255)  :: ZoneTitle
CHARACTER(LEN=255)  :: GroupName
TYPE(tLine),POINTER :: Line
TYPE(tPlane),POINTER:: Plane
REAL,ALLOCATABLE    :: LineData(:,:)
REAL,ALLOCATABLE    :: LineCoord(:) 
REAL,ALLOCATABLE    :: PlaneData(:,:,:)
REAL,ALLOCATABLE    :: PlaneCoord(:,:) 
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE RP TIME AVERAGE DATA TO TECPLOT BINARY FILE..."
!assemble format strings
Format_title(1:255)=' '
offset=1
WRITE(VarString,'(A,A)'),TRIM(CoordNames(2))
Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
Offset=Offset+LEN(TRIM(VarString))

DO iVar=3,nCoords
  WRITE(VarString,'(A,A)')',',TRIM(CoordNames(iVar))
  Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
  Offset=Offset+LEN(TRIM(VarString))
END DO

ALLOCATE(array_0(nCoords),array_1(nCoords),passiveVars(nCoords))
iStat = TECINI112(                  &
        ""//char(0),                    &       ! Title
        TRIM(Format_Title)//char(0),    &       ! Variables
        TRIM(FileString)//char(0),      &       ! FName
        "."//char(0),                   &       ! ScratchDir
        0,                              &       ! FileType
        0,                              &       ! Debug
        1                               )       ! VIsDouble
array_0 = 0
array_1 = 1

! Points 
passiveVars(:)=0
DO iPoint=1,nPoints
  GroupID=Points(iPoint)%GroupID
  GroupName=Groups(GroupID)%Name
  ZoneTitle(1:255)=' '
  WRITE(ZoneTitle,'(A,A,I0.4)')TRIM(GroupName),'_Point',iPoint

  iStat = TECZNE112(          &
          ZoneTitle//char(0), &       ! ZoneTitle
          0,                  &       ! ZoneType
          1,                  &       ! IMxOrNumPts
          1,                  &       ! JMxOrNumElements
          1,                  &       ! KMxOrNumFaces
          0,                  &       ! ICellMax
          0,                  &       ! JCellMax
          0,                  &       ! KCellMax
          0,                  &       ! SolutionTime
          0,                  &       ! StrandID
          0,                  &       ! ParentZone
          1,                  &       ! IsBlock
          0,                  &       ! NumFaceConnections
          0,                  &       ! FaceNeighborMode
          0,                  &       ! TotalNumFaceNodes
          0,                  &       ! NumConnectedBoundaryFaces
          0,                  &       ! TotalNumBoundaryConnections
          passiveVars,        &       ! PassiveVarList
          array_1,            &       ! ValueLocation
          array_0,            &       ! ShareVarFromZone
          0                   )       ! ShareConnectivityFromZone
  !coordinates
  DO iVar=1,3
    iStat = TECDAT112(                  &
            1,                          & ! N
            Points(iPoint)%RP%xF(iVar), & ! Data
            1                           ) ! IsDouble
  END DO !iVar  

END DO ! iPoint

! Lines 
DO iLine=1,nLines
  Line=>Lines(iLine)
  GroupName=Groups(Line%GroupID)%Name
  ZoneTitle(1:255)=' '
  WRITE(ZoneTitle,'(A,A,A)')TRIM(GroupName),'_',TRIM(Line%Name)

  iStat = TECZNE112(          &
          ZoneTitle//char(0), &       ! ZoneTitle
          0,                  &       ! ZoneType
          1                 , &       ! IMxOrNumPts
          Line%nRP          , &       ! JMxOrNumElements
          1,                  &       ! KMxOrNumFaces
          0,                  &       ! ICellMax
          0,                  &       ! JCellMax
          0,                  &       ! KCellMax
          0,                  &       ! SolutionTime
          0,                  &       ! StrandID
          0,                  &       ! ParentZone
          1,                  &       ! IsBlock
          0,                  &       ! NumFaceConnections
          0,                  &       ! FaceNeighborMode
          0,                  &       ! TotalNumFaceNodes
          0,                  &       ! NumConnectedBoundaryFaces
          0,                  &       ! TotalNumBoundaryConnections
          array_0,            &       ! PassiveVarList
          array_1,            &       ! ValueLocation
          array_0,            &       ! ShareVarFromZone
          0                   )       ! ShareConnectivityFromZone

  !coordinates
  ALLOCATE(LineCoord(Line%nRP))
  DO iVar=1,3
    DO iPoint=1,Line%nRP
      LineCoord(iPoint)=Line%RP_ptr(iPoint)%RP%xF(iVar)
    END DO ! iPoint
    iStat = TECDAT112(                  &
            Line%nRP,                   &       ! N
            LineCoord(:),               &       ! Data
            1                           )       ! IsDouble
  END DO !iVar  
  DEALLOCATE(LineCoord)
END DO ! iLine

! Planes 
passiveVars(:)=0
DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  GroupName=Groups(Plane%GroupID)%Name
  ZoneTitle(1:255)=' '
  WRITE(ZoneTitle,'(A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name)

  iStat = TECZNE112(          &
          ZoneTitle//char(0), &       ! ZoneTitle
          0,                  &       ! ZoneType
          1                 , &       ! IMxOrNumPts
          Plane%nRP(1)      , &       ! JMxOrNumElements
          Plane%nRP(2)      , &       ! KMxOrNumFaces
          0,                  &       ! ICellMax
          0,                  &       ! JCellMax
          0,                  &       ! KCellMax
          0,                  &       ! SolutionTime
          0,                  &       ! StrandID
          0,                  &       ! ParentZone
          1,                  &       ! IsBlock
          0,                  &       ! NumFaceConnections
          0,                  &       ! FaceNeighborMode
          0,                  &       ! TotalNumFaceNodes
          0,                  &       ! NumConnectedBoundaryFaces
          0,                  &       ! TotalNumBoundaryConnections
          passiveVars,        &       ! PassiveVarList
          array_1,            &       ! ValueLocation
          array_0,            &       ! ShareVarFromZone
          0                   )       ! ShareConnectivityFromZone

  ALLOCATE(PlaneCoord(Plane%nRP(1),Plane%nRP(2)))
  !coordinates
  DO iVar=1,3
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        PlaneCoord(i,j)=Plane%RP_ptr(i,j)%RP%xF(iVar) 
      END DO ! i
    END DO ! j
    iStat = TECDAT112(                  &
            Plane%nRP(1)*Plane%nRP(2),  &       ! N
            PlaneCoord(:,:),            &       ! Data
            1                           )       ! IsDouble
  END DO !iVar  
  DEALLOCATE(PlaneCoord)
END DO ! iPlane

!close file
iStat = TECEND112()

SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteRPDataToTecplotBinary
#endif

END MODULE MOD_VisuRP
