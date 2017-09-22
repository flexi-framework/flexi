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
#ifdef HASTECPLOT
USE MOD_RPSet_Vars      ,ONLY: nRP_global
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef HASTECPLOT
CHARACTER(LEN=255)            :: FileName,strOutputFile
#endif
!===================================================================================================================================
IF(doVisuRP) THEN
#ifdef HASTECPLOT
  nCoords=4
  ALLOCATE(CoordNames(4)) ! Time,x,y,z 
  CoordNames(1)='Time'
  CoordNames(2)='CoordinateX'
  CoordNames(3)='CoordinateY'
  CoordNames(4)='CoordinateZ'
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RPvisu'
  WRITE(UNIT_StdOut,'(132("-"))')
  strOutputFile=TRIM(FileName)//'.plt'
  WRITE(UNIT_stdOut,'(A,A)')' WRITING RP POSITIONS TO ',strOutputFile
  CALL WriteDataToTecplotBinary(nRP_global,strOutputFile)
#else
  CALL WriteStructuredDataToVTK(ProjectName)
#endif
  WRITE(UNIT_StdOut,'(132("-"))')
END IF
END SUBROUTINE VisuRP


!===================================================================================================================================
!> Subroutine to write 2D or 3D point data to VTK format
!===================================================================================================================================
SUBROUTINE WriteStructuredDataToVTK(ProjectName)
! MODULES
USE MOD_Globals
USE MOD_Parameters
USE MOD_RPSet_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: Projectname  !< Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: ivtk=44
INTEGER                     :: nBytes,Offset
REAL(KIND=4)                :: FLOATdummy
CHARACTER(LEN=35)           :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=200)          :: Buffer2
CHARACTER(LEN=255)          :: ZoneTitle
CHARACTER(LEN=255)          :: FileName
CHARACTER(LEN=255)          :: GroupName
CHARACTER(LEN=1)            :: lf
REAL,ALLOCATABLE            :: PlaneCoord(:,:,:),LineCoord(:,:),PointCoord(:,:)
TYPE(tPlane),POINTER        :: Plane
TYPE(tLine),POINTER         :: Line
INTEGER                     :: iVar,i,j,iPlane,iLine,nSets,iSet
CHARACTER(LEN=255),ALLOCATABLE :: ZoneNames(:),FileNamesVTS(:)
!===================================================================================================================================

WRITE(UNIT_stdOut,'(A,I1,A)')" WRITE Structured Plane to VTK ... "

nSets= MERGE( 1, 0, nPoints.GT.0) +nLines+nPlanes
ALLOCATE(FileNamesVTS(nSets))
ALLOCATE(ZoneNames(nSets))
iSet=0

! Points
IF(nPoints.GT.0) THEN
  FileName=TRIM(ProjectName)//'_Points.vts'
  iSet=iSet+1
  FileNamesVTS(iSet)=FileName
  ZoneNames(iSet)='Points'
  WRITE(UNIT_stdOut,'(A,A)')' WRITING POINT RP POSITIONS TO ',FileName

  ! write header of VTK file
  ! Line feed character
  lf = char(10)

  ! Write file
  OPEN(UNIT=ivtk,FILE=TRIM(FileName),ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  Buffer ='  <StructuredGrid WholeExtent="'
  Buffer2='    <Piece Extent="'
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') nPoints-1
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') 0 
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') 0
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) // '">'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) // '">'//lf;WRITE(ivtk) TRIM(Buffer2)

  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </StructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)

  ALLOCATE(PointCoord(1:3,nPoints))
  !coordinates
  DO iVar=1,3
    DO i=1,nPoints
      PointCoord(iVar,i)=Points(i)%RP%xF(iVar)
    END DO ! i
  END DO !iVar  

  nBytes = nPoints*SIZEOF_F(FLOATdummy) * 3
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL(PointCoord(:,:),4)
  DEALLOCATE(PointCoord)


  ! Footer
  lf = char(10)
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
END IF



! Lines
DO iLine=1,nLines
  Line=>Lines(iLine)
  ZoneTitle(1:255)=' '
  GroupName=Groups(Line%GroupID)%Name
  WRITE(ZoneTitle,'(A,A,A)')TRIM(GroupName),'_',TRIM(Line%Name)
  FileName=TRIM(ProjectName)//'_'//TRIM(ZoneTitle)//'.vts'
  iSet=iSet+1
  FileNamesVTS(iSet)=FileName
  ZoneNames(iSet)=TRIM(ZoneTitle)
  WRITE(UNIT_stdOut,'(A,A)')' WRITING LINE RP POSITIONS TO ',FileName

! write header of VTK file
  ! Line feed character
  lf = char(10)

  ! Write file
  OPEN(UNIT=ivtk,FILE=TRIM(FileName),ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  Buffer ='  <StructuredGrid WholeExtent="'
  Buffer2='    <Piece Extent="'
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') Line%nRP-1
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') 0 
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') 0
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) // '">'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) // '">'//lf;WRITE(ivtk) TRIM(Buffer2)

  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </StructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)

  ALLOCATE(LineCoord(1:3,Line%nRP))
  !coordinates
  DO iVar=1,3
    DO i=1,Line%nRP
      LineCoord(iVar,i)=Line%RP_ptr(i)%RP%xF(iVar) 
    END DO ! i
  END DO !iVar  

  nBytes = Line%nRP*SIZEOF_F(FLOATdummy) * 3
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL(LineCoord(:,:),4)
  DEALLOCATE(LineCoord)

  ! Footer
  lf = char(10)
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
END DO


! Planes
DO iPlane=1, nPlanes
  Plane=>Planes(iPlane)
  ZoneTitle(1:255)=' '
  GroupName=Groups(Plane%GroupID)%Name
  WRITE(ZoneTitle,'(A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name)
  FileName=TRIM(ProjectName)//'_'//TRIM(ZoneTitle)//'.vts'
  iSet=iSet+1
  FileNamesVTS(iSet)=FileName
  ZoneNames(iSet)=TRIM(ZoneTitle)
  WRITE(UNIT_stdOut,'(A,A)')' WRITING PLANE RP POSITIONS TO ',FileName

  ! write header of VTK file
  ! Line feed character
  lf = char(10)

  ! Write file
  OPEN(UNIT=ivtk,FILE=TRIM(FileName),ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  Buffer ='  <StructuredGrid WholeExtent="'
  Buffer2='    <Piece Extent="'
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') Plane%nRP(1)-1
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') Plane%nRP(2)-1
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) 
  WRITE(TempStr1,'(I16)') 0
  WRITE(TempStr2,'(I16)') 0
  Buffer =TRIM(Buffer)  // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) // '">'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer2=TRIM(Buffer2) // ' ' // TRIM(ADJUSTL(TempStr1)) // ' ' // TRIM(ADJUSTL(TempStr2)) // '">'//lf;WRITE(ivtk) TRIM(Buffer2)

  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </StructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)


  ALLOCATE(PlaneCoord(1:3,Plane%nRP(1),Plane%nRP(2)))
  !coordinates
  DO iVar=1,3
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        PlaneCoord(iVar,i,j)=Plane%RP_ptr(i,j)%RP%xF(iVar) 
      END DO ! i
    END DO ! j
  END DO !iVar  

  nBytes = Plane%nRP(1)*Plane%nRP(2)*SIZEOF_F(FLOATdummy) * 3
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL(PlaneCoord(:,:,:),4)
  DEALLOCATE(PlaneCoord)


  ! Footer
  lf = char(10)
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
END DO

CALL WriteVTKMultiBlockDataSetRP(ProjectName,nSets,FileNamesVTS,ZoneNames)
DEALLOCATE(FileNamesVTS,ZoneNames)

SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteStructuredDataToVTK



!===================================================================================================================================
!> Links structured VTK data files together
!===================================================================================================================================
SUBROUTINE WriteVTKMultiBlockDataSetRP(ProjectName,nSets,FileNamesVTS,ZoneNames)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: ProjectName          !< Projectname
INTEGER, INTENT(IN)         :: nSets                !< Number of VTS files to link
CHARACTER(LEN=*),INTENT(IN) :: FileNamesVTS(nSets)  !< Filenames of structured datasets 
CHARACTER(LEN=*),INTENT(IN) :: ZoneNames(nSets)     !< Zone names of structured datasets 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ivtk=44
INTEGER            :: iSet
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
CHARACTER(LEN=35)  :: TempStr1
CHARACTER(LEN=255) :: FileStringOut 
!===================================================================================================================================
FileStringOut=TRIM(ProjectName)//'_RPVisu.vtm'
! write multiblock file
OPEN(UNIT=ivtk,FILE=TRIM(FileStringOut),ACCESS='STREAM')
! Line feed character
lf = char(10)
Buffer='<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf
WRITE(ivtk) TRIM(BUFFER)
Buffer='  <vtkMultiBlockDataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)

DO iSet=1,nSets
  WRITE(TempStr1,'(I16)') iSet
  Buffer='    <DataSet index="' // TRIM(ADJUSTL(TempStr1)) // '" name="' // TRIM(ZoneNames(iSet)) // '" file="'&
            //TRIM(FileNamesVTS(iSet))// '">'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    </DataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
END DO

Buffer='  </vtkMultiBlockDataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(BUFFER)
CLOSE(ivtk)
END SUBROUTINE WriteVTKMultiBlockDataSetRP



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
