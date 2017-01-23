#ifdef WITHTECPLOT
#include "flexi.h"

MODULE MOD_Tecplot
!===================================================================================================================================
! Module for generic data output in Tecplot fromat
!
! WARNING: WriteDataToTecplot works only for POSTPROCESSING
!
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDataToTecplotBinary
  MODULE PROCEDURE WriteDataToTecplotBinary
END INTERFACE

INTERFACE WriteTimeAvgDataToTecplotBinary
  MODULE PROCEDURE WriteTimeAvgDataToTecplotBinary
END INTERFACE

INTERFACE WriteBLPropsToTecplotBinary
  MODULE PROCEDURE WriteBLPropsToTecplotBinary
END INTERFACE

INTERFACE WriteDataToTecplot
  MODULE PROCEDURE WriteDataToTecplot
END INTERFACE

PUBLIC::WriteDataToTecplotBinary
PUBLIC::WriteTimeAvgDataToTecplotBinary
PUBLIC::WriteBLPropsToTecplotBinary
PUBLIC::WriteDataToTecplot
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteDataToTecplotBinary(nSamples,nRP,nVal,VarNames,Time,Value,FileString)
!===================================================================================================================================
! Subroutine to write point data to Tecplot format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Parameters      ,ONLY:Line_LocalCoords,Plane_LocalCoords
USE MOD_Parameters      ,ONLY:OutputPlanes,OutputLines,OutputPoints
USE MOD_RPSet_Vars          ,ONLY:GroupNames 
USE MOD_RPSet_Vars          ,ONLY:OutputGroup 
USE MOD_RPSet_Vars          ,ONLY:nPoints,Points_IDlist,Points_GroupIDlist
USE MOD_RPSet_Vars          ,ONLY:nLines,Lines,tLine
USE MOD_RPSet_Vars          ,ONLY:nPlanes,Planes,tPlane
USE MOD_RPSet_Vars          ,ONLY:xF_RP
USE MOD_OutputRPVisu_Vars         ,ONLY:nCoords,CoordNames
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nSamples                      ! Number of Samples 
INTEGER,INTENT(IN)            :: nRP                           ! Number of RP to be visualized 
INTEGER,INTENT(IN)            :: nVal                          ! Number of nodal output variables
CHARACTER(LEN=255),INTENT(IN) :: VarNames(nVal)                ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Value(1:nVal,nRP,nSamples)    ! Statevector 
REAL,INTENT(IN)               :: Time(nSamples)                ! Time 
CHARACTER(LEN=*),INTENT(IN)   :: FileString                    ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iVar
INTEGER              :: TECINI112,TECZNE112,TECDAT112,TECEND112
INTEGER              :: iStat
INTEGER              :: offset
INTEGER,ALLOCATABLE  :: array_0(:),array_1(:),passiveVars(:)
CHARACTER(LEN=255)   :: Format_Title
CHARACTER(LEN=35)    :: VarString
INTEGER              :: iPoint,iLine,iSample,iPlane,i,j
INTEGER              :: GroupID
CHARACTER(LEN=255)   :: ZoneTitle
CHARACTER(LEN=255)   :: GroupName
REAL                 :: PointData(1:nVal,nSamples)
REAL                 :: Coord(nSamples)
TYPE(tLine),POINTER  :: Line
TYPE(tPlane),POINTER :: Plane
REAL,ALLOCATABLE     :: LineData(:,:,:)
REAL,ALLOCATABLE     :: LineCoord(:,:) 
REAL,ALLOCATABLE     :: PlaneData(:,:,:,:)
REAL,ALLOCATABLE     :: PlaneCoord(:,:,:) 
!===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE RP DATA TO TECPLOT BINARY FILE..."
  !assemble format strings
  Format_title(1:255)=' '
  IF(nVal.EQ.0) THEN
    WRITE(*,*) 'No Variables are given for visualization, no solution file written!!'
    RETURN
  END IF
  offset=1
  WRITE(VarString,'(A,A)')TRIM(CoordNames(1))
  Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
  Offset=Offset+LEN(TRIM(VarString))
  DO iVar=2,nCoords
    WRITE(VarString,'(A,A)')',',TRIM(CoordNames(iVar))
    Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
    Offset=Offset+LEN(TRIM(VarString))
  END DO
  DO iVar=1,nVal
    WRITE(VarString,'(A,A)')',',TRIM(VarNames(iVar))
    Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
    Offset=Offset+LEN(TRIM(VarString))
  END DO
  ALLOCATE(array_0(nVal+nCoords),array_1(nVal+nCoords),passiveVars(nVal+nCoords))
  iStat = TECINI112(                      &
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
IF(OutputPoints)THEN
  passiveVars(:)=0
  offset=0
  IF(Line_LocalCoords) THEN
  ! points dont have local line coord, so set var passive
    passiveVars(2)=1
    offset=2
  END IF
  IF(Plane_LocalCoords) THEN
  ! points dont have local line coord, so set var passive
    passiveVars(offset+1:offset+2)=1
  END IF
  DO iPoint=1,nPoints
    GroupID=Points_GroupIDlist(iPoint)
    IF(.NOT.OutputGroup(GroupID)) CYCLE
    GroupName=GroupNames(GroupID)
    ZoneTitle(1:255)=' '
    WRITE(ZoneTitle,'(A,A,I0.4)')TRIM(GroupName),'_Point',iPoint
  
    iStat = TECZNE112(          &
            ZoneTitle//char(0), &       ! ZoneTitle
            0,                  &       ! ZoneType
            nSamples          , &       ! IMxOrNumPts
            1                 , &       ! JMxOrNumElements
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
    !time  
      iStat = TECDAT112(                 &
              nSamples,                  &       ! N
              Time(:),                   &       ! Data
              1                           )      ! IsDouble
    !coordinates
    DO iVar=1,3
      DO iSample=1,nSamples
        Coord(iSample)=xF_RP(iVar,Points_IDlist(iPoint))
      END DO ! iSample
        iStat = TECDAT112(                  &
                nSamples,                   &       ! N
                Coord(:),                   &       ! Data
                1                           )       ! IsDouble
    END DO !iVar  
    !values
    PointData(:,:)=Value(:,Points_IDlist(iPoint),:)
    DO iVar=1,nVal
        iStat = TECDAT112(                 &
                nSamples,                  &       ! N
                PointData(iVar,:),         &       ! Data
                1                          )      ! IsDouble
    END DO       ! iVar
  END DO ! iPoint
END IF!OutputPoint

! Lines 
IF(OutputLines)THEN
  PassiveVars=0
  offset=0
  IF(Line_LocalCoords) offset=2
  IF(Plane_LocalCoords) THEN
    passiveVars(offset+1:offset+2)=1
  END IF
  DO iLine=1,nLines
    Line=>Lines(iLine)
    IF(.NOT.OutputGroup(Line%GroupID)) CYCLE
    GroupName=GroupNames(Line%GroupID)
    ZoneTitle(1:255)=' '
    WRITE(ZoneTitle,'(A,A,A)')TRIM(GroupName),'_',TRIM(Line%Name)
  
    iStat = TECZNE112(          &
            ZoneTitle//char(0), &       ! ZoneTitle
            0,                  &       ! ZoneType
            nSamples          , &       ! IMxOrNumPts
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
            passiveVars,        &       ! PassiveVarList
            array_1,            &       ! ValueLocation
            array_0,            &       ! ShareVarFromZone
            0                   )       ! ShareConnectivityFromZone
  
    ALLOCATE(LineData(1:nVal,nSamples,Line%nRP))
  
    !time  
    DO iPoint=1,Line%nRP
      LineData(1,:,iPoint)=Time(:)
    END DO ! iPoint
      iStat = TECDAT112(                 &
              nSamples*Line%nRP,         &       ! N
              LineData(1,:,:),           &       ! Data
              1                           )      ! IsDouble
  
    !coordinates
    ALLOCATE(LineCoord(nSamples,Line%nRP))
    ! local coord if required
    IF(Line_LocalCoords) THEN
      DO iSample=1,nSamples
        LineCoord(iSample,:)=Line%LocalCoord(:)             
      END DO ! iSample
      iStat = TECDAT112(                  &
              nSamples*Line%nRP,          &       ! N
              LineCoord(:,:),             &       ! Data
              1                           )       ! IsDouble
    END IF
    ! global xyz coordinates
    DO iVar=1,3
      DO iSample=1,nSamples
        DO iPoint=1,Line%nRP
        LineCoord(iSample,iPoint)=xF_RP(iVar,Line%IDlist(iPoint))
        END DO ! iPoint
      END DO ! iSample
      iStat = TECDAT112(                  &
              nSamples*Line%nRP,          &       ! N
              LineCoord(:,:),             &       ! Data
              1                           )       ! IsDouble
    END DO !iVar
    DEALLOCATE(LineCoord)
  
    !values
    DO iPoint=1,Line%nRP
      LineData(:,:,iPoint)=Value(:,Line%IDlist(iPoint),:)
    END DO ! iPoint
    DO iVar=1,nVal
        iStat = TECDAT112(                 &
                nSamples*Line%nRP,         &       ! N
                LineData(iVar,:,:),        &       ! Data
                1                          )      ! IsDouble
    END DO       ! iVar
    DEALLOCATE(LineData)
  END DO ! iLine
END IF !OutputLines

! Planes 
IF(OutputPlanes)THEN
  passiveVars(:)=0
  offset=0
  IF(Line_LocalCoords) THEN
  ! planes dont have local line coord, so set var passive
    passiveVars(2)=1
    offset=2
  END IF
  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    IF(.NOT.OutputGroup(Plane%GroupID)) CYCLE
    GroupName=GroupNames(Plane%GroupID)
    ZoneTitle(1:255)=' '
    WRITE(ZoneTitle,'(A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name)
    IF(Plane_localCoords.AND..NOT.ALLOCATED(Plane%LocalCoord)) THEN
      ! for planes without localcoords, set them passive
      passiveVars(offset+1:offset+2)=1
    ELSE
      passiveVars(offset+1:offset+2)=0
    END IF
  
    iStat = TECZNE112(          &
            ZoneTitle//char(0), &       ! ZoneTitle
            0,                  &       ! ZoneType
            nSamples          , &       ! IMxOrNumPts
            Plane%nRP(1)     ,  &       ! JMxOrNumElements
            Plane%nRP(2)     ,  &       ! KMxOrNumFaces
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
  
    ALLOCATE(PlaneData(1:nVal,nSamples,Plane%nRP(1),Plane%nRP(2)))
    !time  
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        PlaneData(1,:,i,j)=Time(:)
      END DO ! i
    END DO ! j
      iStat = TECDAT112(                           &
              nSamples*Plane%nRP(1)*Plane%nRP(2),  &       ! N
              PlaneData(1,:,:,:),                  &       ! Data
              1                                     )      ! IsDouble
  
    !coordinates
    ALLOCATE(PlaneCoord(nSamples,Plane%nRP(1),Plane%nRP(2)))
    ! local coord if required
    IF(ALLOCATED(Plane%LocalCoord)) THEN
      DO iVar=1,2
        DO iSample=1,nSamples
          PlaneCoord(iSample,:,:)=Plane%LocalCoord(iVar,:,:)
        END DO ! iSample
        iStat = TECDAT112(                           &
                nSamples*Plane%nRP(1)*Plane%nRP(2),  &       ! N
                PlaneCoord(:,:,:),                   &       ! Data
                1                                    )       ! IsDouble
      END DO !iVar
    END IF
    ! global xyz coordinates
    DO iVar=1,3
      DO iSample=1,nSamples
        DO j=1,Plane%nRP(2)
          DO i=1,Plane%nRP(1)
            PlaneCoord(iSample,i,j)=xF_RP(iVar,Plane%IDlist(i,j))
          END DO ! i
        END DO ! j
      END DO ! iSample
      iStat = TECDAT112(                           &
              nSamples*Plane%nRP(1)*Plane%nRP(2),  &       ! N
              PlaneCoord(:,:,:),                   &       ! Data
              1                                    )       ! IsDouble
    END DO !iVar
    DEALLOCATE(PlaneCoord)
  
    !values
     DO j=1,Plane%nRP(2)
       DO i=1,Plane%nRP(1)
         PlaneData(:,:,i,j)=Value(:,Plane%IDlist(i,j),:)
       END DO ! i
     END DO ! j
    DO iVar=1,nVal
        iStat = TECDAT112(                           &
                nSamples*Plane%nRP(1)*Plane%nRP(2),  &       ! N
                PlaneData(iVar,:,:,:),               &       ! Data
                1                                     )      ! IsDouble
    END DO       ! iVar
    DEALLOCATE(PlaneData)
  END DO ! iPlane
END IF!OutputPlanes

!close file
iStat = TECEND112()

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToTecplotBinary



SUBROUTINE WriteTimeAvgDataToTecplotBinary(nRP,nVal,VarNames,Value,FileString)
!===================================================================================================================================
! Subroutine to write point data to Tecplot format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Parameters      ,ONLY:Line_LocalCoords,Plane_LocalCoords
USE MOD_Parameters      ,ONLY:OutputPlanes,OutputLines,OutputPoints
USE MOD_RPSet_Vars          ,ONLY:GroupNames 
USE MOD_RPSet_Vars          ,ONLY:nPoints,Points_IDlist,Points_GroupIDlist
USE MOD_RPSet_Vars          ,ONLY:nLines,Lines,tLine
USE MOD_RPSet_Vars          ,ONLY:nPlanes,Planes,tPlane
USE MOD_RPSet_Vars          ,ONLY:OutputGroup
USE MOD_RPSet_Vars          ,ONLY:xF_RP
USE MOD_OutputRPVisu_Vars         ,ONLY:nCoords,CoordNames
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nRP                         ! Number of RP to be visualized 
INTEGER,INTENT(IN)            :: nVal                        ! Number of nodal output variables
CHARACTER(LEN=255),INTENT(IN) :: VarNames(nVal)              ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Value(1:nVal,nRP)           ! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString                  ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar
INTEGER             :: TECINI112,TECZNE112,TECDAT112,TECEND112
INTEGER             :: iStat
INTEGER             :: offset
INTEGER             :: nCoords_loc
INTEGER,ALLOCATABLE :: array_0(:),array_1(:),passiveVars(:)
CHARACTER(LEN=255)  :: Format_Title
CHARACTER(LEN=35)   :: VarString
INTEGER             :: iPoint,iLine,i,j,iPlane,GroupID
CHARACTER(LEN=255)  :: ZoneTitle
CHARACTER(LEN=255)  :: GroupName
TYPE(tLine),POINTER :: Line
TYPE(tPlane),POINTER :: Plane
REAL,ALLOCATABLE    :: LineData(:,:)
REAL,ALLOCATABLE    :: LineCoord(:) 
REAL,ALLOCATABLE    :: PlaneData(:,:,:)
REAL,ALLOCATABLE    :: PlaneCoord(:,:) 
!===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE RP TIME AVERAGE DATA TO TECPLOT BINARY FILE..."
  !assemble format strings
  Format_title(1:255)=' '
  IF(nVal.EQ.0) THEN
    WRITE(*,*) 'No Variables are given for visualization, no solution file written!!'
    RETURN
  END IF
  offset=1
    WRITE(VarString,'(A,A)'),TRIM(CoordNames(2))
    Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
    Offset=Offset+LEN(TRIM(VarString))
  nCoords_loc=nCoords-1 ! since we dont use the time 
  DO iVar=3,nCoords
    WRITE(VarString,'(A,A)')',',TRIM(CoordNames(iVar))
    Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
    Offset=Offset+LEN(TRIM(VarString))
  END DO
  DO iVar=1,nVal
    WRITE(VarString,'(A,A)')',',TRIM(VarNames(iVar))
    Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
    Offset=Offset+LEN(TRIM(VarString))
  END DO
  ALLOCATE(array_0(nVal+nCoords_loc),array_1(nVal+nCoords_loc),passiveVars(nVal+nCoords_loc))
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
IF(OutputPoints) THEN
  passiveVars(:)=0
  offset=0
  IF(Line_LocalCoords) THEN
  ! points dont have local line coord, so set var passive
    passiveVars(2)=1
    offset=2
  END IF
  IF(Plane_LocalCoords) THEN
  ! points dont have local line coord, so set var passive
    passiveVars(offset+1:offset+2)=1
  END IF
  DO iPoint=1,nPoints
    GroupID=Points_GroupIDlist(iPoint)
    IF(.NOT.OutputGroup(GroupID)) CYCLE
    GroupName=GroupNames(GroupID)
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
        iStat = TECDAT112(                        &
                1,                                &       ! N
                xF_RP(iVar,Points_IDlist(iPoint)),&       ! Data
                1                                 )       ! IsDouble
    END DO !iVar  
  
    !values
    DO iVar=1,nVal
        iStat = TECDAT112(                          &
                1,                                  &       ! N
                Value(iVar,Points_IDlist(iPoint)),  &       ! Data
                1                                   )       ! IsDouble
    END DO       ! iVar
  END DO ! iPoint
END IF! OutputPoints

! Lines 
IF(OutputLines) THEN
  PassiveVars=0
  offset=0
  IF(Line_LocalCoords) offset=2
  IF(Plane_LocalCoords) THEN
  ! points dont have local line coord, so set var passive
    passiveVars(offset+1:offset+2)=1
  END IF
  DO iLine=1,nLines
    Line=>Lines(iLine)
    IF(.NOT.OutputGroup(Line%GroupID)) CYCLE
    GroupName=GroupNames(Line%GroupID)
    ZoneTitle(1:255)=' '
    WRITE(ZoneTitle,'(A,A,A)')TRIM(GroupName),'_',TRIM(Line%Name)
  
    iStat = TECZNE112(          &
            ZoneTitle//char(0), &       ! ZoneTitle
            0,                  &       ! ZoneType
            1                 , &       ! IMxOrNumPts
            Line%nRP         , &       ! JMxOrNumElements
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
  
    ALLOCATE(LineData(1:nVal,Line%nRP))
    !coordinates
    ! local coord if required
    ALLOCATE(LineCoord(Line%nRP))
    IF(Line_LocalCoords) THEN
      LineCoord(:)=Line%LocalCoord(:)             
      iStat = TECDAT112(                &
              Line%nRP,                &       ! N
              LineCoord(:),             &       ! Data
              1                         )       ! IsDouble
    END IF
    DO iVar=1,3
      DO iPoint=1,Line%nRP
        LineCoord(iPoint)=xF_RP(iVar,Line%IDlist(iPoint))
      END DO ! iPoint
      iStat = TECDAT112(                  &
              Line%nRP,                  &       ! N
              LineCoord(:),               &       ! Data
              1                           )       ! IsDouble
    END DO !iVar  
    DEALLOCATE(LineCoord)
  
    !values
    DO iPoint=1,Line%nRP
      LineData(:,iPoint)=Value(:,Line%IDlist(iPoint))
    END DO ! iPoint
    DO iVar=1,nVal
        iStat = TECDAT112(                 &
                Line%nRP,                 &       ! N
                LineData(iVar,:),          &       ! Data
                1                          )       ! IsDouble
    END DO       ! iVar
    DEALLOCATE(LineData)
  END DO ! iLine
END IF !OutputLines

! Planes 
IF(OutputPlanes) THEN
  passiveVars(:)=0
  IF(Line_LocalCoords) THEN
  ! planes dont have local line coord, so set var passive
    passiveVars(1)=1
  END IF
  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    IF(.NOT.OutputGroup(Plane%GroupID)) CYCLE
    GroupName=GroupNames(Plane%GroupID)
    ZoneTitle(1:255)=' '
    WRITE(ZoneTitle,'(A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name)
    IF(Plane_localCoords.AND..NOT.ALLOCATED(Plane%LocalCoord)) THEN
      ! for planes without localcoords, set them passive
      passiveVars(offset+1:offset+2)=1
    ELSE
      passiveVars(offset+1:offset+2)=0
    END IF
  
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
  
    !coordinates
    ALLOCATE(PlaneCoord(Plane%nRP(1),Plane%nRP(2)))
    ! local coord if required
    IF(ALLOCATED(Plane%LocalCoord)) THEN
      DO iVar=1,2
        PlaneCoord(:,:)=Plane%LocalCoord(iVar,:,:)
        iStat = TECDAT112(                           &
                Plane%nRP(1)*Plane%nRP(2),           &       ! N
                PlaneCoord(:,:),                     &       ! Data
                1                                    )       ! IsDouble
      END DO !iVar
    END IF
    DO iVar=1,3
      DO j=1,Plane%nRP(2)
        DO i=1,Plane%nRP(1)
          PlaneCoord(i,j)=xF_RP(iVar,Plane%IDlist(i,j))
        END DO ! i
      END DO ! j
      iStat = TECDAT112(                  &
              Plane%nRP(1)*Plane%nRP(2),  &       ! N
              PlaneCoord(:,:),            &       ! Data
              1                           )       ! IsDouble
    END DO !iVar  
    DEALLOCATE(PlaneCoord)
  
    !values
    ALLOCATE(PlaneData(1:nVal,Plane%nRP(1),Plane%nRP(2)))
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        PlaneData(:,i,j)=Value(:,Plane%IDlist(i,j))
      END DO ! i
    END DO ! j
    DO iVar=1,nVal
        iStat = TECDAT112(                 &
                Plane%nRP(1)*Plane%nRP(2), &       ! N
                PlaneData(iVar,:,:),       &       ! Data
                1                          )       ! IsDouble
    END DO       ! iVar
    DEALLOCATE(PlaneData)
    
  END DO ! iPlane
END IF!OutputPlanes

!close file
iStat = TECEND112()

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteTimeAvgDataToTecplotBinary



SUBROUTINE WriteBLPropsToTecplotBinary(FileString)
!===================================================================================================================================
! Subroutine to write BL properties (line data associated to planes) to tecplot
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RPSet_Vars          ,ONLY:nPlanes,Planes,tPlane
USE MOD_RPSet_Vars          ,ONLY:OutputGroup,GroupNames
USE MOD_RPSet_Vars          ,ONLY:xF_RP
USE MOD_OutputRPVisu_Vars         ,ONLY:nCoords,CoordNames
USE MOD_Equation_Vars       ,ONLY:nBLProps,VarNames_BLProps
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileString                  ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar
INTEGER             :: TECINI112,TECZNE112,TECDAT112,TECEND112
INTEGER             :: iStat
INTEGER             :: offset
INTEGER             :: nCoords_loc
INTEGER,ALLOCATABLE :: array_0(:),array_1(:),passiveVars(:)
CHARACTER(LEN=255)  :: Format_Title
CHARACTER(LEN=255)  :: VarNames(nBLProps)
CHARACTER(LEN=35)   :: VarString
INTEGER             :: iPoint,iLine,i,j,iPlane,GroupID
CHARACTER(LEN=255)  :: ZoneTitle
CHARACTER(LEN=255)  :: GroupName
TYPE(tPlane),POINTER:: Plane
REAL,ALLOCATABLE    :: LineCoord(:) 
!===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE BOUNDARY LAYER PROPERTY DATA TO TECPLOT BINARY FILE..."
  !assemble format strings
  Format_title(1:255)=' '
  offset=1
    WRITE(VarString,'(A,A)'),TRIM(CoordNames(2))
    Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
    Offset=Offset+LEN(TRIM(VarString))
  nCoords_loc=nCoords-1 ! since we dont use the time 
  DO iVar=3,nCoords
    WRITE(VarString,'(A,A)')',',TRIM(CoordNames(iVar))
    Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
    Offset=Offset+LEN(TRIM(VarString))
  END DO
  DO iVar=1,nBLProps
    WRITE(VarString,'(A,A)')',',TRIM(VarNames_BLProps(iVar))
    Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
    Offset=Offset+LEN(TRIM(VarString))
  END DO
  ALLOCATE(array_0(nBLprops+nCoords_loc),array_1(nBLprops+nCoords_loc),passiveVars(nBLprops+nCoords_loc))
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
  passiveVars=0

  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    IF(Plane%Type.EQ.2) THEN ! BLPlane
      IF(.NOT.OutputGroup(Plane%GroupID)) CYCLE
      GroupName=GroupNames(Plane%GroupID)
      ZoneTitle(1:255)=' '
      WRITE(ZoneTitle,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name),'_BLProps'
    
      iStat = TECZNE112(          &
              ZoneTitle//char(0), &       ! ZoneTitle
              0,                  &       ! ZoneType
              1                 , &       ! IMxOrNumPts
              Plane%nRP(1)      , &       ! JMxOrNumElements
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
      ALLOCATE(LineCoord(Plane%nRP(1)))
      ! local coord always active
      DO iVar=1,2
        LineCoord(:)=Plane%LocalCoord(iVar,:,1)
        iStat = TECDAT112(                &
                Plane%nRP(1),             &       ! N
                LineCoord(:),             &       ! Data
                1                         )       ! IsDouble
      END DO !iVar
      DO iVar=1,3
        DO iPoint=1,Plane%nRP(1)
          LineCoord(iPoint)=xF_RP(iVar,Plane%IDlist(iPoint,1))
        END DO ! iPoint
        iStat = TECDAT112(                  &
                Plane%nRP(1),               &       ! N
                LineCoord(:),               &       ! Data
                1                           )       ! IsDouble
      END DO !iVar  
      DEALLOCATE(LineCoord)
    
      !values
      DO iVar=1,nBLProps
          iStat = TECDAT112(                &
                  Plane%nRP(1),             &       ! N
                  Plane%BLProps(iVar,:),    &       ! Data
                  1                          )       ! IsDouble
      END DO       ! iVar
    END IF!(Plane%Type.EQ.2.AND.Plane_doBLProps) BLPlane
  END DO ! iPlane
!close file
iStat = TECEND112()

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteBLPropsToTecplotBinary




SUBROUTINE WriteDataToTecplot(nSamples,nVal,VarNames,Value,FileString)
!===================================================================================================================================
! Subroutine to write point data to Tecplot format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RPSet_Vars          ,ONLY: nRP_global
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nSamples                ! Number of Samples 
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
CHARACTER(LEN=255),INTENT(IN) :: VarNames(nVal)          ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Value(1:nVal,nRP_global,nSamples)  ! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString              ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER            :: i,j,k,iVal,Offset
!CHARACTER(LEN=255) :: Format_nVal
!CHARACTER(LEN=255) :: Format_Title
!CHARACTER(LEN=35)  :: VarString
!INTEGER            :: iElem
!INTEGER            :: NodeIDElem,nPlot_p1_2,nPlot_p1_3
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE DATA TO TECPLOT ASCII FILE NOT YET IMPLEMENTED!"
!WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE DATA TO TECPLOT ASCII FILE..."
!!assemble format strings
!WRITE(Format_nVal,'(A1,I2,A17)')'(',nVal,'(1X,E21.10))'
!Format_Title(1:51)='VARIABLES="CoordinateX","CoordinateY","CoordinateZ"'
!Offset = 52
!DO iVal=1,nVal
!  WRITE(VarString,'(A2,A,A1)')',"',TRIM(VarNames(iVal)),'"'
!  Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
!  Offset=Offset+LEN(TRIM(VarString))
!END DO
!
!!connected 3D FEM data
!NPlot_p1_2=(NPlot+1)**2
!NPlot_p1_3=(NPlot+1)**3
!
!OPEN(44,FILE=TRIM(FileString),Status="REPLACE")
!WRITE(44,'(A)')Format_Title(1:Offset-1)
!WRITE(44,'(A,I8,A,I8,A)')'ZONE T="",DATAPACKING=POINT,&
! NODES=',nElems*NPlot_p1_3,', ELEMENTS=',nElems*(NPlot)**3,',ZONETYPE=FEBRICK'
!DO iElem=1,nElems
!  DO k=0,NPlot
!    DO j=0,NPlot
!      DO i=0,NPlot
!        WRITE(44,Format_nVal)&
!        Coord(:,i,j,k,iElem),Value(:,i,j,k,iElem)
!      END DO 
!    END DO 
!  END DO 
!END DO 
!!element connectivity
!NodeIDElem=0
!DO iElem=1,nElems
!  DO k=1,NPlot
!    DO j=1,NPlot
!      DO i=1,NPlot
!        !visuHexaElem  
!        WRITE(44,'(8(I8,1X))')&
!          NodeIDElem+i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2,      & !P1(CGNS=tecplot standard)
!          NodeIDElem+i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2,      & !P2
!          NodeIDElem+i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2,      & !P3     
!          NodeIDElem+i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2,      & !P4
!          NodeIDElem+i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2,      & !P5
!          NodeIDElem+i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2,      & !P6
!          NodeIDElem+i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2,      & !P7     
!          NodeIDElem+i+   j   *(NPlot+1)+ k   *NPlot_p1_2         !P8
!      END DO 
!    END DO 
!  END DO 
!  NodeIDElem=NodeIDElem+NPlot_p1_3
!END DO 
!
!CLOSE(44)
!WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToTecplot

END MODULE MOD_Tecplot
#endif WITHTECPLOT
