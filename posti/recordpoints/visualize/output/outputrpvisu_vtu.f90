#include "flexi.h"

!===================================================================================================================================
!> \brief Module for visualize recordpoints data output in ParaView format
!>
!> Since the vtk file system does not allow for more than one time step in a file, we will write a seperate file for each timestep.
!===================================================================================================================================
MODULE MOD_OutputRPVisu_VTK
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDataToVTK
  MODULE PROCEDURE WriteDataToVTK
END INTERFACE

#ifdef WITHBLPROPS
INTERFACE WriteBLPropsToVTK
  MODULE PROCEDURE WriteBLPropsToVTK
END INTERFACE
PUBLIC::WriteBLPropsToVTK
#endif

PUBLIC::WriteDataToVTK
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write point data to VTK file 
!===================================================================================================================================
SUBROUTINE WriteDataToVTK(nSamples,nRP,nVal,VarNames,Time,Value)
! MODULES
USE MOD_Globals
USE MOD_ParametersVisu      ,ONLY:ProjectName
USE MOD_ParametersVisu      ,ONLY:Line_LocalCoords,Plane_LocalCoords
USE MOD_ParametersVisu      ,ONLY:OutputPlanes,OutputLines,OutputPoints
USE MOD_RPSetVisuVisu_Vars  ,ONLY:GroupNames 
USE MOD_RPSetVisuVisu_Vars  ,ONLY:OutputGroup 
USE MOD_RPSetVisuVisu_Vars  ,ONLY:nPoints,Points_IDlist,Points_GroupIDlist
USE MOD_RPSetVisuVisu_Vars  ,ONLY:nLines,Lines,tLine
USE MOD_RPSetVisuVisu_Vars  ,ONLY:nPlanes,Planes,tPlane
USE MOD_RPSetVisuVisu_Vars  ,ONLY:xF_RP
USE MOD_VTKStructuredOutput
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nSamples                      !< Number of Samples 
INTEGER,INTENT(IN)            :: nRP                           !< Number of RP to be visualized 
INTEGER,INTENT(IN)            :: nVal                          !< Number of nodal output variables
CHARACTER(LEN=255),INTENT(IN) :: VarNames(nVal)                !< Names of all variables that will be written out
REAL,INTENT(IN)               :: Value(1:nVal,nRP,nSamples)    !< Statevector 
REAL,INTENT(IN)               :: Time(nSamples)                !< Time 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iSample
INTEGER                   :: iPoint,iLine,iPlane,i,j
INTEGER                   :: GroupID
CHARACTER(LEN=255)        :: ZoneTitle
CHARACTER(LEN=255)        :: GroupName
TYPE(tLine),POINTER       :: Line
TYPE(tPlane),POINTER      :: Plane
REAL,ALLOCATABLE          :: PlaneData(:,:,:)
REAL,ALLOCATABLE          :: PlaneCoord(:,:,:) 
TYPE(RPPoint)             :: RPPoints
TYPE(RPLine),ALLOCATABLE  :: RPLines(:)
TYPE(RPPlane),ALLOCATABLE :: RPPlanes(:)
INTEGER                   :: nPointsOutput,nLinesOutput,nPlanesOutput
INTEGER                   :: iPointsOutput,iLinesOutput,iPlanesOutput
!===================================================================================================================================
! Count the number of points, lines and planes for output. Allocate the data types used for the output.
! Points
IF (OutputPoints) THEN
  nPointsOutput = 0
  DO iPoint=1,nPoints
    GroupID=Points_GroupIDlist(iPoint)
    IF(.NOT.OutputGroup(GroupID)) CYCLE
    nPointsOutput = nPointsOutput + 1
  END DO ! iPoint
  RPPoints%nRPs = nPointsOutput
  IF (nPointsOutput.GT.0) THEN
    ALLOCATE(RPPoints%Coords(3   ,nPointsOutput))
    ALLOCATE(RPPoints%Val(   nVal,nPointsOutput))
  END IF
ELSE
  RPPoints%nRPs = 0
END IF


! Lines
IF (OutputLines) THEN
  nLinesOutput = 0
  DO iLine=1,nLines
    Line=>Lines(iLine)
    IF(.NOT.OutputGroup(Line%GroupID)) CYCLE
    nLinesOutput = nLinesOutput + 1
  END DO ! iLine
  IF (nLinesOutput.GT.0) THEN
    ALLOCATE(RPLines(nLinesOutput))
  END IF
  iLinesOutput = 0
  DO iLine=1,nLines
    Line=>Lines(iLine)
    IF(.NOT.OutputGroup(Line%GroupID)) CYCLE
    iLinesOutput = iLinesOutput + 1
    RPLines(iLinesOutput)%nRPs = Line%nRP
    ALLOCATE(RPLines(iLinesOutput)%Coords(3   ,Line%nRP))
    ALLOCATE(RPLines(iLinesOutput)%Val(   nVal,Line%nRP))
    IF (Line_LocalCoords) RPLines(iLinesOutput)%Coords(2:3,:) = 0.
  END DO ! iLine
ELSE
  nLinesOutput = 0
END IF

! Planes
IF (OutputPlanes) THEN
  nPlanesOutput = 0
  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    IF(.NOT.OutputGroup(Plane%GroupID)) CYCLE
    nPlanesOutput = nPlanesOutput + 1
  END DO ! iPlane
  IF (nPlanesOutput.GT.0) THEN
    ALLOCATE(RPPlanes(nPlanesOutput))
  END IF
  iPlanesOutput = 0
  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    IF(.NOT.OutputGroup(Plane%GroupID)) CYCLE
    iPlanesOutput = iPlanesOutput + 1
    RPPlanes(iPlanesOutput)%nRPs = Plane%nRP
    ALLOCATE(RPPlanes(iPlanesOutput)%Coords(3   ,Plane%nRP(1),Plane%nRP(2)))
    ALLOCATE(RPPlanes(iPlanesOutput)%Val(   nVal,Plane%nRP(1),Plane%nRP(2)))
    IF (Plane_LocalCoords) RPPlanes(iPlanesOutput)%Coords(3,:,:) = 0.
  END DO ! iPlane
ELSE
  nPlanesOutput = 0
END IF

! Loop over the time samples and perform the actual output
DO iSample=1,nSamples
  ! Points
  IF (OutputPoints) THEN
    iPointsOutput = 0
    DO iPoint=1,nPoints
      GroupID=Points_GroupIDlist(iPoint)
      IF(.NOT.OutputGroup(GroupID)) CYCLE
      iPointsOutput = iPointsOutput + 1
      ! coordinates
      RPPoints%Coords(:,iPointsOutput) = xF_RP(:,Points_IDlist(iPoint))
      ! values
      RPPoints%Val(   :,iPointsOutput) = Value(:,Points_IDlist(iPoint),iSample)
    END DO
  END IF
  
  ! Lines
  IF (OutputLines) THEN
    iLinesOutput = 0
    DO iLine=1,nLines
      Line=>Lines(iLine)
      IF(.NOT.OutputGroup(Line%GroupID)) CYCLE
      iLinesOutput = iLinesOutput + 1
      ! coordinates
      IF (Line_LocalCoords) THEN
        RPLines(iLinesOutput)%Coords(1,:) = Line%LocalCoord(:)
      ELSE
        RPLines(iLinesOutput)%Coords(:,:) = xF_RP(:,Line%IDlist(:))
      END IF
      ! values
      RPLines(iLinesOutput)%Val(:,:) = Value(:,Line%IDlist(:),iSample)
      GroupName=GroupNames(Line%GroupID)
      RPPlanes(iLinesOutput)%name = TRIM(GroupName)//'_'//TRIM(Line%Name)
    END DO ! iLine
  END IF
  
  ! Planes
  IF (OutputPlanes) THEN
    iPlanesOutput = 0
    DO iPlane=1,nPlanes
      Plane=>Planes(iPlane)
      IF(.NOT.OutputGroup(Plane%GroupID)) CYCLE
      iPlanesOutput = iPlanesOutput + 1
      ! coordinates
      IF (Plane_LocalCoords) THEN
        RPPlanes(iPlanesOutput)%Coords(1:2,:,:) = Plane%LocalCoord(:,:,:)
      ELSE
        ! global xyz coordinates
        ALLOCATE(PlaneCoord(3,Plane%nRP(1),Plane%nRP(2)))
        DO j=1,Plane%nRP(2)
          DO i=1,Plane%nRP(1)
            PlaneCoord(:,i,j)=xF_RP(:,Plane%IDlist(i,j))
          END DO ! i
        END DO ! j
        RPPlanes(iPlanesOutput)%Coords(:,:,:) = PlaneCoord
        DEALLOCATE(PlaneCoord)
      END IF
      ! values
      ALLOCATE(PlaneData(1:nVal,Plane%nRP(1),Plane%nRP(2)))
      DO j=1,Plane%nRP(2)
        DO i=1,Plane%nRP(1)
          PlaneData(:,i,j)=Value(:,Plane%IDlist(i,j),iSample)
        END DO ! i
      END DO ! j
      RPPlanes(iPlanesOutput)%Val(:,:,:) = PlaneData
      DEALLOCATE(PlaneData)
      ! name
      GroupName=GroupNames(Plane%GroupID)
      RPPlanes(iPlanesOutput)%name = TRIM(GroupName)//'_'//TRIM(Plane%Name)
    END DO ! iPlane
  END IF

  ! Write to VTK
  ZoneTitle = TIMESTAMP(ProjectName,Time(iSample))
  CALL WriteStructuredDataToVTK(ZoneTitle,nLinesOutput,nPlanesOutput,RPPoints,RPLines,RPPlanes,.TRUE.,nVal,VarNames)
END DO

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK

#ifdef WITHBLPROPS
!===================================================================================================================================
!> Subroutine to write point data to VTU file 
!===================================================================================================================================
SUBROUTINE WriteBLPropsToVTU(FileString)
! MODULES
USE MOD_Globals
USE VTU
USE MOD_IO_VTU
USE MOD_VTU_Output
USE MOD_Equation_Vars   ,ONLY:nBLProps,VarNames_BLProps
USE MOD_ParametersVisu      ,ONLY:ProjectName
USE MOD_RPSetVisuVisu_Vars      ,ONLY:GroupNames 
USE MOD_RPSetVisuVisu_Vars      ,ONLY:nPlanes,Planes,tPlane
USE MOD_RPSetVisuVisu_Vars      ,ONLY:xF_RP
USE MOD_OutputRPVisu_Vars     ,ONLY:nCoords,CoordNames
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileString !< Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iVar
INTEGER              :: nCoords_loc
INTEGER              :: iPoint,iPlane
INTEGER              :: size_offsetdim,offsetVar
CHARACTER(LEN=255)   :: ZoneTitle
CHARACTER(LEN=255)   :: GroupName
CHARACTER(LEN=255)   :: CoordNames_loc(nCoords-1)
TYPE(tPlane),POINTER :: Plane
REAL,ALLOCATABLE     :: LineCoord(:)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')" WRITE BOUNDARY LAYER PROPERTY DATA TO VTU FILE '",TRIM(FileString),"'..."
CALL OpenDataFile(Filestring,create=.TRUE.,single=.FALSE.,readOnly=.FALSE.)

! Write dataset attributes
CALL WriteAttribute(File_ID,'File_Type',1,StrScalar='RP_Output')
CALL WriteAttribute(File_ID,'ProjectName',1,StrScalar=TRIM(ProjectName))
CALL WriteAttribute(File_ID,'VarNames',nBLProps,StrArray=VarNames_BLProps)
nCoords_loc=nCoords-1
CoordNames_loc=CoordNames(2:nCoords)
CALL WriteAttribute(File_ID,'CoordNames',nCoords_loc,StrArray=CoordNames_loc)

! local coord always active for BL props
size_offsetdim=5
! Planes 
DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF(Plane%Type.EQ.2) THEN ! BLPlane
    GroupName=GroupNames(Plane%GroupID)
  
    !coordinates
    ZoneTitle(1:255)=' '
    WRITE(ZoneTitle,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name),'_BLProps_X'
    ALLOCATE(LineCoord(Plane%nRP(1)))
    OffsetVar=0
    !local coordinates
    DO iVar=1,2
      LineCoord(:)=Plane%LocalCoord(iVar,:,1)
      CALL WriteArray(TRIM(ZoneTitle),size_offsetdim,2,(/1,Plane%nRP(1)/),OffsetVar,1,RealArray=LineCoord(:))
      OffsetVar=OffsetVar+1
    END DO !iVar
    !global coordinates
    DO iPoint=1,Plane%nRP(1)
      LineCoord(iPoint)=xF_RP(iVar,Plane%IDlist(iPoint,1))
    END DO ! iPoint
    CALL WriteArray(TRIM(ZoneTitle),size_offsetdim,2,(/3,Plane%nRP(1)/),OffsetVar,1,RealArray=LineCoord(:))
    DEALLOCATE(LineCoord)
  
    !values
    WRITE(ZoneTitle,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name),'_BLProps'
    size_offsetdim=nBLprops
    OffsetVar=0
    CALL WriteArray(TRIM(ZoneTitle),size_offsetdim,2,(/nBLProps,Plane%nRP(1)/),&
                    OffsetVar,1,RealArray=Plane%BLProps(:,:))
  END IF!(Plane%Type.EQ.2) THEN ! BLPlane
END DO ! iPlane

! Close the file.
CALL CloseDataFile()

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteBLPropsToVTU
#endif

END MODULE MOD_OutputRPVisu_VTK
