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

INTERFACE WriteTimeAvgDataToVTK
  MODULE PROCEDURE WriteTimeAvgDataToVTK
END INTERFACE

INTERFACE WriteBLPropsToVTK
  MODULE PROCEDURE WriteBLPropsToVTK
END INTERFACE

PUBLIC::WriteDataToVTK,WriteTimeAvgDataToVTK,WriteBLPropsToVTK
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write time-accurate data to VTK file. This routine will loop over all points, lines and planes and collect the ones
!> that should be visualized in the respective data types. Each timestep will be written to .vts files by calling the structured
!> output routine. The files will be collected in a subdirectory 'timeseries', and .pvd will be written in the main directory.
!> Those contain links to all the single timestep files and can be openend e.g. in ParaView to have them all together and with
!> correct time information.
!===================================================================================================================================
SUBROUTINE WriteDataToVTK(nSamples,nRP,nVal,VarNames,Time,Value,FileName)
! MODULES
USE MOD_Globals
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
CHARACTER(LEN=255),INTENT(IN) :: FileName                      !< First part of the file name
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
CHARACTER(LEN=255)        :: FileNamePVD
CHARACTER(LEN=255)        :: TimestepString
INTEGER                   :: ivtk=44
CHARACTER(LEN=1)          :: lf
CHARACTER(LEN=200)        :: Buffer
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

! Create a subdirectory that contains all the output files - since there are A LOT of them for timeseries
CALL SYSTEM('mkdir timeseries')

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
  ZoneTitle = 'timeseries/'//TRIM(TIMESTAMP(FileName,Time(iSample)))
  CALL WriteStructuredDataToVTK(ZoneTitle,nLinesOutput,nPlanesOutput,RPPoints,RPLines,RPPlanes,.TRUE.,nVal,VarNames)
END DO

! Write .pvd collections to easily open the timeseries

! Line feed character
lf = char(10)

! Points
IF (nPointsOutput.GT.0) THEN
  FileNamePVD=TRIM(FileName)//'_Points.pvd'
  OPEN(UNIT=ivtk,FILE=TRIM(FileNamePVD),ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  <Collection>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Loop and write the single collection entries
  DO iSample = 1,nSamples
    WRITE(TimestepString,'(F17.9)') Time(iSample)
    DO i=1,LEN(TRIM(TimestepString))
      IF(TimestepString(i:i).EQ.' ') TimestepString(i:i)='0'
    END DO
    ZoneTitle = 'timeseries/'//TRIM(TIMESTAMP(FileName,Time(iSample)))
    Buffer='    <DataSet timestep="'//TRIM(TimestepString)//'" part="0" file="'//TRIM(ZoneTitle)//&
                 '_Points.vts"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  END DO
  ! Write Footer
  Buffer='  </Collection>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
END IF

! Lines
IF (nLinesOutput.GT.0) THEN
  iLinesOutput = 0
  DO iLine=1,nLines
  Line=>Lines(iLine)
  IF(.NOT.OutputGroup(Line%GroupID)) CYCLE
    iLinesOutput = iLinesOutput + 1
    FileNamePVD=TRIM(FileName)//'_'//TRIM(RPLines(iLinesOutput)%name)//'.pvd'
    OPEN(UNIT=ivtk,FILE=TRIM(FileNamePVD),ACCESS='STREAM')
    ! Write header
    Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
    Buffer='<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
    Buffer='  <Collection>'//lf;WRITE(ivtk) TRIM(Buffer)
    ! Loop and write the single collection entries
    DO iSample = 1,nSamples
      WRITE(TimestepString,'(F17.9)') Time(iSample)
      DO i=1,LEN(TRIM(TimestepString))
        IF(TimestepString(i:i).EQ.' ') TimestepString(i:i)='0'
      END DO
      ZoneTitle = 'timeseries/'//TRIM(TIMESTAMP(FileName,Time(iSample)))
      Buffer='    <DataSet timestep="'//TRIM(TimestepString)//'" part="0" file="'//TRIM(ZoneTitle)//'_'//&
                   TRIM(RPLines(iLinesOutput)%name)//'.vts"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    END DO
    ! Write Footer
    Buffer='  </Collection>'//lf;WRITE(ivtk) TRIM(Buffer)
    Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
    CLOSE(ivtk)
  END DO
END IF

! Planes
IF (nPlanesOutput.GT.0) THEN
  iPlanesOutput = 0
  DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF(.NOT.OutputGroup(Plane%GroupID)) CYCLE
    iPlanesOutput = iPlanesOutput + 1
    FileNamePVD=TRIM(FileName)//'_'//TRIM(RPPlanes(iPlanesOutput)%name)//'.pvd'
    OPEN(UNIT=ivtk,FILE=TRIM(FileNamePVD),ACCESS='STREAM')
    ! Write header
    Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
    Buffer='<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
    Buffer='  <Collection>'//lf;WRITE(ivtk) TRIM(Buffer)
    ! Loop and write the single collection entries
    DO iSample = 1,nSamples
      WRITE(TimestepString,'(F17.9)') Time(iSample)
      DO i=1,LEN(TRIM(TimestepString))
        IF(TimestepString(i:i).EQ.' ') TimestepString(i:i)='0'
      END DO
      ZoneTitle = 'timeseries/'//TRIM(TIMESTAMP(FileName,Time(iSample)))
      Buffer='    <DataSet timestep="'//TRIM(TimestepString)//'" part="0" file="'//TRIM(ZoneTitle)//'_'//&
                   TRIM(RPPlanes(iPlanesOutput)%name)//'.vts"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    END DO
    ! Write Footer
    Buffer='  </Collection>'//lf;WRITE(ivtk) TRIM(Buffer)
    Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
    CLOSE(ivtk)
  END DO
END IF

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK

!===================================================================================================================================
!> Subroutine to write time-averaged data to VTK file. This routine will loop over all points, lines and planes and collect the ones
!> that should be visualized in the respective data types. The average will be written to .vts files by calling the structured
!> output routine.
!===================================================================================================================================
SUBROUTINE WriteTimeAvgDataToVTK(nRP,nVal,VarNames,Value,FileName)
! MODULES
USE MOD_Globals
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
INTEGER,INTENT(IN)            :: nRP                           !< Number of RP to be visualized 
INTEGER,INTENT(IN)            :: nVal                          !< Number of nodal output variables
CHARACTER(LEN=255),INTENT(IN) :: VarNames(nVal)                !< Names of all variables that will be written out
REAL,INTENT(IN)               :: Value(1:nVal,nRP)             !< Statevector 
CHARACTER(LEN=255),INTENT(IN) :: FileName                      !< First part of the file name
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iPoint,iLine,iPlane,i,j
INTEGER                   :: GroupID
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

! Collect the coordinates and values
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
    RPPoints%Val(   :,iPointsOutput) = Value(:,Points_IDlist(iPoint))
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
    RPLines(iLinesOutput)%Val(:,:) = Value(:,Line%IDlist(:))
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
        PlaneData(:,i,j)=Value(:,Plane%IDlist(i,j))
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
CALL WriteStructuredDataToVTK(FileName,nLinesOutput,nPlanesOutput,RPPoints,RPLines,RPPlanes,.TRUE.,nVal,VarNames)

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteTimeAvgDataToVTK

!===================================================================================================================================
!> Subroutine to write the boundary layer properties to the VTK format.
!> The boundary layer properties are calculated for the boundary layer planes, and the output is done on single line - the one
!> at the wall.
!===================================================================================================================================
SUBROUTINE WriteBLPropsToVTK(FileString)
! MODULES
USE MOD_Globals
USE MOD_VTKStructuredOutput
USE MOD_EquationRP_Vars    ,ONLY: nBLProps,VarNames_BLProps
USE MOD_ParametersVisu     ,ONLY: Plane_LocalCoords
USE MOD_RPSetVisuVisu_Vars ,ONLY: GroupNames 
USE MOD_RPSetVisuVisu_Vars ,ONLY: nPlanes,Planes,tPlane
USE MOD_RPSetVisuVisu_Vars ,ONLY: OutputGroup 
USE MOD_RPSetVisuVisu_Vars ,ONLY: xF_RP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileString !< Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iPlane,i
CHARACTER(LEN=255)        :: GroupName
TYPE(tPlane),POINTER      :: Plane
TYPE(RPPlane),ALLOCATABLE :: RPPlanes(:)
TYPE(RPLine),ALLOCATABLE  :: RPLines(:)
TYPE(RPPoint)             :: RPPoints
INTEGER                   :: nPlanesOutput
INTEGER                   :: iPlanesOutput
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')" WRITE BOUNDARY LAYER PROPERTY DATA TO VTU FILE '",TRIM(FileString),"'.vtu ..."

! Count the number of boundary layer planes
! Planes
nPlanesOutput = 0
DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF((.NOT.OutputGroup(Plane%GroupID)).OR.(Plane%Type.NE.2)) CYCLE
  nPlanesOutput = nPlanesOutput + 1
END DO ! iPlane
! The output of the boundary layer properties is done as a line (the bottom of the plane)
IF (nPlanesOutput.GT.0) THEN
  ALLOCATE(RPLines(nPlanesOutput))
END IF
iPlanesOutput = 0
DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF((.NOT.OutputGroup(Plane%GroupID)).OR.(Plane%Type.NE.2)) CYCLE
  iPlanesOutput = iPlanesOutput + 1
  RPLines(iPlanesOutput)%nRPs = Plane%nRP(1)
  ALLOCATE(RPLines(iPlanesOutput)%Coords(3    ,Plane%nRP(1)))
  ALLOCATE(RPLines(iPlanesOutput)%Val(nBLProps,Plane%nRP(1)))
  IF (Plane_LocalCoords) RPLines(iPlanesOutput)%Coords(3,:) = 0.
END DO ! iPlane

! Collect the coordinates and values
! Planes
iPlanesOutput = 0
DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF((.NOT.OutputGroup(Plane%GroupID)).OR.(Plane%Type.NE.2)) CYCLE
  iPlanesOutput = iPlanesOutput + 1
  ! coordinates
  IF (Plane_LocalCoords) THEN
    RPLines(iPlanesOutput)%Coords(1:2,:) = Plane%LocalCoord(:,:,1)
  ELSE
    ! global xyz coordinates
    DO i=1,Plane%nRP(1)
      RPLines(iPlanesOutput)%Coords(:,i) = xF_RP(:,Plane%IDlist(i,1))
    END DO ! i
  END IF
  ! values
  RPLines(iPlanesOutput)%Val(:,:) = Plane%BLProps(:,:)
  ! name
  GroupName=GroupNames(Plane%GroupID)
  RPLines(iPlanesOutput)%name = TRIM(GroupName)//'_'//TRIM(Plane%Name)
END DO ! iPlane

! Write to VTK
CALL WriteStructuredDataToVTK(FileString,0,nPlanesOutput,RPPoints,RPLines,RPPlanes,.TRUE.,nBLProps,VarNames_BLProps)

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"

END SUBROUTINE WriteBLPropsToVTK

END MODULE MOD_OutputRPVisu_VTK
