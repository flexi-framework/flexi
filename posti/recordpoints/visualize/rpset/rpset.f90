#include "flexi.h"

!===================================================================================================================================
!> Module to handle the Recordpoints
!===================================================================================================================================
MODULE MOD_RPSetVisu
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitRPSet
  MODULE PROCEDURE InitRPSet
END INTERFACE

INTERFACE ChangeRPSet
  MODULE PROCEDURE ChangeRPSet
END INTERFACE

INTERFACE FinalizeRPSet
  MODULE PROCEDURE FinalizeRPSet
END INTERFACE

PUBLIC :: InitRPSet
PUBLIC :: ChangeRPSet
PUBLIC :: FinalizeRPSet
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize all necessary information to perform filtering
!===================================================================================================================================
SUBROUTINE InitRPSet(RP_DefFile_in)
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input
USE MOD_ParametersVisu   ,ONLY: Line_LocalCoords,Line_LocalVel,Plane_LocalCoords
USE MOD_ParametersVisu   ,ONLY: nGroups_visu,GroupNames_visu
USE MOD_RPSetVisuVisu_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: RP_DefFile_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tLine),POINTER       :: aLine
TYPE(tLine),POINTER       :: Lines_tmp(:)
TYPE(tPlane),POINTER      :: Plane
TYPE(tPlane),POINTER      :: Planes_tmp(:)
INTEGER                   :: iLine,iLine2,iPlane,iPlane2,iGr1,iGr2,iRP,i,j,iPoint
INTEGER                   :: nRP_output
INTEGER                   :: nLines_tmp,nPlanes_tmp,nPoints_tmp
LOGICAL                   :: DSexists
LOGICAL                   :: found(nGroups_visu)
LOGICAL                   :: LinesInFile=.FALSE.,PlanesInFile=.FALSE.,PointsInFile=.FALSE.  
REAL,ALLOCATABLE          :: xF_tmp(:,:)
INTEGER,ALLOCATABLE       :: Points_IDlist_tmp(:),Points_GroupIDlist_tmp(:) 
CHARACTER(LEN=255)        :: tmp255,PlaneType
!===================================================================================================================================
IF(RPSetInitIsDone)THEN
   CALL CollectiveStop(__STAMP__, &
        'InitRPSet not ready to be called or already called.')
END IF
#ifdef MPI
 IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
      'This tool is designed only for single execution!',nProcessors)
#endif /*MPI*/
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT RECORDPOINT SET...'
WRITE(UNIT_stdOut,'(A)')' Read recordpoint definitions from data file "'//TRIM(RP_DefFile_in)//'" ...'
 
! Open data file
CALL OpenDataFile(RP_DefFile_in,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Readin Groups   
CALL GetDataSize(File_ID,'GroupNames',nDims,HSize)
nGroups=HSize(1) !number of groups
DEALLOCATE(HSize)
ALLOCATE(GroupNames(1:nGroups)) 
CALL ReadArray(TRIM('GroupNames'),1,(/nGroups/),0,1,StrArray=GroupNames)


! generate output map for groups
ALLOCATE(OutputGroup(nGroups))
IF(nGroups_visu.LT.1) THEN
  OutputGroup=.TRUE.
ELSE
  found(:)=.FALSE.
  OutputGroup=.FALSE.
  DO iGr1=1,nGroups
    DO iGr2=1,nGroups_visu
      IF(TRIM(GroupNames(iGr1)).EQ.TRIM(GroupNames_visu(iGr2)))THEN
        OutputGroup(iGr1)=.TRUE.
        found(iGr2)=.TRUE.
      END IF
    END DO
  END DO !iGr1
  IF(.NOT.ALL(found)) THEN
    WRITE(UNIT_stdOut,'(A)') 'One or more of the required Groups are not in this RPSet file!'; STOP
  END IF
END IF! (nGroups_visu.LT.1)
nRP_output=0

! Readin coordinates
CALL GetDataSize(File_ID,'xF_RP',nDims,HSize)
nRP_HDF5=HSize(2) !global number of RecordPoints
DEALLOCATE(HSize)
ALLOCATE(xF_RP(3,nRP_HDF5)) 
CALL ReadArray('xF_RP',2,(/3,nRP_HDF5/),0,2,RealArray=xF_RP)

! Readin Lines   
CALL DatasetExists(File_ID,'LineNames',DSexists)
nLines_tmp=0
nLines=0
IF(DSexists) THEN
  LinesInFile=.TRUE.
  CALL GetDataSize(File_ID,'LineNames',nDims,HSize)
  nLines_tmp=HSize(1) !number of lines
  DEALLOCATE(HSize)
  ALLOCATE(Lines_tmp(1:nLines_tmp)) 
  CALL ReadArray('LineNames',1,(/nLines_tmp/),0,1,StrArray=Lines_tmp(:)%Name)
  DO iLine=1,nLines_tmp
    aLine=>Lines_tmp(iLine)
    CALL ReadAttribute(File_ID,'GroupID',1,DatasetName=TRIM(aLine%Name),IntegerScalar=aLine%GroupID)
    ! if this line is for output, get its recordpoints
    IF(OutputGroup(aLine%GroupID)) THEN
      CALL GetDataSize(File_ID,TRIM(aLine%Name),nDims,HSize)
      aLine%nRP=HSize(1) !number of recordpoints on line
      DEALLOCATE(HSize)
      ALLOCATE(aLine%IDlist(aLine%nRP))
      CALL ReadArray(TRIM(aLine%Name),1,(/aLine%nRP/),0,1,IntegerArray=aLine%IDlist)
      nLines=nLines+1
      nRP_output=nRP_output + aLine%nRP
    ! if this line is not for output, deallocate it
    END IF
  END DO
  ! now build the line list with only those lines which are in an output group
  ALLOCATE(Lines(1:nLines)) 
  iLine2=0
  DO iLine=1,nLines_tmp
    aLine=>Lines_tmp(iLine)
    IF(OutputGroup(aLine%GroupID)) THEN
      iLine2=iLine2+1
      Lines(iLine2)=Lines_tmp(iLine)
    END IF
  END DO
  DEALLOCATE(Lines_tmp)
  ! if reqired, calculate local line coordinates
  IF(Line_LocalCoords)  CALL CalcLine_LocalCoords()
  IF(Line_LocalVel)     CALL CalcLine_LocalVelTransform()
END IF!DSexists


! Readin Points
CALL DatasetExists(File_ID,'Points_IDlist',DSexists)
nPoints=0
IF(DSexists) THEN
  PointsInFile=.TRUE.
  CALL GetDataSize(File_ID,'Points_IDlist',nDims,HSize)
  nPoints_tmp=HSize(1) !number of points on file
  DEALLOCATE(HSize)
  ! first read in all points from file
  ALLOCATE(Points_IDlist_tmp(nPoints_tmp))
  CALL ReadArray('Points_IDlist',1,(/nPoints_tmp/),0,1,IntegerArray=Points_IDlist_tmp(:))
  ALLOCATE(Points_GroupIDlist_tmp(nPoints_tmp))
  CALL ReadArray('Points_GroupIDlist',1,(/nPoints_tmp/),0,1,IntegerArray=Points_GroupIDlist_tmp(:))
  !check if group is for output
  DO iRP=1,nPoints_tmp
    IF(OutputGroup(Points_GroupIDlist_tmp(iRP))) THEN
      nPoints=nPoints+1
    END IF
  END DO !iRP
  ! now filter out all points that are not in an output group
  ALLOCATE(Points_IDlist(nPoints))
  ALLOCATE(Points_GroupIDlist(nPoints))
  iPoint=0
  DO iRP=1,nPoints_tmp
    IF(OutputGroup(Points_GroupIDlist_tmp(iRP))) THEN
      iPoint=iPoint+1
      Points_IDlist(iPoint)=Points_IDlist_tmp(iRP)
      Points_GroupIDlist(iPoint)=Points_GroupIDlist_tmp(iRP)
    END IF
  END DO !iRP
  
  nRP_output=nRP_output + nPoints
END IF!DSexists

! Readin Planes   
CALL DatasetExists(File_ID,'PlaneNames',DSexists)
nPlanes=0
IF(DSexists) THEN
  PlanesInFile=.TRUE.
  CALL GetDataSize(File_ID,'PlaneNames',nDims,HSize)
  nPlanes_tmp=HSize(1) !number of Planes
  DEALLOCATE(HSize)
  ALLOCATE(Planes_tmp(1:nPlanes_tmp)) 
  CALL ReadArray('PlaneNames',1,(/nPlanes_tmp/),0,1,StrArray=Planes_tmp(:)%Name)
  DO iPlane=1,nPlanes_tmp 
    Plane=>Planes_tmp(iPlane)
    CALL ReadAttribute(File_ID,'GroupID',1,DatasetName=TRIM(Plane%Name),IntegerScalar=Plane%GroupID)
    !check if group is for output
    IF(OutputGroup(Plane%GroupID)) THEN
      nPlanes=nPlanes+1
      CALL GetDataSize(File_ID,TRIM(Plane%Name),nDims,HSize)
      Plane%nRP(1)=HSize(1) !i number of recordpoints on Plane
      Plane%nRP(2)=HSize(2) !j number of recordpoints on Plane
      nRP_output=nRP_output+Plane%nRP(1)*Plane%nRP(2)
      DEALLOCATE(HSize)
      ALLOCATE(Plane%IDlist(Plane%nRP(1),Plane%nRP(2)))
      CALL ReadArray(TRIM(Plane%Name),2,(/Plane%nRP(1),Plane%nRP(2)/),0,1,IntegerArray=Plane%IDlist)
      ! readin norm and tangential vectors if suitable
      PlaneType=TRIM(Plane%Name(1:5))
      IF(PlaneType.EQ.TRIM("Spher")) THEN
        Plane%Type=1
      ELSEIF(PlaneType.EQ.TRIM("BLPla")) THEN
        Plane%Type=2
        ALLOCATE(Plane%NormVec(3,Plane%nRP(1))) 
        WRITE(tmp255,'(A,A)')TRIM(Plane%Name),'_NormVec'
        CALL ReadArray(tmp255,2,(/3,Plane%nRP(1)/),0,2,RealArray=Plane%NormVec)
        ALLOCATE(Plane%TangVec(3,Plane%nRP(1))) 
        WRITE(tmp255,'(A,A)')TRIM(Plane%Name),'_TangVec'
        CALL ReadArray(tmp255,2,(/3,Plane%nRP(1)/),0,2,RealArray=Plane%TangVec)
      END IF
    END IF
  END DO
  ! now build the plane list for those planes in an output group
  ALLOCATE(Planes(1:nPlanes)) 
  iPlane2=0
  DO iPlane=1,nPlanes_tmp
    Plane=>Planes_tmp(iPlane)
    IF(OutputGroup(Plane%GroupID)) THEN
      iPlane2=iPlane2+1
      Planes(iPlane2)=Planes_tmp(iPlane)
    END IF
  END DO
  DEALLOCATE(Planes_tmp)
  ! if reqired, calculate local (fitted) plane coordinates
  IF(Plane_LocalCoords)  CALL CalcPlane_LocalCoords()
END IF!DSexists

! build mapping from output RPs (1:nRP_output) to RPs on RPSet and RPData files (1:nRP_HDF5).
ALLOCATE(RPOutMap(1:nRP_output))
iRP=0
! lines
IF(LinesInFile) THEN
  DO iLine=1,nLines
    aLine=>Lines(iLine)
    RPOutMap(iRP+1:iRP+aLine%nRP)=aLine%IDlist(:)
! redefine line mapping to newly generated global RP list
!    aLine%IDlist(:) = (/ (l, l=1,aLine%nRP,1) /) +iRP
    DO i=1,aLine%nRP
      aLine%IDlist(i)=iRP+i
    END DO
    iRP=iRP+aLine%nRP
  END DO
END IF
! points
IF(PointsInFile) THEN
  RPOutMap(iRP+1:iRP+nPoints)=Points_IDlist(:)
  ! redefine points mapping to newly generated global RP list
  DO i=1,nPoints
    Points_IDlist(i)=iRP+i
  END DO
  iRP=iRP+nPoints
END IF
!Planes
IF(PlanesInFile) THEN
  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    DO j=1,Plane%nRP(2)  
      DO i=1,Plane%nRP(1)  
        RPOutMap(iRP+1)=Plane%IDlist(i,j)
        ! redefine plane mapping to newly generated global RP list
        Plane%IDlist(i,j)=iRP+1
        iRP=iRP+1
      END DO
    END DO
  END DO
END IF

ALLOCATE(xF_tmp(3,nRP_HDF5)) 
xF_tmp=xF_RP
DEALLOCATE(xF_RP)
ALLOCATE(xF_RP(3,nRP_output)) 
xF_RP(1:3,:)=xF_tmp(1:3,RPOutMap(:))
nRP_global=nRP_output
WRITE(UNIT_stdOut,*)' ',nRP_global,' recordpoints to process..'
CALL CloseDataFile() 

RPSetInitIsDone = .TRUE.
WRITE(UNIT_stdOut,'(A)')' INIT RECORDPOINTS SET DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRPSet



!===================================================================================================================================
!> Change the current RP set to the new from RP_DefFile_in.
!> We check if the RPs from the old set are in the new one and generate the necessary mappings
!===================================================================================================================================
SUBROUTINE ChangeRPSet(RP_DefFile_in)
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input
USE MOD_RPSetVisuVisu_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: RP_DefFile_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iRP1,iRP2
REAL                      :: dist,distRP(nRP_global)
REAL,ALLOCATABLE          :: xF_newset(:,:)
!===================================================================================================================================
IF(.NOT.RPSetInitIsDone)THEN
   CALL abort(__STAMP__, &
        'InitRPSet not called yet!!')
END IF
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' CHANGING RECORDPOINT SET...'
WRITE(UNIT_stdOut,'(A)')' Read recordpoint definitions from data file "'//TRIM(RP_DefFile_in)//'" ...'
 
! Open data file
CALL OpenDataFile(RP_DefFile_in,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Readin coordinates
CALL GetDataSize(File_ID,'xF_RP',nDims,HSize)
nRP_HDF5=HSize(2) !global number of RecordPoints
DEALLOCATE(HSize)
ALLOCATE(xF_newset(3,nRP_HDF5)) 
CALL ReadArray('xF_RP',2,(/3,nRP_HDF5/),0,2,RealArray=xF_newset)
CALL CloseDataFile()

! compare coordinates to find the mapping between the new RP set and the RPs we want to process
distRP=HUGE(1.)
DO iRP1=1,nRP_HDF5
  DO iRP2=1,nRP_global
    dist=NORM2(xF_newset(1:3,iRP1)-xF_RP(1:3,iRP2))
    IF(dist.LT.distRP(iRP2))THEN
      RPOutMap(iRP2)=iRP1
      distRP(iRP2)=dist
    END IF
  END DO! iRP2=1,nRP_global
END DO! iRP1=1,nRP_HDF5
DEALLOCATE(xF_newset)

! check if the found RPs are within tolerance
dist=0.
DO iRP2=1,nRP_global
  dist=MAX(distRP(iRP2),dist)
  IF(distRP(iRP2).GT.1e-9) THEN 
    CALL abort(__STAMP__, &
         'Not all RPs can be found in the new RP set!!')
  END IF
END DO! iRP2=1,nRP_global
WRITE(UNIT_stdOut,*)' All RPs found. Max. deviation in RP coordinates is ',dist 

WRITE(UNIT_stdOut,'(A)')' CHANGE RECORDPOINTS SET DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE ChangeRPSet



!===================================================================================================================================
!> Calculate the local coordinate from start to end point of the line for each RP on the line
!===================================================================================================================================
SUBROUTINE CalcLine_LocalCoords()
! MODULES
USE MOD_Globals
USE MOD_RPSetVisuVisu_Vars, ONLY: nLines,tLine,xF_RP, Lines
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(len=5)    :: LineType
INTEGER             :: iLine,iPoint
REAL                :: NormLineVec(3),LocalVec(3)
REAL                :: circCenter(3),circAxis(3),circDir(3),RotMat(3,3),x(3)
TYPE(tLine),POINTER :: aLine
!===================================================================================================================================
DO iLine=1,nLines
  aLine=>Lines(iLine)
  ALLOCATE(aLine%LocalCoord(aLine%nRP))
  ! Line parallel vector
  NormLineVec = xF_RP(:,aLine%IDlist(aLine%nRP)) - xF_RP(:,aLine%IDlist(1))
  LineType=TRIM(aLine%Name(1:5))
  IF(LineType.EQ.TRIM("Line_")) THEN  ! linear line
    NormLineVec = NormLineVec/NORM2(NormLineVec)
    aLine%LocalCoord(1)=0.
    DO iPoint=2,aLine%nRP
      LocalVec=xF_RP(:,aLine%IDlist(iPoint)) - xF_RP(:,aLine%IDlist(1))
      aLine%LocalCoord(iPoint) = SUM(LocalVec(1:3)*NormLineVec(1:3))
    END DO !iPoint  
  ELSEIF(LineType.EQ.TRIM("Circl")) THEN  ! circle
    ! get circle center
    circCenter(1)=SUM(xF_RP(1,aLine%IDlist(:)))
    circCenter(2)=SUM(xF_RP(2,aLine%IDlist(:)))
    circCenter(3)=SUM(xF_RP(3,aLine%IDlist(:)))
    circCenter   = circCenter/REAL(aLine%nRP)
    ! vector from center to first point defines phi=0
    circDir   = xF_RP(:,aLine%IDlist(1))-circCenter
    circDir   = circDir/NORM2(circDir)
    ! get circle axis
    circAxis=CROSS(xF_RP(:,aLine%IDlist(1))-circCenter, &
                   xF_RP(:,aLine%IDlist(INT(0.5*aLine%nRP)))-circCenter)
    circAxis  = circAxis/NORM2(circAxis)
    RotMat(1,:)=circDir(:)
    RotMat(2,:)=CROSS(circAxis,circDir)                  ! right hand system
    RotMat(3,:)=circAxis(:)
    aLine%LocalCoord(1)=0.
    DO iPoint=2,aLine%nRP
      x = xF_RP(:,aLine%IDlist(iPoint)) -circCenter
      x = MATMUL(RotMat,x)
      aLine%LocalCoord(iPoint)=ATAN2(x(2),x(1))
    END DO
  ELSEIF(LineType.EQ.TRIM("Custo")) THEN  ! custom
    aLine%LocalCoord(1)=0.
    DO iPoint=2,aLine%nRP
      aLine%LocalCoord(iPoint)= aLine%LocalCoord(iPoint-1) &
                               +NORM2(xF_RP(:,aLine%IDlist(iPoint)) - xF_RP(:,aLine%IDlist(iPoint-1)))
    END DO !iPoint  
  ELSE
    WRITE(UNIT_StdOut,'(A,A,A)')' The type of Line "',LineType,'" is not known!'; STOP
  END IF
END DO !iLine
END SUBROUTINE CalcLine_LocalCoords



!===================================================================================================================================
!> Calculate the local coordinate from start to end point of the line for each RP on the line
!===================================================================================================================================
SUBROUTINE CalcLine_LocalVelTransform()
! MODULES
USE MOD_Globals
USE MOD_RPSetVisuVisu_Vars    ,ONLY:nLines,tLine,xF_RP, Lines
USE MOD_ParametersVisu,ONLY:Line_LocalVel_vec
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iLine
REAL                :: NormLineVec(3)
TYPE(tLine),POINTER :: aLine
!===================================================================================================================================
DO iLine=1,nLines
  aLine=>Lines(iLine)
  ALLOCATE(aLine%Tmat(3,3))
  NormLineVec = xF_RP(:,aLine%IDlist(aLine%nRP)) - xF_RP(:,aLine%IDlist(1))
  NormLineVec = NormLineVec/SQRT(SUM(NormLineVec(1:3)*NormLineVec(1:3)))
! T1= orthogonalized main velocity vector
  aLine%Tmat(1,1:3) = Line_localVel_vec-SUM(Line_localVel_vec(1:3)*NormLineVec(1:3))*NormLineVec 
! T2= normalized Line Vector
  aLine%Tmat(2,1:3) = NormLineVec      
! T3= crossprod to guarantee right hand system
  aLine%Tmat(3,1:3) = CROSS(aLine%Tmat(1,1:3),aLine%Tmat(2,1:3))
END DO !iLine
END SUBROUTINE CalcLine_LocalVelTransform



!===================================================================================================================================
!> Calculate the local coordinate for BLPlanes through a spline representation
!===================================================================================================================================
SUBROUTINE CalcPlane_LocalCoords()
! MODULES
USE MOD_Globals
USE MOD_RPSetVisuVisu_Vars, ONLY: nPlanes,tPlane,xF_RP, Planes
USE MOD_Spline,     ONLY: GetSpline
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER              :: iPlane,i,j,iS
INTEGER              :: Nx,Ny,nSuper
REAL                 :: x_loc(3),x_loc_old(3),t_loc
REAL,ALLOCATABLE     :: xPlane_tmp(:,:,:),localCoordX(:),localCoordY(:),t_Nx(:),coeff(:,:,:)
TYPE(tPlane),POINTER :: Plane
!===================================================================================================================================
DO iPlane=1,nPlanes
  Plane=>Planes(iPlane)
  IF(Plane%Type.NE.2) CYCLE
  Nx=Plane%nRP(1)
  Ny=Plane%nRP(2)
  ALLOCATE(Plane%LocalCoord(2,Nx,Ny))
  ALLOCATE(xPlane_tmp(3,Nx,Ny))
  ALLOCATE(localCoordX(Nx))
  ALLOCATE(localCoordY(Ny))
  DO i=1,Nx; DO j=1,Ny
    xPlane_tmp(:,i,j)=xF_RP(:,Plane%IDlist(i,j))
  END DO; END DO

  ! Plane parallel direction -------------------------------------------------------
  ALLOCATE(t_Nx(Nx))
  ALLOCATE(coeff(3,4,Nx-1))
  CALL GetSpline(3,Nx,xPlane_tmp(:,:,1),coeff,t_Nx)! get spline through the first tangential layer
  localCoordX=0.
  ! calculate arclength on supersampled points
  nSuper=10
  x_loc=xPlane_tmp(:,1,1)
  DO i=2,Nx
    localCoordX(i)=localCoordX(i-1)
    DO iS=2,nSuper
      t_loc=t_Nx(i-1)+REAL(iS-1)/REAL(nSuper-1)*(t_Nx(i)-t_Nx(i-1))
      x_loc_old=x_loc
      x_loc(:)=coeff(:,1,i-1)+coeff(:,2,i-1)*(t_loc-t_Nx(i-1)) &
              +coeff(:,3,i-1)*(t_loc-t_Nx(i-1))**2 + coeff(:,4,i-1)*(t_loc-t_Nx(i-1))**3
      localCoordX(i)=localCoordX(i)+NORM2(x_loc-x_loc_old)
    END DO
  END DO!i
  DO j=1,Ny
    Plane%LocalCoord(1,:,j)=localCoordX(:)
  END DO

  ! Plane normal direction -------------------------------------------------------
  localCoordY=0.
  DO i=1,Nx
    DO j=2,Ny
      localCoordY(j)=localCoordY(j-1)+NORM2(xPlane_tmp(:,i,j)-xPlane_tmp(:,i,j-1))
    END DO! j=1,Ny
    Plane%LocalCoord(2,i,:)=localCoordY(:)
  END DO! i=1,Nx
  DEALLOCATE(xPlane_tmp,localCoordX,localCoordY,t_Nx,coeff)
END DO !iPlane
END SUBROUTINE CalcPlane_LocalCoords



!===================================================================================================================================
!> Deallocate global variable for Recordpoints
!===================================================================================================================================
SUBROUTINE FinalizeRPSet()
! MODULES
USE MOD_RPSetVisuVisu_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: iLine,iPlane
!===================================================================================================================================
SDEALLOCATE(GroupNames)
SDEALLOCATE(RPOutMap)
SDEALLOCATE(Points_IDlist)
SDEALLOCATE(Points_GroupIDlist)
SDEALLOCATE(xF_RP)
SDEALLOCATE(x_RP)
SDEALLOCATE(OutputGroup)

IF(nLines.GT.1)THEN
  DO iLine=1,nLines
    SDEALLOCATE(Lines(iLine)%IDlist)
    SDEALLOCATE(Lines(iLine)%LocalCoord)  
    SDEALLOCATE(Lines(iLine)%Tmat)
  END DO! iLine=1,nLines
  DEALLOCATE(Lines)
END IF

IF(nPlanes.GT.1)THEN
  DO iPlane=1,nPlanes
    SDEALLOCATE(Planes(iPlane)%IDlist)
    SDEALLOCATE(Planes(iPlane)%NormVec)
    SDEALLOCATE(Planes(iPlane)%TangVec)
  END DO! iPlane=1,nPlanes
  DEALLOCATE(Planes)
END IF

RPSetInitIsDone = .FALSE.
END SUBROUTINE FinalizeRPSet

END MODULE MOD_RPSetVisu
