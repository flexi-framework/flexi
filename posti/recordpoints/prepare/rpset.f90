#include "flexi.h"

!===================================================================================================================================
!> Module to handle the Recordpoints
!===================================================================================================================================
MODULE MOD_RPSet
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersRPSet
  MODULE PROCEDURE DefineParametersRPSet
END INTERFACE

INTERFACE InitRPSet
  MODULE PROCEDURE InitRPSet
END INTERFACE

PUBLIC :: DefineParametersRPSet,InitRPSet
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize the record point structure by reading the desired types from the parameter file
!===================================================================================================================================
SUBROUTINE DefineParametersRPSet()
! MODULES
USE MOD_Parameters
USE MOD_ReadInTools ,ONLY: prms
!===================================================================================================================================
CALL prms%SetSection("Prepare Record Points: RPSet definition")
CALL prms%CreateStringOption(   'GroupName'         ,"TODO",multiple=.TRUE.)
!Line
CALL prms%CreateIntOption(      'Line_GroupID'      ,"TODO",multiple=.TRUE.)
CALL prms%CreateIntOption(      'Line_nRP'          ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Line_xstart'       ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Line_xend'         ,"TODO",multiple=.TRUE.)
!Circle
CALL prms%CreateIntOption(      'Circle_GroupID'    ,"TODO",multiple=.TRUE.)
CALL prms%CreateIntOption(      'Circle_nRP'        ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Circle_Center'     ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Circle_Axis'       ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Circle_Dir'        ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Circle_Radius'     ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Circle_Angle'      ,"TODO",multiple=.TRUE.)
!Custom lines
CALL prms%CreateIntOption(      'CustomLine_GroupID',"TODO",multiple=.TRUE.)
CALL prms%CreateIntOption(      'CustomLine_nRP'    ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('CustomLine_x'      ,"TODO",multiple=.TRUE.)
!Points
CALL prms%CreateIntOption(      'Point_GroupID'     ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Point_x'           ,"TODO",multiple=.TRUE.)
!Plane
CALL prms%CreateIntOption(      'Plane_GroupID'     ,"TODO",multiple=.TRUE.)
CALL prms%CreateIntArrayOption( 'Plane_nRP'         ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Plane_CornerX'     ,"TODO",multiple=.TRUE.)
!Sphere
CALL prms%CreateIntOption(      'Sphere_GroupID'    ,"TODO",multiple=.TRUE.)
CALL prms%CreateIntArrayOption( 'Sphere_nRP'        ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Sphere_Center'     ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Sphere_Axis'       ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Sphere_Dir'        ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Sphere_Radius'     ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Sphere_Angle'      ,"TODO",multiple=.TRUE.)
!Boundary Layer Plane
CALL prms%CreateIntOption(      'BLPlane_GroupID'   ,"TODO",multiple=.TRUE.)
CALL prms%CreateIntArrayOption( 'BLPlane_nRP'       ,"TODO",multiple=.TRUE.)
CALL prms%CreateIntOption(      'BLPlane_nCP'       ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('BLPlane_CP'        ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealOption(     'BLPlane_fac'       ,"TODO",multiple=.TRUE.)
CALL prms%CreateRealOption(     'BLPlane_height'    ,"TODO",multiple=.TRUE.)
END SUBROUTINE DefineParametersRPSet

!===================================================================================================================================
!> Initialize the record point structure by reading the desired types from the parameter file
!===================================================================================================================================
SUBROUTINE InitRPSet()
! MODULES
USE MOD_Globals
USE MOD_Readintools   ,ONLY:GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY,GETINTARRAY,CountOption
USE MOD_BLProjection  ,ONLY:GetBLPlane
USE MOD_RPSet_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                      :: anythingThere
INTEGER                      :: iGr,iLine,iP,iRP,iRP_gr,iPlane,i,j
REAL                         :: x(3),x_dummy(12)
REAL                         :: xi,eta
INTEGER                      :: nlinLines,nCircles,nCustomLines
INTEGER                      :: nflatPlanes,nSphericPlanes,nBLPlanes
REAL                         :: Circle_Center(3),Circle_Radius,Circle_Angle,Circle_Axis(3),Circle_dir(3),RotMat(3,3),dphi,phi,pi
REAL                         :: Sphere_Center(3),Sphere_Radius,Sphere_Angle,Sphere_Axis(3),Sphere_dir(3),theta,dtheta
INTEGER                      :: iCP,nCP
REAL                         :: fac
REAL,ALLOCATABLE             :: xCP(:,:),height(:)
TYPE(tRP),POINTER            :: aRP
TYPE(tLine),POINTER          :: aLine
TYPE(tPlane),POINTER         :: Plane
!===================================================================================================================================
IF(RPSetInitIsDone)THEN
   CALL abort(__STAMP__, &
        'InitRPSet not ready to be called or already called.')
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RECORDPOINT SET...'
 
SWRITE(UNIT_stdOut,'(A)')' Read recordpoint definitions from parameter file...'

pi=ACOS(-1.)

! Groups
nGroups=CountOption('GroupName')
ALLOCATE(Groups(1:nGroups))
DO iGr=1,nGroups
  Groups(iGr)%Name=GETSTR('GroupName')
  Groups(iGr)%nRP=0
END DO ! iGr

nRP_global=0

anythingThere=.FALSE.
! ----------------------------------------------------------------------------------------------------
! Lines and Circles
! ----------------------------------------------------------------------------------------------------
nlinLines    = CountOption('Line_GroupID')
nCircles     = CountOption('Circle_GroupID')
nCustomLines = CountOption('CustomLine_GroupID')
nLines       = nlinLines+nCircles+nCustomLines
IF(nLines.GT.0) THEN
  anythingThere=.TRUE.
  ALLOCATE(Lines(1:nLines))
! linear lines
  DO iLine=1,nlinLines 
    aLine=>Lines(iLine)
    aLine%GroupID=GETINT('Line_GroupID')
    WRITE(aLine%Name,'(A5,I6.6)')'Line_',iLine
    aLine%nRP   =GETINT('Line_nRP')
    ALLOCATE(aLine%RP_ptr(1:aLine%nRP))
    aLine%xStart=GETREALARRAY('Line_xstart',3) 
    aLine%xEnd  =GETREALARRAY('Line_xend',3) 
    DO iP=1,aLine%nRP
      x=aLine%xStart + (aLine%xEnd-aLine%xStart)*(iP-1)/(aLine%nRP-1)
      CALL GetNewRP(aLine%RP_ptr(iP)%RP,aLine%GroupID,x)
    END DO ! iP
  END DO ! iLine
! circle "lines"
  DO iLine=nlinLines+1,nlinLines+nCircles
    aLine=>Lines(iLine)
    aLine%GroupID=GETINT('Circle_GroupID')
    WRITE(aLine%Name,'(A7,I6.6)')'Circle_',iLine
    aLine%nRP   =GETINT('Circle_nRP')
    ALLOCATE(aLine%RP_ptr(1:aLine%nRP))
    Circle_Center = GETREALARRAY('Circle_Center',3) 
    Circle_Radius = GETREAL('Circle_Radius')
    ! angle in degrees. 360 is full circle
    Circle_Angle  = GETREAL('Circle_Angle')
    ! rotation axis of the circle
    Circle_Axis   = GETREALARRAY('Circle_Axis',3)
    Circle_Axis   = Circle_Axis/NORM2(Circle_Axis)
    ! 0 vector to define circumferential coordinate phi=0, i.e. first point
    Circle_dir    = GETREALARRAY('Circle_dir',3)    
    Circle_dir    = Circle_dir - SUM(Circle_dir*Circle_Axis)*Circle_Axis ! ensure orthogonality between dir and axis
    IF(NORM2(Circle_dir).LT.1e-9) THEN
      SWRITE(UNIT_stdOut,'(A)') 'Check definitions: Circle_dir seems to be parallel to Circle_Axis!!!'; STOP
    END IF
    Circle_dir   = Circle_dir/NORM2(Circle_dir)
    RotMat(:,1)=Circle_dir(:)
    RotMat(:,2)=CROSS(Circle_Axis,Circle_dir)                  ! right hand system
    RotMat(:,3)=Circle_Axis(:)
    IF(Circle_Angle.EQ.360.) THEN
      dphi=2.*pi/REAL(aLine%nRP)
    ELSE
      dphi=Circle_Angle*pi/180./REAL(aLine%nRP-1)
    END IF
    DO iP=1,aLine%nRP
      ! circle local coordinates
      phi=REAL(iP-1)*dphi
      x=(/COS(phi),SIN(phi),0./)
      ! coordinate transform to physical space
      x=Circle_Center +Circle_Radius*MATMUL(RotMat,x)
      CALL GetNewRP(aLine%RP_ptr(iP)%RP,aLine%GroupID,x)
    END DO ! iP
  END DO ! iLine
! custom Lines 
  DO iLine=nlinLines+nCircles+1,nLines
    aLine=>Lines(iLine)
    aLine%GroupID=GETINT('CustomLine_GroupID')
    WRITE(aLine%Name,'(A11,I6.6)')'CustomLine_',iLine
    aLine%nRP   =GETINT('CustomLine_nRP')
    ALLOCATE(aLine%RP_ptr(1:aLine%nRP))
    DO iP=1,aLine%nRP
      x=GETREALARRAY('CustomLine_x',3)
      CALL GetNewRP(aLine%RP_ptr(iP)%RP,aLine%GroupID,x)
    END DO ! iP
  END DO ! iLine
END IF
  
! ----------------------------------------------------------------------------------------------------
! Points
! ----------------------------------------------------------------------------------------------------
nPoints =CountOption('Point_GroupID')
IF(nPoints.GT.0) THEN
  anythingThere=.TRUE.
  ALLOCATE(Points(1:nPoints))
  DO iP=1,nPoints 
    Points(iP)%GroupID=GETINT('Point_GroupID')
    x     =GETREALARRAY('Point_x',3) 
    CALL GetNewRP(Points(iP)%RP,Points(iP)%GroupID,x)
  END DO ! iP  
END IF

! ----------------------------------------------------------------------------------------------------
! Planes
! ----------------------------------------------------------------------------------------------------
nFlatPlanes    =CountOption('Plane_GroupID')
nSphericPlanes =CountOption('Sphere_GroupID')
nBLPlanes      =CountOption('BLPlane_GroupID')
nPlanes=nflatPlanes+nSphericPlanes+nBLPlanes
IF(nPlanes.GT.0) THEN
  anythingThere=.TRUE.
  ALLOCATE(Planes(1:nPlanes))
! flat planes
  DO iPlane=1,nflatPlanes
    Plane=>Planes(iPlane)
    Plane%GroupID=GETINT('Plane_GroupID')
    WRITE(Plane%Name,'(A5,I6.6)')'Plane_',iPlane
    x_dummy(1:12) = GETREALARRAY('Plane_CornerX',12)
    DO iP=1,4
      Plane%x(1:3,iP)=x_dummy(1+3*(iP-1):3+3*(iP-1))
    END DO ! iPoint
    Plane%nRP(1:2)   =GETINTARRAY('Plane_nRP',2) 
    ALLOCATE(Plane%RP_ptr(1:Plane%nRP(1),1:Plane%nRP(2)))
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        ! bilinear mapping
        xi = REAL(i-1)/REAL(Plane%nRP(1)-1)
        eta= REAL(j-1)/REAL(Plane%nRP(2)-1)
        x(1:3)=   Plane%x(1:3,1) * (1.-xi) * (1.-eta) &
               +  Plane%x(1:3,2) * (   xi) * (1.-eta) &
               +  Plane%x(1:3,3) * (   xi) * (   eta) & 
               +  Plane%x(1:3,4) * (1.-xi) * (   eta) 
        CALL GetNewRP(Plane%RP_ptr(i,j)%RP,Plane%GroupID,x)
      END DO ! i
    END DO ! j
  END DO ! iPlane

! spherical planes
  DO iPlane=nFlatPlanes+1,nFlatplanes+nSphericPlanes
    Plane=>Planes(iPlane)
    Plane%GroupID=GETINT('Sphere_GroupID')
    WRITE(Plane%Name,'(A5,I6.6)')'Sphere_',iPlane
    Sphere_Center = GETREALARRAY('Sphere_Center',3)
    Sphere_Radius = GETREAL('Sphere_Radius')
    ! Azimuth angle theta in degrees. 360 is full Sphere
    Sphere_Angle= GETREAL('Sphere_Angle')
    ! rotation axis of the Sphere
    Sphere_Axis   = GETREALARRAY('Sphere_Axis',3)
    Sphere_Axis   = Sphere_Axis/NORM2(Sphere_Axis)
    ! 0 vector to define circumferential coordinate phi=0, i.e. first point
    Sphere_dir    = GETREALARRAY('Sphere_dir',3)    
    Sphere_dir    = Sphere_dir - SUM(Sphere_dir*Sphere_Axis)*Sphere_Axis ! ensure orthogonality between dir and axis
    IF(NORM2(Sphere_dir).LT.1e-9) THEN
      SWRITE(UNIT_stdOut,'(A)') 'Check definitions: Sphere_dir seems to be parallel to Sphere_Axis!!!'; STOP
    END IF
    Sphere_dir   = Sphere_dir/NORM2(Sphere_dir)
    RotMat(:,1)=Sphere_dir(:)
    RotMat(:,2)=CROSS(Sphere_Axis,Sphere_dir)                  ! right hand system
    RotMat(:,3)=Sphere_Axis(:)
    Plane%nRP(1:2)   =GETINTARRAY('Sphere_nRP',2) 
    IF(Sphere_Angle.EQ.360.) THEN
      dphi=2.*pi/REAL(Plane%nRP(1))
    ELSE
      dphi=Sphere_Angle*pi/180./REAL(Plane%nRP(1)-1)
    END IF
    dtheta=0.5*pi/REAL(Plane%nRP(2)+1)     ! we dont have the points in the singularities,
                                           ! but for proper parametrization we have to count them in 
    ALLOCATE(Plane%RP_ptr(1:Plane%nRP(1),1:Plane%nRP(2)))
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        phi  = REAL(i-1)*dphi
        theta= REAL(j)*dtheta
        x=(/COS(phi)*SIN(theta),SIN(phi)*SIN(theta),COS(theta)/)
        x=Sphere_Center +Sphere_Radius*MATMUL(RotMat,x)
        CALL GetNewRP(Plane%RP_ptr(i,j)%RP,Plane%GroupID,x)
        ! coordinate transform to physical space
      END DO ! i
    END DO ! j
  END DO ! iPlane

! BL planes
  DO iPlane=nPlanes-nBLPlanes+1,nPlanes
    Plane=>Planes(iPlane)
    Plane%GroupID=GETINT('BLPlane_GroupID')
    WRITE(Plane%Name,'(A5,I6.6)')'BLPlane_',iPlane
    Plane%nRP(1:2)   =GETINTARRAY('BLPlane_nRP',2) 
    ALLOCATE(Plane%RP_ptr(1:Plane%nRP(1),1:Plane%nRP(2)))
    nCP   =GETINT('BLPlane_nCP') ! points to define spline, at least two
    fac              =GETREAL('BLPlane_fac','1.')  ! growth factor of the BL mesh
    ALLOCATE(xCP(3,nCP))
    ALLOCATE(height(nCP))
    DO iCP=1,nCP
       xCP(:,iCP) =GETREALARRAY('BLPlane_CP',3)
       height(iCP)=GETREAL('BLPlane_height') ! height of the BL mesh
    END DO! iCP=1,nCP
    ! get coordinates of the rps and allocate pointers
    CALL GetBLPlane(Plane,nCP,height,fac,xCP) 
    DEALLOCATE(xCP,height)
  END DO! iPlane
END IF

IF(.NOT.anythingThere) THEN
  SWRITE(UNIT_StdOut,*) 'No RP infos specified in parameter file, exiting...'
  CALL abort(__STAMP__,'Code stopped!')
END IF 

! Create global RP array
ALLOCATE(RPlist(nRP_global))
iRP=0
! fill Planes
IF(nPlanes.GT.0) THEN
  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
      iRP=iRP+1
      RPlist(iRP)%RP=>Planes(iPlane)%RP_ptr(i,j)%RP
      END DO ! i
    END DO ! j
  END DO !iPlane
END IF

! fill line 
IF(nLines.GT.0) THEN
  DO iLine=1,nLines
    DO iP=1,Lines(iLine)%nRP
      iRP=iRP+1
      RPlist(iRP)%RP=>Lines(iLine)%RP_ptr(iP)%RP
    END DO !iP
  END DO !iLine
END IF
! fill points
IF(nPoints.GT.0) THEN
  DO iP=1,nPoints
    iRP=iRP+1
    RPlist(iRP)%RP=>Points(iP)%RP
  END DO !iP
END IF

! Create pointers from groups to RPs
DO iGr=1,nGroups
  ALLOCATE(Groups(iGr)%RP_ptr(Groups(iGr)%nRP))
  iRP_gr=0
  DO iRP=1,nRP_global
    aRP=>RPlist(iRP)%RP
    IF(aRP%GroupID.EQ.iGr) THEN
      iRP_gr=iRP_gr+1
      Groups(iGr)%RP_ptr(iRP_gr)%RP=>aRP
    END IF
  END DO !iRP  
END DO !iGr
RPSetInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RECORDPOINTS SET DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRPSet

END MODULE MOD_RPSet
