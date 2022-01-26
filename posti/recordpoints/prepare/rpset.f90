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
!> Module to handle the creation of the RP set. Will read in the definitions of the groups and calculate the physical coordinates
!> of the resulting points.
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
!> Define the available parameters for the RP definition.
!===================================================================================================================================
SUBROUTINE DefineParametersRPSet()
! MODULES
USE MOD_Parameters
USE MOD_ReadInTools ,ONLY: prms
!===================================================================================================================================
CALL prms%SetSection("Prepare Record Points: RPSet definition")
CALL prms%CreateStringOption(   'GroupName'         ,"Name of the RP group (one for each group!)",multiple=.TRUE.)
!Line
CALL prms%CreateIntOption(      'Line_GroupID'      ,"ID of a straight line group, defined by start and end coordinates and&
                                                      & the number of points along that line, used to allocate the definition to a&
                                                      & specific group",multiple=.TRUE.)
CALL prms%CreateIntOption(      'Line_nRP'          ,"Number of RPs on line",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Line_xstart'       ,"Coordinates of start of line",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Line_xend'         ,"Coordinates of end of line",multiple=.TRUE.)
!Circle
CALL prms%CreateIntOption(      'Circle_GroupID'    ,"ID of a circular group, used to allocate the definition to a specific group",&
                                                      multiple=.TRUE.)
CALL prms%CreateIntOption(      'Circle_nRP'        ,"Number of RPs along circle",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Circle_Center'     ,"Coordinates of circle center",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Circle_Axis'       ,"Axis vector of circle",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Circle_Dir'        ,"Vector defining the start point on the circle",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Circle_Radius'     ,"Radius of the circle",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Circle_Angle'      ,"Angle from the start point, 360° is a full circle",multiple=.TRUE.)
!Custom lines
CALL prms%CreateIntOption(      'CustomLine_GroupID',"ID of a custom line, defined by an arbitrary number of RPs&
                                                      &, used to allocate the definition to a specific group",multiple=.TRUE.)
CALL prms%CreateIntOption(      'CustomLine_nRP'    ,"Number of points on the custom line",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('CustomLine_x'      ,"Coordinates of the points on the custom line",multiple=.TRUE.)
!Points
CALL prms%CreateIntOption(      'Point_GroupID'     ,"ID of a point group, used to allocate the definition to a specific group",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Point_x'           ,"Coordinates of the single point",multiple=.TRUE.)
!Plane
CALL prms%CreateIntOption(      'Plane_GroupID'     ,"ID of a plane group, defined by the corner points and the number of points in&
                                                      & both directions, used to allocate the definition to a specific group",multiple=.TRUE.)
CALL prms%CreateIntArrayOption( 'Plane_nRP'         ,"Number of points in the plane",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Plane_CornerX'     ,"Coordinates of the 4 corner points (x1,y1,z1,x2,y2,z2,...)",multiple=.TRUE.)
!Box
CALL prms%CreateIntOption(      'Box_GroupID'       ,"ID of a box group, defined by the corner points and the number of points in&
                                                      & both directions, used to allocate the definition to a specific group",multiple=.TRUE.)
CALL prms%CreateIntArrayOption( 'Box_nRP'           ,"Number of points in the box",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Box_CornerX'       ,"Coordinates of the 8 corner points (x1,y1,z1,x2,y2,z2,...)",multiple=.TRUE.)
!Sphere
CALL prms%CreateIntOption(      'Sphere_GroupID'    ,"ID of a spherical group, with points on the circumference, used to allocate&
                                                     & the definition to a specific group",multiple=.TRUE.)
CALL prms%CreateIntArrayOption( 'Sphere_nRP'        ,"Number of points on the spere in phi and theta direction",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Sphere_Center'     ,"Coordinates of sphere center",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Sphere_Axis'       ,"Axis vector of sphere",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('Sphere_Dir'        ,"Vector defining the start point on the sphere",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Sphere_Radius'     ,"Radius of the sphere",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Sphere_Angle'      ,"Phi angle of the sphere (360° is a full sphere)",multiple=.TRUE.)
!Boundary Layer Plane
CALL prms%CreateIntOption(      'BLPlane_GroupID'   ,"ID of a boundary layer group - works like a plane group, but the plane is&
                                                      & created by projecting the points of a spline to the nearest boundary and&
                                                      & extruding the plane along the normal with a stretching factor&
                                                      &, used to allocate the definition to a specific group",multiple=.TRUE.)
CALL prms%CreateIntArrayOption( 'BLPlane_nRP'       ,"Number of RPs along and normal to the boundary",multiple=.TRUE.)
CALL prms%CreateIntOption(      'BLPlane_nCP'       ,"Number of control points defining the spline (at least two)",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('BLPlane_CP'        ,"Coordinates of the spline control points",multiple=.TRUE.)
CALL prms%CreateRealOption(     'BLPlane_fac'       ,"Factor of geometrical stretching in wall-normal direction",multiple=.TRUE.)
CALL prms%CreateRealOption(     'BLPlane_height'    ,"Wall-normal extend of the plane for each control point",multiple=.TRUE.)
!Boundary Layer Box
CALL prms%CreateIntOption(      'BLBox_GroupID'     ,"ID of a boundary layer group - works like a box group, but the box is&
                                                     & created by projecting the points of a spline to the nearest boundary and&
                                                     & extruding the box along the normal with a stretching factor&
                                                     &, used to allocate the definition to a specific group",multiple=.TRUE.)
CALL prms%CreateIntArrayOption( 'BLBox_nRP'         ,"Number of RPs along and normal to the boundary",multiple=.TRUE.)
CALL prms%CreateIntOption(      'BLBox_nCP'         ,"Number of control points defining the spline (at least two)",multiple=.TRUE.)
CALL prms%CreateIntOption(      'BLBox_nSP'         ,"Number of splines in z direction (at least one)",multiple=.TRUE.)
CALL prms%CreateRealArrayOption('BLBox_CP'          ,"Coordinates of the spline control points",multiple=.TRUE.)
CALL prms%CreateRealOption(     'BLBox_fac'         ,"Factor of geometrical stretching in wall-normal direction",multiple=.TRUE.)
CALL prms%CreateRealOption(     'BLBox_height'      ,"Wall-normal extend of the box for each control point",multiple=.TRUE.)
END SUBROUTINE DefineParametersRPSet

!===================================================================================================================================
!> Initialize the record point structure by reading the desired types from the parameter file and calculating the physical
!> coordinates of the RPs.
!===================================================================================================================================
SUBROUTINE InitRPSet()
! MODULES
USE MOD_Globals
USE MOD_Mathtools     ,ONLY:CROSS
USE MOD_Readintools   ,ONLY:GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY,GETINTARRAY,CountOption
USE MOD_BLProjection  ,ONLY:GetBLPlane,GetBLBox
USE MOD_RPSet_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                      :: anythingThere
INTEGER                      :: iGr,iLine,iP,iRP,iRP_gr,iPlane,iBox,i,j,k
REAL                         :: x(3),x_dummy(24)
REAL                         :: xi,eta,zeta
INTEGER                      :: nlinLines,nCircles,nCustomLines
INTEGER                      :: nflatPlanes,nSphericPlanes,nBLPlanes
INTEGER                      :: nFlatBoxes,nBLBoxes
REAL                         :: Circle_Center(3),Circle_Radius,Circle_Angle,Circle_Axis(3),Circle_dir(3),RotMat(3,3),dphi,phi,pi
REAL                         :: Sphere_Center(3),Sphere_Radius,Sphere_Angle,Sphere_Axis(3),Sphere_dir(3),theta,dtheta
INTEGER                      :: iCP,nCP,iSP,nSP
REAL                         :: fac
REAL,ALLOCATABLE             :: xCP(:,:,:),height(:,:)
TYPE(tRP),POINTER            :: aRP
TYPE(tLine),POINTER          :: aLine
TYPE(tPlane),POINTER         :: Plane
TYPE(tBox),POINTER           :: Box
!===================================================================================================================================
IF(RPSetInitIsDone)THEN
   CALL Abort(__STAMP__, &
        'InitRPSet not ready to be called or already called.')
   RETURN
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
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
    Circle_Axis   = Circle_Axis/SQRT(DOT_PRODUCT(Circle_Axis,Circle_Axis))
    ! 0 vector to define circumferential coordinate phi=0, i.e. first point
    Circle_dir    = GETREALARRAY('Circle_dir',3)
    Circle_dir    = Circle_dir - SUM(Circle_dir*Circle_Axis)*Circle_Axis ! ensure orthogonality between dir and axis
    IF(SQRT(DOT_PRODUCT(Circle_dir,Circle_dir)).LT.1e-9) THEN
      SWRITE(UNIT_stdOut,'(A)') 'Check definitions: Circle_dir seems to be parallel to Circle_Axis!!!'; STOP
    END IF
    Circle_dir   = Circle_dir/SQRT(DOT_PRODUCT(Circle_dir,Circle_dir))
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
  DO iLine=nlinLines+nCircles+1,nlinLines+nCircles+nCustomLines
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
    Sphere_Axis   = Sphere_Axis/SQRT(DOT_PRODUCT(Sphere_Axis,Sphere_Axis))
    ! 0 vector to define circumferential coordinate phi=0, i.e. first point
    Sphere_dir    = GETREALARRAY('Sphere_dir',3)
    Sphere_dir    = Sphere_dir - SUM(Sphere_dir*Sphere_Axis)*Sphere_Axis ! ensure orthogonality between dir and axis
    IF(SQRT(DOT_PRODUCT(Sphere_dir,Sphere_dir)).LT.1e-9) THEN
      SWRITE(UNIT_stdOut,'(A)') 'Check definitions: Sphere_dir seems to be parallel to Sphere_Axis!!!'; STOP
    END IF
    Sphere_dir   = Sphere_dir/SQRT(DOT_PRODUCT(Sphere_dir,Sphere_dir))
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
  DO iPlane=nFlatplanes+nSphericPlanes+1,nFlatplanes+nSphericPlanes+nBLPlanes
    Plane=>Planes(iPlane)
    Plane%GroupID=GETINT('BLPlane_GroupID')
    WRITE(Plane%Name,'(A5,I6.6)')'BLPlane_',iPlane
    Plane%nRP(1:2)   =GETINTARRAY('BLPlane_nRP',2)
    ALLOCATE(Plane%RP_ptr(1:Plane%nRP(1),1:Plane%nRP(2)))
    nCP   =GETINT('BLPlane_nCP') ! points to define spline, at least two
    fac   =GETREAL('BLPlane_fac','1.')  ! growth factor of the BL mesh
    ALLOCATE(xCP(3,nCP,1))
    ALLOCATE(height(nCP,1))
    DO iCP=1,nCP
       xCP(:,iCP,1) =GETREALARRAY('BLPlane_CP',3)
       height(iCP,1)=GETREAL('BLPlane_height') ! height of the BL mesh
    END DO! iCP=1,nCP
    ! get coordinates of the rps and allocate pointers
    CALL GetBLPlane(Plane,nCP,height(:,1),fac,xCP(:,:,1))
    DEALLOCATE(xCP,height)
  END DO! iPlane
END IF

! ----------------------------------------------------------------------------------------------------
! Boxes
! ----------------------------------------------------------------------------------------------------
nFlatBoxes    =CountOption('Box_GroupID')
nBLBoxes      =CountOption('BLBox_GroupID')
nBoxes=nFlatBoxes+nBLBoxes
IF(nBoxes.GT.0) THEN
  anythingThere=.TRUE.
  ALLOCATE(Boxes(1:nBoxes))
! flat boxes
  DO iBox=1,nFlatBoxes
    Box=>Boxes(iBox)
    Box%GroupID=GETINT('Box_GroupID')
    WRITE(Box%Name,'(A4,I6.6)')'Box_',iBox
    x_dummy(1:24) = GETREALARRAY('Box_CornerX',24)
    DO iP=1,8
      Box%x(1:3,iP)=x_dummy(1+3*(iP-1):3+3*(iP-1))
    END DO ! iPoint
    Box%nRP(1:3)   =GETINTARRAY('Box_nRP',3)
    ALLOCATE(Box%RP_ptr(1:Box%nRP(1),1:Box%nRP(2),1:Box%nRP(3)))
    DO k=1,Box%nRP(3)
      DO j=1,Box%nRP(2)
        DO i=1,Box%nRP(1)
          ! bilinear mapping
          xi  = REAL(i-1)/REAL(Box%nRP(1)-1)
          eta = REAL(j-1)/REAL(Box%nRP(2)-1)
          zeta= REAL(k-1)/REAL(Box%nRP(3)-1)
          x(1:3)=   Box%x(1:3,1) * (1.-xi) * (1.-eta) * (1.-zeta)&
                 +  Box%x(1:3,2) * (   xi) * (1.-eta) * (1.-zeta)&
                 +  Box%x(1:3,3) * (   xi) * (   eta) * (1.-zeta)&
                 +  Box%x(1:3,4) * (1.-xi) * (   eta) * (1.-zeta)&
                 +  Box%x(1:3,5) * (1.-xi) * (1.-eta) * (   zeta)&
                 +  Box%x(1:3,6) * (   xi) * (1.-eta) * (   zeta)&
                 +  Box%x(1:3,7) * (   xi) * (   eta) * (   zeta)&
                 +  Box%x(1:3,8) * (1.-xi) * (   eta) * (   zeta)
          CALL GetNewRP(Box%RP_ptr(i,j,k)%RP,Box%GroupID,x)
        END DO ! i
      END DO ! j
    END DO
  END DO ! iBox

! BL Box
  DO iBox=nFlatBoxes+1,nFlatBoxes+nBLBoxes
    Box=>Boxes(iBox)
    Box%GroupID=GETINT('BLBox_GroupID')
    WRITE(Box%Name,'(A5,I6.6)')'BLBox_',iBox
    Box%nRP(1:3)   =GETINTARRAY('BLBox_nRP',3)
    ALLOCATE(Box%RP_ptr(1:Box%nRP(1),1:Box%nRP(2),1:Box%nRP(3)))
    nCP   =GETINT('BLBox_nCP') ! points to define spline, at least two
    nSP   =GETINT('BLBox_nSP') ! number of splines in z
    fac   =GETREAL('BLBox_fac','1.')  ! growth factor of the BL mesh
    ALLOCATE(xCP(3,nCP,nSP))
    ALLOCATE(height(nCP,nSP))
    DO iSP=1,nSP
      DO iCP=1,nCP
        xCP(:,iCP,iSP) =GETREALARRAY('BLBox_CP',3)
        height(iCP,iSP)=GETREAL('BLBox_height') ! height of the BL mesh
      END DO! iCP=1,nCP
    END DO! iSP=1,nSP
    ! get coordinates of the rps and allocate pointers
    CALL GetBLBox(Box,nCP,nSP,height,fac,xCP)
    DEALLOCATE(xCP,height)
  END DO! iBox

END IF


IF(.NOT.anythingThere) THEN
  SWRITE(UNIT_stdOut,*) 'No RP infos specified in parameter file, exiting...'
  CALL Abort(__STAMP__,'Code stopped!')
END IF

! Create global RP array
ALLOCATE(RPlist(nRP_global))
iRP=0

! fill Boxes
IF(nBoxes.GT.0) THEN
  DO iBox=1,nBoxes
    Box=>Boxes(iBox)
    DO k=1,Box%nRP(3)
      DO j=1,Box%nRP(2)
        DO i=1,Box%nRP(1)
        iRP=iRP+1
        RPlist(iRP)%RP=>Boxes(iBox)%RP_ptr(i,j,k)%RP
        END DO ! i
      END DO ! j
    END DO ! k
  END DO !iBox
END IF

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
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitRPSet

END MODULE MOD_RPSet
