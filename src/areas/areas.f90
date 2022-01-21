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

!==================================================================================================================================
!> Subroutines for creating arbitrary areas in FLEXI and flaging elements inside the area
!==================================================================================================================================
MODULE MOD_Areas
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitArea
  MODULE PROCEDURE InitArea
END INTERFACE

INTERFACE FinalizeArea
  MODULE PROCEDURE FinalizeArea
END INTERFACE

PUBLIC::InitArea,FinalizeArea,PointInPoly
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> \brief Init the area
!===================================================================================================================================
SUBROUTINE InitArea(AreaStr,locArea,AreaShape)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Areas_Vars
USE MOD_ReadInTools
USE MOD_Mesh_Vars,           ONLY: nElems,Elem_xGP
USE MOD_Interpolation_Vars,  ONLY: NodeType
#if USE_MPI
USE MOD_Mesh_Readin,         ONLY: ELEMIPROC
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=*),INTENT(IN)                 :: AreaStr         !< Area name
TYPE(tArea),INTENT(INOUT),TARGET            :: locArea         !< Area type, will be filled in the routine
INTEGER,INTENT(IN)                          :: AreaShape
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tAreaList),ALLOCATABLE             :: AreasTmp(:)   ! Temporary list of all active areas
TYPE(tShape),POINTER                    :: locShape
LOGICAL                                 :: applyArea(nElems)
CHARACTER(LEN=255)                      :: StrTmp
INTEGER                                 :: i,j,k,iAreaElem,iElem,iVertex
REAL,DIMENSION(  0:PP_N,0:PP_N,0:PP_NZ) :: x_star
REAL                                    :: r_vec(PP_dim),c_pt(PP_dim),dist
#if PP_dim==3
REAL                                    :: t_min
#endif
LOGICAL                                 :: applyPolygonArea
!===================================================================================================================================
! A new area is created, increase the number of considered areas and add to array
IF (postiMode) nAreas = 0
ALLOCATE(AreasTmp(nAreas+1))
IF (nAreas.GT.0) THEN
  AreasTmp(1:nAreas) = Areas(:)
END IF
nAreas = nAreas + 1
AreasTmp(nAreas)%pArea => locArea
SDEALLOCATE(Areas)
ALLOCATE(Areas(nAreas))
Areas = AreasTmp
DEALLOCATE(AreasTmp)

! set name for area
locArea%AreaStr = AreaStr
SWRITE(UNIT_stdOut,'(A)')' Building '//TRIM(AreaStr)//' area!'

!==================== Get information on the Shape and flag elements ===================!

! set shape
locArea%AreaShape = AreaShape
locShape => locArea%Shape

! Readin of geometrical parameters for different Area shapes
SELECT CASE(locArea%AreaShape)
  CASE(SHAPE_REGION) ! ramp aligned with a vector
    locShape%Vec(:)    = GETREALARRAY(TRIM(locArea%AreaStr)//'Dir',3,'(/1.,0.,0./)')
    locShape%xStart(:) = GETREALARRAY(TRIM(locArea%AreaStr)//'XStart',3,'(/0.,0.,0./)')
    ! FLEXI normally stores REAL with KIND8. However, we must be able to perform arithmetic operations on locShape%xEnd, so only use KIND4
    WRITE(StrTmp,'(A,E15.7,A,E15.7,A,E15.7,A)') '(/',HUGE(REAL(1.,KIND=4)),',',HUGE(REAL(1.,KIND=4)),',',HUGE(REAL(1.,KIND=4)),'/)'
    locShape%xEnd(:)   = GETREALARRAY(TRIM(locArea%AreaStr)//'XEnd',3,StrTmp)
#if PP_dim==2
    IF(locShape%Vec(3).NE.0) &
      CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set '//TRIM(locArea%AreaStr)//'Dir'//'(3) = 0!')
#endif
    locShape%Vec(:) = locShape%Vec(:)/SQRT(DOT_PRODUCT(locShape%Vec(:),locShape%Vec(:))) ! Normalize locShape%Vec

  CASE(SHAPE_CUBOID_CARTESIAN)
    locShape%xStart(:)  = GETREALARRAY(TRIM(locArea%AreaStr)//'XStart',3,'(/0.,0.,0./)')
    locShape%xEnd (:)   = GETREALARRAY(TRIM(locArea%AreaStr)//'XEnd',  3,'(/0.,0.,0./)')
    locShape%xCenter(:) = (locShape%xStart(1:PP_dim)+locShape%xEnd(1:PP_dim))/2.

  CASE(SHAPE_CUBE_CARTESIAN) ! cube
    locShape%xCenter(:) = GETREALARRAY(TRIM(locArea%AreaStr)//'XCenter',3,'(/0.,0.,0./)')
    locShape%Radius     = GETREAL(TRIM(locArea%AreaStr)//'Radius')
    locShape%xStart(:)  = locShape%xCenter(:)-locShape%Radius
    locShape%xEnd(:)    = locShape%xCenter(:)+locShape%Radius

  CASE(SHAPE_CYLINDRICAL_OUTER) ! circular Area
    IF (locShape%Radius.EQ.0.) &
      locShape%Radius     = GETREAL(TRIM(locArea%AreaStr)//'Radius')
#if PP_dim==3
    locShape%xStart(:)  = GETREALARRAY(TRIM(locArea%AreaStr)//'XStart',3,'(/0.,0.,0./)')
    locShape%Axis(:)    = GETREALARRAY(TRIM(locArea%AreaStr)//'Axis',3,'(/0.,0.,1./)')
    locShape%xEnd(:)    = locShape%xStart(:)+locShape%Axis(:)
    locShape%xCenter(:) = (locShape%xStart(:)+locShape%xEnd(:))/2.
#else
    locShape%xCenter(:) = GETREALARRAY(TRIM(locArea%AreaStr)//'XCenter',3,'(/0.,0.,0./)')
#endif

  CASE(SHAPE_CYLINDRICAL_INNER) ! circular Area
    IF (locShape%Radius.EQ.0.) &
      locShape%Radius     = GETREAL(TRIM(locArea%AreaStr)//'Radius')
#if PP_dim==3
    locShape%xStart(:)  = GETREALARRAY(TRIM(locArea%AreaStr)//'XStart',3,'(/0.,0.,0./)')
    locShape%xEnd(:)    = GETREALARRAY(TRIM(locArea%AreaStr)//'XEnd',3,'(/0.,0.,0./)')
    locShape%Axis(:)    = locShape%xEnd(:)-locShape%xStart(:)
    locShape%xCenter(:) = (locShape%xStart(:)+locShape%xEnd(:))/2.
#else
    locShape% Center(:) = GETREALARRAY(TRIM(locArea%AreaStr)//'XCenter',3,'(/0.,0.,0./)')
#endif

  CASE(SHAPE_SPHERE) ! sphere
    IF (locShape%Radius.EQ.0.) &
      locShape%Radius     = GETREAL(TRIM(locArea%AreaStr)//'Radius')
    locShape%xCenter(:) = GETREALARRAY(TRIM(locArea%AreaStr)//'XCenter',3,'(/0.,0.,0./)')

  CASE(SHAPE_POLYGON) ! Polygone Area
    locShape%nAreaVertices = GETINT('n'//TRIM(locArea%AreaStr)//'Vertices')
END SELECT

ALLOCATE(locShape%AreaVertex(locShape%nAreaVertices,3))
locShape%AreaVertex = 0.
SELECT CASE(locArea%AreaShape)
  CASE(SHAPE_POLYGON) ! Polygone Area
    DO iVertex=1,locShape%nAreaVertices
      locShape%AreaVertex(iVertex,:) = GETREALARRAY(TRIM(locArea%AreaStr)//'Vertex',3,'(/0.,0.,0./)')
    END DO
END SELECT

applyArea = .FALSE.
DO iElem=1,nElems
  x_star=0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    SELECT CASE(locArea%AreaShape)
      CASE(SHAPE_REGION) ! ramp aligned with a vector
        ! Region between xStart and xEnd
        IF (SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-locShape%xStart  (1:PP_dim))   *locShape%Vec(1:PP_dim)).GE.0 .AND. &
            SUM((locShape%xEnd    (1:PP_dim)   -Elem_xGP(1:PP_dim,i,j,k,iElem))*locShape%Vec(1:PP_dim)).GE.0) THEN
              x_star(i,j,k) = SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-locShape%xStart(1:PP_dim))*locShape%Vec(1:PP_dim))
        END IF

      CASE(SHAPE_CUBOID_CARTESIAN) ! cuboid cartesian aligned defined by two points
        IF(ABS(Elem_xGP(1,i,j,k,iElem)-locShape%xCenter(1)).LT.ABS(locShape%xStart(1)-locShape%xCenter(1))) THEN
          IF(ABS(Elem_xGP(2,i,j,k,iElem)-locShape%xCenter(2)).LT.ABS(locShape%xStart(2)-locShape%xCenter(2))) THEN
            IF(ABS(Elem_xGP(3,i,j,k,iElem)-locShape%xCenter(3)).LT.ABS(locShape%xStart(3)-locShape%xCenter(3))) THEN
              x_star(i,j,k) = 1.
            END IF
          END IF
        END IF

      CASE(SHAPE_CUBE_CARTESIAN) ! cartesian cube
        IF(ABS(Elem_xGP(1,i,j,k,iElem)-locShape%xCenter(1)).LT.ABS(locShape%Radius-locShape%xCenter(1))) THEN
          IF(ABS(Elem_xGP(2,i,j,k,iElem)-locShape%xCenter(2)).LT.ABS(locShape%Radius-locShape%xCenter(2))) THEN
            IF(ABS(Elem_xGP(3,i,j,k,iElem)-locShape%xCenter(3)).LT.ABS(locShape%Radius-locShape%xCenter(3))) THEN
              x_star(i,j,k) = 1.
            END IF
          END IF
        END IF

!       CASE(SHAPE_CYLINDRICAL_OUTER) ! cylindrical outer area
!         r_vec(:) = Elem_xGP(:,i,j,k,iElem)-locShape%xStart(1:PP_dim)
! #if(PP_dim==3)
!         r_vec    = r_vec - SUM((Elem_xGP(:,i,j,k,iElem)-locShape%xStart(:))*locShape%Axis(:))*locShape%Axis(:)
! #endif
!         x_star(i,j,k) = (SQRT(SUM(r_vec*r_vec))-locShape%Radius)

      CASE(SHAPE_CYLINDRICAL_INNER,SHAPE_CYLINDRICAL_OUTER) ! cylindrical inner area
#if(PP_dim==3)
        r_vec = locShape%xEnd(1:PP_dim)-locShape%xStart(1:PP_dim)
        t_min = (DOT_PRODUCT(Elem_xGP(1:PP_dim,i,j,k,iElem),r_vec)-DOT_PRODUCT(locShape%xStart(1:PP_dim),r_vec))/DOT_PRODUCT(r_vec,r_vec)
        c_pt  = locShape%xStart(1:PP_dim) + t_min*r_vec
#else
        c_pt  = locShape%xCenter(1:PP_dim)
#endif
        dist  = SQRT(DOT_PRODUCT(c_pt-Elem_xGP(1:PP_dim,i,j,k,iElem),c_pt-Elem_xGP(1:PP_dim,i,j,k,iElem)))
        SELECT CASE(SHAPE_REGION)
          CASE(SHAPE_CYLINDRICAL_INNER)
            IF (dist.LT.locShape%Radius) THEN
              x_star(i,j,k) = 1.
            END IF
          CASE(SHAPE_CYLINDRICAL_OUTER)
            IF (dist.GT.locShape%Radius) THEN
              x_star(i,j,k) = 1.
            END IF
        END SELECT

      CASE(SHAPE_SPHERE) ! sphere
        IF(((Elem_xGP(1,i,j,k,iElem)-locShape%xCenter(1))**2.+(Elem_xGP(2,i,j,k,iElem)-locShape%xCenter(2))**2. + &
            (Elem_xGP(3,i,j,k,iElem)-locShape%xCenter(3))**2.).LT.locShape%Radius**2.) THEN
          x_star(i,j,k) = 1.
        END IF

      CASE(SHAPE_POLYGON)
        applyPolygonArea = .FALSE.
        CALL PointInPoly(Elem_xGP(1,i,j,k,iElem),Elem_xGP(2,i,j,k,iElem),locShape%AreaVertex(1:locShape%nAreaVertices,1), &
                         locShape%AreaVertex(1:locShape%nAreaVertices,2),locShape%nAreaVertices,applyPolygonArea)
        IF (applyPolygonArea) x_star(i,j,k) = 1.
    END SELECT
  END DO; END DO; END DO

  IF (ANY(x_star.GT.0.)) THEN
    applyArea(iElem) = .TRUE.
    CYCLE
  END IF
END DO ! iElem = 1,nElems

! Get Area count and build Area mapping
locArea%nAreaElems = COUNT(applyArea)
ALLOCATE(locArea%AreaMap(locArea%nAreaElems))
iAreaElem = 0
DO iElem = 1,nElems
  IF (applyArea(iElem)) THEN
    iAreaElem = iAreaElem+1
    locArea%AreaMap(iAreaElem) = iElem
  END IF
END DO

SWRITE(UNIT_stdOut,'(A)')' Building '//TRIM(AreaStr)//' area done!'

END SUBROUTINE InitArea

!----------------------------------------------------------------------------------------------------------------------------------!
! Check if a point is inside a polygon(PolyX,PolyY)
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE PointInPoly(PointX,PointY,PolyX,PolyY,PolyN,Inside)
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)         :: PointX,PointY,PolyX(PolyN),PolyY(PolyN)
INTEGER,INTENT(IN)      :: PolyN
LOGICAL,INTENT(INOUT)   :: Inside
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,InOrOut
REAL                    :: xj,yj,xi,yi
LOGICAL                 :: ix , iy , jx , jy , EOR
!===================================================================================================================================

! EXCLUSIVE OR STATEMENT FUNCTION.
EOR(ix,iy) = (ix .OR. iy) .AND. .NOT.(ix .AND. iy)

Inside = .FALSE.
InOrOut = -1

DO i=1,PolyN
  xi = PolyX(i) - PointX
  yi = PolyY(i) - PointY

  ! CHECK WHETHER THE POINT IN QUESTION IS AT THIS VERTEX.
  IF ( xi.EQ.0.0 .AND. yi.EQ.0.0 ) THEN
     InOrOut = 0
     RETURN
  ENDIF

  ! j IS NEXT VERTEX NUMBER OF POLYGON.
  j = 1 + MOD(i,PolyN) !Check if last point in polygon
  xj = PolyX(j) - PointX
  yj = PolyY(j) - PointY

  ! IS THIS LINE OF 0 LENGTH ?
  IF ( xi.EQ.xj .AND. yi.EQ.yj ) CYCLE
  ix = (xi.GE.0.0)
  iy = (yi.GE.0.0)
  jx = (xj.GE.0.0)
  jy = (yj.GE.0.0)

  ! CHECK WHETHER (PointX,PointY) IS ON VERTICAL SIDE OF POLYGON.
  IF ( xi.EQ.0.0 .AND. xj.EQ.0.0 .AND. EOR(iy,jy) ) THEN
    InOrOut = 0
    RETURN
  ENDIF
  ! CHECK WHETHER (PointX,PointY) IS ON HORIZONTAL SIDE OF POLYGON.
  IF ( yi.EQ.0.0 .AND. yj.EQ.0.0 .AND. EOR(ix,jx) ) THEN
    InOrOut = 0
    RETURN
  ENDIF

  ! CHECK WHETHER BOTH ENDS OF THIS SIDE ARE COMPLETELY 1) TO RIGHT
  ! OF, 2) TO LEFT OF, OR 3) BELOW (PointX,PointY).
  IF ( .NOT.((iy .OR. jy) .AND. EOR(ix,jx)) ) CYCLE

  ! DOES THIS SIDE OBVIOUSLY CROSS LINE RISING VERTICALLY FROM (PointX,PointY)
  IF ( .NOT.(iy .AND. jy .AND. EOR(ix,jx)) ) THEN
    IF ( (yi*xj-xi*yj)/(xj-xi).LT.0.0 ) THEN
      CYCLE
    ELSEIF ( (yi*xj-xi*yj)/(xj-xi).EQ.0.0 ) THEN
      InOrOut = 0
      RETURN
    ELSE
      InOrOut = -InOrOut
    ENDIF
  ELSE
    InOrOut = -InOrOut
  ENDIF
END DO
IF (InOrOut .GE. 0) Inside = .TRUE.

END SUBROUTINE PointInPoly

!===================================================================================================================================
!> \brief Finalize the area
!>
!===================================================================================================================================
SUBROUTINE FinalizeArea(locArea)
! MODULES
USE MOD_Globals
USE MOD_Areas_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(tArea),INTENT(INOUT),TARGET            :: locArea          !< area type, will be filled in the routine
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(locArea%AreaMap)
SDEALLOCATE(locArea%Shape%AreaVertex)
SDEALLOCATE(Areas)

END SUBROUTINE FinalizeArea

END MODULE MOD_Areas
