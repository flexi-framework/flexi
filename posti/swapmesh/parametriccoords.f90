!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
!> Module to evaluate parametric coordinates of interpolation points (IP) of new state in old mesh
!===================================================================================================================================
MODULE MOD_SMParametricCoordinates
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC:: GetParametricCoordinates

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Routine to evaluate parametric coordinates of interpolation points (IP) of new state in old mesh.
!> In a first step, the centroids and radii of old and new mesh elements are computed. We only consider elements for the search
!> when they overlap (within a certain tolerance).
!> Then a check for equal elements is performed, they can simply be copied from the old mesh.
!> For non-equal elements the procedure is as follows:
!> Employ a Newton algorithm to find the parametric coordinates xi of an interpolation point in an element defined by CL points on
!> NGeo. We try to solve F(xi) = x(xi) - x_InterpolationPoint = 0 for the parametric coordinates.
!> Newton iteration: xi_(n+1) = xi_n - (J(xi_n))^(-1)*F(xi_n), the Jacobian is the derivative of the mesh coordinates w.r.t. the
!> parametric coordinates.
!> At the end, a check is performed if all new interpolation points have been found. A certain range outside of [-1,1] is allowed.
!===================================================================================================================================
SUBROUTINE GetParametricCoordinates()
! MODULES
USE MOD_Globals
USE MOD_Basis,                 ONLY: ChebyGaussLobNodesAndWeights,BarycentricWeights
USE MOD_ChangeBasisByDim,      ONLY: ChangeBasisVolume
USE MOD_Interpolation,         ONLY: GetVandermonde,GetDerivativeMatrix
USE MOD_Interpolation_Vars,    ONLY: NodeTypeCL
USE MOD_Newton,                ONLY: Newton
USE MOD_Output,                ONLY: PrintPercentage
USE MOD_SwapMesh_Vars
#if USE_OPENMP
USE OMP_Lib
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k     ! CL,NSuper
INTEGER            :: ii,jj,kk  ! IP
INTEGER            :: iElemNew,iElemOld
INTEGER            :: l
LOGICAL            :: equal
INTEGER            :: ElemCounter,nEqualElems,nNotfound
! Element search
LOGICAL            :: ElemDone(   nElemsNew)
LOGICAL            :: ElemDoneOld(nElemsOld)
REAL               :: xCOld(PP_dim,nElemsOld),radOld(nElemsOld),radSqOld(nElemsOld)
REAL               :: xCNew(PP_dim,nElemsNew),radNew(nElemsNew),radSqNew(nElemsNew)
! Point search
REAL               :: dist,maxDist,best
REAL               :: DCL_NGeo(    0:NGeoOld,0:NGeoOld)
REAL               :: Xi_CLNGeo(   0:NGeoOld)
REAL               :: wBary_CLNGeo(0:NGeoOld)
REAL               :: dXCL_NGeo(1:PP_dim,1:PP_dim,0:NGeoOld,0:NGeoOld,0:ZDIM(NGeoOld))
! Interpolation/supersampling
LOGICAL            :: IPOverlaps(       0:NInter,0:NInter,0:ZDIM(NInter))
REAL               :: X_NSuper(1:PP_dim,0:NSuper,0:NSuper,0:ZDIM(NSuper))
REAL               :: Xi_NSuper(0:NSuper)
REAL               :: xInter(PP_dim)
REAL               :: xi(1:PP_dim)
! Progress indicator
REAL               :: StartT,EndT,percent
!===================================================================================================================================
StartT = OMP_FLEXITIME()
SWRITE(UNIT_stdOut,'(A)') ' EVALUATING PARAMETRIC COORDINATES...'

! Prepare CL basis evaluation
CALL ChebyGaussLobNodesAndWeights(NGeoOld,Xi_CLNGeo)
CALL BarycentricWeights(NGeoOld,Xi_CLNGeo,wBary_CLNGeo)
CALL GetDerivativeMatrix(NGeoOld,NodeTypeCL,DCL_NGeo)

! Equidistant points (including the edges) for supersampling
DO i=0,NSuper
  Xi_NSuper(i) = 2./REAL(NSuper) * REAL(i) - 1.
END DO

! Compute centroids and radii for the old mesh on CL points and the interpolation mesh (new mesh on NInter)
DO iElemOld=1,nElemsOld
 CALL getCentroidAndRadius(xCLOld(:,:,:,:,iElemOld),NGeoOld,xCOld(:,iElemOld),radOld(iElemOld))
END DO
DO iElemNew=1,nElemsNew
 CALL getCentroidAndRadius(xCLInter(:,:,:,:,iElemNew),NInter, xCNew(:,iElemNew),radNew(iElemNew))
END DO
! add 5 % tolerance to old
radOld   = radOld*1.05
! Precompute square of radii
radSqOld = radOld*radOld
radSqNew = radNew*radNew

xiInter     = HUGE(1.)
InterToElem = -999

IPDone      = .FALSE.       ! Mark single interpolation point as done
ElemDone    = .FALSE.       ! Mark new element as done
ElemDoneOld = .FALSE.       ! Mark old element as done
ElemCounter = 0
nEqualElems = 0
equalElem   = -999

! look for identical elements and mark them
maxDist = 1e-9
! only for same Ngeo so far
IF(NgeoOld.EQ.NGeoNew) THEN
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(iElemOld,equal,iElemNew,i,j,k,dist) SCHEDULE(DYNAMIC,100)
  DO iElemOld = 1,nElemsOld
    DO iElemNew = 1,nElemsNew
      IF (ElemDone(iElemNew)) CYCLE

      equal = .FALSE.
      ! check each dimension first only before getting 2Norm
      ! we compare the distance of the centroids here for a first check!
      IF(ABS(  xCOld(1,iElemOld)-xCNew(1,iElemNew)).GT.maxDist) CYCLE
      IF(ABS(  xCOld(2,iElemOld)-xCNew(2,iElemNew)).GT.maxDist) CYCLE
#if PP_dim == 3
      IF(ABS(  xCOld(3,iElemOld)-xCNew(3,iElemNew)).GT.maxDist) CYCLE
#endif
      dist = SUM((xCOld(1:PP_dim,iElemOld)-xCNew(1:PP_dim,iElemNew))**2)

      ! Element definitely out of range
      IF(dist.GT.maxDist**2) CYCLE

      ! Now the coordinates of the mesh nodes themselves are compared
      equal = .TRUE.
      DO k = 0,ZDIM(NGeoOld); DO j = 0,NGeoOld; DO i = 0,NGeoOld
        dist = SUM((xCLOld(1:PP_dim,i,j,k,iElemOld) - xCLNew(1:PP_dim,i,j,k,iElemNew))**2)
        IF (dist.GT.maxDist**2) equal = .FALSE.
      END DO; END DO; END DO

      IF (equal) THEN
!$OMP CRITICAL
        equalElem(  iElemNew) = iElemOld
        ElemDone(   iElemNew) = .TRUE.
        ElemDoneOld(iElemOld) = .TRUE.
        nEqualElems = nEqualElems+1
!$OMP END CRITICAL
      END IF ! equal
    END DO! iElemNew=1,nElemsNew
  END DO! iElemOld=1,nElemsOld
!$OMP END DO
!$OMP END PARALLEL
  SWRITE(UNIT_stdOut,'(I12,A,I12,A)') nEqualElems,' equal elements, ',nElemsOld-nEqualElems,' remaining.'
END IF ! (NgeoOld.EQ.NGeoNew)


! If extrusion is done, only elements with k=1 are checked here. All others are marked as already done
IF (ExtrudeTo3D) THEN
  DO iElemNew=1,nElemsNew
    IF (Elem_IJK(3,iElemNew).NE.ExtrudeK) THEN
      ElemDone(iElemNew) =.TRUE.
      IPDone(:,:,:,iElemNew) = .TRUE.
    END IF
  END DO
END IF

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(iElemOld,iElemNew,X_NSuper,dXCL_NGeo,IPOverlaps), &
!$OMP& PRIVATE(ii,jj,kk,xInter,best,i,j,k,l,dist,xi,maxDist) SCHEDULE(DYNAMIC,100)
DO iElemOld=1,nElemsOld
  ! Already found matching new element
  IF (ElemDoneOld(iElemOld)) CYCLE

  ! Find initial guess for Newton, get supersampled element
  CALL ChangeBasisVolume(PP_dim,NGeoOld,NSuper,Vdm_CLNGeo_EquiNSuper,xCLOld(1:PP_dim,:,:,:,iElemOld),X_NSuper)

  ! Compute Jacobian of element mapping for each CL point
  dXCL_NGeo=0.
  DO k=0,ZDIM(NGeoOld); DO j=0,NGeoOld; DO i=0,NGeoOld
    DO l=0,NGeoOld
      ! Matrix-vector multiplication
      dXCL_NGeo(1:PP_dim,1,i,j,k)=dXCL_NGeo(1:PP_dim,1,i,j,k) + DCL_NGeo(i,l)*xCLOld(1:PP_dim,l,j,k,iElemOld)
      dXCL_NGeo(1:PP_dim,2,i,j,k)=dXCL_NGeo(1:PP_dim,2,i,j,k) + DCL_NGeo(j,l)*xCLOld(1:PP_dim,i,l,k,iElemOld)
#if PP_dim == 3
      dXCL_NGeo(:,3,i,j,k)=dXCL_NGeo(:,3,i,j,k) + DCL_NGeo(k,l)*xCLOld(:,i,j,l,iElemOld)
#endif
    END DO
  END DO; END DO; END DO

  DO iElemNew=1,nElemsNew
    ! All IP already found
    IF (ElemDone(iElemNew)) CYCLE

    ! Check if the two elements could be overlapping by comparing the distance between the centroids with the sum of the radii
    maxDist=radNew(iElemNew)+radOld(iElemOld)
    IF(SUM((xCOld(1:PP_dim,iElemOld)-xCNew(1:PP_dim,iElemNew))**2).GT.maxDist**2) CYCLE

    ! Check which IP of new elem overlaps with old elem
    IPOverlaps=.FALSE.
    DO kk=0,ZDIM(NInter); DO jj=0,NInter; DO ii=0,NInter
      ! Already found matching interpolation point
      IF(IPDone(ii,jj,kk,iElemNew)) CYCLE

      dist=SUM((xCOld(1:PP_dim,iElemOld)-xCLInter(1:PP_dim,ii,jj,kk,iElemNew))**2)
      IPOverlaps(ii,jj,kk)=(dist.LE.radSqOld(iElemOld))
    END DO; END DO; END DO

    ! Get smallest distance to supersampled points for starting Newton
    DO kk=0,ZDIM(NInter); DO jj=0,NInter; DO ii=0,NInter
      IF(.NOT.IPOverlaps(ii,jj,kk)) CYCLE

      ! Get initial guess for Newton by finding nearest supersampled point
      xInter = xCLInter(1:PP_dim,ii,jj,kk,iElemNew)
      best=HUGE(1.)
      DO i=0,NSuper; DO j=0,NSuper; DO k=0,ZDIM(NSuper)
        dist=SUM((xInter-X_NSuper(:,i,j,k))**2)
        IF (dist.LT.best) THEN
          best=dist
#if PP_dim == 3
          xi=(/Xi_NSuper(i),Xi_NSuper(j),Xi_NSuper(k)/)
#else
          xi=(/Xi_NSuper(i),Xi_NSuper(j)/)
#endif
        END IF
      END DO; END DO; END DO

      ! Find coordinates in reference space with Newton starting from initial guess
#if PP_dim == 3
      CALL Newton(NGeoOld,XInter,dXCL_NGeo(1:PP_dim,1:PP_dim,:,:,:),Xi_CLNGeo,wBary_CLNGeo,xCLOld(1:PP_dim,:,:,:,iElemOld),xi)
#else
      CALL Newton(NGeoOld,XInter,dXCL_NGeo(1:PP_dim,1:PP_dim,:,:,0),Xi_CLNGeo,wBary_CLNGeo,xCLOld(1:PP_dim,:,:,0,iElemOld),xi)
#endif /*PP_dim == 3*/

      ! Check if result is better than previous result
      ! (Result might be outside of element but still within user-specified tolerance)
!$OMP CRITICAL
      IF(MAXVAL(ABS(Xi)).LT.MAXVAL(ABS(xiInter(:,ii,jj,kk,iElemNew)))) THEN
        ! If point is inside element, stop searching
        IF (MAXVAL(ABS(Xi)).LE.1.) IPDone(ii,jj,kk,iElemNew) = .TRUE.
        xiInter(1:PP_dim,ii,jj,kk,iElemNew)=Xi
        InterToElem(ii,jj,kk,iElemNew) = iElemOld
      END IF
!$OMP END CRITICAL
    END DO; END DO; END DO ! ii,jj,kk (IP loop)
    ! Mark Element as done, if all its interpolation points are found
    IF(ALL(IPDone(:,:,:,iElemNew))) ElemDone(iElemNew)=.TRUE.
  END DO ! iElem
!$OMP CRITICAL
  ElemCounter=ElemCounter+1

  ! Print progress
  percent = REAL(ElemCounter)/REAL(nElemsNew)*100.
  CALL PrintPercentage(NameOpt='Searching coordinates',percent=percent)

!$OMP END CRITICAL
END DO ! iElemOld
!$OMP END DO
!$OMP END PARALLEL

! check if all interpolation points have been found
nNotFound = 0
DO iElemNew = 1,nElemsNew
  IF (equalElem(iElemNew).GT.0) CYCLE

  DO ii=0,NInter; DO jj=0,NInter; DO kk=0,ZDIM(NInter)
    IF (IPDone(ii,jj,kk,iElemNew)) CYCLE

    ! Only mark as invalid if greater then max tolerance
    IF (MAXVAL(ABS(xiInter(:,ii,jj,kk,iElemNew))).GT.maxTol) THEN
      ! RP has not been found
      IF(printTroublemakers)THEN
        SWRITE(UNIT_stdOut,*) 'IP with ID:',       ii,jj,kk,iElemNew,'is a troublemaker!'
        SWRITE(UNIT_stdOut,*) 'Coords:',xCLInter(:,ii,jj,kk,iElemNew)
        SWRITE(UNIT_stdOut,*) 'Xi:'    ,xiInter( :,ii,jj,kk,iElemNew)
      END IF

      ! Abort if severly broken and no refstate present (with refstate present abortTol=HUGE)
      IF(MAXVAL(ABS(xiInter(1:PP_dim,ii,jj,kk,iElemNew))).GT.abortTol)&
        CALL Abort(__STAMP__, 'IP not found.')
      nNotFound = nNotFound+1
    ELSE
      IPDone(ii,jj,kk,iElemNew) = .TRUE.
    END IF
  END DO; END DO; END DO ! ii,jj,kk (IP loop)
END DO
SWRITE(Unit_stdOut,'(I12,A)') nNotFound,' nodes not found.'

EndT = OMP_FLEXITIME()
CALL DisplayMessageAndTime(EndT-StartT,'DONE!',DisplayLine=.TRUE.)

END SUBROUTINE GetParametricCoordinates


!=================================================================================================================================
!> Computes the centroid and the radius of an element
!=================================================================================================================================
SUBROUTINE getCentroidAndRadius(elem,NGeo,xC,radius)
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: NGeo                               !> Polynomial degree of mesh representation
REAL,INTENT(IN)    :: elem(3,0:NGeo,0:NGeo,0:ZDIM(NGeo)) !> Coordinates of single element
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: xC(PP_dim)                   !> Coordinates of centroid
REAL,INTENT(OUT)   :: radius                       !> Radius of element
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,nNodes
LOGICAL            :: onSide(3)
!=================================================================================================================================
! Compute centroid of element
nNodes= (NGeo+1)**PP_dim
xC(1) = SUM(elem(1,:,:,:))/nNodes
xC(2) = SUM(elem(2,:,:,:))/nNodes
#if PP_dim == 3
xC(3) = SUM(elem(3,:,:,:))/nNodes
#endif

! Compute max distance from bary to surface nodes and return as radius
radius=0.
DO k=0,ZDIM(NGeo)
  onSide(3)=((k.EQ.0).OR.(k.EQ.NGeo))
  DO j=0,NGeo
    onSide(2)=((j.EQ.0).OR.(j.EQ.NGeo))
    DO i=0,NGeo
      onSide(1)=((i.EQ.0).OR.(i.EQ.NGeo))
      IF(.NOT.ANY(onSide(1:PP_dim))) CYCLE
      radius=MAX(radius,NORM2(elem(1:PP_dim,i,j,k)-xC(1:PP_dim)))
    END DO
  END DO
END DO

END SUBROUTINE getCentroidAndRadius

END MODULE MOD_SMParametricCoordinates
