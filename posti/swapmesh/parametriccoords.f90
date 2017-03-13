!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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
!> Module to evaluate parametric coordinates of interpolation points (IP) of new state in old mesh
!===================================================================================================================================
MODULE MOD_SMParametricCoordinates
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE GetParametricCoordinates
  MODULE PROCEDURE GetParametricCoordinates
END INTERFACE

PUBLIC :: GetParametricCoordinates

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
USE MOD_SwapMesh_Vars
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys,ChebyGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Interpolation,         ONLY: GetVandermonde,GetDerivativeMatrix
USE MOD_Interpolation_Vars,    ONLY: NodeTypeCL
USE MOD_ChangeBasis,           ONLY: ChangeBasis3D
USE MOD_Mathtools,             ONLY: INVERSE2,INVERSE3,DET3
use omp_lib
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
INTEGER            :: l,iter
REAL               :: xCOld(3,nElemsOld),radOld(nElemsOld),radSqOld(nElemsOld)
REAL               :: xCNew(3,nElemsNew),radNew(nElemsNew),radSqNew(nElemsNew)
LOGICAL            :: IPOverlaps(0:NInter,0:NInter,0:NInter)
LOGICAL            :: ElemDone(   nElemsNew)
LOGICAL            :: ElemDoneOld(nElemsOld)
LOGICAL            :: equal
REAL               :: X_NSuper(1:3,0:NSuper,0:NSuper,0:NSuper)
REAL               :: dist,maxDist,best
REAL               :: DCL_NGeo(0:NGeoOld,0:NGeoOld)
REAL               :: dXCL_NGeo(1:3,1:3,0:NGeoOld,0:NGeoOld,0:NGeoOld)
REAL               :: Xi_CLNGeo(0:NGeoOld),wBary_CLNGeo(0:NGeoOld)
REAL               :: Xi_NSuper(0:NSuper)
REAL               :: xInter(3)
REAL               :: LagXi(0:NGeoOld),LagEta(0:NGeoOld),LagZeta(0:NGeoOld),LagVol(0:NGeoOld,0:NGeoOld,0:NGeoOld)
REAL               :: F(1:3),eps_F,xi(1:3),Jac(1:3,1:3),sdetJac,sJac(1:3,1:3)
INTEGER            :: ElemCounter,nEqualElems
INTEGER            :: nNotfound
REAL               :: Time
!===================================================================================================================================
Time=FLEXITIME()
! Prepare CL basis evaluation
CALL ChebyGaussLobNodesAndWeights(NGeoOld,Xi_CLNGeo)
CALL BarycentricWeights(NGeoOld,Xi_CLNGeo,wBary_CLNGeo)
CALL GetDerivativeMatrix(NGeoOld,NodeTypeCL,DCL_NGeo)

! Equdistant points (including the edges) for supersampling
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
radOld=radOld*1.05
radNew=radNew
! Precompute square of radii
radSqOld=radOld*radOld
radSqNew=radNew*radNew

xiInter=HUGE(1.)
InterToElem=-999

IPDone=.FALSE.            ! Mark single interpolation point as done
ElemDone=.FALSE.          ! Mark new element as done
ElemDoneOld=.FALSE.       ! Mark old element as done
StartTime=FLEXITIME()
ElemCounter=0
nEqualElems=0
equalElem=-999
! look for identical elements and mark them 
! only for same Ngeo so far
IF(NgeoOld.EQ.NGeoNew) THEN
  maxDist=1e-9
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(iElemOld,equal,iElemNew,i,j,k,dist) SCHEDULE(DYNAMIC,100)
  DO iElemOld=1,nElemsOld
    !print*, OMP_get_thread_num()
    DO iElemNew=1,nElemsNew
      IF(ElemDone(iElemNew)) CYCLE
      equal=.FALSE.
      ! check each dimension first only before getting 2Norm
      ! we compare the distance of the centroids here for a first check!
      IF(ABS(  xCOld(1,iElemOld)-xCNew(1,iElemNew)).GT.maxDist) CYCLE
      IF(ABS(  xCOld(2,iElemOld)-xCNew(2,iElemNew)).GT.maxDist) CYCLE
      IF(ABS(  xCOld(3,iElemOld)-xCNew(3,iElemNew)).GT.maxDist) CYCLE
      dist=SUM((xCOld(:,iElemOld)-xCNew(:,iElemNew))**2)
      IF(dist.GT.maxDist**2) CYCLE
      ! now the coordinates of the mesh nodes themselves are compared
      equal=.TRUE.
      DO k=0,NGeoOld; DO j=0,NGeoOld; DO i=0,NGeoOld
        dist=SUM((xCLOld(:,i,j,k,iElemOld)-xCLNew(:,i,j,k,iElemNew))**2)
        IF(dist.GT.maxDist**2) THEN
          equal=.FALSE.
        END IF
      END DO; END DO; END DO
      IF(equal) THEN
!$OMP CRITICAL 
        equalElem(iElemNew)=iElemOld
        ElemDone(iElemNew)=.TRUE.
        ElemDoneOld(iElemOld)=.TRUE.
        nEqualElems=nEqualElems+1
!$OMP END CRITICAL 
      END IF ! equal
    END DO! iElemNew=1,nElemsNew
  END DO! iElemOld=1,nElemsOld
!$OMP END DO
!$OMP END PARALLEL
  WRITE(UNIT_StdOut,*)nEqualElems,' equal elements, ',nElemsOld-nEqualElems,' remaining.'
END IF ! (NgeoOld.EQ.NGeoNew)

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(iElemOld,iElemNew,X_NSuper,dXCL_NGeo,IPOverlaps), &
!$OMP& PRIVATE(ii,jj,kk,xInter,best,i,j,k,dist,xi,LagXi,LagEta,LagZeta,LagVol,F,iter,Jac,sdetJac,sJac) SCHEDULE(DYNAMIC,100)
DO iElemOld=1,nElemsOld
  IF(ElemDoneOld(iElemOld)) CYCLE ! already found matching new element

  ! Find initial guess for Newton, get supersampled element
  CALL ChangeBasis3D(3,NGeoOld,NSuper,Vdm_CLNGeo_EquiNSuper,xCLOld(:,:,:,:,iElemOld),X_NSuper)
  ! Compute Jacobian of element mapping for each CL point
  dXCL_NGeo=0.
  DO k=0,NGeoOld; DO j=0,NGeoOld; DO i=0,NGeoOld
    DO l=0,NGeoOld
      ! Matrix-vector multiplication
      dXCL_NGeo(:,1,i,j,k)=dXCL_NGeo(:,1,i,j,k) + DCL_NGeo(i,l)*xCLOld(:,l,j,k,iElemOld)
      dXCL_NGeo(:,2,i,j,k)=dXCL_NGeo(:,2,i,j,k) + DCL_NGeo(j,l)*xCLOld(:,i,l,k,iElemOld)
      dXCL_NGeo(:,3,i,j,k)=dXCL_NGeo(:,3,i,j,k) + DCL_NGeo(k,l)*xCLOld(:,i,j,l,iElemOld)
    END DO
  END DO; END DO; END DO

  DO iElemNew=1,nElemsNew
    IF(ElemDone(iElemNew)) CYCLE ! all IP already found
    ! Check if the two elements could be overlapping by comparing the distance between the centroids with the sum of the radii
    maxDist=radNew(iElemNew)+radOld(iElemOld)
    IF(SUM((xCOld(:,iElemOld)-xCNew(:,iElemNew))**2).GT.maxDist**2) CYCLE

    ! Check which IP of new elem overlaps with old elem
    IPOverlaps=.FALSE.
    DO kk=0,NInter; DO jj=0,NInter; DO ii=0,NInter
      IF(IPDone(ii,jj,kk,iElemNew)) CYCLE

      dist=SUM((xCOld(:,iElemOld)-xCLInter(:,ii,jj,kk,iElemNew))**2)
      IPOverlaps(ii,jj,kk)=(dist.LE.radSqOld(iElemOld))
    END DO; END DO; END DO

    ! TODO: Move to own routine
    ! Get smallest distance to supersampled points for starting Newton
    DO kk=0,NInter; DO jj=0,NInter; DO ii=0,NInter
      IF(.NOT.IPOverlaps(ii,jj,kk)) CYCLE

      xInter = xCLInter(:,ii,jj,kk,iElemNew)
      best=HUGE(1.)
      DO i=0,NSuper; DO j=0,NSuper; DO k=0,NSuper
        dist=SUM((xInter-X_NSuper(:,i,j,k))**2)
        IF (dist.LT.best) THEN
          best=dist
          xi=(/Xi_NSuper(i),Xi_NSuper(j),Xi_NSuper(k)/)
        END IF
      END DO; END DO; END DO

      CALL LagrangeInterpolationPolys(Xi(1),NGeoOld,Xi_CLNGeo,wBary_CLNGeo,LagXi)
      CALL LagrangeInterpolationPolys(Xi(2),NGeoOld,Xi_CLNGeo,wBary_CLNGeo,LagEta)
      CALL LagrangeInterpolationPolys(Xi(3),NGeoOld,Xi_CLNGeo,wBary_CLNGeo,LagZeta)

      ! F(xi) = x(xi) - xInter
      F=-xInter
      DO k=0,NGeoOld; DO j=0,NGeoOld; DO i=0,NGeoOld
        LagVol(i,j,k)=LagXi(i)*LagEta(j)*LagZeta(k)
        F=F+xCLOld(:,i,j,k,iElemOld)*LagVol(i,j,k)
      END DO; END DO; END DO

      !eps_F=1.E-16
      eps_F=1.E-8*SUM(F*F) ! relative error to initial guess
      iter=0
      DO WHILE ((SUM(F*F).GT.eps_F).AND.(iter.LT.100))
        iter=iter+1
        ! Compute F Jacobian dx/dXi
        Jac=0.
        DO k=0,NGeoOld; DO j=0,NGeoOld; DO i=0,NGeoOld
          Jac=Jac+dXCL_NGeo(:,:,i,j,k)*LagVol(i,j,k)
        END DO; END DO; END DO

        ! Compute inverse of Jacobian
        sdetJac=DET3(Jac)
        IF(sdetJac.NE.0.) THEN
         sdetJac=1./sdetJac
        ELSE
          ! Newton has not converged !?!?
          ! allow Newton to fail without aborting, may happen when far outside of reference space [-1,1]
!          WRITE(UNIT_stdOut,'(A)')' Newton has not converged! skipping...'        
!          WRITE(UNIT_stdOut,*)' Xi,Eta,Zeta = ', Xi
!         CALL abort(__STAMP__, &
!              'Newton method is singular!')
          EXIT
        ENDIF
        sJac=INVERSE3(Jac,sdetJac)

        ! Iterate Xi using Newton step
        Xi = Xi - MATMUL(sJac,F)
        ! if Newton gets outside reference space range [-1,1], exit. 
        ! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
        IF((iter.GT.3).AND.(ANY(ABS(Xi).GT.1.2))) EXIT

        ! Compute function value
        CALL LagrangeInterpolationPolys(Xi(1),NGeoOld,Xi_CLNGeo,wBary_CLNGeo,LagXi)
        CALL LagrangeInterpolationPolys(Xi(2),NGeoOld,Xi_CLNGeo,wBary_CLNGeo,LagEta)
        CALL LagrangeInterpolationPolys(Xi(3),NGeoOld,Xi_CLNGeo,wBary_CLNGeo,LagZeta)
        ! F(xi) = x(xi) - xInter
        F=-xInter
        DO k=0,NGeoOld; DO j=0,NGeoOld; DO i=0,NGeoOld
          LagVol(i,j,k)=LagXi(i)*LagEta(j)*LagZeta(k)
          F=F+xCLOld(:,i,j,k,iElemOld)*LagVol(i,j,k)
        END DO; END DO; END DO
      END DO !newton

      ! check if result is better than previous result
!$OMP CRITICAL 
      IF(MAXVAL(ABS(Xi)).LT.MAXVAL(ABS(xiInter(:,ii,jj,kk,iElemNew)))) THEN
        IF(MAXVAL(ABS(Xi)).LE.1.) IPDone(ii,jj,kk,iElemNew) = .TRUE. ! if point is inside element, stop searching
        xiInter(:,ii,jj,kk,iElemNew)=Xi
        InterToElem(ii,jj,kk,iElemNew) = iElemOld
      END IF
!$OMP END CRITICAL 
    END DO; END DO; END DO ! ii,jj,kk (IP loop)
    IF(ALL(IPDone(:,:,:,iElemNew))) ElemDone(iElemNew)=.TRUE.
  END DO ! iElem
!$OMP CRITICAL 
  ElemCounter=ElemCounter+1
  IF(MOD(ElemCounter,1000).EQ.0)THEN
    WRITE(*,*) ElemCounter,'elements of',nElemsOld-nEqualElems,' processed. Time= ', FLEXITIME()-StartTime
    StartTime=FLEXITIME()
  END IF
!$OMP END CRITICAL 
END DO ! iElemOld
!$OMP END DO
!$OMP END PARALLEL

! check if all interpolation points have been found
nNotFound=0
DO iElemNew=1,nElemsNew
  IF(equalElem(iElemNew).GT.0) CYCLE
  DO ii=0,NInter; DO jj=0,NInter; DO kk=0,NInter
    IF(.NOT.IPDone(ii,jj,kk,iElemNew))THEN
      ! Only mark as invalid if greater then max tolerance
      IF(MAXVAL(ABS(xiInter(:,ii,jj,kk,iElemNew))).GT.maxTol)THEN
        ! RP has not been found
        IF(printTroublemakers)THEN
          WRITE(*,*) 'IP with ID:',ii,jj,kk,iElemNew,'is a troublemaker!'
          WRITE(*,*) 'Coords:',xCLInter(:,ii,jj,kk,iElemNew)
          WRITE(*,*) 'Xi:',xiInter(:,ii,jj,kk,iElemNew)
        END IF
        !abort if severly broken and no refstate present (with refstate present abortTol=HUGE)
        IF(MAXVAL(ABS(xiInter(:,ii,jj,kk,iElemNew))).GT.abortTol)&
          CALL abort(__STAMP__, 'IP not found.')
        nNotFound=nNotFound+1
      ELSE
        IPDone(ii,jj,kk,iElemNew)=.TRUE.
      END IF
    END IF
  END DO; END DO; END DO ! ii,jj,kk (IP loop)
END DO
WRITE(*,*) nNotFound,' nodes not found.'
Time=FLEXITIME() -Time
WRITE(UNIT_stdOut,'(A,F0.3,A)')' DONE  [',Time,'s]'
END SUBROUTINE GetParametricCoordinates

!=================================================================================================================================
!> Computes the centroid and the radius of an element
!=================================================================================================================================
SUBROUTINE getCentroidAndRadius(elem,NGeo,xC,radius)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: NGeo                         !> Polynomial degree of mesh representation
REAL,INTENT(IN)    :: elem(3,0:NGeo,0:NGeo,0:NGeo) !> Coordinates of single element
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: xC(3)                        !> Coordinates of centroid
REAL,INTENT(OUT)   :: radius                       !> Radius of element
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,nNodes
LOGICAL            :: onSide(3)
!=================================================================================================================================
! Compute centroid of element
nNodes= (NGeo+1)**3
xC(1) = SUM(elem(1,:,:,:))/nNodes
xC(2) = SUM(elem(2,:,:,:))/nNodes
xC(3) = SUM(elem(3,:,:,:))/nNodes

! Compute max distance from bary to surface nodes and return as radius
radius=0.
DO k=0,NGeo
  onSide(3)=((k.EQ.0).OR.(k.EQ.NGeo))
  DO j=0,NGeo
    onSide(2)=((j.EQ.0).OR.(j.EQ.NGeo))
    DO i=0,NGeo
      onSide(1)=((i.EQ.0).OR.(i.EQ.NGeo))
      IF(.NOT.ANY(onSide)) CYCLE
      radius=MAX(radius,NORM2(elem(:,i,j,k)-xC))
    END DO
  END DO
END DO

END SUBROUTINE getCentroidAndRadius

END MODULE MOD_SMParametricCoordinates
