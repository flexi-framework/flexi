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
!> Module to search for the reference coordinates of the RPs by inverting the mapping using Newton's method and to sort them
!> according to the element number for later parallel read-in.
!===================================================================================================================================
MODULE MOD_RPParametricCoords
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE GetRecordPoints
  MODULE PROCEDURE GetRecordPoints
END INTERFACE

PUBLIC :: GetRecordPoints
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Wrapper routine to call the search and sorting algorithms
!===================================================================================================================================
SUBROUTINE GetRecordPoints()
! MODULES
USE MOD_Globals
USE MOD_RPSet_Vars        ,ONLY: RPSetInitIsDone
USE MOD_Mesh_Vars,         ONLY: MeshInitIsDone
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys,ChebyGaussLobNodesAndWeights,BarycentricWeights,InitializeVandermonde
USE MOD_ReadIntools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.RPSetInitIsDone.OR.(.NOT.MeshInitIsDone))THEN
  CALL Abort(__STAMP__, &
       'GetRecordPoints not ready to be called or already called.')
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' FIND RECORDPOINTS IN MESH...'
CALL GetParametricCoordinates()
CALL SortRP()

SWRITE(UNIT_stdOut,'(A)')' FINDING RECORDPOINTS DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE GetRecordPoints



!===================================================================================================================================
!> Computes parametric coordinates and element of recordpoints given the physical coordinates iteratively
!===================================================================================================================================
SUBROUTINE GetParametricCoordinates()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys,ChebyGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis,             ONLY: PolynomialDerivativeMatrix
USE MOD_ChangeBasis,       ONLY: ChangeBasis3D,ChangeBasis2D
USE MOD_Interpolation,     ONLY: GetVandermonde,GetNodesAndWeights,GetDerivativeMatrix
USE MOD_Interpolation_Vars
USE MOD_Mesh_Vars,         ONLY: NGeo,Elem_xGP,SideToElem,nBCSides,Face_xGP,NormVec,nElems
USE MOD_Newton,            ONLY: Newton
USE MOD_Parameters,        ONLY: NSuper,maxTol
USE MOD_RPSet_Vars,        ONLY: RPlist,nRP_global,tRP
USE MOD_RPSet_Vars,        ONLY: nRP_global
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: DCL_NGeo(0:Ngeo,0:Ngeo)
REAL                  :: XCL_Ngeo(3,0:Ngeo,0:Ngeo,0:Ngeo)          !< mapping X(xi) P\in Ngeo
REAL                  :: hmax2
REAL                  :: X_NSuper(1:3,0:NSuper,0:NSuper,0:NSuper)
REAL                  :: Winner_Dist2,Dist2
REAL                  :: dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo)
REAL                  :: Xi_CLNGeo(0:NGeo),wBary_CLNGeo(0:NGeo)
REAL                  :: Xi_NSuper(0:NSuper)
REAL                  :: xRP(3),xC(3)
REAL                  :: Vdm_N_CLNGeo(0:NGeo,0:PP_N)
REAL                  :: Vdm_CLNGeo_EquiNSuper(0:NSuper,0:NGeo)
REAL                  :: Lag(1:3,0:NGeo)
REAL                  :: F(1:3),Xi(1:3)
INTEGER               :: i,j,k,l,iElem,nNodes
INTEGER               :: iRP
CHARACTER(LEN=255)    :: NodeType_Super
TYPE(tRP),POINTER     :: aRP
LOGICAL               :: changeBasisDone,calcJacobianDone
LOGICAL               :: RPFound(1:nRP_Global),onSide(3)
! boundary projection
INTEGER               :: SideID,locSideID,iWinner,jWinner,nRP_notFound,iRP2
INTEGER,ALLOCATABLE   :: mapRP(:)
REAL,ALLOCATABLE      :: dist2RP(:)
REAL                  :: Vdm_N_EquiNSuper(0:NSuper,0:PP_N),xBC_NSuper(3,0:NSuper,0:NSuper)
REAL                  :: NormVec_NSuper(3,0:NSuper,0:NSuper)
!newton method
REAL                  :: wBary_NSuper(0:NSuper), D_NSuper(0:NSuper,0:NSuper), Lag_NSuper(1:2,0:NSuper)
REAL                  :: dxBC_NSuper(3,2,0:NSuper,0:NSuper)
REAL                  :: Gmat(2,0:NSuper,0:NSuper),dGmat(2,2,0:NSuper,0:NSuper)
REAL                  :: Xi2(2),xWinner(3),NormVecWinner(3)
!===================================================================================================================================
! Prepare CL basis evaluation
CALL GetNodesAndWeights( NGeo, NodeTypeCL, Xi_CLNGeo,wIPBary=wBary_CLNGeo)
CALL GetDerivativeMatrix(Ngeo, NodeTypeCL, DCL_Ngeo)

DO i=0,NSuper
  Xi_NSuper(i) = 2./REAL(NSuper) * REAL(i) - 1.
END DO

NodeType_Super='VISU'
CALL GetVanderMonde(PP_N,NodeType  ,NGeo  ,NodeTypeCL    ,Vdm_N_CLNGeo)
CALL GetVanderMonde(NGeo,NodeTypeCL,NSuper,NodeType_Super,Vdm_CLNGeo_EquiNSuper)
! Define maximum mesh size (square) hmax2
nNodes=(NGeo+1)**3
RPFound=.FALSE.

DO iElem=1,nElems
  ! We have the Elem_xGP on N, compute XCL_NGeo for minimal changes
  CALL ChangeBasis3D(3,PP_N,NGeo,Vdm_N_CLNGeo,Elem_xGP(:,:,:,:,iElem),XCL_NGeo)

  ! Compute centroid of element
  xC(1) = SUM(XCL_NGeo(1,:,:,:))/nNodes
  xC(2) = SUM(XCL_NGeo(2,:,:,:))/nNodes
  xC(3) = SUM(XCL_NGeo(3,:,:,:))/nNodes

  ! Compute max distance from bary to surface nodes
  hmax2=0.
  DO k=0,NGeo
    onSide(3)=((k.EQ.0).OR.(k.EQ.NGeo))
    DO j=0,NGeo
      onSide(2)=((j.EQ.0).OR.(j.EQ.NGeo))
      DO i=0,NGeo
        onSide(1)=((i.EQ.0).OR.(i.EQ.NGeo))
        IF(.NOT.ANY(onSide)) CYCLE
        hmax2=MAX(hmax2,SUM((XCL_NGeo(:,i,j,k)-xC)*(XCL_NGeo(:,i,j,k)-xC)))
      END DO
    END DO
  END DO
  hmax2=hmax2*1.10 ! 10% tolerance

  changeBasisDone=.FALSE.
  calcJacobianDone=.FALSE.
  DO iRP=1,nRP_global
    IF(RPFound(iRP)) CYCLE

    aRP=>RPlist(iRP)%RP
    xRP=aRP%x
    ! Coarse search of possible elems
    IF(SUM((xC-xRP)*((xC-xRP))).GT.hmax2) CYCLE

    ! Find initial guess for Newton, get supersampled element for first recordpoint
    IF(.NOT.changeBasisDone) &
      CALL ChangeBasis3D(3,NGeo,NSuper,Vdm_CLNGeo_EquiNSuper,XCL_NGeo,X_NSuper)
    changeBasisDone=.TRUE.

    Winner_Dist2=HUGE(1.)
    DO i=0,NSuper; DO j=0,NSuper; DO k=0,NSuper
      Dist2=SUM((xRP-X_NSuper(:,i,j,k))*(xRP-X_NSuper(:,i,j,k)))
      IF (Dist2.LT.Winner_Dist2) THEN
        Winner_Dist2=Dist2
        Xi=(/Xi_NSuper(i),Xi_NSuper(j),Xi_NSuper(k)/)
      END IF
    END DO; END DO; END DO

    ! Compute Jacobian of Winner element Mapping for each CL point
    IF(.NOT.calcJacobianDone)THEN
      dXCL_NGeo=0.
      DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
        ! Matrix-vector multiplication
        DO l=0,NGeo
          dXCL_NGeo(:,1,i,j,k)=dXCL_NGeo(:,1,i,j,k) + DCL_NGeo(i,l)*XCL_NGeo(:,l,j,k)
          dXCL_NGeo(:,2,i,j,k)=dXCL_NGeo(:,2,i,j,k) + DCL_NGeo(j,l)*XCL_NGeo(:,i,l,k)
          dXCL_NGeo(:,3,i,j,k)=dXCL_NGeo(:,3,i,j,k) + DCL_NGeo(k,l)*XCL_NGeo(:,i,j,l)
        END DO !l=0,NGeo
      END DO; END DO; END DO
      calcJacobianDone=.TRUE.
    END IF

    CALL Newton(NGeo,XRP,dXCL_NGeo,Xi_CLNGeo,wBary_CLNGeo,XCL_NGeo(:,:,:,:),Xi,LagOut=Lag,FOut=F)

    ! check if result is better then previous result
    IF(MAXVAL(ABS(Xi)).LT.MAXVAL(ABS(aRP%xi))) THEN
      IF(MAXVAL(ABS(Xi)).LE.1.) RPFound(iRP) = .TRUE. ! if point is inside element, stop searching
      aRP%xi=Xi
      aRP%ElemID = iElem
      F=0.
      DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
        F=F+XCL_NGeo(:,i,j,k)*Lag(1,i)*Lag(2,j)*Lag(3,k)
      END DO; END DO; END DO
      aRP%xF=F
    END IF

  END DO
END DO

DO iRP=1,nRP_Global
  IF(.NOT.RPFound(iRP))THEN
    aRP=>RPlist(iRP)%RP
    ! mark as valid if greater than max tolerance
    IF(MAXVAL(ABS(aRP%Xi)).LT.maxTol) RPFound(iRP)=.TRUE.
  END IF
END DO

! Remaining points: If the point is close to a boundary, project the point on the boundary
IF(ANY(.NOT.RPFound)) THEN
  nRP_notfound=nRP_global-COUNT(RPFound)
  SWRITE(UNIT_stdOut,'(A,I4,A,I4,A)')' ',nRP_notfound,' of ',nRP_global,' RPs have not been found inside the mesh.'
  SWRITE(UNIT_stdOut,'(A)')' Attempting to project them on the closest boundary...'
  ALLOCATE(dist2RP(nRP_notfound))
  ALLOCATE(mapRP(nRP_global))
  iRP2=0
  mapRP=-999
  DO iRP=1,nRP_global
    IF(RPFound(iRP)) CYCLE
    iRP2=iRP2+1
    mapRP(iRP)=iRP2
  END DO! iRP=1,nRP_global
  dist2RP=HUGE(1.)
  IF(NSuper.LT.Ngeo*2) &
    WRITE(*,*)'Warning: NSuper<2Ngeo, derivative may be wrong.'
  ! Prepare equidistant basis
  CALL BarycentricWeights(NSuper,Xi_NSuper,wBary_NSuper)
  CALL GetVandermonde(PP_N,NodeType,NSuper,NodeType_Super,Vdm_N_EquiNSuper)
  CALL PolynomialDerivativeMatrix(NSuper,Xi_NSuper,D_NSuper)
  nNodes=(NSuper+1)**2
  DO SideID=1,nBCSides
    calcJacobianDone=.FALSE.
    ! Supersampling of the side to equidistant grid for search
    CALL ChangeBasis2D(3,PP_N,NSuper,Vdm_N_EquiNSuper,Face_xGP(:,:,:,0,SideID),xBC_NSuper)
    CALL ChangeBasis2D(3,PP_N,NSuper,Vdm_N_EquiNSuper,NormVec( :,:,:,0,SideID),NormVec_NSuper)
    ! get the BCSide centroid
    xC(1) = SUM(xBC_NSuper(1,:,:))/nNodes
    xC(2) = SUM(xBC_NSuper(2,:,:))/nNodes
    xC(3) = SUM(xBC_NSuper(3,:,:))/nNodes
    ! calculate the max. distance within the side to the centroid
    hmax2=0.
    DO j=0,NSuper
      DO i=0,NSuper
        hmax2=MAX(hmax2,SUM((xBC_NSuper(:,i,j)-xC)*(xBC_NSuper(:,i,j)-xC)))
      END DO
    END DO
    hmax2=hmax2*1.20 ! 20% tolerance
    DO iRP=1,nRP_Global
      IF(RPFound(iRP)) CYCLE

      aRP=>RPlist(iRP)%RP
      xRP=aRP%x
      ! Check if the point is close to the boundary side
      ! Coarse search of possible elems
      IF(SUM((xC-xRP)*((xC-xRP))).GT.hmax2) CYCLE
      ! Get closest point on the supersampled side as starting value
      locSideID= SideToElem(S2E_LOC_SIDE_ID,SideID)
      iElem    = SideToElem(S2E_ELEM_ID,SideID)
      Winner_Dist2=HUGE(1.)
      DO i=0,NSuper; DO j=0,NSuper
        Dist2=SUM((xRP-xBC_NSuper(:,i,j))*(xRP-xBC_NSuper(:,i,j)))
        IF (Dist2.LT.Winner_Dist2) THEN
          Winner_Dist2=Dist2
          iWinner=i
          jWinner=j
        END IF
      END DO; END DO
      Xi2=(/Xi_NSuper(iWinner),Xi_NSuper(jWinner)/)
      xWinner=xBC_NSuper(:,iWinner,jWinner)
      NormVecWinner=NormVec_NSuper(:,iWinner,jWinner)
      F=xRP-xWinner

      ! Newton to find the minimum distance
      ! Calculate the surface jacobian
      IF(.NOT.calcJacobianDone)THEN
        dxBC_NSuper=0.
        DO j=0,NSuper
          DO i=0,NSuper
          ! Matrix-vector multiplication
            DO l=0,NSuper
              dxBC_NSuper(:,1,i,j)=dxBC_NSuper(:,1,i,j) + D_NSuper(i,l)*xBC_NSuper(:,l,j)
              dxBC_NSuper(:,2,i,j)=dxBC_NSuper(:,2,i,j) + D_NSuper(j,l)*xBC_NSuper(:,i,l)
            END DO !l=0,NSuper
          END DO !i=0,NSuper
        END DO !j=0,NSuper
        calcJacobianDone=.TRUE.
      END IF

      ! for Newton we first need the function Gmat(:,i,j) and its gradient in parameter space
      ! G= d/dXi((xBC-xRP)²)=0, degree of G 2NGeo
      Gmat=0.
      DO j=0,NSuper
        DO i=0,NSuper
          Gmat(:,i,j)=Gmat(:,i,j)+dxBC_nSuper(1,:,i,j)*2*(xBC_NSuper(1,i,j)-xRP(1)) &
                                 +dxBC_nSuper(2,:,i,j)*2*(xBC_NSuper(2,i,j)-xRP(2)) &
                                 +dxBC_nSuper(3,:,i,j)*2*(xBC_NSuper(3,i,j)-xRP(3))
        END DO! i=0,NSuper
      END DO! j=0,NSuper

      dGmat=0.
      DO j=0,NSuper
        DO i=0,NSuper
          ! Matrix-vector multiplication
          DO l=0,NSuper
            dGmat(:,1,i,j)=dGmat(:,1,i,j) + D_NSuper(i,l)*Gmat(:,l,j)
            dGmat(:,2,i,j)=dGmat(:,2,i,j) + D_NSuper(j,l)*Gmat(:,i,l)
          END DO !l=0,NSuper
        END DO! i=0,NSuper
      END DO! j=0,NSuper
      ! get initial value of the functional G

      CALL Newton(NSuper,(/0.,0./),dGmat,Xi_NSuper,wBary_NSuper,Gmat,Xi2,LagOut=Lag_NSuper)

      ! use Newton result if minimum is within parameter range, else see if supersampled
      ! initial guess is better than previous result
      IF(MAXVAL(ABS(Xi2)).LE.1.) THEN ! use newton result
        ! calculate new distance
        xWinner=0.
        NormVecWinner=0.
        DO j=0,NSuper
          DO i=0,NSuper
            xWinner=xWinner+xBC_NSuper(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
            NormVecWinner=NormVecWinner+NormVec_NSuper(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
          END DO! i=0,NSuper
        END DO! j=0,NSuper
        Winner_Dist2=SUM((xWinner-xRP)*(xWinner-xRP))
      END IF

      NormVecWinner=NormVecWinner/NORM2(NormVecWinner)
      F=(xRP-xWinner)/NORM2(xRP-xWinner)
      IF((Winner_Dist2.LE.dist2RP(mapRP(iRP))))THEN
        SELECT CASE(locSideID)
        CASE(XI_MINUS)
          ! switch to right hand system
          aRP%xi=(/-1.,Xi2(2),Xi2(1)/)
        CASE(ETA_MINUS)
          aRP%xi=(/Xi2(1),-1.,Xi2(2)/)
        CASE(ZETA_MINUS)
          ! switch to right hand system
          aRP%xi=(/Xi2(2),Xi2(1),-1./)
        CASE(XI_PLUS)
          aRP%xi=(/1.,Xi2(1),Xi2(2)/)
        CASE(ETA_PLUS)
          ! switch to right hand system
          aRP%xi=(/-Xi2(1),1.,Xi2(2)/)
        CASE(ZETA_PLUS)
          aRP%xi=(/Xi2(1),Xi2(2),1./)
        END SELECT
        aRP%ElemID = iElem
        aRP%xF=xWinner
        ! keep the RP in the not found list, find the side with the best angle and minimum distance
        dist2RP(mapRP(iRP))=Winner_Dist2
      END IF
    END DO! iRP=1,nRP_Global
  END DO! SideID=1,nBCSides
  SWRITE(UNIT_stdOut,'(A)')' done.'
  SWRITE(UNIT_stdOut,'(A,F15.8)')'  Max. distance: ',SQRT(MAXVAL(dist2RP))

  DEALLOCATE(mapRP,dist2RP)
END IF!(.NOT.ANY(RPFound)

DO iRP=1,nRP_Global
  IF(.NOT.RPFound(iRP))THEN
    aRP=>RPlist(iRP)%RP
    ! Only mark as invalid if greater then max tolerance
    IF(MAXVAL(ABS(aRP%Xi)).GT.maxTol)THEN
      ! RP has not been found
      WRITE(*,*) 'Record Point with ID :',iRP,' and Coordinates ',aRP%x, ' is a troublemaker!'
      CALL Abort(__STAMP__, &
           'Newton has reached 50 Iter, Point not found')
    END IF
  END IF
END DO

END SUBROUTINE GetParametricCoordinates


!===================================================================================================================================
!> sort points according to element numbering and prepare Offset Array - needed for parallel processing
!===================================================================================================================================
SUBROUTINE SortRP()
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars         ,ONLY: nElems
USE MOD_RPSet_Vars        ,ONLY: RPlist,nRP_global,tRP,tRPlist
USE MOD_RPSet_Vars        ,ONLY: OffsetRP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem
INTEGER               :: iRP,iRP_new
TYPE(tRP),POINTER     :: aRP
TYPE(tRPlist),POINTER :: RPlist_tmp(:)
!===================================================================================================================================

! sort points according to element numbering
ALLOCATE(OffsetRP(2,nElems))
OffsetRP=0
ALLOCATE(RPlist_tmp(nRP_global))
DO iRP=1,nRP_global
  RPlist_tmp(iRP)%RP=>RPlist(iRP)%RP
END DO !iRP

iRP_new=0
DO iElem=1,nElems
  IF(iElem .GT. 1) THEN
    OffsetRP(1,iElem)=OffsetRP(2,iElem-1)
    OffsetRP(2,iElem)=OffsetRP(2,iElem-1)
  END IF
  DO iRP=1,nRP_global
    aRP=>RPlist_tmp(iRP)%RP
    IF(aRP%ElemID.EQ.iElem) THEN
      iRP_new=iRP_new+1
      RPlist(iRP_new)%RP=>aRP
      RPlist(iRP_new)%RP%ID=iRP_new
      OffsetRP(2,iElem)=iRP_new
    END IF
  END DO !iRP
END DO !iElem

END SUBROUTINE SortRP

END MODULE MOD_RPParametricCoords

