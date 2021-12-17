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
!> Module that includes the routines necessary to generate boundary layer planes, which are created by the projection of a spline
!> onto the nearest boundary and extrusion of the resulting line along the wall normal direction
!===================================================================================================================================
MODULE MOD_BLProjection
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE GetBLPlane
  MODULE PROCEDURE GetBLPlane
END INTERFACE

INTERFACE GetBLBox
  MODULE PROCEDURE GetBLBox
END INTERFACE

PUBLIC :: GetBLPlane,GetBLBox
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Define a spline from control points, create equidistant parametrization along that spline, project it on
!> the closest boundary and create a plane along the boundary normals
!===================================================================================================================================
SUBROUTINE GetBLPlane(Plane,nCP,height,fac,xCP)
! MODULES
USE MOD_Globals
USE MOD_RPSet_Vars,ONLY:tPlane,tRPlist
USE MOD_RPSet_Vars,ONLY:GetNewRP
USE MOD_Spline    ,ONLY:GetSpline,GetEquiPoints,EvalSpline,EvalSplineDeriv,EvalEquiError
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tPlane),POINTER            :: Plane
INTEGER,INTENT(IN)              :: nCP
REAL,INTENT(IN)                 :: height(nCP),fac
REAL,INTENT(IN)                 :: xCP(3,nCP)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nRP(2),i,j,iCP,iter
REAL                            :: x_loc(3),s_loc,dx_loc(3)
TYPE(tRPlist),POINTER           :: RPlist_tmp(:)
REAL,ALLOCATABLE                :: xRP(:,:),xRP_tmp(:,:),NormVecRP(:,:),TangVecRP(:,:),dh(:)
REAL,ALLOCATABLE                :: s(:),s_mod(:),s_equi(:),coeff(:,:,:)
REAL                            :: h,t,height_loc,EquiErr
!===================================================================================================================================
nRP=Plane%nRP(1:2)
ALLOCATE(RPlist_tmp(nRP(1)))
ALLOCATE(NormVecRP(3,nRP(1)))
ALLOCATE(TangVecRP(3,nRP(1)))
ALLOCATE(xRP(3,nRP(1)))
ALLOCATE(xRP_tmp(3,nRP(1)))
ALLOCATE(s_equi(nRP(1)))
ALLOCATE(coeff(3,4,nCP-1),s(nCP))

CALL GetSpline(3,nCP,xCP,coeff,s) ! get spline coefficients and arclength variable s

! Sampling of the spline to nRP(1) and allocate first row of RPs
DO i=1,nRP(1)
  s_loc=s(nCP)*(i-1)/(nRP(1)-1)
  CALL EvalSpline(3,nCP,s_loc,s,coeff,x_loc)
  CALL GetNewRP(Plane%RP_ptr(i,1)%RP,Plane%GroupID,x_loc)
  RPlist_tmp(i)%RP=>Plane%RP_ptr(i,1)%RP
END DO

! projection of the supersampled points on the closest boundary
CALL ProjectRPtoBC(nRP(1),RPlist_tmp,NormVecRP)
DO i=1,nRP(1)
  RPlist_tmp(i)%RP%x=RPlist_tmp(i)%RP%xF
  xRP(:,i)=RPlist_tmp(i)%RP%xF
END DO

!iterate between projection and equidistant partitioning
DO iter=1,10
  ! get equidistant distribution along projected spline
  xRP_tmp = xRP
  CALL GetEquiPoints(3,nRP(1),nRP(1),xRP_tmp,xRP,s_equi)
  DO i=1,nRP(1)
    RPlist_tmp(i)%RP%x=xRP(:,i)
  END DO
  ! evaluate deviation from equidistancy
  CALL EvalEquiError(3,nRP(1),xRP,EquiErr)
  WRITE(*,*) 'After partitioning',EquiErr
  ! project the equidistant points again to get their normal vector
  CALL ProjectRPtoBC(nRP(1),RPlist_tmp,NormVecRP)
  DO i=1,nRP(1)
    RPlist_tmp(i)%RP%x=RPlist_tmp(i)%RP%xF
    xRP(:,i)=RPlist_tmp(i)%RP%xF
  END DO
  ! evaluate deviation from equidistancy
  CALL EvalEquiError(3,nRP(1),xRP,EquiErr)
  WRITE(*,*) 'After projection',iter, EquiErr
END DO

! calculate the tangent vector
DEALLOCATE(coeff)
ALLOCATE(s_mod(nRP(1)),coeff(3,4,nRP(1)-1))
CALL GetSpline(3,nRP(1),xRP,coeff,s_mod) ! get the spline through the projected points
! get the tangent vector in each point
DO i=1,nRP(1)
  !derivative of the spline
  CALL EvalSplineDeriv(3,nRP(1),s_mod(i),s_mod,coeff,dx_loc)
  !project it on the local surface
  tangVecRP(:,i)=dx_loc - SUM(dx_loc(1:3)*NormVecRP(1:3,i))*NormVecRP(:,i)
  tangVecRP(:,i)=tangVecRP(:,i)/SQRT(DOT_PRODUCT(tangVecRP(:,i),tangVecRP(:,i)))
END DO

! extrapolation of the BL mesh along the boundary normals to height
ALLOCATE(dh(2:nRP(2)))
dh(2)=1.
DO i=3,nRP(2)
  dh(i)=dh(i-1)*fac
END DO
dh(:)=dh(:)/SUM(dh(:))
iCP=2
DO i=1,nRP(1)
  s_loc=s(nCP)*(i-1)/(nRP(1)-1)
  ! linear interpolation of the height between control points
  DO WHILE(s(iCP).LT.s_loc)
    iCP=iCP+1
  END DO
  iCP=MIN(iCP,nCP)
!  IF(s_loc.GE.s(iCP)) iCP=MAX(iCP+1,nCP)
  t=(s_loc-s(iCP-1))/(s(iCP)-s(iCP-1))
  height_loc=height(iCP-1)*(1-t)+height(iCP)*t
  h=0.
  DO j=2,nRP(2)
    h=h+dh(j)
    x_loc=Plane%RP_ptr(i,1)%RP%x+NormVecRP(:,i)*h*height_loc
    CALL GetNewRP(Plane%RP_ptr(i,j)%RP,Plane%GroupID,x_loc)
  END DO
END DO

ALLOCATE(Plane%NormVec(3,nRP(1)))
ALLOCATE(Plane%TangVec(3,nRP(1)))
Plane%NormVec=NormVecRP
Plane%TangVec=TangVecRP


DEALLOCATE(RPlist_tmp,NormVecRP,TangVecRP,coeff,s)
END SUBROUTINE GetBLPlane

!===================================================================================================================================
!> Define a spline from control points, create equidistant parametrization along that spline, project it on
!> the closest boundary and create a box along the boundary normals
!===================================================================================================================================
SUBROUTINE GetBLBox(Box,nCP,nSP,height,fac,xCP)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_RPSet_Vars,ONLY:tBox,tRPlist
USE MOD_RPSet_Vars,ONLY:GetNewRP
USE MOD_Spline    ,ONLY:GetSpline,GetEquiPoints,EvalSpline,EvalSplineDeriv,EvalEquiError
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tBox),POINTER              :: Box
INTEGER,INTENT(IN)              :: nCP
INTEGER,INTENT(IN)              :: nSP
REAL,INTENT(IN)                 :: height(nCP,nSP),fac
REAL,INTENT(IN)                 :: xCP(3,nCP,nSP)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nRP(3),i,j,k,iCP,iSP,iter
REAL                            :: x_loc(3),s_loc,dx_loc(3)
TYPE(tRPlist),POINTER           :: RPlist_tmp(:,:)
REAL,ALLOCATABLE                :: xRP(:,:,:),xRP_tmp(:,:,:),NormVecRP(:,:,:),TangVecRP(:,:,:),dh(:)
REAL,ALLOCATABLE                :: s(:,:),s_mod(:),s_equi(:,:),coeff(:,:,:,:)
REAL                            :: h,t,height_loc,EquiErr,EquiErrSum
INTEGER,ALLOCATABLE             :: Points(:)
INTEGER                         :: iBase
REAL,ALLOCATABLE                :: xBase(:,:,:),heightBase(:,:),depth(:)
!===================================================================================================================================

! get resolution
nRP=Box%nRP(1:3)

ALLOCATE(RPlist_tmp(nRP(1),nRP(3)))
ALLOCATE(NormVecRP(3,nRP(1),nRP(3)))
ALLOCATE(TangVecRP(3,nRP(1),nRP(3)))
ALLOCATE(xRP(3,nRP(1),nRP(3)))
ALLOCATE(xRP_tmp(3,nRP(1),nRP(3)))
ALLOCATE(s_equi(nRP(1),nRP(3)))
ALLOCATE(coeff(3,4,nCP-1,nRP(3)),s(nCP,nRP(3)))

! linear interpolation of data
ALLOCATE(depth(nSP-1))
ALLOCATE(Points(nSP-1))
DO iSP=2,nSP
  depth(iSP-1) = (xCP( 3,1,iSP)-xCP( 3,1,iSP-1))/(xCP( 3,1,nSP)-xCP( 3,1,1))
END DO
Points=FLOOR((REAL(nRP(3))-REAL(nSP))*depth)
iSP = 1
DO WHILE(SUM(Points).LT.(nRP(3)-nSP))
  Points(iSP) = Points(iSP)+1
  iSP = iSP+1
  IF (iSP.EQ.nSP) iSP = 1
END DO
ALLOCATE(xBase(3,nCP,nRP(3)))
ALLOCATE(heightBase(nCP,nRP(3)))
iBase = 0
iSP = 1
DO k=1,nRP(3)
  IF (iBase.EQ.(Points(iSP)+1)) THEN
    iBase = 0
    iSP = iSP+1
  END IF
  IF (iBase.EQ.0) THEN ! reuse provided data
    xBase(:,:,k)    = xCP( :,:,iSP)
    heightBase(:,k) = height(:,iSP)
  ELSE ! linear interpolation between lines
      ! in x
      xBase(1,:,k)    = xCP(1,:,iSP)+(xCP(1,:,iSP+1)-xCP(1,:,iSP))/(Points(iSP)+1)*iBase
      ! in y
      xBase(2,:,k)    = xCP(2,:,iSP)+(xCP(2,:,iSP+1)-xCP(2,:,iSP))/(Points(iSP)+1)*iBase
      ! in z
      xBase(3,:,k)    = xCP(3,:,iSP)+(xCP(3,:,iSP+1)-xCP(3,:,iSP))/(Points(iSP)+1)*iBase
      ! height
      heightBase(:,k) = height(:,iSP)+(height(:,iSP+1)-height(:,iSP))/(Points(iSP)+1)*iBase
  END IF
  iBase = iBase+1
END DO

DO k=1,nRP(3) ! Iterate along each point in z direction
  CALL GetSpline(3,nCP,xBase(:,:,k),coeff(:,:,:,k),s(:,k)) ! get spline coefficients and arclength variable s

  ! Sampling of the spline to nRP(1) and allocate first row of RPs
  DO i=1,nRP(1)
    s_loc=s(nCP,k)*(i-1)/(nRP(1)-1)
    CALL EvalSpline(3,nCP,s_loc,s(:,k),coeff(:,:,:,k),x_loc)
    CALL GetNewRP(Box%RP_ptr(i,1,k)%RP,Box%GroupID,x_loc)
    RPlist_tmp(i,k)%RP=>Box%RP_ptr(i,1,k)%RP
  END DO ! i

  ! projection of the supersampled points on the closest boundary
  CALL ProjectRPtoBC(nRP(1),RPlist_tmp(:,k),NormVecRP(:,:,k))
  DO i=1,nRP(1)
    RPlist_tmp(i,k)%RP%x=RPlist_tmp(i,k)%RP%xF
    xRP(:,i,k)=RPlist_tmp(i,k)%RP%xF
  END DO
END DO

!iterate between projection and equidistant partitioning
DO iter=1,10
  ! get equidistant distribution along projected spline
  xRP_tmp = xRP
  EquiErrSum = 0.
  DO k=1,nRP(3)
    CALL GetEquiPoints(3,nRP(1),nRP(1),xRP_tmp(:,:,k),xRP(:,:,k),s_equi(:,k))
    DO i=1,nRP(1)
      RPlist_tmp(i,k)%RP%x=xRP(:,i,k)
    END DO
    ! evaluate deviation from equidistancy
    CALL EvalEquiError(3,nRP(1),xRP(:,:,k),EquiErr)
    EquiErrSum = EquiErrSum+EquiErr
  END DO
  WRITE(*,*) 'After partitioning ',EquiErrSum/nRP(3)
  ! project the equidistant points again to get their normal vector
  EquiErrSum = 0.
  DO k=1,nRP(3)
    CALL ProjectRPtoBC(nRP(1),RPlist_tmp,NormVecRP)
    DO i=1,nRP(1)
      RPlist_tmp(i,k)%RP%x=RPlist_tmp(i,k)%RP%xF
      xRP(:,i,k)=RPlist_tmp(i,k)%RP%xF
    END DO
    ! evaluate deviation from equidistancy
    CALL EvalEquiError(3,nRP(1),xRP(:,:,k),EquiErr)
    EquiErrSum = EquiErrSum + EquiErr
  END DO
  WRITE(*,*) 'After projection ',iter, EquiErrSum/nRP(3)
END DO

! calculate the tangent vector
DEALLOCATE(coeff)
ALLOCATE(s_mod(nRP(1)),coeff(3,4,nRP(1)-1,nRP(3)))
DO k=1,nRP(3) ! Iterate along each point in z direction
  CALL GetSpline(3,nRP(1),xRP(:,:,k),coeff(:,:,:,k),s_mod) ! get the spline through the projected points

  ! get the tangent vector in each point
  DO i=1,nRP(1)
    !derivative of the spline
    CALL EvalSplineDeriv(3,nRP(1),s_mod(i),s_mod,coeff(:,:,:,k),dx_loc)
    !project it on the local surface
    tangVecRP(:,i,k)=dx_loc - SUM(dx_loc(1:3)*NormVecRP(1:3,i,k))*NormVecRP(:,i,k)
    tangVecRP(:,i,k)=tangVecRP(:,i,k)/SQRT(DOT_PRODUCT(tangVecRP(:,i,k),tangVecRP(:,i,k)))
  END DO
END DO

! extrapolation of the BL mesh along the boundary normals to height
ALLOCATE(dh(2:nRP(2)))
dh(2)=1.
DO i=3,nRP(2)
  dh(i)=dh(i-1)*fac
END DO
dh(:)=dh(:)/SUM(dh(:))
iCP=2
DO k=1,nRP(3)
  DO i=1,nRP(1)
    s_loc=s(nCP,k)*(i-1)/(nRP(1)-1)
    ! linear interpolation of the height between control points
    DO WHILE(s(iCP,k).LT.(s_loc-100.*PP_RealTolerance))
      iCP=iCP+1
    END DO
    iCP=MIN(iCP,nCP)
  !  IF(s_loc.GE.s(iCP)) iCP=MAX(iCP+1,nCP)
    t=(s_loc-s(iCP-1,k))/(s(iCP,k)-s(iCP-1,k))
    height_loc=heightBase(iCP-1,k)*(1-t)+heightBase(iCP,k)*t
    h=0.
    DO j=2,nRP(2)
      h=h+dh(j)
      x_loc=Box%RP_ptr(i,1,k)%RP%x+NormVecRP(:,i,k)*h*height_loc
      CALL GetNewRP(Box%RP_ptr(i,j,k)%RP,Box%GroupID,x_loc)
    END DO
  END DO
END DO

ALLOCATE(Box%NormVec(3,nRP(1),nRP(3)))
ALLOCATE(Box%TangVec(3,nRP(1),nRP(3)))
Box%NormVec=NormVecRP
Box%TangVec=TangVecRP

DEALLOCATE(RPlist_tmp,NormVecRP,TangVecRP,coeff,s,depth,xBase,heightBase,Points)
END SUBROUTINE GetBLBox


!===================================================================================================================================
!> Project a list of RP to the closest boundary
!===================================================================================================================================
SUBROUTINE ProjectRPtoBC(nRP,RPlist_in,NormVecRP)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Parameters        ,ONLY: NSuper,maxTol
USE MOD_Interpolation     ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars,ONLY: NodeType
USE MOD_RPSet_Vars        ,ONLY: tRP,tRPlist
USE MOD_Mesh_Vars,         ONLY: SideToElem,nBCSides,Face_xGP,NormVec,NGeo
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys,ChebyGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis,             ONLY: PolynomialDerivativeMatrix
USE MOD_ChangeBasis,       ONLY: ChangeBasis2D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: nRP
TYPE(tRPlist)                   :: RPlist_in(nRP)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: NormVecRP(3,nRP)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: Winner_Dist2,Dist2
REAL                  :: Xi_NSuper(0:NSuper)
REAL                  :: xRP(3)
INTEGER               :: i,j,l,iElem,nNodes
INTEGER               :: iRP,NewtonIter
TYPE(tRP),POINTER     :: aRP
LOGICAL               :: calcJacobianDone
INTEGER               :: SideID,locSideID,iWinner,jWinner,iRP2
REAL                  :: dist2RP(nRP)
REAL                  :: Vdm_GP_EquiNSuper(0:NSuper,0:PP_N),xBC_NSuper(3,0:NSuper,0:NSuper)
REAL                  :: NormVec_NSuper(3,0:NSuper,0:NSuper)
REAL                  :: wBary_NSuper(0:NSuper), D_NSuper(0:NSuper,0:NSuper), Lag_NSuper(1:2,0:NSuper)
REAL                  :: dxBC_NSuper(3,2,0:NSuper,0:NSuper)
REAL                  :: Gmat(2,0:NSuper,0:NSuper),dGmat(2,2,0:NSuper,0:NSuper)
REAL                  :: G(2),Xi2(2),Jac2(2,2),sJac2(2,2),xWinner(3),NormVecWinner(3)
REAL                  :: F(1:3),eps_F

!===================================================================================================================================
SWRITE(UNIT_StdOut,'(A,I4,A)')' Project ',nRP,' RPs on the closest boundary...'
iRP2=0
dist2RP=HUGE(1.)
NormVecRP=0.
! Prepare basis
IF(NSuper.LT.Ngeo*2) &
  WRITE(*,*)'Warning: NSuper<2Ngeo, derivative may be wrong.'

DO i=0,NSuper
  Xi_NSuper(i) = 2./REAL(NSuper) * REAL(i) - 1.
END DO
CALL BarycentricWeights(NSuper,Xi_NSuper,wBary_NSuper)
CALL GetVandermonde(PP_N,NodeType,NSuper,'VISU',Vdm_GP_EquiNSuper)
CALL PolynomialDerivativeMatrix(NSuper,Xi_NSuper,D_NSuper)
nNodes=(NSuper+1)**2

DO SideID=1,nBCSides
  calcJacobianDone=.FALSE.
  ! Supersampling of the side to equidistant grid for search
  CALL ChangeBasis2D(3,PP_N,NSuper,Vdm_GP_EquiNSuper,Face_xGP(:,:,:,0,SideID),xBC_NSuper)
  CALL ChangeBasis2D(3,PP_N,NSuper,Vdm_GP_EquiNSuper,NormVec( :,:,:,0,SideID),NormVec_NSuper)
  DO iRP=1,nRP
    aRP=>RPlist_in(iRP)%RP
    xRP=aRP%x
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
    xWinner(:)=xBC_NSuper(:,iWinner,jWinner)
    NormVecWinner(:)=NormVec_NSuper(:,iWinner,jWinner)
    F=xRP-xWinner
    NormVecWinner=NormVecWinner/SQRT(DOT_PRODUCT(NormVecWinner,NormVecWinner))
    F=(xRP-xWinner)/SQRT(DOT_PRODUCT(xRP-xWinner,xRP-xWinner))


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

    ! Get initial value of the functional G
    CALL LagrangeInterpolationPolys(Xi2(1),NSuper,Xi_NSuper,wBary_NSuper,Lag_NSuper(1,:))
    CALL LagrangeInterpolationPolys(Xi2(2),NSuper,Xi_NSuper,wBary_NSuper,Lag_NSuper(2,:))
    G=0.
    DO j=0,NSuper
      DO i=0,NSuper
        G=G+Gmat(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
      END DO! i=0,NSuper
    END DO! j=0,NSuper

    ! Start Newton
    eps_F=1.E-10*(SUM(G*G))
    NewtonIter=0
    DO WHILE ((SUM(G*G).GT.eps_F).AND.(NewtonIter.LT.50))
      NewtonIter=NewtonIter+1

      ! Compute G Jacobian dG/dXi
      Jac2=0.
      DO j=0,NSuper
        DO i=0,NSuper
          Jac2=Jac2 + dGmat(:,:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
        END DO !l=0,NSuper
      END DO !i=0,NSuper

      ! Compute inverse of Jacobian
      sJac2=getInv2(Jac2)

      ! Iterate Xi using Newton step
      Xi2 = Xi2 - MATMUL(sJac2,G)
      ! if Newton gets outside reference space range [-1,1], exit.
      ! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
      IF((NewtonIter.GT.4).AND.(ANY(ABS(Xi2).GT.1.2))) EXIT

      ! Compute function value
      CALL LagrangeInterpolationPolys(Xi2(1),NSuper,Xi_NSuper,wBary_NSuper,Lag_NSuper(1,:))
      CALL LagrangeInterpolationPolys(Xi2(2),NSuper,Xi_NSuper,wBary_NSuper,Lag_NSuper(2,:))
      ! Exit if we are far enough outside of [-1,1] for the basis to reach 'Infinity' overflow
      IF (ANY(ABS(Lag_NSuper(:,:)).GT.(HUGE(1.)))) EXIT
      G=0.
      DO j=0,NSuper
       DO i=0,NSuper
         G=G+Gmat(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
       END DO! i=0,NSuper
     END DO! j=0,NSuper
    END DO !newton
    ! use Newton result if minimum is within parameter range, else see if supersampled
    ! initial guess is better than previous result
    IF(MAXVAL(ABS(Xi2)).LE.maxTol) THEN ! use newton result
      ! calculate new distance and normal vector
      xWinner=0.
      NormVecWinner=0.
      DO j=0,NSuper
        DO i=0,NSuper
          xWinner(:)=xWinner(:)+xBC_NSuper(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
          NormVecWinner(:)=NormVecWinner(:)+NormVec_NSuper(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
        END DO! i=0,NSuper
      END DO! j=0,NSuper
      Winner_Dist2=SUM((xWinner-xRP)*(xWinner-xRP))
    END IF

    !! Check the angle, should approach the normal of the face
    !NormVecWinner=NormVecWinner/NORM2(NormVecWinner)
    !F=(xRP-xWinner)/NORM2(xRP-xWinner)
    !ang=ACOS(ABS(SUM(F*NormVecWinner)))
    IF((Winner_Dist2.LE.dist2RP(iRP)))THEN
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
      ! find the side with minimum distance
      dist2RP(iRP)=Winner_Dist2
      NormVecRP(:,iRP)=-NormVecWinner
!       SWRITE(UNIT_StdOut,'(A,I4,A,3F8.4,A,3F8.4,A)')' Projected RP ',iRP,' with coordinates ',aRP%x,' to ',ARP%xF,'.'
    END IF
  END DO! iRP=1,nRP
END DO! SideID=1,nBCSides

SWRITE(UNIT_StdOut,'(A)')' done.'
SWRITE(UNIT_StdOut,'(A,F15.8)')'  Max. distance: ',SQRT(MAXVAL(dist2RP))
END SUBROUTINE ProjectRPtoBC



!=================================================================================================================================
!> Computes the inverse of a 2x2 matrix
!=================================================================================================================================
FUNCTION getInv2(Mat)
! MODULES
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(2,2)
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getInv2(2,2)
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sdet
!=================================================================================================================================
sdet=1./(Mat(1,1) * Mat(2,2) - Mat(1,2)*Mat(2,1))
getInv2(1,1) = (  Mat(2,2) ) * sdet
getInv2(1,2) = (- Mat(1,2) ) * sdet
getInv2(2,1) = (- Mat(2,1) ) * sdet
getInv2(2,2) = (  Mat(1,1) ) * sdet
END FUNCTION getInv2

END MODULE MOD_BLProjection
