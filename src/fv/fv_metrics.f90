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
#if FV_ENABLED
#include "flexi.h"

!==================================================================================================================================
!> Calculate metric terms for FV subcells. (see also metrics.f90: a lot of FV metric terms are computed there)
!==================================================================================================================================
MODULE MOD_FV_Metrics
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitFV_Metrics
  MODULE PROCEDURE InitFV_Metrics
END INTERFACE

INTERFACE FinalizeFV_Metrics
  MODULE PROCEDURE FinalizeFV_Metrics
END INTERFACE

PUBLIC::InitFV_Metrics
PUBLIC::FinalizeFV_Metrics
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute the remaining metric terms for FV subcells, that are not computed in metrics.f90. 
!> Normal, tangential vectors, SurfElems, ... for FV subcells.
!==================================================================================================================================
SUBROUTINE InitFV_Metrics()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars
USE MOD_FV_Basis
USE MOD_Mesh_Vars          ,ONLY: nElems,nSides,firstMPISide_YOUR,lastMPISide_YOUR
USE MOD_Mesh_Vars          ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,sJ
USE MOD_Mesh_Vars          ,ONLY: NGeoRef,DetJac_Ref,MortarType,MortarInfo
USE MOD_Mesh_Vars          ,ONLY: NormalDirs,TangDirs,NormalSigns,SideToElem
USE MOD_Mesh_Vars          ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP,Ja_Face
USE MOD_ChangeBasis        ,ONLY: ChangeBasis1D,ChangeBasis2D,ChangeBasis3D
USE MOD_Metrics            ,ONLY: SurfMetricsFromJa
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
#if FV_RECONSTRUCT
USE MOD_Mesh_Vars          ,ONLY: firstBCSide,lastBCSide,firstInnerSide,lastMPISide_MINE
USE MOD_Mesh_Vars          ,ONLY: S2V2,ElemToSide,dXCL_N
USE MOD_Interpolation      ,ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars ,ONLY: NodeTypeG,NodeTypeCL,xGP
#if USE_MPI
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_MPI_Vars           ,ONLY: nNbProcs
#endif
#endif
USE MOD_FillMortar1        ,ONLY: U_Mortar1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                :: i,j,k,p,q,l,d,iSide,iElem,iLocSide
INTEGER                                :: dd,NormalDir,TangDir
REAL                                   :: NormalSign
REAL                                   :: Vdm_Gauss_FVboundary(0:PP_N+1,0:PP_N)
REAL                                   :: JaVol(3,3,0:PP_N+1,0:PP_N+1,0:PP_N+1)
REAL                                   :: FV_Ja_Face(3,3,0:PP_N,0:PP_N)
REAL,DIMENSION(0:PP_N,0:NgeoRef)       :: Vdm_NgeoRef_N,Vdm_NgeoRef_FV
REAL                                   :: FV_DetJac(1,0:PP_N,0:PP_N,0:PP_N)
INTEGER                                :: flip, SideID, iMortar
#if FV_RECONSTRUCT
INTEGER                                :: ijk(2), locSideID
REAL                                   :: FV_dx_Face(0:PP_N,0:PP_N,1:3)
REAL                                   :: DG_dx_slave (1,0:PP_N,0:PP_N,1:nSides)
REAL                                   :: DG_dx_master(1,0:PP_N,0:PP_N,1:nSides)

REAL                                   :: tmp2(3,0:PP_N)
REAL,DIMENSION(0:PP_N)                 :: xG,wG,wBaryG
REAL,DIMENSION(0:PP_N,0:PP_N)          :: Vdm_CLN_FV, Vdm_CLN_GaussN,length
REAL,DIMENSION(3,0:PP_N,0:PP_N,0:PP_N) :: FV_Path_XI, FV_Path_ETA, FV_Path_ZETA
REAL                                   :: x0, xN
REAL,POINTER                           :: FV_dx_P(:,:)
#if USE_MPI
INTEGER                                :: MPIRequest(nNbProcs,2)
#endif
#endif
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  Build Metrics ...'

#if FV_RECONSTRUCT
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! All the following distances are computed in PHYSICAL space, even so the variables are named in reference space notation.
! since they are used in a '1D tensor-product'-way (first index is xi-, second eta-, third zeta-direction)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ---------------------------------
! |       |       |       |       | 
! |<->.   |<->.   |<->.   |<->.   |    FV_dx_XI_L: left distance between Face and center of FV subcell 
! |       |       |       |       |                (see Attention to storage order below!)
! ---------------------------------
! |       |       |       |       |
! |   .<->|   .<->|   .<->|   .<->|    FV_dx_XI_R: right distance between subcell-Faces and center of FV subcell
! |       |       |       |       |                (see Attention to storage order below!)
! ---------------------------------
! |       |       |       |       |
! |   .<----->.<----->.<----->.   |    FV_dx_XI: distance between centers of FV subcells 
! |       |       |       |       |              computed as FV_dx_XI_R(i-1) + FV_dx_XI_L(i)
! ---------------------------------
! |       |       |       |       |    DG_dx_master/slave: distance from DG interface to first Gauss point  
! |<->.   |   .   |   .   |   .<->|    FV_dx_master/slave: distance from DG interface to first center of FV subell 
! |       |       |       |       |                        stored in Face-based coordinates
! ---------------------------------
! |       |       |       |       |
! |   .   |   .   |   .   |   .<------>.  FV_dx_Face: -both FV:        FV_dx_master + FV_dx_slave
! |       |       |       |       |                   -mixed DG/FV:    FV_dx_master + DG_dx_slave 
! ---------------------------------                                 or DG_dx_master + FV_dx_slave
! |       ^       ^       ^       |
! |   .   |   .   |   .   |   .   |    FV_SurfElemXi_sw: SurfElem of !inner! subcell faces (see Attention to storage order below!)
! |       v       v       v       | 
! ---------------------------------
! |       |       |       |       |
! |   .   x   .   x   .   x   .   |    FV_NormVecXi/TangVec.Xi: normal/tangent vectors at !inner! subcell faces (positions x)
! |       |       |       |       |                             (see Attention to storage order below!)
! ---------------------------------

ALLOCATE(FV_sdx_XI    (0:PP_N,0:PP_N,1:PP_N,nElems)) ! 1. / FV_dx_XI    Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_sdx_ETA   (0:PP_N,0:PP_N,1:PP_N,nElems)) ! 1. / FV_dx_ETA   Attention: storage order is (i,k,j,iElem)
ALLOCATE(FV_sdx_ZETA  (0:PP_N,0:PP_N,1:PP_N,nElems)) ! 1. / FV_dx_ZETA  Attention: storage order is (i,j,k,iElem)

ALLOCATE(FV_sdx_Face (0:PP_N,0:PP_N,1:3,1:lastMPISide_MINE)) ! 1. / FV_dx_Face

ALLOCATE(FV_dx_XI_L  (0:PP_N,0:PP_N,0:PP_N,nElems))  ! Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_dx_XI_R  (0:PP_N,0:PP_N,0:PP_N,nElems))  ! Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_dx_ETA_L (0:PP_N,0:PP_N,0:PP_N,nElems))  ! Attention: storage order is (i,k,j,iElem)
ALLOCATE(FV_dx_ETA_R (0:PP_N,0:PP_N,0:PP_N,nElems))  ! Attention: storage order is (i,k,j,iElem)
ALLOCATE(FV_dx_ZETA_L(0:PP_N,0:PP_N,0:PP_N,nElems))  ! Attention: storage order is (i,j,k,iElem)
ALLOCATE(FV_dx_ZETA_R(0:PP_N,0:PP_N,0:PP_N,nElems))  ! Attention: storage order is (i,j,k,iElem)

ALLOCATE(FV_dx_slave (1,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(FV_dx_master(1,0:PP_N,0:PP_N,1:nSides))
#endif
ALLOCATE(FV_SurfElemXi_sw  (0:PP_N,0:PP_N,1:PP_N,nElems)) ! Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_SurfElemEta_sw (0:PP_N,0:PP_N,1:PP_N,nElems)) ! Attention: storage order is (i,k,j,iElem)
ALLOCATE(FV_SurfElemZeta_sw(0:PP_N,0:PP_N,1:PP_N,nElems)) ! Attention: storage order is (i,j,k,iElem)

ALLOCATE(FV_NormVecXi   (3,0:PP_N,0:PP_N,1:PP_N,nElems))  ! Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_TangVec1Xi  (3,0:PP_N,0:PP_N,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_TangVec2Xi  (3,0:PP_N,0:PP_N,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_NormVecEta  (3,0:PP_N,0:PP_N,1:PP_N,nElems))  ! Attention: storage order is (i,k,j,iElem)
ALLOCATE(FV_TangVec1Eta (3,0:PP_N,0:PP_N,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_TangVec2Eta (3,0:PP_N,0:PP_N,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_NormVecZeta (3,0:PP_N,0:PP_N,1:PP_N,nElems))  ! Attention: storage order is (i,j,k,iElem)
ALLOCATE(FV_TangVec1Zeta(3,0:PP_N,0:PP_N,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_TangVec2Zeta(3,0:PP_N,0:PP_N,1:PP_N,nElems))  !  -"-

#if PARABOLIC
ALLOCATE(FV_Metrics_fTilde_sJ(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(FV_Metrics_gTilde_sJ(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(FV_Metrics_hTilde_sJ(3,0:PP_N,0:PP_N,0:PP_N,nElems))
#endif

ALLOCATE(FV_Elems_master(1:nSides)) ! Moved from InitFV to here, since needed in U_Mortar below.

! compute FV NormVec, TangVec,.. on boundary of DG-cells
DO iSide=1,nSides
  IF(iSide.GE.firstMPISide_YOUR.AND.iSide.LE.lastMPISide_YOUR) CYCLE
  CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,Face_xGP(:,:,:,0,iSide),Face_xGP(:,:,:,1,iSide))

  iLocSide=SideToElem(S2E_LOC_SIDE_ID,iSide)
  IF(iLocSide.LT.1) CYCLE
  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide)

  CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,Ja_Face( 1,:,:,:,iSide),FV_Ja_Face(1,:,:,:))
  CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,Ja_Face( 2,:,:,:,iSide),FV_Ja_Face(2,:,:,:))
  CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,Ja_Face( 3,:,:,:,iSide),FV_Ja_Face(3,:,:,:))
  CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
                         NormVec( :,:,:,1,iSide),TangVec1(:,:,:,1,iSide),&
                         TangVec2(:,:,:,1,iSide),SurfElem(  :,:,1,iSide))

  IF(MortarType(1,iSide).LE.0) CYCLE ! no mortars
  DO iMortar=1,MERGE(4,2,MortarType(1,iSide).EQ.1)
    SideID=MortarInfo(MI_SIDEID,iMortar,MortarType(2,iSide))
    Flip  =MortarInfo(MI_FLIP,iMortar,MortarType(2,iSide))
    IF(flip.NE.0) CYCLE ! for MPI sides some sides are built from the inside and for type 2/3 there are only 2 neighbours
    CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,Ja_Face( 1,:,:,:,SideID),FV_Ja_Face(1,:,:,:))
    CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,Ja_Face( 2,:,:,:,SideID),FV_Ja_Face(2,:,:,:))
    CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,Ja_Face( 3,:,:,:,SideID),FV_Ja_Face(3,:,:,:))
    CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
                           NormVec( :,:,:,1,SideID),TangVec1(:,:,:,1,SideID),&
                           TangVec2(:,:,:,1,SideID),SurfElem(  :,:,1,SideID))
  END DO
END DO

CALL FV_Build_Vdm_Gauss_FVboundary(PP_N, Vdm_Gauss_FVboundary)
CALL GetVandermonde(NgeoRef, NodeType, PP_N, NodeType, Vdm_NgeoRef_N, modal=.TRUE.)
Vdm_NgeoRef_FV = MATMUL(FV_Vdm, Vdm_NgeoRef_N)
DO iElem=1,nElems
  ! compute Jacobian
  CALL ChangeBasis3D(1,NGeoRef,PP_N,Vdm_NgeoRef_FV,DetJac_Ref(:,:,:,:,iElem),FV_DetJac)
  sJ(:,:,:,iElem,1) = 1./FV_DetJac(1,:,:,:)
  CALL ChangeBasis3D(3,PP_N,PP_N,FV_Vdm,Metrics_fTilde(:,:,:,:,iElem,0),Metrics_fTilde(:,:,:,:,iElem,1))
  CALL ChangeBasis3D(3,PP_N,PP_N,FV_Vdm,Metrics_gTilde(:,:,:,:,iElem,0),Metrics_gTilde(:,:,:,:,iElem,1))
  CALL ChangeBasis3D(3,PP_N,PP_N,FV_Vdm,Metrics_hTilde(:,:,:,:,iElem,0),Metrics_hTilde(:,:,:,:,iElem,1))

  !================================================================
  ! Compute FV NormVec,TangVec,.. at inner cell boundaries...START
  !================================================================

  ! Xi direction
  DO k=0,PP_N; DO j=0,PP_N
    ! interpolate Metrics to boundaries of FV subcells in XI direction, other directions stay DG
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_fTilde(:,:,j,k,iElem,0),JaVol(1,:,:,j,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_gTilde(:,:,j,k,iElem,0),JaVol(2,:,:,j,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_hTilde(:,:,j,k,iElem,0),JaVol(3,:,:,j,k))
  END DO; END DO ! j,k=0,PP_N
  DO l=1,PP_N
    ! at every inner interface/slice between FV subcells in XI direction: 
    ! convert metrics in the other directions from DG to FV subcells
    DO dd=1,3
      CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,JaVol(dd,1:3,l,0:PP_N,0:PP_N),FV_Ja_Face(dd,:,:,:))
    END DO
    ! use metrics to build normal/tangential vectors and surelem at the inner interfaces/slices
    NormalDir=1; TangDir=2; NormalSign=1.
    CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
        FV_NormVecXi  (:,:,:,l,iElem),&
        FV_TangVec1Xi (:,:,:,l,iElem),&
        FV_TangVec2Xi (:,:,:,l,iElem),&
        FV_SurfElemXi_sw(:,:,l,iElem))
    ! multiply FV_SurfElemXi with 1/FV_w
    FV_SurfElemXi_sw(:,:,l,iElem) = FV_SurfElemXi_sw(:,:,l,iElem)*FV_w_inv
  END DO

  ! Eta direction
  DO k=0,PP_N; DO i=0,PP_N
    ! interpolate Metrics to boundaries of FV subcells in ETA direction, other directions stay DG
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_fTilde(:,i,:,k,iElem,0),JaVol(1,:,i,:,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_gTilde(:,i,:,k,iElem,0),JaVol(2,:,i,:,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_hTilde(:,i,:,k,iElem,0),JaVol(3,:,i,:,k))
  END DO; END DO ! i,k=0,PP_N
  DO l=1,PP_N
    ! at every inner interface/slice between FV subcells in ETA direction: 
    ! convert metrics in the other directions from DG to FV subcells
    DO dd=1,3
      CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,JaVol(dd,1:3,0:PP_N,l,0:PP_N),FV_Ja_Face(dd,:,:,:))
    END DO
    ! use metrics to build normal/tangential vectors and surelem at the inner interfaces/slices
    NormalDir=2; TangDir=3; NormalSign=1.
    CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
        FV_NormVecEta  (:,:,:,l,iElem),&
        FV_TangVec1Eta (:,:,:,l,iElem),&
        FV_TangVec2Eta (:,:,:,l,iElem),&
        FV_SurfElemEta_sw(:,:,l,iElem))
    ! multiply FV_SurfElemEta with 1/FV_w
    FV_SurfElemEta_sw(:,:,l,iElem) = FV_SurfElemEta_sw(:,:,l,iElem)*FV_w_inv
  END DO

  ! Zeta direction
  DO j=0,PP_N; DO i=0,PP_N
    ! interpolate Metrics to boundaries of FV subcells in ZETA direction, other directions stay DG
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_fTilde(:,i,j,:,iElem,0),JaVol(1,:,i,j,:))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_gTilde(:,i,j,:,iElem,0),JaVol(2,:,i,j,:))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_hTilde(:,i,j,:,iElem,0),JaVol(3,:,i,j,:))
  END DO; END DO ! i,j=0,PP_N
  DO l=1,PP_N
    ! at every inner interface/slice between FV subcells in ZETA direction: 
    ! convert metrics in the other directions from DG to FV subcells
    DO dd=1,3
      CALL ChangeBasis2D(3,PP_N,PP_N,FV_Vdm,JaVol(dd,1:3,0:PP_N,0:PP_N,l),FV_Ja_Face(dd,:,:,:))
    END DO
    ! use metrics to build normal/tangential vectors and surelem at the inner interfaces/slices
    NormalDir=3; TangDir=1; NormalSign=1.
    CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
        FV_NormVecZeta  (:,:,:,l,iElem),&
        FV_TangVec1Zeta (:,:,:,l,iElem),&
        FV_TangVec2Zeta (:,:,:,l,iElem),&
        FV_SurfElemZeta_sw(:,:,l,iElem))
    ! multiply FV_SurfElemZeta with 1/FV_w
    FV_SurfElemZeta_sw(:,:,l,iElem) = FV_SurfElemZeta_sw(:,:,l,iElem)*FV_w_inv
  END DO
  !==================================================================
  ! Compute FV NormVec,TangVec,.. at inner cell boundaries...DONE
  !==================================================================
END DO 


#if FV_RECONSTRUCT
!=====================================================================================
! Compute distances between FV subcells and between first Gauss point and interface
! Integrate path given by dXCL_N, to ensure free stream preservation for N<Ngeo
!=====================================================================================
CALL GetVandermonde(PP_N,NodeTypeCL,PP_N,NodeTypeG, Vdm_CLN_GaussN, modal=.FALSE.)
Vdm_CLN_FV = MATMUL(FV_Vdm, Vdm_CLN_GaussN)
CALL GetNodesAndWeights(PP_N,NodeTypeG,xG,wG,wBaryG)

DO iElem=1,nElems
  DO l=0,PP_N
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_FV, dXCL_N(1,:,l,:,:,iElem), FV_Path_XI  (:,l,:,:))
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_FV, dXCL_N(2,:,:,l,:,iElem), FV_Path_ETA (:,l,:,:))
    CALL ChangeBasis2D(3,PP_N,PP_N,Vdm_CLN_FV, dXCL_N(3,:,:,:,l,iElem), FV_Path_ZETA(:,l,:,:))
  END DO ! i=0,PP_N
  DO q=0,PP_N; DO p=0,PP_N
    tmp2 = FV_Path_XI(:,:,p,q)
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN, tmp2, FV_Path_XI(:,:,p,q))
    tmp2 = FV_Path_ETA(:,:,p,q)
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN, tmp2, FV_Path_ETA(:,:,p,q))
    tmp2 = FV_Path_ZETA(:,:,p,q)
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN, tmp2, FV_Path_ZETA(:,:,p,q))
  END DO; END DO! p,q=0,PP_N

  ! Calculate distances between FV subcells
  DO l=0,PP_N
    ! left
    x0 = FV_BdryX(l)
    xN = x0 + (FV_BdryX(l+1) - FV_BdryX(l)) * 0.5
    CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, x0,xN, FV_Path_XI  , FV_dx_XI_L  (:,:,l,iElem))
    CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, x0,xN, FV_Path_ETA , FV_dx_ETA_L (:,:,l,iElem))
    CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, x0,xN, FV_Path_ZETA, FV_dx_ZETA_L(:,:,l,iElem))
    
    ! right
    xN = FV_BdryX(l+1)
    x0 = xN - (FV_BdryX(l+1) - FV_BdryX(l)) * 0.5
    CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, x0,xN, FV_Path_XI  , FV_dx_XI_R  (:,:,l,iElem))
    CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, x0,xN, FV_Path_ETA , FV_dx_ETA_R (:,:,l,iElem))
    CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, x0,xN, FV_Path_ZETA, FV_dx_ZETA_R(:,:,l,iElem))
  END DO ! l=0,PP_N

  ! build inverse distances (volumes)
  FV_sdx_XI  (:,:,:,iElem) = 1. / (FV_dx_XI_R  (:,:,0:PP_N-1,iElem) +  FV_dx_XI_L  (:,:,1:PP_N,iElem)) ! 1 / FV_dx_XI
  FV_sdx_ETA (:,:,:,iElem) = 1. / (FV_dx_ETA_R (:,:,0:PP_N-1,iElem) +  FV_dx_ETA_L (:,:,1:PP_N,iElem)) ! 1 / FV_dx_ETA
  FV_sdx_ZETA(:,:,:,iElem) = 1. / (FV_dx_ZETA_R(:,:,0:PP_N-1,iElem) +  FV_dx_ZETA_L(:,:,1:PP_N,iElem)) ! 1 / FV_dx_ZETA

  ! Calculate distance between first GaussPoint and interface
  DO locSideID=1,6
    length=0.
    SideID = ElemToSide(E2S_SIDE_ID,locSideID,iElem)
    flip   = ElemToSide(E2S_FLIP,   locSideID,iElem)
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, -1., xGP(0),   FV_Path_XI  , length)
      FV_dx_P(0:PP_N,0:PP_N) => FV_dx_XI_L(:,:,0,iElem)
    CASE(ETA_MINUS)
      CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, -1., xGP(0),   FV_Path_ETA , length)
      FV_dx_P(0:PP_N,0:PP_N) => FV_dx_ETA_L(:,:,0,iElem)
    CASE(ZETA_MINUS)
      CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, -1., xGP(0),   FV_Path_ZETA, length)
      FV_dx_P(0:PP_N,0:PP_N) => FV_dx_ZETA_L(:,:,0,iElem)
    CASE(XI_PLUS)
      CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, xGP(PP_N), 1., FV_Path_XI  , length)
      FV_dx_P(0:PP_N,0:PP_N) => FV_dx_XI_R(:,:,PP_N,iElem)
    CASE(ETA_PLUS)
      CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, xGP(PP_N), 1., FV_Path_ETA , length)
      FV_dx_P(0:PP_N,0:PP_N) => FV_dx_ETA_R(:,:,PP_N,iElem)
    CASE(ZETA_PLUS)
      CALL Integrate_Path(PP_N,PP_N, xG,wG,wBaryG, xGP(PP_N), 1., FV_Path_ZETA, length)
      FV_dx_P(0:PP_N,0:PP_N) => FV_dx_ZETA_R(:,:,PP_N,iElem)
    CASE DEFAULT
      STOP 'Local side index out of range (1...6).'
    END SELECT

    IF (flip.EQ.0) THEN ! master side
      DO q=0,PP_N; DO p=0,PP_N
        ijk = S2V2(:,p,q,flip,locSideID)
        DG_dx_master(1,p,q,SideID) = length(ijk(1),ijk(2))
        FV_dx_master(1,p,q,SideID) = FV_dx_P(ijk(1),ijk(2))
      END DO; END DO
    ELSE ! slave side
      DO q=0,PP_N; DO p=0,PP_N
        ijk = S2V2(:,p,q,flip,locSideID)
        DG_dx_slave(1,p,q,SideID) = length(ijk(1),ijk(2))
        FV_dx_slave(1,p,q,SideID) = FV_dx_P(ijk(1),ijk(2))
      END DO; END DO
    END IF
  END DO
END DO

! distances at big mortar interfaces must be distributed to the smaller sides
FV_Elems_master = 1 ! Force use of FV mortar matrices in U_Mortar routine
#if USE_MPI
MPIRequest=0
! distances at MPI slave sides must be transmitted to master sides 
CALL U_Mortar1(FV_dx_master,FV_dx_slave,doMPISides=.TRUE.)
CALL StartReceiveMPIData(FV_dx_slave, (PP_N+1)**2, 1,nSides,MPIRequest(:,SEND),SendID=2)
CALL StartSendMPIData(   FV_dx_slave, (PP_N+1)**2, 1,nSides,MPIRequest(:,RECV),SendID=2)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest) !Send MINE -receive YOUR
#endif
CALL U_Mortar1(FV_dx_master,FV_dx_slave,doMPISides=.FALSE.)

FV_Elems_master = 0 ! Force use of DG mortar matrices in U_Mortar routine
#if USE_MPI
CALL U_Mortar1(DG_dx_master,DG_dx_slave,doMPISides=.TRUE.)
CALL StartReceiveMPIData(DG_dx_slave, (PP_N+1)**2, 1,nSides,MPIRequest(:,SEND),SendID=2)
CALL StartSendMPIData(   DG_dx_slave, (PP_N+1)**2, 1,nSides,MPIRequest(:,RECV),SendID=2)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest) !Send MINE -receive YOUR
#endif
CALL U_Mortar1(DG_dx_master,DG_dx_slave,doMPISides=.FALSE.)

! calculate distances of boundary sides
DO SideID=firstBCSide,lastBCSide
  FV_sdx_Face(:,:,1,SideID) = -99. ! dummy, should never be used
  FV_sdx_Face(:,:,2,SideID) = -99. ! dummy, should never be used
  FV_sdx_Face(:,:,3,SideID) = 1. / FV_dx_master(1,:,:,SideID)
END DO

! calculate distances of inner and MPI_MINE sides
DO SideID=firstInnerSide,lastMPISide_MINE
  ! master=FV, slave=DG 
  FV_dx_Face(:,:,1) = DG_dx_slave(1,:,:,SideID) + FV_dx_master(1,:,:,SideID)
  ! master=DG, sla                                                                                
  FV_dx_Face(:,:,2) = FV_dx_slave(1,:,:,SideID) + DG_dx_master(1,:,:,SideID)
  ! master=FV, sla                                                                                
  FV_dx_Face(:,:,3) = FV_dx_slave(1,:,:,SideID) + FV_dx_master(1,:,:,SideID)
  ! precompute inverse
  FV_sdx_Face(:,:,:,SideID) = 1. / FV_dx_Face
END DO

! scale metrics for equidistant subcells
#if PARABOLIC
DO iElem=1,nElems
  DO d=1,3
    DO i=0,PP_N
      FV_Metrics_fTilde_sJ(d,i,:,:,iElem)=FV_w_inv*Metrics_fTilde(d,i,:,:,iElem,1)*&
          (FV_dx_XI_L  (:,:,i,iElem)+FV_dx_XI_R  (:,:,i,iElem))
    END DO ! i=0,PP_N
    DO j=0,PP_N
      FV_Metrics_gTilde_sJ(d,:,j,:,iElem)=FV_w_inv*Metrics_gTilde(d,:,j,:,iElem,1)*&
          (FV_dx_ETA_L (:,:,j,iElem)+FV_dx_ETA_R (:,:,j,iElem))
    END DO ! j=0,PP_N
    FV_Metrics_hTilde_sJ(d,:,:,:,iElem)=FV_w_inv*Metrics_hTilde(d,:,:,:,iElem,1)*&
        (FV_dx_ZETA_L(:,:,:,iElem)+FV_dx_ZETA_R(:,:,:,iElem))
    FV_Metrics_fTilde_sJ(d,:,:,:,iElem)=FV_Metrics_fTilde_sJ(d,:,:,:,iElem)*sJ(:,:,:,iElem,1)
    FV_Metrics_gTilde_sJ(d,:,:,:,iElem)=FV_Metrics_gTilde_sJ(d,:,:,:,iElem)*sJ(:,:,:,iElem,1)
    FV_Metrics_hTilde_sJ(d,:,:,:,iElem)=FV_Metrics_hTilde_sJ(d,:,:,:,iElem)*sJ(:,:,:,iElem,1)
  END DO
END DO
#endif /* PARABOLIC */
#endif /* FV_RECONSTRUCT */

SWRITE(UNIT_stdOut,'(A)')' Done !'

END SUBROUTINE InitFV_Metrics


#if FV_RECONSTRUCT
!==================================================================================================================================
!> Computes the distance between two points along a path given in 1D reference coordinates
!==================================================================================================================================
SUBROUTINE Integrate_Path(Nloc,Nloc2,xGP,wGP,wBary,x0,xN,FV_Path_1D,FV_Length)
! MODULES
USE MOD_Basis       ,ONLY: InitializeVandermonde
USE MOD_ChangeBasis ,ONLY: ChangeBasis1D
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN) :: Nloc                                 !< degree of path polynomial
INTEGER,INTENT(IN) :: Nloc2                                !< number of points to compute (Nloc2+1)**2
REAL,INTENT(IN)    :: xGP(  0:Nloc)                        !< parametric coords
REAL,INTENT(IN)    :: wGP(  0:Nloc)                        !< integration weights
REAL,INTENT(IN)    :: wBary(0:Nloc)                        !< interpolations weights
REAL,INTENT(IN)    :: x0                                   !< start point
REAL,INTENT(IN)    :: xN                                   !< end point
REAL,INTENT(IN)    :: FV_Path_1D(3,0:Nloc,0:Nloc2,0:Nloc2) !< path polynomial
REAL,INTENT(OUT)   :: FV_Length(          0:Nloc2,0:Nloc2) !< distance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: VDM(0:Nloc,0:Nloc)
REAL               :: SubxGP(1,0:Nloc)
REAL               :: FV_Path_Cut(3,0:Nloc)
INTEGER            :: q,p,l
!===================================================================================================================================
subxGP(1,:) = x0 + (xGP + 1.)/2. * (xN-x0)
CALL InitializeVandermonde(Nloc,Nloc,wBary,xGP,subxGP(1,:),Vdm)

FV_Length=0.
DO q=0,Nloc2; DO p=0,Nloc2
  ! path to integrate in ref space [-1,1]
  CALL ChangeBasis1D(3,Nloc,Nloc,Vdm,FV_Path_1D(:,:,p,q), FV_Path_Cut)
  ! integrate path
  DO l=0,Nloc
    FV_Length(p,q) = FV_Length(p,q) + NORM2(FV_Path_Cut(:,l)) * wGP(l)
  END DO      
END DO; END DO! p,q=0,Nloc2
FV_Length=FV_Length*0.5*(xN-x0) ! *0.5 since reference element has width=2
END SUBROUTINE Integrate_Path
#endif


!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeFV_Metrics()
! MODULES
USE MOD_FV_Vars
IMPLICIT NONE
!==================================================================================================================================
#if FV_RECONSTRUCT
SDEALLOCATE(FV_sdx_XI)
SDEALLOCATE(FV_sdx_ETA)
SDEALLOCATE(FV_sdx_ZETA)

SDEALLOCATE(FV_sdx_Face)

SDEALLOCATE(FV_dx_XI_L)
SDEALLOCATE(FV_dx_XI_R)
SDEALLOCATE(FV_dx_ETA_L)
SDEALLOCATE(FV_dx_ETA_R)
SDEALLOCATE(FV_dx_ZETA_L)
SDEALLOCATE(FV_dx_ZETA_R)

SDEALLOCATE(FV_dx_slave)
SDEALLOCATE(FV_dx_master)
#endif

SDEALLOCATE(FV_SurfElemXi_sw)
SDEALLOCATE(FV_SurfElemEta_sw)
SDEALLOCATE(FV_SurfElemZeta_sw)

SDEALLOCATE(FV_NormVecXi)
SDEALLOCATE(FV_TangVec1Xi)
SDEALLOCATE(FV_TangVec2Xi)
SDEALLOCATE(FV_NormVecEta)
SDEALLOCATE(FV_TangVec1Eta)
SDEALLOCATE(FV_TangVec2Eta)
SDEALLOCATE(FV_NormVecZeta)
SDEALLOCATE(FV_TangVec1Zeta)
SDEALLOCATE(FV_TangVec2Zeta)

#if PARABOLIC
SDEALLOCATE(FV_Metrics_fTilde_sJ)
SDEALLOCATE(FV_Metrics_gTilde_sJ)
SDEALLOCATE(FV_Metrics_hTilde_sJ)
#endif

SDEALLOCATE(FV_Elems_master) 
END SUBROUTINE FinalizeFV_Metrics


END MODULE MOD_FV_Metrics
#endif /* FV_ENABLED */

