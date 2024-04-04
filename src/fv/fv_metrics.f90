!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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

INTERFACE FV_CalcMetrics
  MODULE PROCEDURE FV_CalcMetrics
END INTERFACE

INTERFACE FinalizeFV_Metrics
  MODULE PROCEDURE FinalizeFV_Metrics
END INTERFACE

PUBLIC::InitFV_Metrics
PUBLIC::FV_CalcMetrics
PUBLIC::FinalizeFV_Metrics
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute the remaining metric terms for FV subcells, that are not computed in metrics.f90.
!> Normal, tangential vectors, SurfElems, ... for FV subcells.
!==================================================================================================================================
SUBROUTINE InitFV_Metrics()
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars
USE MOD_Mesh_Vars          ,ONLY: nElems,nSides,lastMPISide_MINE
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
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
! |<->.   |   .   |   .   |   .<->|    FV_dx_master/slave: distance from DG interface to first center of FV subcell
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

ALLOCATE(FV_sdx_XI    (0:PP_N,0:PP_NZ,1:PP_N,nElems)) ! 1. / FV_dx_XI    Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_sdx_ETA   (0:PP_N,0:PP_NZ,1:PP_N,nElems)) ! 1. / FV_dx_ETA   Attention: storage order is (i,k,j,iElem)
#if (PP_dim == 3)
ALLOCATE(FV_sdx_ZETA  (0:PP_N,0:PP_N ,1:PP_N,nElems)) ! 1. / FV_dx_ZETA  Attention: storage order is (i,j,k,iElem)
#endif

ALLOCATE(FV_sdx_Face (0:PP_N,0:PP_NZ,1:3,1:lastMPISide_MINE)) ! 1. / FV_dx_Face

ALLOCATE(FV_dx_XI_L  (0:PP_N,0:PP_NZ,0:PP_N,nElems))  ! Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_dx_XI_R  (0:PP_N,0:PP_NZ,0:PP_N,nElems))  ! Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_dx_ETA_L (0:PP_N,0:PP_NZ,0:PP_N,nElems))  ! Attention: storage order is (i,k,j,iElem)
ALLOCATE(FV_dx_ETA_R (0:PP_N,0:PP_NZ,0:PP_N,nElems))  ! Attention: storage order is (i,k,j,iElem)
#if (PP_dim == 3)
ALLOCATE(FV_dx_ZETA_L(0:PP_N,0:PP_N ,0:PP_N,nElems))  ! Attention: storage order is (i,j,k,iElem)
ALLOCATE(FV_dx_ZETA_R(0:PP_N,0:PP_N ,0:PP_N,nElems))  ! Attention: storage order is (i,j,k,iElem)
#endif

ALLOCATE(FV_dx_slave (1,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(FV_dx_master(1,0:PP_N,0:PP_NZ,1:nSides))
#endif /* FV_RECONSTRUCT */
ALLOCATE(FV_SurfElemXi_sw  (0:PP_N,0:PP_NZ,1:PP_N,nElems)) ! Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_SurfElemEta_sw (0:PP_N,0:PP_NZ,1:PP_N,nElems)) ! Attention: storage order is (i,k,j,iElem)
#if (PP_dim == 3)
ALLOCATE(FV_SurfElemZeta_sw(0:PP_N,0:PP_N ,1:PP_N,nElems)) ! Attention: storage order is (i,j,k,iElem)
#endif

ALLOCATE(FV_NormVecXi   (3,0:PP_N,0:PP_NZ,1:PP_N,nElems))  ! Attention: storage order is (j,k,i,iElem)
ALLOCATE(FV_TangVec1Xi  (3,0:PP_N,0:PP_NZ,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_TangVec2Xi  (3,0:PP_N,0:PP_NZ,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_NormVecEta  (3,0:PP_N,0:PP_NZ,1:PP_N,nElems))  ! Attention: storage order is (i,k,j,iElem)
ALLOCATE(FV_TangVec1Eta (3,0:PP_N,0:PP_NZ,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_TangVec2Eta (3,0:PP_N,0:PP_NZ,1:PP_N,nElems))  !  -"-
#if (PP_dim == 3)
ALLOCATE(FV_NormVecZeta (3,0:PP_N,0:PP_N ,1:PP_N,nElems))  ! Attention: storage order is (i,j,k,iElem)
ALLOCATE(FV_TangVec1Zeta(3,0:PP_N,0:PP_N ,1:PP_N,nElems))  !  -"-
ALLOCATE(FV_TangVec2Zeta(3,0:PP_N,0:PP_N ,1:PP_N,nElems))  !  -"-
#endif

#if FV_RECONSTRUCT
#if VOLINT_VISC
ALLOCATE(FV_Metrics_fTilde_sJ(3,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(FV_Metrics_gTilde_sJ(3,0:PP_N,0:PP_N,0:PP_NZ,nElems))
#if (PP_dim == 3)
ALLOCATE(FV_Metrics_hTilde_sJ(3,0:PP_N,0:PP_N,0:PP_NZ,nElems))
#endif

! Metrics for rotation from slave to master coordinate system
ALLOCATE(FV_Metrics_NormVec_master (3,0:PP_N,0:PP_NZ,nSides))
ALLOCATE(FV_Metrics_TangVec1_master(3,0:PP_N,0:PP_NZ,nSides))
ALLOCATE(FV_Metrics_NormVec_slave (3,0:PP_N,0:PP_NZ,nSides))
ALLOCATE(FV_Metrics_TangVec1_slave(3,0:PP_N,0:PP_NZ,nSides))
FV_Metrics_NormVec_master =0.
FV_Metrics_TangVec1_master=0.
FV_Metrics_NormVec_slave  =0.
FV_Metrics_TangVec1_slave =0.
#if PP_dim==3
ALLOCATE(FV_Metrics_TangVec2_master(3,0:PP_N,0:PP_NZ,nSides))
ALLOCATE(FV_Metrics_TangVec2_slave(3,0:PP_N,0:PP_NZ,nSides))
FV_Metrics_TangVec2_slave=0.
FV_Metrics_TangVec2_master=0.
#endif /* PP_dim==3 */

ALLOCATE(FV_Metrics_fTilde_sJ_xi  (3,0:PP_N-1,0:PP_N  ,0:PP_NZ  ,nElems))
ALLOCATE(FV_Metrics_gTilde_sJ_xi  (3,0:PP_N-1,0:PP_N  ,0:PP_NZ  ,nElems))
#if (PP_dim == 3)
ALLOCATE(FV_Metrics_hTilde_sJ_xi  (3,0:PP_N-1,0:PP_N  ,0:PP_NZ  ,nElems))
#endif
ALLOCATE(FV_Metrics_fTilde_sJ_eta (3,0:PP_N  ,0:PP_N-1,0:PP_NZ  ,nElems))
ALLOCATE(FV_Metrics_gTilde_sJ_eta (3,0:PP_N  ,0:PP_N-1,0:PP_NZ  ,nElems))
#if (PP_dim == 3)
ALLOCATE(FV_Metrics_hTilde_sJ_eta (3,0:PP_N  ,0:PP_N-1,0:PP_NZ  ,nElems))
#endif
ALLOCATE(FV_Metrics_fTilde_sJ_zeta(3,0:PP_N  ,0:PP_N  ,0:PP_NZ-1,nElems))
ALLOCATE(FV_Metrics_gTilde_sJ_zeta(3,0:PP_N  ,0:PP_N  ,0:PP_NZ-1,nElems))
#if (PP_dim == 3)
ALLOCATE(FV_Metrics_hTilde_sJ_zeta(3,0:PP_N  ,0:PP_N  ,0:PP_NZ-1,nElems))
#endif
#endif /* VOLINT_VISC */
#endif /* FV_RECONSTRUCT */

ALLOCATE(FV_Elems_master(1:nSides)) ! Moved from InitFV to here, since needed in U_Mortar below.

END SUBROUTINE InitFV_Metrics

!==================================================================================================================================
!> Compute the remaining metric terms for FV subcells, that are not computed in metrics.f90.
!> Normal, tangential vectors, SurfElems, ... for FV subcells.
!==================================================================================================================================
SUBROUTINE FV_CalcMetrics()
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
USE MOD_Mesh_Vars          ,ONLY: sJ_master,sJ_slave
USE MOD_ChangeBasis        ,ONLY: ChangeBasis1D,ChangeBasis2D,ChangeBasis3D
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf,ChangeBasisVolume
USE MOD_Metrics            ,ONLY: SurfMetricsFromJa
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
#if FV_RECONSTRUCT
USE MOD_Mesh_Vars          ,ONLY: firstBCSide,lastBCSide,firstInnerSide,lastMPISide_MINE
USE MOD_Mesh_Vars          ,ONLY: S2V2,ElemToSide,dXCL_N
USE MOD_Interpolation      ,ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars ,ONLY: NodeTypeCL,xGP,wGP,wBary
#if VOLINT_VISC
USE MOD_Mesh_Vars          ,ONLY: SideToElem,S2V
#endif /* VOLINT_VISC */
#endif
#if USE_MPI
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_MPI_Vars           ,ONLY: nNbProcs
USE MOD_Mesh_Vars          ,ONLY: firstMPISide_MINE
#endif
USE MOD_2D
USE MOD_FillMortar1        ,ONLY: U_Mortar1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: i,j,k,l,iSide,iElem,iLocSide
INTEGER                                :: dd,NormalDir,TangDir
REAL                                   :: NormalSign
REAL                                   :: Vdm_Gauss_FVboundary(0:PP_N+1,0:PP_N)
REAL                                   :: JaVol(3,3,0:PP_N+1,0:PP_N+1,0:PP_NZ+1)
REAL                                   :: FV_Ja_Face(3,3,0:PP_N,0:PP_NZ)
REAL,DIMENSION(0:PP_N,0:NgeoRef)       :: Vdm_NgeoRef_N,Vdm_NgeoRef_FV
REAL                                   :: FV_DetJac(1,0:PP_N,0:PP_N,0:PP_NZ)
INTEGER                                :: flip, SideID, iMortar
#if FV_RECONSTRUCT
INTEGER                                :: ijk(2),p,q,locSideID
REAL                                   :: FV_dx_Face    (0:PP_N,0:PP_NZ,1:3)
REAL                                   :: DG_dx_slave (1,0:PP_N,0:PP_NZ,1:nSides)
REAL                                   :: DG_dx_master(1,0:PP_N,0:PP_NZ,1:nSides)
REAL                                   :: tmp2(3,0:PP_N)
REAL,DIMENSION(0:PP_N,0:PP_N)          :: Vdm_CLN_FV, Vdm_CLN_GaussN,length
REAL,DIMENSION(0:PP_N+1,0:PP_N)        :: Vdm_CLN_FVboundary
REAL,DIMENSION(3,0:PP_N,0:PP_N,0:PP_NZ):: FV_Path_XI, FV_Path_ETA, FV_Path_ZETA
REAL                                   :: x0, xN
REAL,POINTER                           :: FV_dx_P(:,:)
#if USE_MPI
INTEGER                                :: MPIRequest(nNbProcs,2)
#endif
#if VOLINT_VISC
INTEGER                                :: d,ijk2(3),ElemID,NB_ElemID
#endif
#endif
#if USE_MPI
INTEGER                                :: MPIRequest_Geo(nNbProcs,2)
REAL,ALLOCATABLE                       :: Geo(:,:,:,:)
#endif
#if VOLINT_VISC
REAL                                   :: dXFace  (3,3,0:PP_N  ,0:PP_NZ           )
REAL                                   :: dXVol   (3,3,0:PP_N+1,0:PP_N+1,0:PP_NZ+1)
REAL                                   :: Metrics_fTilde_xi(3,0:PP_N-1,0:PP_N,0:PP_NZ)
REAL                                   :: Metrics_gTilde_xi(3,0:PP_N-1,0:PP_N,0:PP_NZ)
REAL                                   :: Metrics_hTilde_xi(3,0:PP_N-1,0:PP_N,0:PP_NZ)
REAL                                   :: Metrics_fTilde_eta(3,0:PP_N,0:PP_N-1,0:PP_NZ)
REAL                                   :: Metrics_gTilde_eta(3,0:PP_N,0:PP_N-1,0:PP_NZ)
REAL                                   :: Metrics_hTilde_eta(3,0:PP_N,0:PP_N-1,0:PP_NZ)
REAL                                   :: Metrics_fTilde_zeta(3,0:PP_N,0:PP_N,0:PP_NZ-1)
REAL                                   :: Metrics_gTilde_zeta(3,0:PP_N,0:PP_N,0:PP_NZ-1)
REAL                                   :: Metrics_hTilde_zeta(3,0:PP_N,0:PP_N,0:PP_NZ-1)
REAL                                   :: dX_xi   (3,3,0:PP_N-1,0:PP_N  ,0:PP_NZ  )
REAL                                   :: dX_eta  (3,3,0:PP_N  ,0:PP_N-1,0:PP_NZ  )
#if PP_dim==3
REAL                                   :: dX_zeta (3,3,0:PP_N  ,0:PP_N  ,0:PP_NZ-1)
#endif
REAL                                   :: sJ_xi  (0:PP_N-1,0:PP_N  ,0:PP_NZ  )
REAL                                   :: sJ_eta (0:PP_N  ,0:PP_N-1,0:PP_NZ  )
#if PP_dim==3
REAL                                   :: sJ_zeta(0:PP_N  ,0:PP_N  ,0:PP_NZ-1)
#endif
REAL                                   :: FV_dX_1_xi_LG_eta_FV_zeta_BC_L(0:PP_N  ,0:PP_NZ+1,0:PP_N)  ! Attention: storage order is (j,k,i,iElem)
REAL                                   :: FV_dX_1_xi_LG_eta_BC_zeta_FV_L(0:PP_N+1,0:PP_NZ  ,0:PP_N)  ! Attention: storage order is (j,k,i,iElem)
#if (PP_dim == 3)
REAL                                   :: FV_dX_3_xi_FV_eta_BC_zeta_LG_L(0:PP_N  ,0:PP_NZ+1,0:PP_N)  ! Attention: storage order is (i,j,k,iElem)
#endif
REAL                                   :: FV_dX_2_xi_FV_eta_LG_zeta_BC_L(0:PP_N  ,0:PP_NZ+1,0:PP_N)  ! Attention: storage order is (i,k,j,iElem)
REAL                                   :: FV_dX_2_xi_BC_eta_LG_zeta_FV_L(0:PP_N+1,0:PP_NZ  ,0:PP_N)  ! Attention: storage order is (i,k,j,iElem)
#if (PP_dim == 3)
REAL                                   :: FV_dX_3_xi_BC_eta_FV_zeta_LG_L(0:PP_N+1,0:PP_NZ  ,0:PP_N)  ! Attention: storage order is (i,j,k,iElem)
#endif
REAL                                   :: FV_dX_1_xi_LG_eta_FV_zeta_BC_R(0:PP_N  ,0:PP_NZ+1,0:PP_N)  ! Attention: storage order is (j,k,i,iElem)
REAL                                   :: FV_dX_1_xi_LG_eta_BC_zeta_FV_R(0:PP_N+1,0:PP_NZ  ,0:PP_N)  ! Attention: storage order is (j,k,i,iElem)
#if (PP_dim == 3)
REAL                                   :: FV_dX_3_xi_FV_eta_BC_zeta_LG_R(0:PP_N  ,0:PP_NZ+1,0:PP_N)  ! Attention: storage order is (i,j,k,iElem)
#endif
REAL                                   :: FV_dX_2_xi_FV_eta_LG_zeta_BC_R(0:PP_N  ,0:PP_NZ+1,0:PP_N)  ! Attention: storage order is (i,k,j,iElem)
REAL                                   :: FV_dX_2_xi_BC_eta_LG_zeta_FV_R(0:PP_N+1,0:PP_NZ  ,0:PP_N)  ! Attention: storage order is (i,k,j,iElem)
#if (PP_dim == 3)
REAL                                   :: FV_dX_3_xi_BC_eta_FV_zeta_LG_R(0:PP_N+1,0:PP_NZ  ,0:PP_N)  ! Attention: storage order is (i,j,k,iElem)
#endif
REAL                                   :: dX_1_xi_LG_eta_CL_zeta_CL_Path       (3,0:PP_N  ,0:PP_N   ,0:PP_NZ  )
REAL                                   :: dX_2_xi_CL_eta_LG_zeta_CL_Path       (3,0:PP_N  ,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_1_xi_LG_eta_FV_zeta_CL_Path    (3,0:PP_N  ,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_1_xi_LG_eta_BC_zeta_CL_Path    (3,0:PP_N  ,0:PP_N+1 ,0:PP_NZ  )
REAL                                   :: FV_dX_2_xi_FV_eta_LG_zeta_CL_Path    (3,0:PP_N  ,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_2_xi_BC_eta_LG_zeta_CL_Path    (3,0:PP_N+1,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_1_xi_LG_eta_FV_zeta_BC_Path    (3,0:PP_N  ,0:PP_N   ,0:PP_NZ+1)
REAL                                   :: FV_dX_1_xi_LG_eta_FV_zeta_BC_Path_pq (3,0:PP_N  ,0:PP_N   ,0:PP_NZ+1)
REAL                                   :: FV_dX_1_xi_LG_eta_BC_zeta_FV_Path    (3,0:PP_N  ,0:PP_N+1 ,0:PP_NZ  )
REAL                                   :: FV_dX_1_xi_LG_eta_BC_zeta_FV_Path_pq (3,0:PP_N  ,0:PP_N+1 ,0:PP_NZ  )
REAL                                   :: FV_dX_2_xi_FV_eta_LG_zeta_BC_Path    (3,0:PP_N  ,0:PP_N   ,0:PP_NZ+1)
REAL                                   :: FV_dX_2_xi_FV_eta_LG_zeta_BC_Path_pq (3,0:PP_N  ,0:PP_N   ,0:PP_NZ+1)
REAL                                   :: FV_dX_2_xi_BC_eta_LG_zeta_FV_Path    (3,0:PP_N+1,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_2_xi_BC_eta_LG_zeta_FV_Path_pq (3,0:PP_N  ,0:PP_N+1 ,0:PP_NZ  )
#if (PP_dim == 3)
REAL                                   :: dX_3_xi_CL_eta_CL_zeta_LG_Path       (3,0:PP_N  ,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_3_xi_FV_eta_CL_zeta_LG_Path    (3,0:PP_N  ,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_3_xi_BC_eta_CL_zeta_LG_Path    (3,0:PP_N+1,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_3_xi_FV_eta_BC_zeta_LG_Path    (3,0:PP_N  ,0:PP_N+1 ,0:PP_NZ  )
REAL                                   :: FV_dX_3_xi_FV_eta_BC_zeta_LG_Path_pq (3,0:PP_N  ,0:PP_N   ,0:PP_NZ+1)
REAL                                   :: FV_dX_3_xi_BC_eta_FV_zeta_LG_Path    (3,0:PP_N+1,0:PP_N   ,0:PP_NZ  )
REAL                                   :: FV_dX_3_xi_BC_eta_FV_zeta_LG_Path_pq (3,0:PP_N  ,0:PP_N+1 ,0:PP_NZ  )
#endif /*(PP_dim == 3)*/
#endif /*VOLINT_VISC*/
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  Build FV-Metrics ...'

! compute FV NormVec, TangVec,.. on boundary of DG-cells
DO iSide=1,nSides
  IF(iSide.GE.firstMPISide_YOUR.AND.iSide.LE.lastMPISide_YOUR) CYCLE
  CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,Face_xGP(:,:,:,0,iSide),Face_xGP(:,:,:,1,iSide))

  iLocSide=SideToElem(S2E_LOC_SIDE_ID,iSide)
  IF(iLocSide.LT.1) CYCLE
  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide)

  CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,Ja_Face( 1,:,:,:,iSide),FV_Ja_Face(1,:,:,:))
  CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,Ja_Face( 2,:,:,:,iSide),FV_Ja_Face(2,:,:,:))
  CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,Ja_Face( 3,:,:,:,iSide),FV_Ja_Face(3,:,:,:))
  CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
                         NormVec( :,:,:,1,iSide),TangVec1(:,:,:,1,iSide),&
                         TangVec2(:,:,:,1,iSide),SurfElem(  :,:,1,iSide))
  CALL ChangeBasisSurf(1,PP_N,PP_N,FV_Vdm,sJ_master(:,0:PP_N,0:PP_NZ,iSide,0),sJ_master(:,0:PP_N,0:PP_NZ,iSide,1))
  CALL ChangeBasisSurf(1,PP_N,PP_N,FV_Vdm,sJ_slave (:,0:PP_N,0:PP_NZ,iSide,0),sJ_slave (:,0:PP_N,0:PP_NZ,iSide,1))

  IF(MortarType(1,iSide).LE.0) CYCLE ! no mortars
  DO iMortar=1,MERGE(4,2,MortarType(1,iSide).EQ.1)
    SideID=MortarInfo(MI_SIDEID,iMortar,MortarType(2,iSide))
    Flip  =MortarInfo(MI_FLIP,iMortar,MortarType(2,iSide))
    IF(flip.NE.0) CYCLE ! for MPI sides some sides are built from the inside and for type 2/3 there are only 2 neighbours
    CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,Ja_Face( 1,:,:,:,SideID),FV_Ja_Face(1,:,:,:))
    CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,Ja_Face( 2,:,:,:,SideID),FV_Ja_Face(2,:,:,:))
    CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,Ja_Face( 3,:,:,:,SideID),FV_Ja_Face(3,:,:,:))
    CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
                           NormVec( :,:,:,1,SideID),TangVec1(:,:,:,1,SideID),&
                           TangVec2(:,:,:,1,SideID),SurfElem(  :,:,1,SideID))
    CALL ChangeBasisSurf(1,PP_N,PP_N,FV_Vdm,sJ_master(:,0:PP_N,0:PP_NZ,SideID,0),sJ_master(:,0:PP_N,0:PP_NZ,SideID,1))
    CALL ChangeBasisSurf(1,PP_N,PP_N,FV_Vdm,sJ_slave (:,0:PP_N,0:PP_NZ,SideID,0),sJ_slave (:,0:PP_N,0:PP_NZ,SideID,1))
  END DO
END DO

#if USE_MPI
! Send surface geomtry informations from mpi master to mpi slave
ALLOCATE(Geo(10,0:PP_N,0:PP_NZ,firstMPISide_MINE:nSides))
Geo=0.
Geo(1,:,:,:)   =SurfElem(  :,0:PP_NZ,1,firstMPISide_MINE:nSides)
Geo(2:4,:,:,:) =NormVec (:,:,0:PP_NZ,1,firstMPISide_MINE:nSides)
Geo(5:7,:,:,:) =TangVec1(:,:,0:PP_NZ,1,firstMPISide_MINE:nSides)
Geo(8:10,:,:,:)=TangVec2(:,:,0:PP_NZ,1,firstMPISide_MINE:nSides)
MPIRequest_Geo=MPI_REQUEST_NULL
CALL StartReceiveMPIData(Geo,10*(PP_N+1)**(PP_dim-1),firstMPISide_MINE,nSides,MPIRequest_Geo(:,RECV),SendID=1) ! Receive YOUR / Geo: master -> slave
CALL StartSendMPIData(   Geo,10*(PP_N+1)**(PP_dim-1),firstMPISide_MINE,nSides,MPIRequest_Geo(:,SEND),SendID=1) ! SEND MINE / Geo: master -> slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Geo)
SurfElem  (:,0:PP_NZ,1,firstMPISide_YOUR:lastMPISide_YOUR)= Geo(1   ,:,:,firstMPISide_YOUR:lastMPISide_YOUR)
NormVec (:,:,0:PP_NZ,1,firstMPISide_YOUR:lastMPISide_YOUR)= Geo(2:4 ,:,:,firstMPISide_YOUR:lastMPISide_YOUR)
TangVec1(:,:,0:PP_NZ,1,firstMPISide_YOUR:lastMPISide_YOUR)= Geo(5:7 ,:,:,firstMPISide_YOUR:lastMPISide_YOUR)
TangVec2(:,:,0:PP_NZ,1,firstMPISide_YOUR:lastMPISide_YOUR)= Geo(8:10,:,:,firstMPISide_YOUR:lastMPISide_YOUR)
DEALLOCATE(Geo)
#endif /*MPI*/

CALL FV_Build_Vdm_Gauss_FVboundary(PP_N, Vdm_Gauss_FVboundary)
CALL GetVandermonde(NgeoRef, NodeType, PP_N, NodeType, Vdm_NgeoRef_N, modal=.TRUE.)
Vdm_NgeoRef_FV = MATMUL(FV_Vdm, Vdm_NgeoRef_N)
CALL GetVandermonde(PP_N,NodeTypeCL,PP_N,NodeType, Vdm_CLN_GaussN, modal=.FALSE.)
Vdm_CLN_FV = MATMUL(FV_Vdm, Vdm_CLN_GaussN)
Vdm_CLN_FVboundary = MATMUL(Vdm_Gauss_FVboundary, Vdm_CLN_GaussN)
DO iElem=1,nElems
  ! compute Jacobian
  CALL ChangeBasisVolume(1,NGeoRef,PP_N,Vdm_NgeoRef_FV,DetJac_Ref(:,:,:,:,iElem),FV_DetJac)
  sJ(:,:,:,iElem,1) = 1./FV_DetJac(1,:,:,:)
  CALL ChangeBasisVolume(3,PP_N,PP_N,FV_Vdm,Metrics_fTilde(:,:,:,:,iElem,0),Metrics_fTilde(:,:,:,:,iElem,1))
  CALL ChangeBasisVolume(3,PP_N,PP_N,FV_Vdm,Metrics_gTilde(:,:,:,:,iElem,0),Metrics_gTilde(:,:,:,:,iElem,1))
  CALL ChangeBasisVolume(3,PP_N,PP_N,FV_Vdm,Metrics_hTilde(:,:,:,:,iElem,0),Metrics_hTilde(:,:,:,:,iElem,1))

  !================================================================
  ! Compute FV NormVec,TangVec,.. at inner cell boundaries...START
  !================================================================

  ! Xi direction
  DO k=0,PP_NZ; DO j=0,PP_N
    ! interpolate Metrics to boundaries of FV subcells in XI direction, other directions stay DG
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_fTilde(:,:,j,k,iElem,0),JaVol(1,:,:,j,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_gTilde(:,:,j,k,iElem,0),JaVol(2,:,:,j,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_hTilde(:,:,j,k,iElem,0),JaVol(3,:,:,j,k))
#if VOLINT_VISC
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(1,:,:,j,k,iElem),dXVol(1,:,:,j,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(2,:,:,j,k,iElem),dXVol(2,:,:,j,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(3,:,:,j,k,iElem),dXVol(3,:,:,j,k))
#endif /*VOLINT_VISC*/
  END DO; END DO ! j,k=0,PP_N
  DO l=1,PP_N
    ! at every inner interface/slice between FV subcells in XI direction:
    ! convert metrics in the other directions from DG to FV subcells
    DO dd=1,3
      CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,JaVol(dd,1:3,l,0:PP_N,0:PP_NZ),FV_Ja_Face(dd,:,:,:))
#if VOLINT_VISC
      CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,dXVol(dd,1:3,l,0:PP_N,0:PP_NZ),dXFace(dd,:,:,:))
#endif /*VOLINT_VISC*/
    END DO
    ! use metrics to build normal/tangential vectors and surelem at the inner interfaces/slices
    NormalDir=NormalDirs(XI_PLUS); TangDir=TangDirs(XI_PLUS); NormalSign=NormalSigns(XI_PLUS)
    CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
        FV_NormVecXi  (:,:,:,l,iElem),&
        FV_TangVec1Xi (:,:,:,l,iElem),&
        FV_TangVec2Xi (:,:,:,l,iElem),&
        FV_SurfElemXi_sw(:,:,l,iElem))
    ! multiply FV_SurfElemXi with 1/FV_w
    !FV_SurfElemXi_sw(:,:,l,iElem) = FV_SurfElemXi_sw(:,:,l,iElem)*FV_w_inv
#if VOLINT_VISC
    DO dd=1,3
      DO k=0,PP_NZ; DO j=0,PP_N
        Metrics_fTilde_xi(dd,l-1,j,k)=FV_Ja_Face(1,dd,j,k)
        Metrics_gTilde_xi(dd,l-1,j,k)=FV_Ja_Face(2,dd,j,k)
        Metrics_hTilde_xi(dd,l-1,j,k)=FV_Ja_Face(3,dd,j,k)
        dX_xi          (1,dd,l-1,j,k)=dXFace    (1,dd,j,k)
        dX_xi          (2,dd,l-1,j,k)=dXFace    (2,dd,j,k)
        dX_xi          (3,dd,l-1,j,k)=dXFace    (3,dd,j,k)
      END DO; END DO ! j,k=0,PP_N
    END DO
#endif /*VOLINT_VISC*/
  END DO

  ! Eta direction
  DO k=0,PP_NZ; DO i=0,PP_N
    ! interpolate Metrics to boundaries of FV subcells in ETA direction, other directions stay DG
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_fTilde(:,i,:,k,iElem,0),JaVol(1,:,i,:,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_gTilde(:,i,:,k,iElem,0),JaVol(2,:,i,:,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_hTilde(:,i,:,k,iElem,0),JaVol(3,:,i,:,k))
#if VOLINT_VISC
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(1,:,i,:,k,iElem),dXVol(1,:,i,:,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(2,:,i,:,k,iElem),dXVol(2,:,i,:,k))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(3,:,i,:,k,iElem),dXVol(3,:,i,:,k))
#endif /*VOLINT_VISC*/
  END DO; END DO ! i,k=0,PP_N
  DO l=1,PP_N
    ! at every inner interface/slice between FV subcells in ETA direction:
    ! convert metrics in the other directions from DG to FV subcells
    DO dd=1,3
      CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,JaVol(dd,1:3,0:PP_N,l,0:PP_NZ),FV_Ja_Face(dd,:,:,:))
#if VOLINT_VISC
      CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,dXVol(dd,1:3,0:PP_N,l,0:PP_NZ),dXFace(dd,:,:,:))
#endif /*VOLINT_VISC*/
    END DO
    ! use metrics to build normal/tangential vectors and surelem at the inner interfaces/slices
    NormalDir=NormalDirs(ETA_PLUS); TangDir=TangDirs(ETA_PLUS); NormalSign=NormalSigns(ETA_PLUS)
    CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
        FV_NormVecEta  (:,:,:,l,iElem),&
        FV_TangVec1Eta (:,:,:,l,iElem),&
        FV_TangVec2Eta (:,:,:,l,iElem),&
        FV_SurfElemEta_sw(:,:,l,iElem))
    ! multiply FV_SurfElemEta with 1/FV_w
    !FV_SurfElemEta_sw(:,:,l,iElem) = FV_SurfElemEta_sw(:,:,l,iElem)*FV_w_inv
#if VOLINT_VISC
    DO dd=1,3
      DO k=0,PP_NZ; DO i=0,PP_N
        Metrics_fTilde_eta(dd,i,l-1,k)=FV_Ja_Face(1,dd,i,k)
        Metrics_gTilde_eta(dd,i,l-1,k)=FV_Ja_Face(2,dd,i,k)
        Metrics_hTilde_eta(dd,i,l-1,k)=FV_Ja_Face(3,dd,i,k)
        dX_eta          (1,dd,i,l-1,k)=dXFace    (1,dd,i,k)
        dX_eta          (2,dd,i,l-1,k)=dXFace    (2,dd,i,k)
        dX_eta          (3,dd,i,l-1,k)=dXFace    (3,dd,i,k)
      END DO; END DO ! i,k=0,PP_N
    END DO
#endif /*VOLINT_VISC*/
  END DO

#if (PP_dim == 3)
  ! Zeta direction
  DO j=0,PP_N; DO i=0,PP_N
    ! interpolate Metrics to boundaries of FV subcells in ZETA direction, other directions stay DG
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_fTilde(:,i,j,:,iElem,0),JaVol(1,:,i,j,:))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_gTilde(:,i,j,:,iElem,0),JaVol(2,:,i,j,:))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_Gauss_FVboundary,Metrics_hTilde(:,i,j,:,iElem,0),JaVol(3,:,i,j,:))
#if VOLINT_VISC
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(1,:,i,j,:,iElem),dXVol(1,:,i,j,:))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(2,:,i,j,:,iElem),dXVol(2,:,i,j,:))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary,dXCL_N(3,:,i,j,:,iElem),dXVol(3,:,i,j,:))
#endif /*VOLINT_VISC*/
  END DO; END DO ! i,j=0,PP_N
  DO l=1,PP_N
    ! at every inner interface/slice between FV subcells in ZETA direction:
    ! convert metrics in the other directions from DG to FV subcells
    DO dd=1,3
      CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,JaVol(dd,1:3,0:PP_N,0:PP_N,l),FV_Ja_Face(dd,:,:,:))
#if VOLINT_VISC
      CALL ChangeBasisSurf(3,PP_N,PP_N,FV_Vdm,dXVol(dd,1:3,0:PP_N,0:PP_N,l),dXFace(dd,:,:,:))
#endif /*VOLINT_VISC*/
    END DO
    ! use metrics to build normal/tangential vectors and surelem at the inner interfaces/slices
    NormalDir=NormalDirs(ZETA_PLUS); TangDir=TangDirs(ZETA_PLUS); NormalSign=NormalSigns(ZETA_PLUS)
    CALL SurfMetricsFromJa(PP_N,NormalDir,TangDir,NormalSign,FV_Ja_Face,&
        FV_NormVecZeta  (:,:,:,l,iElem),&
        FV_TangVec1Zeta (:,:,:,l,iElem),&
        FV_TangVec2Zeta (:,:,:,l,iElem),&
        FV_SurfElemZeta_sw(:,:,l,iElem))
    ! multiply FV_SurfElemZeta with 1/FV_w
    !FV_SurfElemZeta_sw(:,:,l,iElem) = FV_SurfElemZeta_sw(:,:,l,iElem)*FV_w_inv
#if VOLINT_VISC
    DO dd=1,3
      DO j=0,PP_N; DO i=0,PP_N
        Metrics_fTilde_zeta(dd,i,j,l-1)=FV_Ja_Face(1,dd,i,j)
        Metrics_gTilde_zeta(dd,i,j,l-1)=FV_Ja_Face(2,dd,i,j)
        Metrics_hTilde_zeta(dd,i,j,l-1)=FV_Ja_Face(3,dd,i,j)
        dX_zeta          (1,dd,i,j,l-1)=dXFace    (1,dd,i,j)
        dX_zeta          (2,dd,i,j,l-1)=dXFace    (2,dd,i,j)
        dX_zeta          (3,dd,i,j,l-1)=dXFace    (3,dd,i,j)
      END DO; END DO ! i,k=0,PP_N
    END DO
#endif /*VOLINT_VISC*/
  END DO
#endif /* PP_dim == 3 */

#if VOLINT_VISC
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N-1
#if (PP_dim == 3)
    sJ_xi(i,j,k)= 1./ ( &
      + dX_xi(1,1,i,j,k)*(dX_xi(2,2,i,j,k)*dX_xi(3,3,i,j,k) - dX_xi(3,2,i,j,k)*dX_xi(2,3,i,j,k))  &
      + dX_xi(2,1,i,j,k)*(dX_xi(3,2,i,j,k)*dX_xi(1,3,i,j,k) - dX_xi(1,2,i,j,k)*dX_xi(3,3,i,j,k))  &
      + dX_xi(3,1,i,j,k)*(dX_xi(1,2,i,j,k)*dX_xi(2,3,i,j,k) - dX_xi(2,2,i,j,k)*dX_xi(1,3,i,j,k))  )
#else
    sJ_xi(i,j,k)= 1./(dX_xi(1,1,i,j,k)*dX_xi(2,2,i,j,k) - dX_xi(2,1,i,j,k)*dX_xi(1,2,i,j,k))
#endif
  END DO; END DO; END DO

  DO k=0,PP_NZ; DO j=0,PP_N-1; DO i=0,PP_N
#if (PP_dim == 3)
    sJ_eta(i,j,k)= 1./ ( &
      + dX_eta(1,1,i,j,k)*(dX_eta(2,2,i,j,k)*dX_eta(3,3,i,j,k) - dX_eta(3,2,i,j,k)*dX_eta(2,3,i,j,k))  &
      + dX_eta(2,1,i,j,k)*(dX_eta(3,2,i,j,k)*dX_eta(1,3,i,j,k) - dX_eta(1,2,i,j,k)*dX_eta(3,3,i,j,k))  &
      + dX_eta(3,1,i,j,k)*(dX_eta(1,2,i,j,k)*dX_eta(2,3,i,j,k) - dX_eta(2,2,i,j,k)*dX_eta(1,3,i,j,k))  )
#else
    sJ_eta(i,j,k)= 1./(dX_eta(1,1,i,j,k)*dX_eta(2,2,i,j,k) - dX_eta(2,1,i,j,k)*dX_eta(1,2,i,j,k))
#endif

  END DO; END DO; END DO

#if (PP_dim == 3)
  DO k=0,PP_NZ-1; DO j=0,PP_N; DO i=0,PP_N
    sJ_zeta(i,j,k)= 1./ ( &
      + dX_zeta(1,1,i,j,k)*(dX_zeta(2,2,i,j,k)*dX_zeta(3,3,i,j,k) - dX_zeta(3,2,i,j,k)*dX_zeta(2,3,i,j,k))  &
      + dX_zeta(2,1,i,j,k)*(dX_zeta(3,2,i,j,k)*dX_zeta(1,3,i,j,k) - dX_zeta(1,2,i,j,k)*dX_zeta(3,3,i,j,k))  &
      + dX_zeta(3,1,i,j,k)*(dX_zeta(1,2,i,j,k)*dX_zeta(2,3,i,j,k) - dX_zeta(2,2,i,j,k)*dX_zeta(1,3,i,j,k))  )
  END DO; END DO; END DO
#endif
#endif /*VOLINT_VISC*/
  !==================================================================
  ! Compute FV NormVec,TangVec,.. at inner cell boundaries...DONE
  !==================================================================

#if FV_RECONSTRUCT
  !=====================================================================================
  ! Compute distances between FV subcells and between first Gauss point and interface
  ! Integrate path given by dXCL_N, to ensure free stream preservation for N<Ngeo
  !=====================================================================================
  DO l=0,PP_N
    CALL ChangeBasisSurf(3,PP_N,PP_N,Vdm_CLN_FV, dXCL_N(1,:,l,:,:,iElem), FV_Path_XI  (:,l,:,:))
    CALL ChangeBasisSurf(3,PP_N,PP_N,Vdm_CLN_FV, dXCL_N(2,:,:,l,:,iElem), FV_Path_ETA (:,l,:,:))
#if (PP_dim == 3)
    CALL ChangeBasisSurf(3,PP_N,PP_N,Vdm_CLN_FV, dXCL_N(3,:,:,:,l,iElem), FV_Path_ZETA(:,l,:,:))
#endif
  END DO ! i=0,PP_N
  DO q=0,PP_NZ; DO p=0,PP_N
    tmp2 = FV_Path_XI(:,:,p,q)
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN, tmp2, FV_Path_XI(:,:,p,q))
    tmp2 = FV_Path_ETA(:,:,p,q)
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN, tmp2, FV_Path_ETA(:,:,p,q))
#if (PP_dim == 3)
    tmp2 = FV_Path_ZETA(:,:,p,q)
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN, tmp2, FV_Path_ZETA(:,:,p,q))
#endif
  END DO; END DO! p,q=0,PP_N

  ! Calculate distances between FV subcells
  DO l=0,PP_N
    ! left
    x0 = FV_BdryX(l)
    xN = x0 + (FV_BdryX(l+1) - FV_BdryX(l)) * 0.5
    CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, x0,xN, FV_Path_XI  , FV_dx_XI_L  (:,:,l,iElem))
    CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, x0,xN, FV_Path_ETA , FV_dx_ETA_L (:,:,l,iElem))
#if (PP_dim == 3)
    CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, x0,xN, FV_Path_ZETA, FV_dx_ZETA_L(:,:,l,iElem))
#endif

    ! right
    xN = FV_BdryX(l+1)
    x0 = xN - (FV_BdryX(l+1) - FV_BdryX(l)) * 0.5
    CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, x0,xN, FV_Path_XI  , FV_dx_XI_R  (:,:,l,iElem))
    CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, x0,xN, FV_Path_ETA , FV_dx_ETA_R (:,:,l,iElem))
#if (PP_dim == 3)
    CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, x0,xN, FV_Path_ZETA, FV_dx_ZETA_R(:,:,l,iElem))
#endif
  END DO ! l=0,PP_N

  ! build inverse distances (volumes)
  FV_sdx_XI  (:,:,:,iElem) = 1. / (FV_dx_XI_R  (:,:,0:PP_N-1,iElem) +  FV_dx_XI_L  (:,:,1:PP_N,iElem)) ! 1 / FV_dx_XI
  FV_sdx_ETA (:,:,:,iElem) = 1. / (FV_dx_ETA_R (:,:,0:PP_N-1,iElem) +  FV_dx_ETA_L (:,:,1:PP_N,iElem)) ! 1 / FV_dx_ETA
#if PP_dim == 3
  FV_sdx_ZETA(:,:,:,iElem) = 1. / (FV_dx_ZETA_R(:,:,0:PP_N-1,iElem) +  FV_dx_ZETA_L(:,:,1:PP_N,iElem)) ! 1 / FV_dx_ZETA
#endif

#if VOLINT_VISC
  DO d=1,3
    DO i=0,PP_N
      FV_Metrics_fTilde_sJ(d,i,:,:,iElem)=FV_w_inv(i)*Metrics_fTilde(d,i,:,:,iElem,1)*&
          (FV_dx_XI_L  (:,:,i,iElem)+FV_dx_XI_R  (:,:,i,iElem))
    END DO ! i=0,PP_N
    DO j=0,PP_N
      FV_Metrics_gTilde_sJ(d,:,j,:,iElem)=FV_w_inv(j)*Metrics_gTilde(d,:,j,:,iElem,1)*&
          (FV_dx_ETA_L (:,:,j,iElem)+FV_dx_ETA_R (:,:,j,iElem))
    END DO ! j=0,PP_N
    FV_Metrics_fTilde_sJ(d,:,:,:,iElem)=FV_Metrics_fTilde_sJ(d,:,:,:,iElem)*sJ(:,:,:,iElem,1)
    FV_Metrics_gTilde_sJ(d,:,:,:,iElem)=FV_Metrics_gTilde_sJ(d,:,:,:,iElem)*sJ(:,:,:,iElem,1)
#if (PP_dim == 3)
    DO k=0,PP_N
      FV_Metrics_hTilde_sJ(d,:,:,k,iElem)=FV_w_inv(k)*Metrics_hTilde(d,:,:,k,iElem,1)*&
          (FV_dx_ZETA_L(:,:,k,iElem)+FV_dx_ZETA_R(:,:,k,iElem))
    END DO ! k=0,PP_N
    FV_Metrics_hTilde_sJ(d,:,:,:,iElem)=FV_Metrics_hTilde_sJ(d,:,:,:,iElem)*sJ(:,:,:,iElem,1)
#endif
  END DO
#endif /* VOLINT_VISC */

  ! Calculate distance between first GaussPoint and interface
#if PP_dim == 3
  DO locSideID=1,6
#else
  DO locSideID=2,5
#endif
    length=0.
    SideID = ElemToSide(E2S_SIDE_ID,locSideID,iElem)
    flip   = ElemToSide(E2S_FLIP,   locSideID,iElem)
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, -1., xGP(0),   FV_Path_XI  , length)
      FV_dx_P(0:PP_N,0:PP_NZ) => FV_dx_XI_L(:,:,0,iElem)
    CASE(ETA_MINUS)
      CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, -1., xGP(0),   FV_Path_ETA , length)
      FV_dx_P(0:PP_N,0:PP_NZ) => FV_dx_ETA_L(:,:,0,iElem)
    CASE(ZETA_MINUS)
      CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, -1., xGP(0),   FV_Path_ZETA, length)
      FV_dx_P(0:PP_N,0:PP_NZ) => FV_dx_ZETA_L(:,:,0,iElem)
    CASE(XI_PLUS)
      CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, xGP(PP_N), 1., FV_Path_XI  , length)
      FV_dx_P(0:PP_N,0:PP_NZ) => FV_dx_XI_R(:,:,PP_N,iElem)
    CASE(ETA_PLUS)
      CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, xGP(PP_N), 1., FV_Path_ETA , length)
      FV_dx_P(0:PP_N,0:PP_NZ) => FV_dx_ETA_R(:,:,PP_N,iElem)
    CASE(ZETA_PLUS)
      CALL Integrate_Path(PP_N,PP_N, xGP,wGP,wBary, xGP(PP_N), 1., FV_Path_ZETA, length)
      FV_dx_P(0:PP_N,0:PP_NZ) => FV_dx_ZETA_R(:,:,PP_N,iElem)
    CASE DEFAULT
      STOP 'Local side index out of range (1...6).'
    END SELECT

    IF (flip.EQ.0) THEN ! master side
      DO q=0,PP_NZ; DO p=0,PP_N
        ijk = S2V2(:,p,q,flip,locSideID)
        DG_dx_master(1,p,q,SideID) = length(ijk(1),ijk(2))
        FV_dx_master(1,p,q,SideID) = FV_dx_P(ijk(1),ijk(2))
      END DO; END DO
    ELSE ! slave side
      DO q=0,PP_NZ; DO p=0,PP_N
        ijk = S2V2(:,p,q,flip,locSideID)
        DG_dx_slave(1,p,q,SideID) = length(ijk(1),ijk(2))
        FV_dx_slave(1,p,q,SideID) = FV_dx_P(ijk(1),ijk(2))
      END DO; END DO
    END IF
  END DO

#if VOLINT_VISC
  DO q=0,PP_NZ; DO p=0,PP_N
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN,dXCL_N(1,:,:,p,q,iElem), dX_1_xi_LG_eta_CL_zeta_CL_Path(:,:,p,q))
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN,dXCL_N(2,:,p,:,q,iElem), dX_2_xi_CL_eta_LG_zeta_CL_Path(:,p,:,q))
#if (PP_dim == 3)
    CALL ChangeBasis1D(3,PP_N,PP_N,Vdm_CLN_GaussN,dXCL_N(3,:,p,q,:,iElem), dX_3_xi_CL_eta_CL_zeta_LG_Path(:,p,q,:))
#endif
  END DO; END DO

  DO q=0,PP_NZ; DO p=0,PP_N
    CALL ChangeBasis1D(3,PP_N,PP_N  ,Vdm_CLN_FV        , dX_1_xi_LG_eta_CL_zeta_CL_Path(:,p,:,q), FV_dX_1_xi_LG_eta_FV_zeta_CL_Path(:,p,:,q))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary, dX_1_xi_LG_eta_CL_zeta_CL_Path(:,p,:,q), FV_dX_1_xi_LG_eta_BC_zeta_CL_Path(:,p,:,q))

    CALL ChangeBasis1D(3,PP_N,PP_N  ,Vdm_CLN_FV        , dX_2_xi_CL_eta_LG_zeta_CL_Path(:,:,p,q), FV_dX_2_xi_FV_eta_LG_zeta_CL_Path(:,:,p,q))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary, dX_2_xi_CL_eta_LG_zeta_CL_Path(:,:,p,q), FV_dX_2_xi_BC_eta_LG_zeta_CL_Path(:,:,p,q))

#if (PP_dim == 3)
    CALL ChangeBasis1D(3,PP_N,PP_N  ,Vdm_CLN_FV        , dX_3_xi_CL_eta_CL_zeta_LG_Path(:,:,p,q), FV_dX_3_xi_FV_eta_CL_zeta_LG_Path(:,:,p,q))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary, dX_3_xi_CL_eta_CL_zeta_LG_Path(:,:,p,q), FV_dX_3_xi_BC_eta_CL_zeta_LG_Path(:,:,p,q))
#endif
  END DO; END DO

#if (PP_dim == 3)
  DO q=0,PP_N+1; DO p=0,PP_N
    CALL ChangeBasis1D(3,PP_N,PP_N  ,Vdm_CLN_FV        , FV_dX_1_xi_LG_eta_BC_zeta_CL_Path(:,p,q,:), FV_dX_1_xi_LG_eta_BC_zeta_FV_Path(:,p,q,:))
  END DO; END DO

  DO q=0,PP_N  ; DO p=0,PP_N+1
    CALL ChangeBasis1D(3,PP_N,PP_N  ,Vdm_CLN_FV        , FV_dX_2_xi_BC_eta_LG_zeta_CL_Path(:,p,q,:), FV_dX_2_xi_BC_eta_LG_zeta_FV_Path(:,p,q,:))
    CALL ChangeBasis1D(3,PP_N,PP_N  ,Vdm_CLN_FV        , FV_dX_3_xi_BC_eta_CL_zeta_LG_Path(:,p,:,q), FV_dX_3_xi_BC_eta_FV_zeta_LG_Path(:,p,:,q))
  END DO; END DO

  DO q=0,PP_N  ; DO p=0,PP_N
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary, FV_dX_1_xi_LG_eta_FV_zeta_CL_Path(:,p,q,:), FV_dX_1_xi_LG_eta_FV_zeta_BC_Path(:,p,q,:))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary, FV_dX_2_xi_FV_eta_LG_zeta_CL_Path(:,p,q,:), FV_dX_2_xi_FV_eta_LG_zeta_BC_Path(:,p,q,:))
    CALL ChangeBasis1D(3,PP_N,PP_N+1,Vdm_CLN_FVboundary, FV_dX_3_xi_FV_eta_CL_zeta_LG_Path(:,p,:,q), FV_dX_3_xi_FV_eta_BC_zeta_LG_Path(:,p,:,q))
  END DO; END DO
#else
  DO q=0,PP_N+1; DO p=0,PP_N
    FV_dX_1_xi_LG_eta_BC_zeta_FV_Path(:,p,q,:)=FV_dX_1_xi_LG_eta_BC_zeta_CL_Path(:,p,q,:)
  END DO; END DO

  DO q=0,PP_N  ; DO p=0,PP_N+1
    FV_dX_2_xi_BC_eta_LG_zeta_FV_Path(:,p,q,:)=FV_dX_2_xi_BC_eta_LG_zeta_CL_Path(:,p,q,:)
  END DO; END DO

  DO q=0,PP_N  ; DO p=0,PP_N
    FV_dX_1_xi_LG_eta_FV_zeta_BC_Path(:,p,q,0)=FV_dX_1_xi_LG_eta_FV_zeta_CL_Path(:,p,q,0)
    FV_dX_1_xi_LG_eta_FV_zeta_BC_Path(:,p,q,1)=FV_dX_1_xi_LG_eta_FV_zeta_CL_Path(:,p,q,0)
    FV_dX_2_xi_FV_eta_LG_zeta_BC_Path(:,p,q,0)=FV_dX_2_xi_FV_eta_LG_zeta_CL_Path(:,p,q,0)
    FV_dX_2_xi_FV_eta_LG_zeta_BC_Path(:,p,q,1)=FV_dX_2_xi_FV_eta_LG_zeta_CL_Path(:,p,q,0)
  END DO; END DO
#endif

  DO l=0,PP_N
    FV_dX_1_xi_LG_eta_FV_zeta_BC_Path_pq(:,l,:,:)=FV_dX_1_xi_LG_eta_FV_zeta_BC_Path(:,l,:,:)
    FV_dX_1_xi_LG_eta_BC_zeta_FV_Path_pq(:,l,:,:)=FV_dX_1_xi_LG_eta_BC_zeta_FV_Path(:,l,:,:)

    FV_dX_2_xi_FV_eta_LG_zeta_BC_Path_pq(:,l,:,:)=FV_dX_2_xi_FV_eta_LG_zeta_BC_Path(:,:,l,:)
    FV_dX_2_xi_BC_eta_LG_zeta_FV_Path_pq(:,l,:,:)=FV_dX_2_xi_BC_eta_LG_zeta_FV_Path(:,:,l,:)
#if (PP_dim == 3)
    FV_dX_3_xi_FV_eta_BC_zeta_LG_Path_pq(:,l,:,:)=FV_dX_3_xi_FV_eta_BC_zeta_LG_Path(:,:,:,l)
    FV_dX_3_xi_BC_eta_FV_zeta_LG_Path_pq(:,l,:,:)=FV_dX_3_xi_BC_eta_FV_zeta_LG_Path(:,:,:,l)
#endif
  END DO

  ! Calculate distances between FV subcells
  DO l=0,PP_N
    ! left
    x0 = FV_BdryX(l)
    xN = x0 + (FV_BdryX(l+1) - FV_BdryX(l)) * 0.5
    CALL Integrate_Path_BC_zeta(PP_N,PP_N  ,xGP,wGP,wBary,x0,xN,FV_dX_1_xi_LG_eta_FV_zeta_BC_Path_pq,FV_dX_1_xi_LG_eta_FV_zeta_BC_L(:,:,l))
#if (PP_dim == 2)
    CALL Integrate_Path        (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_1_xi_LG_eta_BC_zeta_FV_Path_pq,FV_dX_1_xi_LG_eta_BC_zeta_FV_L(:,:,l))
#else
    CALL Integrate_Path_BC_eta (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_1_xi_LG_eta_BC_zeta_FV_Path_pq,FV_dX_1_xi_LG_eta_BC_zeta_FV_L(:,:,l))
    CALL Integrate_Path_BC_zeta(PP_N,PP_N  ,xGP,wGP,wBary,x0,xN,FV_dX_3_xi_FV_eta_BC_zeta_LG_Path_pq,FV_dX_3_xi_FV_eta_BC_zeta_LG_L(:,:,l))
#endif
    CALL Integrate_Path_BC_zeta(PP_N,PP_N  ,xGP,wGP,wBary,x0,xN,FV_dX_2_xi_FV_eta_LG_zeta_BC_Path_pq,FV_dX_2_xi_FV_eta_LG_zeta_BC_L(:,:,l))
#if (PP_dim == 2)
    CALL Integrate_Path        (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_2_xi_BC_eta_LG_zeta_FV_Path_pq,FV_dX_2_xi_BC_eta_LG_zeta_FV_L(:,:,l))
#else
    CALL Integrate_Path_BC_eta (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_2_xi_BC_eta_LG_zeta_FV_Path_pq,FV_dX_2_xi_BC_eta_LG_zeta_FV_L(:,:,l))
    CALL Integrate_Path_BC_eta (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_3_xi_BC_eta_FV_zeta_LG_Path_pq,FV_dX_3_xi_BC_eta_FV_zeta_LG_L(:,:,l))
#endif

    ! right
    xN = FV_BdryX(l+1)
    x0 = xN - (FV_BdryX(l+1) - FV_BdryX(l)) * 0.5
    CALL Integrate_Path_BC_zeta(PP_N,PP_N  ,xGP,wGP,wBary,x0,xN,FV_dX_1_xi_LG_eta_FV_zeta_BC_Path_pq,FV_dX_1_xi_LG_eta_FV_zeta_BC_R(:,:,l))
#if (PP_dim == 2)
    CALL Integrate_Path        (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_1_xi_LG_eta_BC_zeta_FV_Path_pq,FV_dX_1_xi_LG_eta_BC_zeta_FV_R(:,:,l))
#else
    CALL Integrate_Path_BC_eta (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_1_xi_LG_eta_BC_zeta_FV_Path_pq,FV_dX_1_xi_LG_eta_BC_zeta_FV_R(:,:,l))
    CALL Integrate_Path_BC_zeta(PP_N,PP_N  ,xGP,wGP,wBary,x0,xN,FV_dX_3_xi_FV_eta_BC_zeta_LG_Path_pq,FV_dX_3_xi_FV_eta_BC_zeta_LG_R(:,:,l))
#endif
    CALL Integrate_Path_BC_zeta(PP_N,PP_N  ,xGP,wGP,wBary,x0,xN,FV_dX_2_xi_FV_eta_LG_zeta_BC_Path_pq,FV_dX_2_xi_FV_eta_LG_zeta_BC_R(:,:,l))
#if (PP_dim == 2)
    CALL Integrate_Path        (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_2_xi_BC_eta_LG_zeta_FV_Path_pq,FV_dX_2_xi_BC_eta_LG_zeta_FV_R(:,:,l))
#else
    CALL Integrate_Path_BC_eta (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_2_xi_BC_eta_LG_zeta_FV_Path_pq,FV_dX_2_xi_BC_eta_LG_zeta_FV_R(:,:,l))
    CALL Integrate_Path_BC_eta (PP_N,PP_N+1,xGP,wGP,wBary,x0,xN,FV_dX_3_xi_BC_eta_FV_zeta_LG_Path_pq,FV_dX_3_xi_BC_eta_FV_zeta_LG_R(:,:,l))
#endif
  END DO ! l=0,PP_N

  ! scale metrics for equidistant subcells
  DO d=1,3
    DO i=0,PP_N-1
      FV_Metrics_fTilde_sJ_xi(d,i,:,:,iElem)=(PP_N+1)*0.5*Metrics_fTilde_xi(d,i,:,:)*&
          (FV_dx_XI_L  (:,:,i+1,iElem)+FV_dx_XI_R  (:,:,i,iElem))
    END DO ! i=0,PP_N
    DO j=0,PP_N
      ! NOTE: naive solution: calculation of metrics on FV faces via arithmetic mean
      !FV_Metrics_gTilde_sJ_xi(d,:,j,:,iElem)=(PP_N+1)*0.5*Metrics_gTilde_xi(d,:,j,:,iElem)*&
          !0.5*((FV_dx_ETA_L (0:PP_N-1,:,j,iElem)+FV_dx_ETA_R (0:PP_N-1,:,j,iElem)) + &
               !(FV_dx_ETA_L (1:PP_N  ,:,j,iElem)+FV_dx_ETA_R (1:PP_N  ,:,j,iElem)))
      FV_Metrics_gTilde_sJ_xi(d,:,j,:,iElem)=(PP_N+1)*0.5*Metrics_gTilde_xi(d,:,j,:)*&
          (FV_dX_2_xi_BC_eta_LG_zeta_FV_R(1:PP_N,:,j)+FV_dX_2_xi_BC_eta_LG_zeta_FV_L(1:PP_N,:,j))
    END DO ! j=0,PP_N

    DO i=0,PP_N
      ! NOTE: naive solution: calculation of metrics on FV faces via arithmetic mean
      !FV_Metrics_fTilde_sJ_eta(d,i,:,:,iElem)=(PP_N+1)*0.5*Metrics_fTilde_eta(d,i,:,:,iElem)*&
          !0.5*((FV_dx_XI_L  (0:PP_N-1,:,i,iElem)+FV_dx_XI_R  (0:PP_N-1,:,i,iElem)) + &
               !(FV_dx_XI_L  (1:PP_N  ,:,i,iElem)+FV_dx_XI_R  (1:PP_N  ,:,i,iElem)))
      FV_Metrics_fTilde_sJ_eta(d,i,:,:,iElem)=(PP_N+1)*0.5*Metrics_fTilde_eta(d,i,:,:)*&
          (FV_dX_1_xi_LG_eta_BC_zeta_FV_R(1:PP_N,:,i)+FV_dX_1_xi_LG_eta_BC_zeta_FV_L(1:PP_N,:,i))
    END DO ! i=0,PP_N
    DO j=0,PP_N-1
      FV_Metrics_gTilde_sJ_eta(d,:,j,:,iElem)=(PP_N+1)*0.5*Metrics_gTilde_eta(d,:,j,:)*&
          (FV_dx_ETA_L (:,:,j+1,iElem)+FV_dx_ETA_R (:,:,j,iElem))
    END DO ! j=0,PP_N

    DO i=0,PP_N
      ! NOTE: naive solution: calculation of metrics on FV faces via arithmetic mean
      !FV_Metrics_fTilde_sJ_zeta(d,i,:,:,iElem)=(PP_N+1)*0.5*Metrics_fTilde_zeta(d,i,:,:,iElem)*&
          !0.5*((FV_dx_XI_L  (:,0:PP_NZ-1,i,iElem)+FV_dx_XI_R  (:,0:PP_NZ-1,i,iElem)) + &
               !(FV_dx_XI_L  (:,1:PP_NZ  ,i,iElem)+FV_dx_XI_R  (:,1:PP_NZ  ,i,iElem)))
      FV_Metrics_fTilde_sJ_zeta(d,i,:,:,iElem)=(PP_N+1)*0.5*Metrics_fTilde_zeta(d,i,:,:)*&
          (FV_dX_1_xi_LG_eta_FV_zeta_BC_R(:,1:PP_NZ,i)+FV_dX_1_xi_LG_eta_FV_zeta_BC_L(:,1:PP_NZ,i))
    END DO ! i=0,PP_N
    DO j=0,PP_N
      ! NOTE: naive solution: calculation of metrics on FV faces via arithmetic mean
      !FV_Metrics_gTilde_sJ_zeta(d,:,j,:,iElem)=(PP_N+1)*0.5*Metrics_gTilde_zeta(d,:,j,:,iElem)*&
          !0.5*((FV_dx_ETA_L (:,0:PP_NZ-1,j,iElem)+FV_dx_ETA_R (:,0:PP_NZ-1,j,iElem)) + &
               !(FV_dx_ETA_L (:,1:PP_NZ  ,j,iElem)+FV_dx_ETA_R (:,1:PP_NZ  ,j,iElem)))
      FV_Metrics_gTilde_sJ_zeta(d,:,j,:,iElem)=(PP_N+1)*0.5*Metrics_gTilde_zeta(d,:,j,:)*&
          (FV_dX_2_xi_FV_eta_LG_zeta_BC_R(:,1:PP_NZ,j)+FV_dX_2_xi_FV_eta_LG_zeta_BC_L(:,1:PP_NZ,j))
    END DO ! j=0,PP_N

    FV_Metrics_fTilde_sJ_xi(d,:,:,:,iElem)   = FV_Metrics_fTilde_sJ_xi(d,:,:,:,iElem)   * sJ_xi (:,:,:)
    FV_Metrics_gTilde_sJ_xi(d,:,:,:,iElem)   = FV_Metrics_gTilde_sJ_xi(d,:,:,:,iElem)   * sJ_xi( :,:,:)
    FV_Metrics_fTilde_sJ_eta(d,:,:,:,iElem)  = FV_Metrics_fTilde_sJ_eta(d,:,:,:,iElem)  * sJ_eta (:,:,:)
    FV_Metrics_gTilde_sJ_eta(d,:,:,:,iElem)  = FV_Metrics_gTilde_sJ_eta(d,:,:,:,iElem)  * sJ_eta( :,:,:)

#if (PP_dim == 3)
    DO k=0,PP_NZ
      ! NOTE: naive solution: calculation of metrics on FV faces via arithmetic mean
      !FV_Metrics_hTilde_sJ_xi(d,:,:,k,iElem)=(PP_N+1)*0.5*Metrics_hTilde_xi(d,:,:,k,iElem)*&
        !0.5*((FV_dx_ZETA_L(0:PP_N-1,:,k,iElem)+FV_dx_ZETA_R(0:PP_N-1,:,k,iElem)) + &
             !(FV_dx_ZETA_L(1:PP_N  ,:,k,iElem)+FV_dx_ZETA_R(1:PP_N  ,:,k,iElem)))
      FV_Metrics_hTilde_sJ_xi(d,:,:,k,iElem)=(PP_N+1)*0.5*Metrics_hTilde_xi(d,:,:,k)*&
          (FV_dX_3_xi_BC_eta_FV_zeta_LG_R(1:PP_N,:,k)+FV_dX_3_xi_BC_eta_FV_zeta_LG_L(1:PP_N,:,k))
    END DO ! j=0,PP_N

    DO k=0,PP_NZ
      ! NOTE: naive solution: calculation of metrics on FV faces via arithmetic mean
      !FV_Metrics_hTilde_sJ_eta(d,:,:,k,iElem)=(PP_N+1)*0.5*Metrics_hTilde_eta(d,:,:,k,iElem)*&
        !0.5*((FV_dx_ZETA_L(:,0:PP_N-1,k,iElem)+FV_dx_ZETA_R(:,0:PP_N-1,k,iElem)) + &
             !(FV_dx_ZETA_L(:,1:PP_N  ,k,iElem)+FV_dx_ZETA_R(:,1:PP_N  ,k,iElem)))
      FV_Metrics_hTilde_sJ_eta(d,:,:,k,iElem)=(PP_N+1)*0.5*Metrics_hTilde_eta(d,:,:,k)*&
          (FV_dX_3_xi_FV_eta_BC_zeta_LG_R(:,1:PP_N,k)+FV_dX_3_xi_FV_eta_BC_zeta_LG_L(:,1:PP_N,k))
    END DO ! j=0,PP_N

    DO k=0,PP_NZ-1
      FV_Metrics_hTilde_sJ_zeta(d,:,:,k,iElem)=(PP_N+1)*0.5*Metrics_hTilde_zeta(d,:,:,k)*&
        (FV_dx_ZETA_L(:,:,k+1,iElem)+FV_dx_ZETA_R(:,:,k,iElem))
    END DO ! j=0,PP_N

    FV_Metrics_hTilde_sJ_xi  (d,:,:,:,iElem)=FV_Metrics_hTilde_sJ_xi  (d,:,:,:,iElem)*sJ_xi  (:,:,:)
    FV_Metrics_hTilde_sJ_eta (d,:,:,:,iElem)=FV_Metrics_hTilde_sJ_eta (d,:,:,:,iElem)*sJ_eta( :,:,:)

    FV_Metrics_fTilde_sJ_zeta(d,:,:,:,iElem)=FV_Metrics_fTilde_sJ_zeta(d,:,:,:,iElem)*sJ_zeta(:,:,:)
    FV_Metrics_gTilde_sJ_zeta(d,:,:,:,iElem)=FV_Metrics_gTilde_sJ_zeta(d,:,:,:,iElem)*sJ_zeta(:,:,:)
    FV_Metrics_hTilde_sJ_zeta(d,:,:,:,iElem)=FV_Metrics_hTilde_sJ_zeta(d,:,:,:,iElem)*sJ_zeta(:,:,:)
#endif
  END DO
#endif /* VOLINT_VISC */
#endif /* FV_RECONSTRUCT */
END DO

#if FV_RECONSTRUCT
#if VOLINT_VISC
! Caluclate transformation matrix from xi/eta/zeta  to normal system of master side
DO iSide = 1, nSides
  ! master
  ElemID=SideToElem(S2E_ELEM_ID    ,iSide)
  IF (ElemID .GT. 0) THEN
    locSideID=SideToElem(S2E_LOC_SIDE_ID,iSide)
    flip     = 0
    DO q=0,PP_NZ; DO p=0,PP_N
      ijk2=S2V(:,0,p,q,flip,locSideID)
      FV_Metrics_NormVec_master (1,p,q,iSide)= DOT_PRODUCT(NormVec (:,p,q,1,iSide),FV_Metrics_fTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))
      FV_Metrics_NormVec_master (2,p,q,iSide)= DOT_PRODUCT(NormVec (:,p,q,1,iSide),FV_Metrics_gTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))
      FV_Metrics_TangVec1_master(1,p,q,iSide)= DOT_PRODUCT(TangVec1(:,p,q,1,iSide),FV_Metrics_fTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))
      FV_Metrics_TangVec1_master(2,p,q,iSide)= DOT_PRODUCT(TangVec1(:,p,q,1,iSide),FV_Metrics_gTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))
#if PP_dim==3
      FV_Metrics_NormVec_master (3,p,q,iSide)= DOT_PRODUCT(NormVec (:,p,q,1,iSide),FV_Metrics_hTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))
      FV_Metrics_TangVec1_master(3,p,q,iSide)= DOT_PRODUCT(TangVec1(:,p,q,1,iSide),FV_Metrics_hTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))

      FV_Metrics_TangVec2_master(1,p,q,iSide)= DOT_PRODUCT(TangVec2(:,p,q,1,iSide),FV_Metrics_fTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))
      FV_Metrics_TangVec2_master(2,p,q,iSide)= DOT_PRODUCT(TangVec2(:,p,q,1,iSide),FV_Metrics_gTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))
      FV_Metrics_TangVec2_master(3,p,q,iSide)= DOT_PRODUCT(TangVec2(:,p,q,1,iSide),FV_Metrics_hTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),ElemID))
#endif /* PP_dim==3 */
    END DO; END DO ! p,q=0,PP_N
  END IF
  ! slave
  NB_ElemID=SideToElem(S2E_NB_ELEM_ID    ,iSide)
  IF (NB_ElemID .GT. 0) THEN
    locSideID=SideToElem(S2E_NB_LOC_SIDE_ID,iSide)
    flip     =SideToElem(S2E_FLIP          ,iSide)
    DO q=0,PP_NZ; DO p=0,PP_N
      ijk2=S2V(:,0,p,q,flip,locSideID)
      FV_Metrics_NormVec_slave (1,p,q,iSide)= DOT_PRODUCT(NormVec (:,p,q,1,iSide),FV_Metrics_fTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))
      FV_Metrics_NormVec_slave (2,p,q,iSide)= DOT_PRODUCT(NormVec (:,p,q,1,iSide),FV_Metrics_gTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))
      FV_Metrics_TangVec1_slave(1,p,q,iSide)= DOT_PRODUCT(TangVec1(:,p,q,1,iSide),FV_Metrics_fTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))
      FV_Metrics_TangVec1_slave(2,p,q,iSide)= DOT_PRODUCT(TangVec1(:,p,q,1,iSide),FV_Metrics_gTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))
#if PP_dim==3
      FV_Metrics_NormVec_slave (3,p,q,iSide)= DOT_PRODUCT(NormVec (:,p,q,1,iSide),FV_Metrics_hTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))
      FV_Metrics_TangVec1_slave(3,p,q,iSide)= DOT_PRODUCT(TangVec1(:,p,q,1,iSide),FV_Metrics_hTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))

      FV_Metrics_TangVec2_slave(1,p,q,iSide)= DOT_PRODUCT(TangVec2(:,p,q,1,iSide),FV_Metrics_fTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))
      FV_Metrics_TangVec2_slave(2,p,q,iSide)= DOT_PRODUCT(TangVec2(:,p,q,1,iSide),FV_Metrics_gTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))
      FV_Metrics_TangVec2_slave(3,p,q,iSide)= DOT_PRODUCT(TangVec2(:,p,q,1,iSide),FV_Metrics_hTilde_sJ(:,ijk2(1),ijk2(2),ijk2(3),NB_ElemID))
#endif /* PP_dim==3 */
    END DO; END DO ! p,q=0,PP_N
  END IF
END DO ! iSide = 1, nSides
#endif /* VOLINT_VISC */

! distances at big mortar interfaces must be distributed to the smaller sides
FV_Elems_master = 1 ! Force use of FV mortar matrices in U_Mortar routine
#if USE_MPI
MPIRequest=0
! distances at MPI slave sides must be transmitted to master sides
CALL U_Mortar1(FV_dx_master,FV_dx_slave,doMPISides=.TRUE.)
CALL StartReceiveMPIData(FV_dx_slave, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest(:,SEND),SendID=2)
CALL StartSendMPIData(   FV_dx_slave, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest(:,RECV),SendID=2)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest) !Send MINE -receive YOUR
#endif
CALL U_Mortar1(FV_dx_master,FV_dx_slave,doMPISides=.FALSE.)
#if USE_MPI
MPIRequest=0
! distances at MPI master sides must be transmitted to slave sides (required in preconditioner)
CALL StartReceiveMPIData(FV_dx_master, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest(:,SEND),SendID=1)
CALL StartSendMPIData(   FV_dx_master, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest(:,RECV),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest) !Send MINE -receive YOUR
#endif

FV_Elems_master = 0 ! Force use of DG mortar matrices in U_Mortar routine
#if USE_MPI
CALL U_Mortar1(DG_dx_master,DG_dx_slave,doMPISides=.TRUE.)
CALL StartReceiveMPIData(DG_dx_slave, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest(:,SEND),SendID=2)
CALL StartSendMPIData(   DG_dx_slave, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest(:,RECV),SendID=2)
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
  ! master=DG, slave=FV
  FV_dx_Face(:,:,2) = FV_dx_slave(1,:,:,SideID) + DG_dx_master(1,:,:,SideID)
  ! master=FV, slave=FV
  FV_dx_Face(:,:,3) = FV_dx_slave(1,:,:,SideID) + FV_dx_master(1,:,:,SideID)
  ! precompute inverse
  FV_sdx_Face(:,:,:,SideID) = 1. / FV_dx_Face
END DO
#endif /* FV_RECONSTRUCT */

SWRITE(UNIT_stdOut,'(A)')' Done !'

END SUBROUTINE FV_CalcMetrics


#if FV_RECONSTRUCT
!==================================================================================================================================
!> Computes the distance between two points along a path given in 1D reference coordinates
!==================================================================================================================================
SUBROUTINE Integrate_Path(Nloc2,Nloc,xGP,wGP,wBary,x0,xN,FV_Path_1D,FV_Length)
! MODULES
USE MOD_Basis       ,ONLY: InitializeVandermonde
USE MOD_ChangeBasis ,ONLY: ChangeBasis1D
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc2                                     !< degree of path polynomial
INTEGER,INTENT(IN) :: Nloc                                      !< number of points to compute (Nloc+1)**2
REAL,INTENT(IN)    :: xGP(  0:Nloc2)                            !< parametric coords
REAL,INTENT(IN)    :: wGP(  0:Nloc2)                            !< integration weights
REAL,INTENT(IN)    :: wBary(0:Nloc2)                            !< interpolations weights
REAL,INTENT(IN)    :: x0                                        !< start point
REAL,INTENT(IN)    :: xN                                        !< end point
REAL,INTENT(INOUT) :: FV_Path_1D(3,0:Nloc2,0:Nloc,0:ZDIM(Nloc)) !< path polynomial
REAL,INTENT(OUT)   :: FV_Length(           0:Nloc,0:ZDIM(Nloc)) !< distance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: VDM(0:Nloc2,0:Nloc2)
REAL               :: SubxGP(1,0:Nloc2)
REAL               :: FV_Path_Cut(3,0:Nloc2)
INTEGER            :: q,p,l
!===================================================================================================================================
subxGP(1,:) = x0 + (xGP + 1.)/2. * (xN-x0)
CALL InitializeVandermonde(Nloc2,Nloc2,wBary,xGP,subxGP(1,:),Vdm)

FV_Length=0.
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  ! path to integrate in ref space [-1,1]
  CALL ChangeBasis1D(3,Nloc2,Nloc2,Vdm,FV_Path_1D(:,:,p,q), FV_Path_Cut)
  ! integrate path
  DO l=0,Nloc2
    FV_Length(p,q) = FV_Length(p,q) + NORM2(FV_Path_Cut(:,l)) * wGP(l)
  END DO
END DO; END DO! p,q=0,Nloc
FV_Length=FV_Length*0.5*(xN-x0) ! *0.5 since reference element has width=2
END SUBROUTINE Integrate_Path

#if VOLINT_VISC
!==================================================================================================================================
!> Computes the distance between two points along a path given in 1D reference coordinates
!==================================================================================================================================
SUBROUTINE Integrate_Path_BC_eta(Nloc2,Nloc,xGP,wGP,wBary,x0,xN,FV_Path_1D,FV_Length)
! MODULES
USE MOD_Basis       ,ONLY: InitializeVandermonde
USE MOD_ChangeBasis ,ONLY: ChangeBasis1D
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc2                                   ! < degree of path polynomial
INTEGER,INTENT(IN) :: Nloc                                    ! < number of points to compute (Nloc+1)**2
REAL,INTENT(IN)    :: xGP(  0:Nloc2)                          ! < parametric coords
REAL,INTENT(IN)    :: wGP(  0:Nloc2)                          ! < integration weights
REAL,INTENT(IN)    :: wBary(0:Nloc2)                          ! < interpolations weights
REAL,INTENT(IN)    :: x0                                      ! < start point
REAL,INTENT(IN)    :: xN                                      ! < end point
REAL,INTENT(INOUT) :: FV_Path_1D(3,0:Nloc2,0:Nloc,0:Nloc2)    ! < path polynomial
REAL,INTENT(OUT)   :: FV_Length(           0:Nloc,0:Nloc2)    ! < distance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: VDM(0:Nloc2,0:Nloc2)
REAL               :: SubxGP(1,0:Nloc2)
REAL               :: FV_Path_Cut(3,0:Nloc2)
INTEGER            :: q,p,l
!===================================================================================================================================
subxGP(1,:) = x0 + (xGP + 1.)/2. * (xN-x0)
CALL InitializeVandermonde(Nloc2,Nloc2,wBary,xGP,subxGP(1,:),Vdm)

FV_Length=0.
DO q=0,Nloc2; DO p=0,Nloc
  ! path to integrate in ref space [-1,1]
  CALL ChangeBasis1D(3,Nloc2,Nloc2,Vdm,FV_Path_1D(:,:,p,q), FV_Path_Cut)
  ! integrate path
  DO l=0,Nloc2
    FV_Length(p,q) = FV_Length(p,q) + NORM2(FV_Path_Cut(:,l)) * wGP(l)
  END DO
END DO; END DO! p,q=0,Nloc
FV_Length=FV_Length*0.5*(xN-x0) ! *0.5 since reference element has width=2
END SUBROUTINE Integrate_Path_BC_eta

SUBROUTINE Integrate_Path_BC_zeta(Nloc2,Nloc,xGP,wGP,wBary,x0,xN,FV_Path_1D,FV_Length)
! MODULES
USE MOD_Basis       ,ONLY: InitializeVandermonde
USE MOD_ChangeBasis ,ONLY: ChangeBasis1D
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc2                                   ! < degree of path polynomial
INTEGER,INTENT(IN) :: Nloc                                    ! < number of points to compute (Nloc+1)**2
REAL,INTENT(IN)    :: xGP(  0:Nloc2)                          ! < parametric coords
REAL,INTENT(IN)    :: wGP(  0:Nloc2)                          ! < integration weights
REAL,INTENT(IN)    :: wBary(0:Nloc2)                          ! < interpolations weights
REAL,INTENT(IN)    :: x0                                      ! < start point
REAL,INTENT(IN)    :: xN                                      ! < end point
REAL,INTENT(INOUT) :: FV_Path_1D(3,0:Nloc2,0:Nloc,0:ZDIM(Nloc)+1) ! < path polynomial
REAL,INTENT(OUT)   :: FV_Length(           0:Nloc,0:ZDIM(Nloc)+1) ! < distance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: VDM(0:Nloc2,0:Nloc2)
REAL               :: SubxGP(1,0:Nloc2)
REAL               :: FV_Path_Cut(3,0:Nloc2)
INTEGER            :: q,p,l
!===================================================================================================================================
subxGP(1,:) = x0 + (xGP + 1.)/2. * (xN-x0)
CALL InitializeVandermonde(Nloc2,Nloc2,wBary,xGP,subxGP(1,:),Vdm)

FV_Length=0.
DO q=0,ZDIM(Nloc)+1; DO p=0,Nloc
  ! path to integrate in ref space [-1,1]
  CALL ChangeBasis1D(3,Nloc2,Nloc2,Vdm,FV_Path_1D(:,:,p,q), FV_Path_Cut)
  ! integrate path
  DO l=0,Nloc2
    FV_Length(p,q) = FV_Length(p,q) + NORM2(FV_Path_Cut(:,l)) * wGP(l)
  END DO
END DO; END DO! p,q=0,Nloc
FV_Length=FV_Length*0.5*(xN-x0) ! *0.5 since reference element has width=2
END SUBROUTINE Integrate_Path_BC_zeta
#endif /*VOLINT_VISC*/
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

#if FV_RECONSTRUCT
#if VOLINT_VISC
SDEALLOCATE(FV_Metrics_fTilde_sJ)
SDEALLOCATE(FV_Metrics_gTilde_sJ)
SDEALLOCATE(FV_Metrics_hTilde_sJ)

SDEALLOCATE(FV_Metrics_NormVec_master)
SDEALLOCATE(FV_Metrics_TangVec1_master)
SDEALLOCATE(FV_Metrics_NormVec_slave)
SDEALLOCATE(FV_Metrics_TangVec1_slave)
#if (PP_dim == 3)
SDEALLOCATE(FV_Metrics_TangVec2_master)
SDEALLOCATE(FV_Metrics_TangVec2_slave)
#endif

SDEALLOCATE(FV_Metrics_fTilde_sJ_xi)
SDEALLOCATE(FV_Metrics_gTilde_sJ_xi)
#if (PP_dim == 3)
SDEALLOCATE(FV_Metrics_hTilde_sJ_xi)
#endif
SDEALLOCATE(FV_Metrics_fTilde_sJ_eta)
SDEALLOCATE(FV_Metrics_gTilde_sJ_eta)
#if (PP_dim == 3)
SDEALLOCATE(FV_Metrics_hTilde_sJ_eta)
#endif
SDEALLOCATE(FV_Metrics_fTilde_sJ_zeta)
SDEALLOCATE(FV_Metrics_gTilde_sJ_zeta)
#if (PP_dim == 3)
SDEALLOCATE(FV_Metrics_hTilde_sJ_zeta)
#endif
#endif /* VOLINT_VISC */
#endif

SDEALLOCATE(FV_Elems_master)
END SUBROUTINE FinalizeFV_Metrics


END MODULE MOD_FV_Metrics
#endif /* FV_ENABLED */

