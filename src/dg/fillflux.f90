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
!> \brief Prepares the integrand for the flux integrals, i.e. computes the common inteface fluxes.
!> This module prepares the computation of the fluxes over the sides by filling the two parts (advective and viscous) 
!> of the common interface flux by calling the associated numerical flux functions.
!> We distinguish between inner sides, sides with boundary conditions and mpi sides.
!==================================================================================================================================
MODULE MOD_FillFlux
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------


INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE


PUBLIC::FillFlux
!==================================================================================================================================



CONTAINS



!==================================================================================================================================
!> Computes the fluxes for inner sides, MPI sides where the local proc is "master"  and boundary conditions.
!> The flux computation is performed seperately for advection and diffusion fluxes in case
!> parabolic terms are considered.
!> If selective overintegration is used, the advection fluxes are evaluated at a higher
!> polynomial degree and then projected down to the polynomial degree of the solution
!==================================================================================================================================
SUBROUTINE FillFlux(t,Flux_master,Flux_slave,U_master,U_slave,UPrim_master,UPrim_slave,doMPISides)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,         ONLY: U_masterO,U_slaveO,UPrim_masterO,UPrim_slaveO,FluxO
USE MOD_Mesh_Vars,       ONLY: NormVec, TangVec1, TangVec2, SurfElem, Face_xGP
USE MOD_Mesh_Vars,       ONLY: NormVecO,TangVec1O,TangVec2O,SurfElemO,Face_xGPO
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars,       ONLY: nSides,firstBCSide
USE MOD_ChangeBasis,     ONLY: ChangeBasis2D,ChangeBasis2D_selective
USE MOD_Riemann,         ONLY: Riemann
USE MOD_GetBoundaryFlux, ONLY: GetBoundaryFlux
USE MOD_Overintegration_Vars, ONLY: OverintegrationType,NOver,VdmNOverToN,VdmNToNOver
USE MOD_EOS,             ONLY: ConsToPrim
USE MOD_Mesh_Vars,       ONLY: nBCSides
#if PARABOLIC
USE MOD_Riemann,         ONLY: ViscousFlux
USE MOD_Lifting_Vars,    ONLY: gradUx_master ,gradUy_master ,gradUz_master ,gradUx_slave,gradUy_slave,gradUz_slave
USE MOD_Lifting_Vars,    ONLY: gradUx_masterO,gradUy_masterO,gradUz_masterO
#endif /*PARABOLIC*/
#ifdef EDDYVISCOSITY
USE MOD_EddyVisc_Vars,   ONLY: DeltaS_master,DeltaS_slave,SGS_Ind_master,SGS_Ind_slave
#endif
#if FV_ENABLED
USE MOD_FV
USE MOD_FV_Vars
#endif
USE MOD_EOS             ,ONLY: PrimToCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  !< = .TRUE. only MINE (where the proc is master)  MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(IN)    :: t           !< physical time required for BC state evaluation in case of time dependent BCs
REAL,INTENT(OUT)   :: Flux_master(1:PP_nVar,0:PP_N,0:PP_NZ,1:nSides)  ! sum of advection and diffusion fluxes across the boundary
REAL,INTENT(OUT)   :: Flux_slave (1:PP_nVar,0:PP_N,0:PP_NZ,1:nSides)  ! sum of advection and diffusion fluxes across the boundary
REAL,INTENT(INOUT) :: U_master(PP_nVar,0:PP_N, 0:PP_NZ, 1:nSides) !< solution on master sides
REAL,INTENT(INOUT) :: U_slave( PP_nVar,0:PP_N, 0:PP_NZ, 1:nSides) !< solution on slave sides
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N, 0:PP_NZ, 1:nSides) !< solution on master sides
REAL,INTENT(IN)    :: UPrim_slave( PP_nVarPrim,0:PP_N, 0:PP_NZ, 1:nSides) !< solution on slave sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID
#if PARABOLIC
REAL    :: FluxV_loc(PP_nVar,0:PP_N, 0:PP_NZ)
#endif /*PARABOLIC*/
INTEGER :: FV_Elems_Max(1:nSides) ! 0 if both sides DG, 1 else
LOGICAL :: addToOutput_loc
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver for advection and viscous terms
! Set the side range according to MPI or no MPI 
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides (where the local proc is master) 
  firstSideID_wo_BC = firstMPISide_MINE
  firstSideID = firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  firstSideID = firstBCSide    ! include BCs for master sides
   lastSideID = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
  FV_Elems_Max(SideID) = MAX(FV_Elems_master(SideID),FV_Elems_slave(SideID))
END DO 

! =============================
! Workflow:
!
! numbers marked with curly brackets {..} are only peformed for selective overintegration
! 
! {0.} interpolate U_master/slave and gradU.._master/slave from N to NOver for all DG/DG interfaces
!  1.  compute flux for non-BC sides
!            no overintegration: store in Flux_master
!         [with overintegration: store in FluxO]
!  1.1) advective flux
!  1.2) viscous flux
!  1.3) add up viscous flux to Flux_master 
!  2.  compute flux for BC sides
!            no overintegration: store in Flux_master
!         [with overintegration: store in FluxO]
!  3.  multiply by SurfElem 
! {4.} project FluxO from NOver to N for all DG/DG interfaces and add to Flux_master
!  5.  copy flux from Flux_master to Flux_slave
!  6.  convert FV flux to DG flux at mixed interfaces
!==============================

! 0. interpolate U_master/slave from N to NOver for all pure DG/DG interfaces
!    interpolate gradUx/y/z_master/slave from N to NOver for all DG BC sides 
IF(OverintegrationType.EQ.SELECTIVE)THEN
  ! for surface overintegration, solution (and gradients for BCs) at boundaries at NOver are required:
  ! Interpolate surface states from N to Nover
#if FV_ENABLED
  CALL ChangeBasis2D_selective(PP_nVar,PP_N,Nover,1,nSides,firstSideID    ,lastSideID,VdmNToNOver,U_master,U_masterO,FV_Elems_Sum,0)
  CALL ChangeBasis2D_selective(PP_nVar,PP_N,Nover,1,nSides,firstSideID_wo_BC,lastSideID,VdmNToNOver,U_slave,U_slaveO,FV_Elems_Sum,0)
#else
  CALL ChangeBasis2D_selective(PP_nVar,PP_N,Nover,1,nSides,firstSideID    ,lastSideID,VdmNToNOver,U_master,U_masterO)
  CALL ChangeBasis2D_selective(PP_nVar,PP_N,Nover,1,nSides,firstSideID_wo_BC,lastSideID,VdmNToNOver,U_slave,U_slaveO)
#endif
  DO SideID=firstSideID,lastSideID
    IF (FV_Elems_Sum(SideID).EQ.0) THEN
      CALL ConsToPrim(NOver,UPrim_masterO(:,:,:,SideID), U_masterO(:,:,:,SideID))
    END IF  
  END DO 
  DO SideID=firstSideID_wo_BC,lastSideID  
    IF (FV_Elems_Sum(SideID).EQ.0) THEN
      CALL ConsToPrim(NOver,UPrim_slaveO(:,:,:,SideID), U_slaveO(:,:,:,SideID))
    END IF  
  END DO 

#if PARABOLIC
  IF(.NOT.doMPISides)THEN
#if FV_ENABLED
    CALL ChangeBasis2D_selective(PP_nVarPrim,PP_N,Nover,1,nBCSides,1,nBCSides,VdmNToNOver,&
        gradUx_master(:,:,:,1:nBCSides),gradUx_masterO,FV_Elems_master,0)
    CALL ChangeBasis2D_selective(PP_nVarPrim,PP_N,Nover,1,nBCSides,1,nBCSides,VdmNToNOver,&
        gradUy_master(:,:,:,1:nBCSides),gradUy_masterO,FV_Elems_master,0)
    CALL ChangeBasis2D_selective(PP_nVarPrim,PP_N,Nover,1,nBCSides,1,nBCSides,VdmNToNOver,&
        gradUz_master(:,:,:,1:nBCSides),gradUz_masterO,FV_Elems_master,0)
#else
    CALL ChangeBasis2D_selective(PP_nVarPrim,PP_N,Nover,1,nBCSides,1,nBCSides,VdmNToNOver,&
        gradUx_master(:,:,:,1:nBCSides),gradUx_masterO)
    CALL ChangeBasis2D_selective(PP_nVarPrim,PP_N,Nover,1,nBCSides,1,nBCSides,VdmNToNOver,&
        gradUy_master(:,:,:,1:nBCSides),gradUy_masterO)
    CALL ChangeBasis2D_selective(PP_nVarPrim,PP_N,Nover,1,nBCSides,1,nBCSides,VdmNToNOver,&
        gradUz_master(:,:,:,1:nBCSides),gradUz_masterO)
#endif
  END IF
#endif
END IF

! 1. compute flux for non-BC sides
DO SideID=firstSideID_wo_BC,lastSideID
  ! 1.1) advective part of flux
  IF (FV_Elems_Sum(SideID).EQ.0) THEN ! DG/DG interfaces
    IF(OverintegrationType.EQ.SELECTIVE)THEN
      ! 1.1 [over_DG]) selective overintegration of advective flux
      CALL Riemann(NOver,FluxO(:,:,:,SideID),   &
          U_masterO    ( :,:,:,SideID),U_slaveO    (  :,:,:,SideID),&
          UPrim_masterO( :,:,:,SideID),UPrim_slaveO(  :,:,:,SideID),&
          NormVecO (:,:,:,0,SideID),&
          TangVec1O(:,:,:,0,SideID),&
          TangVec2O(:,:,:,0,SideID),doBC=.FALSE.)
    ELSE
      ! 1.1 [DG]) no selective overintegration of advective flux
      CALL Riemann(PP_N,Flux_master(:,:,:,SideID),&
          U_master    (:,:,:,SideID),U_slave    (:,:,:,SideID),       &
          UPrim_master(:,:,:,SideID),UPrim_slave(:,:,:,SideID),       &
          NormVec (:,:,:,0,SideID), &
          TangVec1(:,:,:,0,SideID), &
          TangVec2(:,:,:,0,SideID),doBC=.FALSE.)
    END IF
  ELSE ! at least one of the faces FV
    ! 1.1 [FV]) no selective overintegration of advective flux for all sides involving at least one FV face
    CALL Riemann(PP_N,Flux_master(:,:,:,SideID),&
        U_master    (:,:,:,SideID),U_slave    (:,:,:,SideID),       &
        UPrim_master(:,:,:,SideID),UPrim_slave(:,:,:,SideID),       &
        NormVec (:,:,:,1,SideID), &
        TangVec1(:,:,:,1,SideID), &
        TangVec2(:,:,:,1,SideID),doBC=.FALSE.)
  END IF

#if PARABOLIC
  ! 1.2) Fill viscous flux for non-BC sides
  CALL ViscousFlux(PP_N,FluxV_loc, UPrim_master(:,:,:,SideID), UPrim_slave  (:,:,:,SideID), &
      gradUx_master(:,:,:,SideID),gradUy_master(:,:,:,SideID), gradUz_master(:,:,:,SideID),&
      gradUx_slave (:,:,:,SideID),gradUy_slave (:,:,:,SideID), gradUz_slave (:,:,:,SideID),&
      NormVec(:,:,:,FV_Elems_Max(SideID),SideID)&
#ifdef EDDYVISCOSITY
      ,DeltaS_master(SideID),DeltaS_slave(SideID),SGS_Ind_master(1,:,:,SideID),SGS_Ind_slave(1,:,:,SideID),&
          Face_xGP(:,:,:,FV_Elems_Max(SideID),SideID)&
#endif
  )
  ! 1.3) add up viscous flux
  IF ((FV_Elems_Sum(SideID).EQ.0).AND.(OverintegrationType.EQ.SELECTIVE)) THEN ! DG/DG interface and selective overintegration
    Flux_master(:,:,:,SideID) = FluxV_loc(:,:,:) ! here misses the advective flux, which will be added in 5.) from FluxO
  ELSE
    Flux_master(:,:,:,SideID) = Flux_master(:,:,:,SideID) + FluxV_loc(:,:,:)
  END IF
#endif /*PARABOLIC*/
END DO ! SideID


! 2. Compute the fluxes at the boundary conditions: 1..nBCSides
IF(.NOT.doMPISides)THEN
  DO SideID=1,nBCSides
    ! if DG element and selective overintegration 
    IF ((FV_Elems_Sum(SideID).EQ.0).AND.(OverintegrationType.EQ.SELECTIVE)) THEN
     CALL GetBoundaryFlux(SideID,t,NOver,FluxO,UPrim_masterO, &
#if PARABOLIC
         gradUx_masterO,gradUy_masterO,gradUz_masterO, &
#endif
         NormVecO,TangVec1O,TangVec2O,Face_xGPO)
    ELSE 
      ! no selective overintegration or FV element
      CALL GetBoundaryFlux(SideID,t,PP_N,Flux_master,UPrim_master, &
#if PARABOLIC
          gradUx_master,gradUy_master,gradUz_master,        &
#endif
              NormVec,TangVec1,TangVec2,Face_xGP)
    END IF
  END DO
END IF ! .NOT. MPISIDES


! 3. multiply by SurfElem
DO SideID=firstSideID,lastSideID
  ! multiply with SurfElem
  DO q=0,PP_NZ; DO p=0,PP_N
    Flux_master(:,p,q,SideID) = Flux_master(:,p,q,SideID) * SurfElem(p,q,FV_Elems_Max(SideID),SideID)
  END DO; END DO
  ! multiply FluxO with SurfElem0 before projection from NOver to N
  IF ((FV_Elems_Sum(SideID).EQ.0).AND.(OverintegrationType.EQ.SELECTIVE)) THEN
    DO q=0,PP_NOverZ; DO p=0,NOver
      FluxO(:,p,q,SideID) = FluxO(:,p,q,SideID) * SurfElemO(p,q,0,SideID)
    END DO; END DO
  END IF
END DO ! SideID


! 4. project flux from NOver to N 
IF (OverintegrationType.EQ.SELECTIVE) THEN
#if PARABOLIC
  addToOutput_loc = .TRUE.
#else 
  addToOutput_loc = .FALSE.
#endif      
    ! project Flux back on N: For parabolic, add the contribution to the flux (on N), for Euler, this is the full flux already
#if FV_ENABLED
    CALL ChangeBasis2D_selective(PP_nVar,Nover,PP_N,1,nSides,firstSideID,lastSideID,VdmNOverToN,FluxO,Flux_master,FV_Elems_Sum,0,&
        addToOutput=addToOutput_loc)
#else
    CALL ChangeBasis2D_selective(PP_nVar,Nover,PP_N,1,nSides,firstSideID,lastSideID,VdmNOverToN,FluxO,Flux_master,&
        addToOutput=addToOutput_loc)
#endif
END IF


! 5. copy flux from master side to slave side
DO SideID=firstSideID ,lastSideID
  Flux_slave(:,:,:,SideID) = Flux_master(:,:,:,SideID)
END DO !SideID


#if FV_ENABLED
! 6. convert flux on FV points to DG points for all DG faces at mixed interfaces
! only inner sides can be mixed (BC do not require a change basis)
CALL ChangeBasis2D_selective(PP_nVar,PP_N,1,nSides,firstSideID_wo_BC,lastSideID,FV_sVdm,Flux_master,FV_Elems_Sum,2)
CALL ChangeBasis2D_selective(PP_nVar,PP_N,1,nSides,firstSideID_wo_BC,lastSideID,FV_sVdm,Flux_slave ,FV_Elems_Sum,1)
#endif

END SUBROUTINE FillFlux



END MODULE MOD_FillFlux
