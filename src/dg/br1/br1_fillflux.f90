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
#if PARABOLIC
#include "flexi.h"

!==================================================================================================================================
!> Contains routines that fill the inter element and boundary surface fluxes for the BR1 scheme.
!==================================================================================================================================
MODULE MOD_Lifting_FillFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_FillFlux
  MODULE PROCEDURE Lifting_FillFlux
END INTERFACE

INTERFACE Lifting_FillFlux_BC
  MODULE PROCEDURE Lifting_FillFlux_BC
END INTERFACE

PUBLIC::Lifting_FillFlux,Lifting_FillFlux_BC
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> \brief Computes the BR1 surface fluxes in direction "dir" used to compute the gradient in direction "dir",
!> surfelem contribution is considered as well.
!>
!> Fills the interior surface fluxes for the BR1 scheme. Is called on a per direction basis. The numerical flux in the
!> BR1 lifting is simply taken as the arithmetic mean of the solution. 
!> The physical flux is  multiplied by the normal vector and the surface element contribution to transform into reference space,
!> see e.g. "Explicit discontinuous Galerkin methods for unsteady problems" (Hindenlang et al. 2012) for details.
!> For the strong form, in the surface integral the inner solution is substracted form the numerical flux. Since the numerical
!> flux is \f$ \frac{1}{2} (U^+ + U^-) \f$ and the inner solution is simply \f$ U^- \f$ for the master side
!> , the surface flux will become \f$ \frac{1}{2} (U^+ + U^-) - U^- = \frac{1}{2} (U^+ - U^-) \f$ for the strong form.
!> 
!> The flux is filled for the master side, the contribution for the slave side (which is different because the inner solution
!> is equal to \f$ U^+ \f$) is taken into account in the SurfInt routine.
!==================================================================================================================================
SUBROUTINE Lifting_FillFlux(dir,UPrimface_master,UPrimface_slave,Flux,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Lifting_Vars,    ONLY: doWeakLifting
USE MOD_Mesh_Vars,       ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY: nSides
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,   lastInnerSide
USE MOD_Mesh_Vars,       ONLY: firstMPISide_MINE,lastMPISide_MINE
#if FV_ENABLED
USE MOD_FV_Vars         ,ONLY: FV_Elems_Sum,FV_sVdm
USE MOD_ChangeBasis     ,ONLY: ChangeBasis2D_selective
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: dir                                                                 !< direction of gradients (1=x,2=y,3=z)
LOGICAL,INTENT(IN) :: doMPISides                                                          !< = .TRUE. only MINE MPISides are filled,
                                                                                          !< =.FALSE. InnerSides
REAL,INTENT(IN)    :: UPrimface_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< Solution on master sides
REAL,INTENT(IN)    :: UPrimface_slave( PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< Solution on slave sides
REAL,INTENT(OUT)   :: Flux(1:PP_nVarPrim,0:PP_N,0:PP_N,nSides)                            !< Lifting-Flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID,sig
#if FV_ENABLED  
REAL               :: UPrim_glob(1:PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
#endif
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  firstSideID = firstInnerSide
   lastSideID =  lastInnerSide
END IF
! BR1 uses arithmetic mean value of states for the Riemann flux
IF(doWeakLifting)THEN
  sig = 1
ELSE
  ! for strong form subtract solution from the inside
  sig = -1
END IF

#if FV_ENABLED
CALL ChangeBasis2D_selective(PP_nVarPrim,PP_N,PP_N,1,nSides,firstSideID,lastSideID,FV_sVdm,UPrimface_master,UPrim_glob,FV_Elems_Sum,1)
CALL ChangeBasis2D_selective(PP_nVarPrim,PP_N,PP_N,1,nSides,firstSideID,lastSideID,FV_sVdm,UPrimface_slave ,UPrim_glob,FV_Elems_Sum,2)
#endif

DO SideID=firstSideID,lastSideID
#if FV_ENABLED
  SELECT CASE(FV_Elems_Sum(SideID))
  CASE(0) ! both DG
#endif
    Flux(:,:,:,SideID) = sig*UPrimface_master(:,:,:,SideID) + UPrimface_slave(:,:,:,SideID)
#if FV_ENABLED
  CASE(1) ! master=FV, slave=DG
    Flux(:,:,:,SideID) = sig*UPrim_glob     (:,:,:,SideID) + UPrimface_slave(:,:,:,SideID)
  CASE(2) ! master=DG, slave=FV
    Flux(:,:,:,SideID) = sig*UPrimface_master(:,:,:,SideID) + UPrim_glob    (:,:,:,SideID)
  CASE(3) ! both FV
    CYCLE 
  END SELECT
#endif

  ! BR1 uses arithmetic mean value of states for the Riemann flux
  DO q=0,PP_N
    DO p=0,PP_N
      Flux(:,p,q,SideID)=0.5*Flux(:,p,q,SideID)*NormVec(dir,p,q,0,SideID)*SurfElem(p,q,0,SideID)
    END DO
  END DO
END DO ! SideID

END SUBROUTINE Lifting_FillFlux



!==================================================================================================================================
!> \brief Computes the BR1 surface fluxes for boundary conditions, surfelem contribution is considered as well.
!> 
!> Fills the boundary surface fluxes for the BR1 scheme. Is called for all directions at once. The numerical flux in the
!> BR1 lifting is simply taken as the arithmetic mean of the solution. 
!> This routine calls the equation system dependant routine Lifting_GetBoundaryFlux which will fill the flux depending on the
!> boundary condition that has to be applied. The Lifting_GetBoundaryFlux routine will also differentiate between weak and
!> strong form and already multiply the flux by the surface element.
!> 
!> The returned flux is multiplied by the normal vector to transform into reference space (in combination with the 
!> surface element).
!> 
!> The flux is filled for the master side, the contribution for the slave side (which is different because the inner solution
!> is equal to \f$ U^+ \f$) is taken into account in the SurfInt routine.
!==================================================================================================================================
SUBROUTINE Lifting_FillFlux_BC(t,UPrim_master,FluxX,FluxY,FluxZ)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,nSides,nBCSides
USE MOD_GetBoundaryFlux, ONLY: Lifting_GetBoundaryFlux
#if FV_ENABLED  
USE MOD_FV_Vars         ,ONLY: FV_Elems_master
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: t                                                !< Current solution time
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< Primitive solution from the inside
REAL,INTENT(OUT)   :: FluxX(       PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< Lifting boundary flux in x direction
REAL,INTENT(OUT)   :: FluxY(       PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< Lifting boundary flux in y direction
REAL,INTENT(OUT)   :: FluxZ(       PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< Lifting boundary flux in z direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q
!==================================================================================================================================
! fill flux for boundary sides

! only compute boundary flux once use FluxZ as temp storage
CALL Lifting_GetBoundaryFlux(t,UPrim_master,FluxZ)

DO SideID=1,nBCSides
#if FV_ENABLED
  IF (FV_Elems_master(SideID).GT.0) CYCLE
#endif
  DO q=0,PP_N; DO p=0,PP_N
    FluxX(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(1,p,q,0,SideID)
  END DO; END DO
END DO
DO SideID=1,nBCSides
#if FV_ENABLED  
  IF (FV_Elems_master(SideID).GT.0) CYCLE
#endif  
  DO q=0,PP_N; DO p=0,PP_N
    FluxY(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(2,p,q,0,SideID)
  END DO; END DO
END DO
DO SideID=1,nBCSides
#if FV_ENABLED  
  IF (FV_Elems_master(SideID).GT.0) CYCLE
#endif
  DO q=0,PP_N; DO p=0,PP_N
    FluxZ(:,p,q,SideID)=FluxZ(:,p,q,SideID)*NormVec(3,p,q,0,SideID)
  END DO; END DO
END DO

END SUBROUTINE Lifting_FillFlux_BC

END MODULE MOD_Lifting_FillFlux
#endif /*PARABOLIC*/
