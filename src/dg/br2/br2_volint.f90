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
!> Routines for computing the lifting volume integral for the BR2 scheme
!> Computes the volume integral contribution based on the derivative of U and updates gradU
!==================================================================================================================================
MODULE MOD_Lifting_VolInt
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_VolInt
  MODULE PROCEDURE Lifting_VolInt_Conservative
  MODULE PROCEDURE Lifting_VolInt_Nonconservative
END INTERFACE

PUBLIC::Lifting_VolInt
!==================================================================================================================================
CONTAINS


!==================================================================================================================================
!> \brief Computes the volume integral of the BR2 scheme in non-conservative form for all directions
!>
!> Requires lifting in strong form
!> In the non conservative form of the volume integral in BR1 we first differentiate the flux (which is the solution in BR2) and
!> then apply the metric terms. This is the fastest implementation of the volume integral and only available in strong form. 
!==================================================================================================================================
SUBROUTINE Lifting_VolInt_Nonconservative(UPrim,gradUx,gradUy,gradUz)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars            ,ONLY: D_T
USE MOD_Mesh_Vars          ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde   ! metrics
USE MOD_Mesh_Vars          ,ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars            ,ONLY: FV_Elems
USE MOD_FV_Vars            ,ONLY: FV_Metrics_fTilde_sJ,FV_Metrics_gTilde_sJ,FV_Metrics_hTilde_sJ  ! metrics
#if PP_dim == 3
USE MOD_FV_Vars            ,ONLY: gradUxi_central, gradUeta_central, gradUzeta_central
#else
USE MOD_FV_Vars            ,ONLY: gradUxi_central, gradUeta_central
#endif
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)            :: UPrim( PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution
REAL,INTENT(OUT)           :: gradUx(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< gradients in x-direction
REAL,INTENT(OUT)           :: gradUy(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< gradients in y-direction
REAL,INTENT(OUT)           :: gradUz(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< gradients in z-direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVarPrim):: gradUxi,gradUeta,gradUzeta ! gradients in xi/eta/zeta directions
INTEGER                    :: iElem,i,j,k,l
!==================================================================================================================================
! volume integral
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.0) THEN ! DG element
#endif
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    gradUxi     =             D_T(0,i)*UPrim(:,0,j,k,iElem)
    gradUeta    =             D_T(0,j)*UPrim(:,i,0,k,iElem)
#if PP_dim == 3
    gradUzeta   =             D_T(0,k)*UPrim(:,i,j,0,iElem)
#endif 
    DO l=1,PP_N
      gradUxi   = gradUxi   + D_T(l,i)*UPrim(:,l,j,k,iElem)
      gradUeta  = gradUeta  + D_T(l,j)*UPrim(:,i,l,k,iElem)
#if PP_dim == 3
      gradUzeta = gradUzeta + D_T(l,k)*UPrim(:,i,j,l,iElem)
#endif 
    END DO
    gradUx(:,i,j,k,iElem) = Metrics_fTilde(1,i,j,k,iElem,0)*gradUxi   &
                          + Metrics_gTilde(1,i,j,k,iElem,0)*gradUeta  &
#if PP_dim == 3
                          + Metrics_hTilde(1,i,j,k,iElem,0)*gradUzeta
#else
                          +0.
#endif 
    gradUy(:,i,j,k,iElem) = Metrics_fTilde(2,i,j,k,iElem,0)*gradUxi   &
                          + Metrics_gTilde(2,i,j,k,iElem,0)*gradUeta  &
#if PP_dim == 3
                          + Metrics_hTilde(2,i,j,k,iElem,0)*gradUzeta
#else
                          +0.
#endif 
#if PP_dim == 3
    gradUz(:,i,j,k,iElem) = Metrics_fTilde(3,i,j,k,iElem,0)*gradUxi   &
                          + Metrics_gTilde(3,i,j,k,iElem,0)*gradUeta  &
                          + Metrics_hTilde(3,i,j,k,iElem,0)*gradUzeta
#else
    gradUz(:,i,j,k,iElem) = 0.
#endif 
  END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
  ! in FV elements, the central gradients are transformed, i.e. they are multiplied by the metrics terms
  ELSE
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      gradUx(:,i,j,k,iElem) =                       FV_Metrics_fTilde_sJ(1,i,j,k,iElem)*gradUxi_central  (:,i,j,k,iElem)
      gradUx(:,i,j,k,iElem) = gradUx(:,i,j,k,iElem)+FV_Metrics_gTilde_sJ(1,i,j,k,iElem)*gradUeta_central (:,i,j,k,iElem)
#if PP_dim == 3
      gradUx(:,i,j,k,iElem) = gradUx(:,i,j,k,iElem)+FV_Metrics_hTilde_sJ(1,i,j,k,iElem)*gradUzeta_central(:,i,j,k,iElem)
#endif
      gradUy(:,i,j,k,iElem) =                       FV_Metrics_fTilde_sJ(2,i,j,k,iElem)*gradUxi_central  (:,i,j,k,iElem)
      gradUy(:,i,j,k,iElem) = gradUy(:,i,j,k,iElem)+FV_Metrics_gTilde_sJ(2,i,j,k,iElem)*gradUeta_central (:,i,j,k,iElem)
#if PP_dim == 3
      gradUy(:,i,j,k,iElem) = gradUy(:,i,j,k,iElem)+FV_Metrics_hTilde_sJ(2,i,j,k,iElem)*gradUzeta_central(:,i,j,k,iElem)
      gradUz(:,i,j,k,iElem) =                       FV_Metrics_fTilde_sJ(3,i,j,k,iElem)*gradUxi_central  (:,i,j,k,iElem)
      gradUz(:,i,j,k,iElem) = gradUz(:,i,j,k,iElem)+FV_Metrics_gTilde_sJ(3,i,j,k,iElem)*gradUeta_central (:,i,j,k,iElem)
      gradUz(:,i,j,k,iElem) = gradUz(:,i,j,k,iElem)+FV_Metrics_hTilde_sJ(3,i,j,k,iElem)*gradUzeta_central(:,i,j,k,iElem)
#endif
    END DO; END DO; END DO! i,j,k=0,PP_N
  END IF
#endif    
END DO ! iElem=1,nElems
END SUBROUTINE Lifting_VolInt_Nonconservative


!==================================================================================================================================
!> \brief Computes the volume integral of the BR2 scheme in conservative form for one direction at a time (x,y,z)
!>
!> Can be computed only in strong form, due to the formulation of the BR2 scheme
!> In the conservative form, the volume integral is calculated from the transformed solution, i.e. the solution is multiplied by teh
!> metrics terms.
!==================================================================================================================================
SUBROUTINE Lifting_VolInt_Conservative(dir,UPrim,gradU)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: D_T
USE MOD_Mesh_Vars    ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde   ! metrics
USE MOD_Mesh_Vars    ,ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
USE MOD_FV_Vars      ,ONLY: FV_Metrics_fTilde_sJ,FV_Metrics_gTilde_sJ,FV_Metrics_hTilde_sJ  ! metrics
USE MOD_FV_Vars      ,ONLY: gradUxi_central, gradUeta_central, gradUzeta_central
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                           :: dir                                          !< direction (x,y,z)
REAL,INTENT(IN)                              :: UPrim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution
REAL,INTENT(OUT)                             :: gradU(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution gradient in direction dir
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ) :: UE_f,UE_g,UE_h ! transformed gradient flux (i.e. transformed solution)
INTEGER                                          :: iElem,i,j,k,l
!==================================================================================================================================
! volume integral
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.0) THEN ! DG element
#endif
  ! transform the gradient "flux" into the reference element coordinates
  CALL Lifting_Metrics(dir,UPrim(:,:,:,:,iElem),&
                       Metrics_fTilde(:,:,:,:,iElem,0),&
                       Metrics_gTilde(:,:,:,:,iElem,0),&
                       Metrics_hTilde(:,:,:,:,iElem,0),&
                       UE_f,UE_g,UE_h)

  ! calculate the volume integral of the gradient "flux"
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    gradU(:,i,j,k,iElem)   =                      D_T(0,i)*UE_f(:,0,j,k)+&
                                                  D_T(0,j)*UE_g(:,i,0,k)+&
#if PP_dim == 3
                                                  D_T(0,k)*UE_h(:,i,j,0)
#else
                                                  0.
#endif 
    DO l=1,PP_N
      gradU(:,i,j,k,iElem) = gradU(:,i,j,k,iElem)+D_T(l,i)*UE_f(:,l,j,k)+&
                                                  D_T(l,j)*UE_g(:,i,l,k)+&
#if PP_dim == 3
                                                  D_T(l,k)*UE_h(:,i,j,l)
#else
                                                  0.
#endif 
    END DO ! l
  END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
  ! in FV elements, the central gradients are transformed, i.e. they are multiplied by the metrics terms
  ELSE
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      gradU(:,i,j,k,iElem) =                      FV_Metrics_fTilde_sJ(dir,i,j,k,iElem)*gradUxi_central  (:,i,j,k,iElem)
      gradU(:,i,j,k,iElem) = gradU(:,i,j,k,iElem)+FV_Metrics_gTilde_sJ(dir,i,j,k,iElem)*gradUeta_central (:,i,j,k,iElem)
#if PP_dim == 3
      gradU(:,i,j,k,iElem) = gradU(:,i,j,k,iElem)+FV_Metrics_hTilde_sJ(dir,i,j,k,iElem)*gradUzeta_central(:,i,j,k,iElem)
#endif 
    END DO; END DO; END DO! i,j,k=0,PP_N
  END IF
#endif    
END DO ! iElem=1,nElems
END SUBROUTINE Lifting_VolInt_Conservative


!==================================================================================================================================
!> \brief Compute the transformed gradient fluxes
!>
!> Transform the gradient terms by multiplying them with the metrics terms
!==================================================================================================================================
SUBROUTINE Lifting_Metrics(dir,UPrim,Mf,Mg,Mh,UPrim_f,UPrim_g,UPrim_h)
! MODULES
USE MOD_DG_Vars,ONLY:nDOFElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: dir                                 !< direction (x,y,z)
REAL,INTENT(IN)    :: Mf(3,nDOFElem)                      !< metrics in xi
REAL,INTENT(IN)    :: Mg(3,nDOFElem)                      !< metrics in eta
REAL,INTENT(IN)    :: Mh(3,nDOFElem)                      !< metrics in zeta
REAL,INTENT(IN)    :: UPrim(PP_nVarPrim,nDOFElem)         !< solution ("flux")
REAL,INTENT(OUT)   :: UPrim_f(PP_nVarPrim,nDOFElem)       !< gradient flux xi
REAL,INTENT(OUT)   :: UPrim_g(PP_nVarPrim,nDOFElem)       !< gradient flux eta
REAL,INTENT(OUT)   :: UPrim_h(PP_nVarPrim,nDOFElem)       !< gradient flux zeta
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i
!==================================================================================================================================
DO i=1,nDOFElem
  UPrim_f(:,i) = Mf(dir,i)*UPrim(:,i)
  UPrim_g(:,i) = Mg(dir,i)*UPrim(:,i)
#if PP_dim == 3
  UPrim_h(:,i) = Mh(dir,i)*UPrim(:,i)
#endif
END DO ! i
END SUBROUTINE Lifting_Metrics


END MODULE MOD_Lifting_VolInt
#endif /*PARABOLIC*/
