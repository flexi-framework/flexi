!=================================================================================================================================
! Copyright (c) 2010-2017 Prof. Claus-Dieter Munz
! Copyright (c) 2016-2017 Gregor Gassner (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Florian Hindenlang (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Andrew Winters (github.com/project-fluxo/fluxo)
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
!>\brief Computes the DGSEM volume integral
!> The volume integral is computed via the weak form of the DG method
!> Computes the volume integral contribution based on U and updates Ut
!> Volume integral is split into integral of advection and diffusion part
!==================================================================================================================================
#include "flexi.h"
MODULE MOD_VolInt
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE VolInt
#ifndef SPLIT_DG
  MODULE PROCEDURE VolInt_weakForm
#else
  MODULE PROCEDURE VolInt_splitForm
#endif /*SPLIT_DG*/
END INTERFACE

PUBLIC::VolInt
!==================================================================================================================================
CONTAINS


!==================================================================================================================================
!> Computes the advection and viscous part volume integral of the weak DG form according to Kopriva
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is overwritten with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolInt_weakForm(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: D_hat_T,nDOFElem,UPrim,U
USE MOD_Mesh_Vars    ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
USE MOD_Flux         ,ONLY: EvalFlux3D      ! computes volume fluxes in local coordinates
#if PARABOLIC
USE MOD_Flux         ,ONLY: EvalDiffFlux3D  ! computes volume fluxes in local coordinates
USE MOD_Lifting_Vars ,ONLY: gradUx,gradUy,gradUz
#endif
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< Time derivative of the volume integral (viscous part)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: f,g,h     !< Volume advective fluxes at GP
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: fv,gv,hv  !< Volume viscous fluxes at GP
#endif
!==================================================================================================================================
! Diffusive part
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.1) CYCLE ! FV Elem
#endif
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  CALL EvalFlux3D(PP_N,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem),f,g,h)
#if PARABOLIC
  CALL EvalDiffFlux3D( UPrim(:,:,:,:,iElem),&
                      gradUx(:,:,:,:,iElem),&
                      gradUy(:,:,:,:,iElem),&
                      gradUz(:,:,:,:,iElem),&
                      fv,gv,hv,iElem)

  f=f+fv
  g=g+gv
#if PP_dim==3
  h=h+hv
#endif
#endif

  CALL VolInt_Metrics(nDOFElem,f,g,h,Metrics_fTilde(:,:,:,:,iElem,0),&
                                     Metrics_gTilde(:,:,:,:,iElem,0),&
                                     Metrics_hTilde(:,:,:,:,iElem,0))

  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      ! Update the time derivative with the spatial derivatives of the transformed fluxes
      Ut(:,i,j,k,iElem) =                     D_Hat_T(0,i)*f(:,0,j,k) + &
#if PP_dim==3
                                              D_Hat_T(0,k)*h(:,i,j,0) + &
#endif
                                              D_Hat_T(0,j)*g(:,i,0,k)
    DO l=1,PP_N
      ! Update the time derivative with the spatial derivatives of the transformed fluxes
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_Hat_T(l,i)*f(:,l,j,k) + &
#if PP_dim==3
                                              D_Hat_T(l,k)*h(:,i,j,l) + &
#endif
                                              D_Hat_T(l,j)*g(:,i,l,k)
    END DO ! l
  END DO; END DO; END DO !i,j,k
END DO ! iElem
END SUBROUTINE VolInt_weakForm

#ifdef SPLIT_DG
!==================================================================================================================================
!> Computes the advection and viscous part volume integral in SplitDG formulation
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is overwritten with the volume flux derivatives
!> Attention 3: the factor of 2 in front of the derivative matrix entries is incorporated into the split fluxes!
!> For details on the derivation see Gassner, Gregor J., Andrew R. Winters, and David A. Kopriva.
!> "Split form nodal discontinuous Galerkin schemes with summation-by-parts property for the compressible Euler equations."
!> Journal of Computational Physics 327 (2016): 39-66.
!==================================================================================================================================
SUBROUTINE VolInt_splitForm(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: DVolSurf,nDOFElem,UPrim,U
USE MOD_Mesh_Vars    ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
USE MOD_Flux         ,ONLY: EvalFlux3D      ! computes volume fluxes in local coordinates
#if PARABOLIC
USE MOD_DG_Vars      ,ONLY: D_Hat_T
USE MOD_Flux         ,ONLY: EvalDiffFlux3D  ! computes volume fluxes in local coordinates
USE MOD_Lifting_Vars ,ONLY: gradUx,gradUy,gradUz
#endif
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
#endif
USE MOD_SplitFlux    ,ONLY:SplitDGVolume_pointer ! computes volume fluxes in split formulation
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< Time derivative of the volume integral (viscous part)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar                     )  :: Flux         !< temp variable for split flux
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: f_c,g_c,h_c  !< Euler fluxes at GP
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: fv,gv,hv     !< Parabolic fluxes at GP
#endif /*PARABOLIC*/
!==================================================================================================================================
! Diffusive part
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.1) CYCLE ! FV Elem
#endif
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  CALL EvalFlux3D(PP_N,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem),f_c,g_c,h_c)
  ! Add metric terms to fluxes
  CALL VolInt_Metrics(nDOFElem,f_c,g_c,h_c,Metrics_fTilde(:,:,:,:,iElem,0),&
                                           Metrics_gTilde(:,:,:,:,iElem,0),&
                                           Metrics_hTilde(:,:,:,:,iElem,0))
#if PARABOLIC
  CALL EvalDiffFlux3D( UPrim(:,:,:,:,iElem),&
                      gradUx(:,:,:,:,iElem),&
                      gradUy(:,:,:,:,iElem),&
                      gradUz(:,:,:,:,iElem),&
                      fv,gv,hv,iElem)

  CALL VolInt_Metrics(nDOFElem,fv,gv,hv,Metrics_fTilde(:,:,:,:,iElem,0),&
                                        Metrics_gTilde(:,:,:,:,iElem,0),&
                                        Metrics_hTilde(:,:,:,:,iElem,0))
#endif

  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    ! consistency: if both points for flux evaluation are the same the standard
    ! euler fluxes are retained.
#if PARABOLIC
    Ut(:,i,j,k,iElem) = DVolSurf(i,i)*f_c(:,i,j,k) + &
                        DVolSurf(j,j)*g_c(:,i,j,k) + &
#if PP_dim==3
                        DVolSurf(k,k)*h_c(:,i,j,k) + &
#endif /*PP_dim==3*/
                        D_Hat_T(i,i)*fv(:,i,j,k)   + &
                        D_Hat_T(j,j)*gv(:,i,j,k)   + &
#if PP_dim==3
                        D_Hat_T(k,k)*hv(:,i,j,k)
#else
                        0.
#endif /*PP_dim==3*/
#else
    Ut(:,i,j,k,iElem) = DVolSurf(i,i)*f_c(:,i,j,k) + &
#if PP_dim==3
                        DVolSurf(j,j)*h_c(:,i,j,k) + &
#endif /*PP_dim==3*/
                        DVolSurf(k,k)*g_c(:,i,j,k)
#endif /*PARABOLIC*/
  END DO; END DO; END DO !i,j,k


  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO l=i+1,PP_N
       ! compute split flux in x-direction
       CALL SplitDGVolume_pointer(U(:,i,j,k,iElem),UPrim(:,i,j,k,iElem), &
                                  U(:,l,j,k,iElem),UPrim(:,l,j,k,iElem), &
                                  Metrics_fTilde(:,i,j,k,iElem,0),Metrics_fTilde(:,l,j,k,iElem,0),Flux)
#if PARABOLIC
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,i)*Flux(:) + D_Hat_T(l,i)*fv(:,l,j,k)
       !symmetry
       Ut(:,l,j,k,iElem) = Ut(:,l,j,k,iElem) + DVolSurf(i,l)*Flux(:) + D_Hat_T(i,l)*fv(:,i,j,k)
#else
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,i)*Flux(:)
       !symmetry
       Ut(:,l,j,k,iElem) = Ut(:,l,j,k,iElem) + DVolSurf(i,l)*Flux(:)
#endif /*PARABOLIC*/
    END DO ! m

    DO l=j+1,PP_N
       ! compute split flux in y-direction
       CALL SplitDGVolume_pointer(U(:,i,j,k,iElem),UPrim(:,i,j,k,iElem), &
                                  U(:,i,l,k,iElem),UPrim(:,i,l,k,iElem), &
                                  Metrics_gTilde(:,i,j,k,iElem,0),Metrics_gTilde(:,i,l,k,iElem,0),Flux)
#if PARABOLIC
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,j)*Flux(:) + D_Hat_T(l,j)*gv(:,i,l,k)
       !symmetry
       Ut(:,i,l,k,iElem) = Ut(:,i,l,k,iElem) + DVolSurf(j,l)*Flux(:) + D_Hat_T(j,l)*gv(:,i,j,k)
#else
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,j)*Flux(:)
       !symmetry
       Ut(:,i,l,k,iElem) = Ut(:,i,l,k,iElem) + DVolSurf(j,l)*Flux(:)
#endif /*PARABOLIC*/
    END DO ! m

#if PP_dim==3
    DO l=k+1,PP_N
       ! compute split flux in z-direction
       CALL SplitDGVolume_pointer(U(:,i,j,k,iElem),UPrim(:,i,j,k,iElem), &
                                  U(:,i,j,l,iElem),UPrim(:,i,j,l,iElem), &
                                  Metrics_hTilde(:,i,j,k,iElem,0),Metrics_hTilde(:,i,j,l,iElem,0),Flux)
#if PARABOLIC
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,k)*Flux(:) + D_Hat_T(l,k)*hv(:,i,j,l)
       !symmetry
       Ut(:,i,j,l,iElem) = Ut(:,i,j,l,iElem) + DVolSurf(k,l)*Flux(:) + D_Hat_T(k,l)*hv(:,i,j,k)
#else
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,k)*Flux(:)
       !symmetry
       Ut(:,i,j,l,iElem) = Ut(:,i,j,l,iElem) + DVolSurf(k,l)*Flux(:)
#endif /*PARABOLIC*/
    END DO ! l
#endif /*PP_dim==3*/

  END DO; END DO; END DO !i,j,k
END DO ! iElem
END SUBROUTINE VolInt_splitForm
#endif /*SPLIT_DG*/

!==================================================================================================================================
!> Compute the tranformed states for all conservative variables using the metric terms
!==================================================================================================================================
PPURE SUBROUTINE VolInt_Metrics(nDOFs,f,g,h,Mf,Mg,Mh)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                          :: nDOFs    !< Number of DOFs per element
                                                        !> Metrics terms
REAL,DIMENSION(3,nDOFs),INTENT(IN)          :: Mf,Mg,Mh
                                                        !> Volume fluxes at GP to be transformed from physical to reference space
REAL,DIMENSION(PP_nVar,nDOFs),INTENT(INOUT) :: f,g,h
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                     :: i
REAL,DIMENSION(PP_nVar)                     :: fTilde,gTilde !< Auxiliary variables needed to store the fluxes at one GP
#if PP_dim==3
REAL,DIMENSION(PP_nVar)                     :: hTilde !< Auxiliary variables needed to store the fluxes at one GP
#endif
!==================================================================================================================================
DO i=1,nDOFs
  fTilde=f(:,i)
  gTilde=g(:,i)
#if PP_dim==3
  hTilde=h(:,i)

  ! Compute the transformed fluxes with the metric terms
  ! Attention 1: we store the transformed fluxes in f,g,h again
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i) + &
           hTilde*Mf(3,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i) + &
           hTilde*Mg(3,i)
  h(:,i) = fTilde*Mh(1,i) + &
           gTilde*Mh(2,i) + &
           hTilde*Mh(3,i)
#else
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i)
#endif
END DO ! i
END SUBROUTINE VolInt_Metrics



END MODULE MOD_VolInt
