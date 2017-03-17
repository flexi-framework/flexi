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
!>\brief Computes the DGSEM volume integral
!> The volume integral is computed via the weak form of the DG method 
!> Computes the volume integral contribution based on U and updates Ut
!> Volume integral is split into integral of advection and diffusion part
!==================================================================================================================================
MODULE MOD_VolInt
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE VolInt
  MODULE PROCEDURE VolInt_weakForm
END INTERFACE

INTERFACE VolIntAdv
  MODULE PROCEDURE VolIntAdv_weakForm
END INTERFACE

#if PARABOLIC
INTERFACE VolIntVisc
  MODULE PROCEDURE VolIntVisc_weakForm
END INTERFACE
#endif


PUBLIC::VolInt
PUBLIC::VolIntAdv
#if PARABOLIC
PUBLIC::VolIntVisc
#endif
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
    DO l=0,PP_N
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

!==================================================================================================================================
!> Computes the advection part volume integral of the weak DG form according to Kopriva
!> Polynomial degree is either N or NOver (overintegration)
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is overwritten with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolIntAdv_weakForm(Nloc,nDOFElem,D_Hat_T,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,U,UPrim,Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars ,ONLY:nElems
USE MOD_Flux      ,ONLY:EvalFlux3D ! computes volume fluxes in local coordinates
#if FV_ENABLED
USE MOD_FV_Vars   ,ONLY:FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc       !< Polynomial degree either N or NOver in case of overintegration
INTEGER,INTENT(IN) :: nDOFElem   !< Number of DOFs per element
REAL,INTENT(IN)    :: D_Hat_T(0:Nloc,0:Nloc)!< Transpose of differentiation matrix premultiplied by 
                                            !< mass matrix, size [0..Nloc,0..Nloc].
                                                                !> Metric terms in \f$ \xi\ / \eta / \zeta \f$-direction
REAL,INTENT(IN)    :: Metrics_fTilde(1:3,0:Nloc,0:Nloc,0:PP_NlocZ,1:nElems)
REAL,INTENT(IN)    :: Metrics_gTilde(1:3,0:Nloc,0:Nloc,0:PP_NlocZ,1:nElems)
REAL,INTENT(IN)    :: Metrics_hTilde(1:3,0:Nloc,0:Nloc,0:PP_NlocZ,1:nElems)
REAL,INTENT(IN)    :: U    (PP_nVar    ,0:Nloc,0:Nloc,0:PP_NlocZ,1:nElems)  !< conservative solution vector 
REAL,INTENT(IN)    :: UPrim(PP_nVarPrim,0:Nloc,0:Nloc,0:PP_NlocZ,1:nElems)  !< primitive solution vector 
REAL,INTENT(OUT)   :: Ut(PP_nVar,0:Nloc,0:Nloc,0:PP_NlocZ,1:nElems) !< Time derivative of the volume integral (advection part)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc,0:PP_NlocZ)  :: f,g,h !< Advective volume fluxes at GP
!==================================================================================================================================
! Advective part
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.1) THEN
    Ut(:,:,:,:,iElem) = 0.
    CYCLE ! FV Elem
  END IF
#endif
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  CALL EvalFlux3D(Nloc,U(:,:,:,:,iElem),UPrim(:,:,:,:,iElem),f,g,h)
  CALL VolInt_Metrics(nDOFElem,f,g,h,Metrics_fTilde(:,:,:,:,iElem),&
                                     Metrics_gTilde(:,:,:,:,iElem),&
                                     Metrics_hTilde(:,:,:,:,iElem))
  DO k=0,PP_NlocZ
    DO j=0,Nloc
      DO i=0,Nloc
        Ut(:,i,j,k,iElem) = D_Hat_T(0,i)*f(:,0,j,k) + &
                            D_Hat_T(0,j)*g(:,i,0,k) + &
                            D_Hat_T(0,k)*h(:,i,j,0)
        DO l=1,Nloc
          ! Update the time derivative with the spatial derivatives of the transformed fluxes
          Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_Hat_T(l,i)*f(:,l,j,k) + &
#if PP_dim==3
                                                  D_Hat_T(l,k)*h(:,i,j,l) + &
#endif
                                                  D_Hat_T(l,j)*g(:,i,l,k) 
        END DO ! l
      END DO !i
    END DO ! j
  END DO ! k
END DO ! iElem
END SUBROUTINE VolIntAdv_weakForm



#if PARABOLIC
!==================================================================================================================================
!> Computes the viscous part volume integral of the weak DG form according to Kopriva
!> Polynomial degree is either N or NOver (overintegration)
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is overwritten with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolIntVisc_weakForm(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: D_hat_T,nDOFElem,UPrim
USE MOD_Mesh_Vars    ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Mesh_Vars    ,ONLY: nElems
USE MOD_Flux         ,ONLY: EvalDiffFlux3D  ! computes volume fluxes in local coordinates
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
#endif
USE MOD_Lifting_Vars ,ONLY: gradUx,gradUy,gradUz
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< Time derivative of the volume integral (viscous part)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: f,g,h           !< Viscous volume fluxes at GP
!==================================================================================================================================
! Diffusive part
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.1) CYCLE ! FV Elem
#endif
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  CALL EvalDiffFlux3D( UPrim(:,:,:,:,iElem),&
                      gradUx(:,:,:,:,iElem),&
                      gradUy(:,:,:,:,iElem),&
                      gradUz(:,:,:,:,iElem),&
                      f,g,h,iElem)
  CALL VolInt_Metrics(nDOFElem,f,g,h,Metrics_fTilde(:,:,:,:,iElem,0),&
                                     Metrics_gTilde(:,:,:,:,iElem,0),&
                                     Metrics_hTilde(:,:,:,:,iElem,0))

  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      ! Update the time derivative with the spatial derivatives of the transformed fluxes
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_Hat_T(l,i)*f(:,l,j,k) + &
#if PP_dim==3
                                              D_Hat_T(l,k)*h(:,i,j,l) +&
#endif                                                  
                                              D_Hat_T(l,j)*g(:,i,l,k) 
    END DO ! l
  END DO; END DO; END DO !i,j,k
END DO ! iElem
END SUBROUTINE VolIntVisc_weakForm
#endif /* PARABOLIC */



!==================================================================================================================================
!> Compute the tranformed states for all conservative variables using the metric terms
!==================================================================================================================================
SUBROUTINE VolInt_Metrics(nDOFs,f,g,h,Mf,Mg,Mh)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
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
