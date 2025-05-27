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
#include "eos.h"

!==================================================================================================================================
!> Contains routines to compute the riemann (Advection, Diffusion) for a given Face
!==================================================================================================================================
#include "flexi.h"
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: Riemann
PUBLIC:: GetFlux
#if PARABOLIC
PUBLIC:: ViscousFlux
#endif
PUBLIC:: FinalizeRiemann
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes sum of numerical and viscous flux
!==================================================================================================================================
SUBROUTINE GetFlux(Nloc,F,U_L,U_R, &
#if PARABOLIC
                   gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R, &
#endif /* PARABOLIC */
                   nv,t1,t2,doBC)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                       :: Nloc                         !< Polynomial degree
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)     :: U_L                          !< Left state
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)     :: U_R                          !< Right state
#if PARABOLIC
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN) :: gradUx_L                     !< Left gradient in x-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN) :: gradUy_L                     !< Left gradient in y-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN) :: gradUz_L                     !< Left gradient in z-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN) :: gradUx_R                     !< Right gradient in x-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN) :: gradUy_R                     !< Right gradient in y-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN) :: gradUz_R                     !< Right gradient in z-direction
#endif
REAL,INTENT(IN)                                          :: nv(3,0:Nloc,0:ZDIM(Nloc))      !< Normal vector
REAL,INTENT(IN)                                          :: t1(3,0:Nloc,0:ZDIM(Nloc))      !< First tangential vector
REAL,INTENT(IN)                                          :: t2(3,0:Nloc,0:ZDIM(Nloc))      !< Second tangential vector
LOGICAL,INTENT(IN)                                       :: doBC                         !< Switch to do BC sides or not
REAL,INTENT(OUT)                                         :: F(PP_nVar,0:Nloc,0:ZDIM(Nloc)) !< Flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: Fv(PP_nVar,0:Nloc,0:ZDIM(Nloc))
!==================================================================================================================================
CALL Riemann(Nloc,F,U_L,U_R,U_L,U_R,nv,t1,t2,doBC=doBC)
#if PARABOLIC
CALL ViscousFlux(Nloc,Fv,U_L,U_R,gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv)
F=F+Fv
#endif /*PARABOLIC*/

END SUBROUTINE GetFlux


!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!==================================================================================================================================
SUBROUTINE Riemann(Nloc,F,U_L,U_R,dummy_L,dummy_R,nv,t1,t2,doBC)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                          :: Nloc                           !< Polynomial degree
REAL,DIMENSION(PP_nVar,    0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L                            !< Left state
REAL,DIMENSION(PP_nVar,    0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R                            !< Right state
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: dummy_L                        !< primitive state (useless here)
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: dummy_R                        !< primitive state (useless here)
REAL,INTENT(IN)                                             :: nv(3,0:Nloc,0:ZDIM(Nloc))      !< Normal vector
REAL,INTENT(IN)                                             :: t1(3,0:Nloc,0:ZDIM(Nloc))      !< First tangential vector
REAL,INTENT(IN)                                             :: t2(3,0:Nloc,0:ZDIM(Nloc))      !< Second tangential vector
LOGICAL,INTENT(IN)                                          :: doBC                           !< Switch to do BC sides or not
REAL,INTENT(OUT)                                            :: F(PP_nVar,0:Nloc,0:ZDIM(Nloc)) !< Flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F_loc
REAL,DIMENSION(PP_nVar) :: G_L,G_R,H_loc
REAL,DIMENSION(PP_nVar) :: H_L,H_R,G_loc
REAL,DIMENSION(PP_nVar) :: U_LL,U_RR
!==================================================================================================================================
! The 2D/3D Burgers equations are NOT rotational invariant!
DO j=0,ZDIM(Nloc); DO i=0,Nloc
  ! Take the sign of the normal vector of the side into account
  U_LL(1) = U_L(1,i,j)*nv(1,i,j)
  U_LL(2) = U_L(2,i,j)*nv(2,i,j)
  U_LL(3) = U_L(3,i,j)*nv(3,i,j)

  U_RR(1) = U_R(1,i,j)*nv(1,i,j)
  U_RR(2) = U_R(2,i,j)*nv(2,i,j)
  U_RR(3) = U_R(3,i,j)*nv(3,i,j)

  ! Solve the Riemann problem in x-direction
  ! Evaluate the fluxes on the left and the right
  F_L(1) = 0.5*U_L(1,i,j)*U_L(1,i,j)
  F_L(2) =     U_L(1,i,j)*U_L(2,i,j)
  F_L(3) =     U_L(1,i,j)*U_L(3,i,j)

  F_R(1) = 0.5*U_R(1,i,j)*U_R(1,i,j)
  F_R(2) =     U_R(1,i,j)*U_R(2,i,j)
  F_R(3) =     U_R(1,i,j)*U_R(3,i,j)

  ! Solve the Riemann problem
  ! Shock
  IF(U_LL(1).GT.U_RR(1))THEN
    IF(0.5*(U_LL(1)+U_RR(1)).GT.0.)THEN
      F_loc = F_L
    ELSE
      F_loc = F_R
    END IF
  ! Rarefaction
  ELSE
    IF(U_LL(1)*U_RR(1).LT.0.)THEN
      F_loc = 0.
    ELSE
      IF(U_LL(1).GT.0.)THEN
        F_loc = F_L
      ELSE
        F_loc = F_R
      END IF
    END IF
  END IF

  ! Solve the Riemann problem in y-direction
  ! Evaluate the fluxes on the left and the right
  G_L(1) =     U_L(2,i,j)*U_L(1,i,j)
  G_L(2) = 0.5*U_L(2,i,j)*U_L(2,i,j)
  G_L(3) =     U_L(2,i,j)*U_L(3,i,j)

  G_R(1) =     U_R(2,i,j)*U_R(1,i,j)
  G_R(2) = 0.5*U_R(2,i,j)*U_R(2,i,j)
  G_R(3) =     U_R(2,i,j)*U_R(3,i,j)

  ! Solve the Riemann problem
  ! Shock
  IF(U_LL(2).GT.U_RR(2))THEN
    IF(0.5*(U_LL(2)+U_RR(2)).GT.0.)THEN
      G_loc = G_L
    ELSE
      G_loc = G_R
    END IF
  ! Rarefaction
  ELSE
    IF(U_LL(2)*U_RR(2).LT.0.)THEN
      G_loc = 0.
    ELSE
      IF(U_LL(2).GT.0.)THEN
        G_loc = G_L
      ELSE
        G_loc = G_R
      END IF
    END IF
  END IF

#if PP_dim==3
  ! Solve the Riemann problem in y-direction
  ! Evaluate the fluxes on the left and the right
  H_L(1) =     U_L(3,i,j)*U_L(1,i,j)
  H_L(2) =     U_L(3,i,j)*U_L(2,i,j)
  H_L(3) = 0.5*U_L(3,i,j)*U_L(3,i,j)

  H_R(1) =     U_R(3,i,j)*U_R(1,i,j)
  H_R(2) =     U_R(3,i,j)*U_R(2,i,j)
  H_R(3) = 0.5*U_R(3,i,j)*U_R(3,i,j)

  ! Solve the Riemann problem
  ! Shock
  IF(U_LL(3).GT.U_RR(3))THEN
    IF(0.5*(U_LL(3)+U_RR(3)).GT.0.)THEN
      H_loc = H_L
    ELSE
      H_loc = H_R
    END IF
  ! Rarefaction
  ELSE
    IF(U_LL(3)*U_RR(3).LT.0.)THEN
      H_loc = 0.
    ELSE
      IF(U_LL(3).GT.0.)THEN
        H_loc = H_L
      ELSE
        H_loc = H_R
      END IF
    END IF
  END IF
#endif

  ! Back Rotate the normal flux into Cartesian direction
  F(1:3,i,j) =   nv(1,i,j)*F_loc(:)  &
               + nv(2,i,j)*G_loc(:)  &
#if PP_dim==3
               + nv(3,i,j)*H_loc(:)
#else
               + 0.
#endif

END DO; END DO
END SUBROUTINE Riemann


#if PARABOLIC
!==================================================================================================================================
!> Computes the viscous diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux(Nloc,F,U_L,U_R, &
                       gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv)
! MODULES
USE MOD_Flux, ONLY:EvalDiffFlux2D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                                   :: Nloc                           !< Polynomial degree
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)           :: U_L                            !< Left state
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)           :: U_R                            !< Right state
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUx_L                       !< Left gradient in x-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUy_L                       !< Left gradient in y-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUz_L                       !< Left gradient in z-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUx_R                       !< Right gradient in x-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUy_R                       !< Right gradient in y-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUz_R                       !< Right gradient in z-direction
REAL,INTENT(IN)                                                      :: nv(3,0:Nloc,0:ZDIM(Nloc))      !< Normal vector
REAL,INTENT(OUT)                                                     :: F(PP_nVar,0:Nloc,0:ZDIM(Nloc)) !< Flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc)) ::  f_L, f_R, g_L, g_R, h_L, h_R
INTEGER                                     ::  p,q
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute viscous Burgers Diffusion flux
CALL EvalDiffFlux2D(Nloc,f_L,g_L,h_L,U_L,gradUx_L,gradUy_L,gradUz_L)
CALL EvalDiffFlux2D(Nloc,f_R,g_R,h_R,U_R,gradUx_R,gradUy_R,gradUz_R)

! BR1 uses arithmetic mean of the fluxes
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  F(1:3,p,q) = 0.5*(nv(1,p,q)*(f_L(1:3,p,q)+f_R(1:3,p,q)) &
                   +nv(2,p,q)*(g_L(1:3,p,q)+g_R(1:3,p,q)))
#if PP_dim==3
  F(1:3,p,q) = F(1:3,p,q)+0.5*(nv(3,p,q)*(h_L(1:3,p,q)+h_R(1:3,p,q)))
#endif
END DO; END DO
END SUBROUTINE ViscousFlux
#endif /* PARABOLIC */


!==================================================================================================================================
!> Finalize Riemann solver routines
!==================================================================================================================================
SUBROUTINE FinalizeRiemann()
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeRiemann

END MODULE MOD_Riemann
