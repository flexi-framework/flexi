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
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

INTERFACE Riemann_Point
  MODULE PROCEDURE Riemann_Point
END INTERFACE

INTERFACE GetFlux
  MODULE PROCEDURE GetFlux
END INTERFACE

#if PARABOLIC
INTERFACE ViscousFlux
  MODULE PROCEDURE ViscousFlux
END INTERFACE
PUBLIC::ViscousFlux
#endif

INTERFACE FinalizeRiemann
  MODULE PROCEDURE FinalizeRiemann
END INTERFACE

PUBLIC::Riemann
PUBLIC::Riemann_Point
PUBLIC::GetFlux
PUBLIC::FinalizeRiemann
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
SUBROUTINE Riemann(Nloc,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,doBC)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                          :: Nloc       !< local polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L        !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R        !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L    !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_R    !< primitive solution at right side of the interface
!> normal vector and tangential vectors at side
REAL,DIMENSION(          3,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv,t1,t2
LOGICAL,INTENT(IN)                                          :: doBC       !< marker whether side is a BC side
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: FOut       !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j
!==================================================================================================================================
DO j=0,ZDIM(Nloc); DO i=0,Nloc
  CALL Riemann_Point(Fout(:,i,j),U_L(:,i,j),U_R(:,i,j),UPrim_L(:,i,j),UPrim_R(:,i,j),nv(:,i,j),t1(:,i,j),t2(:,i,j),doBC)
END DO; END DO
END SUBROUTINE Riemann

!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!==================================================================================================================================
SUBROUTINE Riemann_Point(F,U_L,U_R,dummy_L,dummy_R,nv,t1,t2,doBC)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_L        !< Left state
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_R        !< Right state
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: dummy_L    !< primitive state (useless here)
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: dummy_R    !< primitive state (useless here)
REAL,INTENT(IN)                         :: nv(3)      !< Normal vector
REAL,INTENT(IN)                         :: t1(3)      !< First tangential vector
REAL,INTENT(IN)                         :: t2(3)      !< Second tangential vector
LOGICAL,INTENT(IN)                      :: doBC       !< Switch to do BC sides or not
REAL,INTENT(OUT)                        :: F(PP_nVar) !< Flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: LambdaMax
!==================================================================================================================================
LambdaMax = AdvVel(1)*nv(1) +  AdvVel(2)*nv(2) + AdvVel(3)*nv(3)
! Compute the classic upwind flux into normal direction for each face GP
F(1) = 0.5*( (LambdaMax + ABS(LambdaMax))*U_L(1) + (LambdaMax-ABS(LambdaMax))*U_R(1))
END SUBROUTINE Riemann_Point


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
INTEGER,INTENT(IN)                                              :: Nloc                         !< Polynomial degree
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)        :: U_L                          !< Left state
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)        :: U_R                          !< Right state
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUx_L                     !< Left gradient in x-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUy_L                     !< Left gradient in y-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUz_L                     !< Left gradient in z-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUx_R                     !< Right gradient in x-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUy_R                     !< Right gradient in y-direction
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN),TARGET :: gradUz_R                     !< Right gradient in z-direction
REAL,INTENT(IN)                                                 :: nv(3,0:Nloc,0:ZDIM(Nloc))      !< Normal vector
REAL,INTENT(OUT)                                                :: F(PP_nVar,0:Nloc,0:ZDIM(Nloc)) !< Flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc)),TARGET ::  f_L, f_R, g_L, g_R, h_L, h_R
REAL,DIMENSION(:,:,:),POINTER                    :: Pf_L,Pf_R,Pg_L,Pg_R,Ph_L,Ph_R
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
CALL EvalDiffFlux2D(Nloc,f_L,g_L,h_L,U_L,gradUx_L,gradUy_L,gradUz_L)
CALL EvalDiffFlux2D(Nloc,f_R,g_R,h_R,U_R,gradUx_R,gradUy_R,gradUz_R)
Pf_L => f_L; Pg_L => g_L; Ph_L => h_L
Pf_R => f_R; Pg_R => g_R; Ph_R => h_R
F(1,:,:)=0.5*(nv(1,:,:)*(Pf_L(1,:,:)+Pf_R(1,:,:)) &
             +nv(2,:,:)*(Pg_L(1,:,:)+Pg_R(1,:,:)))
#if PP_dim==3
F(1,:,:)=F(1,:,:)+0.5*nv(3,:,:)*(Ph_L(1,:,:)+Ph_R(1,:,:))
#endif
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
