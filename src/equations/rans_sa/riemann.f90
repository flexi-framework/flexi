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
#include "eos.h"
!==================================================================================================================================
!> Contains routines to compute the riemann (Advection, Diffusion) for a given Face
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  SUBROUTINE RiemannInt(F_L,F_R,U_LL,U_RR,F)
    REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
    REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
    REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
  END SUBROUTINE
END INTERFACE

PROCEDURE(RiemannInt),POINTER :: Riemann_pointer    !< pointer defining the standard inner Riemann solver
PROCEDURE(RiemannInt),POINTER :: RiemannBC_pointer  !< pointer defining the standard BC    Riemann solver

INTEGER,PARAMETER      :: PRM_RIEMANN_SAME          = -1
INTEGER,PARAMETER      :: PRM_RIEMANN_LF            = 1

INTERFACE InitRiemann
  MODULE PROCEDURE InitRiemann
END INTERFACE

INTERFACE Riemann
  MODULE PROCEDURE Riemann
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


PUBLIC::InitRiemann
PUBLIC::Riemann
PUBLIC::FinalizeRiemann
!==================================================================================================================================

PUBLIC::DefineParametersRiemann
CONTAINS


!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Riemann")
CALL prms%CreateIntFromStringOption('Riemann',   "Riemann solver to be used: only LF for RANS-SA!", "lf")
CALL addStrListEntry('Riemann','lf',           PRM_RIEMANN_LF)
CALL prms%CreateIntFromStringOption('RiemannBC', "Riemann solver used for boundary conditions: Same, LF", "Same")
CALL addStrListEntry('RiemannBC','lf',           PRM_RIEMANN_LF)
CALL addStrListEntry('RiemannBC','same',         PRM_RIEMANN_SAME)
END SUBROUTINE DefineParametersRiemann

!==================================================================================================================================!
!> Initialize Riemann solver routines, read inner and BC Riemann solver parameters and set pointers
!==================================================================================================================================!
SUBROUTINE InitRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: Riemann
!==================================================================================================================================
Riemann = GETINTFROMSTR('Riemann')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_LF)
  Riemann_pointer => Riemann_LF
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'Riemann solver not defined!')
END SELECT

Riemann = GETINTFROMSTR('RiemannBC')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_SAME)
  RiemannBC_pointer => Riemann_pointer
CASE(PRM_RIEMANN_LF)
  RiemannBC_pointer => Riemann_LF
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'RiemannBC solver not defined!')
END SELECT

END SUBROUTINE InitRiemann

!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
SUBROUTINE Riemann(Nloc,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,doBC)
! MODULES
USE MOD_Flux         ,ONLY:EvalEulerFlux1D_fast
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                        :: Nloc       !< local polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L        !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R        !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L    !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_R    !< primitive solution at right side of the interface
!> normal vector and tangential vectors at side
REAL,DIMENSION(          3,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv,t1,t2
LOGICAL,INTENT(IN)                                        :: doBC       !< marker whether side is a BC side
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: FOut       !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F
REAL,DIMENSION(PP_2Var) :: U_LL,U_RR
PROCEDURE(RiemannInt),POINTER :: Riemann_loc !< pointer defining the standard inner Riemann solver
!==================================================================================================================================
IF (doBC) THEN
  Riemann_loc => RiemannBC_pointer
ELSE
  Riemann_loc => Riemann_pointer
END IF

! Momentum has to be rotatet using the normal system individual for each
DO j=0,ZDIM(Nloc); DO i=0,Nloc
  ! left state: U_L
  U_LL(DENS)=U_L(DENS,i,j)
  U_LL(SRHO)=1./U_LL(DENS)
  U_LL(ENER)=U_L(5,i,j)
  U_LL(PRES)=UPrim_L(5,i,j)
  U_LL(MUSA)=U_L(6,i,j)
  ! rotate velocity in normal and tangential direction 
  U_LL(VEL1)=DOT_PRODUCT(UPrim_L(2:4,i,j),nv(:,i,j))
  U_LL(VEL2)=DOT_PRODUCT(UPrim_L(2:4,i,j),t1(:,i,j))
  U_LL(MOM1)=U_LL(DENS)*U_LL(VEL1)
  U_LL(MOM2)=U_LL(DENS)*U_LL(VEL2)
#if PP_dim==3
  U_LL(VEL3)=DOT_PRODUCT(UPrim_L(2:4,i,j),t2(:,i,j))
  U_LL(MOM3)=U_LL(DENS)*U_LL(VEL3)
#else
  U_LL(VEL3)=0.
  U_LL(MOM3)=0.
#endif
  ! right state: U_R
  U_RR(DENS)=U_R(DENS,i,j)
  U_RR(SRHO)=1./U_RR(DENS)
  U_RR(ENER)=U_R(5,i,j)
  U_RR(PRES)=UPrim_R(5,i,j)
  U_RR(MUSA)=U_R(6,i,j)
  ! rotate momentum in normal and tangential direction 
  U_RR(VEL1)=DOT_PRODUCT(UPRIM_R(2:4,i,j),nv(:,i,j))
  U_RR(VEL2)=DOT_PRODUCT(UPRIM_R(2:4,i,j),t1(:,i,j))
  U_RR(MOM1)=U_RR(DENS)*U_RR(VEL1)
  U_RR(MOM2)=U_RR(DENS)*U_RR(VEL2)
#if PP_dim==3
  U_RR(VEL3)=DOT_PRODUCT(UPRIM_R(2:4,i,j),t2(:,i,j))
  U_RR(MOM3)=U_RR(DENS)*U_RR(VEL3)
#else
  U_RR(VEL3)=0.
  U_RR(MOM3)=0.
#endif

# ifndef SPLIT_DG
  CALL EvalEulerFlux1D_fast(U_LL,F_L)
  CALL EvalEulerFlux1D_fast(U_RR,F_R)
#endif /*SPLIT_DG*/

  CALL Riemann_loc(F_L,F_R,U_LL,U_RR,F)

  ! Back Rotate the normal flux into Cartesian direction
  Fout(DENS,i,j)=F(DENS)
  Fout(MOMV,i,j)=nv(:,i,j)*F(MOM1)     &
                  + t1(:,i,j)*F(MOM2)  &
#if PP_dim==3
                  + t2(:,i,j)*F(MOM3) 
#else
                  + 0.
#endif
  Fout(ENER,i,j)=F(ENER)
  Fout(MUSA,i,j)=F(MUSA)
END DO; END DO

END SUBROUTINE Riemann



#if PARABOLIC
!==================================================================================================================================
!> Computes the viscous RANS SA diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux(Nloc,F,UPrim_L,UPrim_R, &
                       gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv)
! MODULES
USE MOD_Flux,ONLY: EvalDiffFlux3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                         :: Nloc     !< local polynomial degree
                                                           !> solution in primitive variables at left/right side of the interface 
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)   :: UPrim_L,UPrim_R
                                                           !> solution gradients in x/y/z-direction left/right of the interface 
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)   :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,INTENT(IN)                                            :: nv(3,0:Nloc,0:ZDIM(Nloc)) !< normal vector
REAL,INTENT(OUT)                                           :: F(PP_nVar,0:Nloc,0:ZDIM(Nloc)) !< viscous flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                              :: p,q
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc))            :: diffFluxX_L,diffFluxY_L,diffFluxZ_L
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc))            :: diffFluxX_R,diffFluxY_R,diffFluxZ_R
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
  CALL EvalDiffFlux3D(Nloc,UPrim_L,gradUx_L,   gradUy_L,   gradUz_L, &
                                diffFluxX_L,diffFluxY_L,diffFluxZ_L)
  CALL EvalDiffFlux3D(Nloc,UPrim_R,gradUx_R,   gradUy_R,   gradUz_R, &
                                diffFluxX_R,diffFluxY_R,diffFluxZ_R)
! BR1 uses arithmetic mean of the fluxes
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  F(:,p,q)=0.5*(nv(1,p,q)*(diffFluxX_L(:,p,q)+diffFluxX_R(:,p,q)) &
               +nv(2,p,q)*(diffFluxY_L(:,p,q)+diffFluxY_R(:,p,q)) &
               +nv(3,p,q)*(diffFluxZ_L(:,p,q)+diffFluxZ_R(:,p,q)))
END DO; END DO
END SUBROUTINE ViscousFlux
#endif /* PARABOLIC */





!==================================================================================================================================
!> Local Lax-Friedrichs (Rusanov) Riemann solver
!==================================================================================================================================
SUBROUTINE Riemann_LF(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                                !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                                !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: LambdaMax
!==================================================================================================================================
! Lax-Friedrichs
LambdaMax = MAX( ABS(U_RR(VEL1)),ABS(U_LL(VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )
#ifndef SPLIT_DG
F = 0.5*((F_L+F_R) - LambdaMax*(U_RR(CONS) - U_LL(CONS)))
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! compute surface flux
F = F - 0.5*LambdaMax*(U_RR(CONS) - U_LL(CONS))
#endif /*SPLIT_DG*/

END SUBROUTINE Riemann_LF

!==================================================================================================================================
!> Finalize Riemann solver routines
!==================================================================================================================================
SUBROUTINE FinalizeRiemann()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeRiemann


END MODULE MOD_Riemann
