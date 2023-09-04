!=================================================================================================================================
! Copyright (c) 2010-2024 Prof. Claus-Dieter Munz
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
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEENTROPYFIX = 33

INTERFACE InitRiemann
  MODULE PROCEDURE InitRiemann
END INTERFACE

INTERFACE Riemann
  MODULE PROCEDURE Riemann_Point
  MODULE PROCEDURE Riemann_Side
END INTERFACE

#if PARABOLIC
INTERFACE ViscousFlux
  MODULE PROCEDURE ViscousFlux_Point
  MODULE PROCEDURE ViscousFlux_Side
END INTERFACE
#endif

INTERFACE FinalizeRiemann
  MODULE PROCEDURE FinalizeRiemann
END INTERFACE

PUBLIC::DefineParametersRiemann
PUBLIC::InitRiemann
PUBLIC::Riemann
PUBLIC::FinalizeRiemann
#if PARABOLIC
PUBLIC::ViscousFlux
#endif
!==================================================================================================================================

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
CALL addStrListEntry('Riemann','roeentropyfix',PRM_RIEMANN_ROEENTROPYFIX)
CALL prms%CreateIntFromStringOption('RiemannBC', "Riemann solver used for boundary conditions: Same, LF", "Same")
CALL addStrListEntry('RiemannBC','lf',           PRM_RIEMANN_LF)
CALL addStrListEntry('RiemannBC','roeentropyfix',PRM_RIEMANN_ROEENTROPYFIX)
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
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  Riemann_pointer => Riemann_RoeEntropyFix
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
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  RiemannBC_pointer => Riemann_RoeEntropyFix
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
SUBROUTINE Riemann_Side(Nloc,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,doBC)
! MODULES
USE MOD_Flux         ,ONLY:EvalEulerFlux1D_fast
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                          :: Nloc       !< local polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L        !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R        !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L    !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_R    !< primitive solution at right side of the interface
REAL,DIMENSION(          3,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv,t1,t2   !> normal vector and tangential vectors at side
LOGICAL,INTENT(IN)                                          :: doBC       !< marker whether side is a BC side
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

! Momentum has to be rotated using the normal system individual for each
DO j=0,ZDIM(Nloc); DO i=0,Nloc
  ! left state: U_L
  U_LL(EXT_DENS)=U_L(DENS,i,j)
  U_LL(EXT_SRHO)=1./U_LL(EXT_DENS)
  U_LL(EXT_ENER)=U_L(ENER,i,j)
  U_LL(EXT_PRES)=UPrim_L(PRES,i,j)
  U_LL(EXT_MUSA)=U_L(MUSA,i,j)
  ! rotate velocity in normal and tangential direction
  U_LL(EXT_VEL1)=DOT_PRODUCT(UPrim_L(VELV,i,j),nv(:,i,j))
  U_LL(EXT_VEL2)=DOT_PRODUCT(UPrim_L(VELV,i,j),t1(:,i,j))
  U_LL(EXT_MOM1)=U_LL(EXT_DENS)*U_LL(EXT_VEL1)
  U_LL(EXT_MOM2)=U_LL(EXT_DENS)*U_LL(EXT_VEL2)
#if PP_dim==3
  U_LL(EXT_VEL3)=DOT_PRODUCT(UPrim_L(VELV,i,j),t2(:,i,j))
  U_LL(EXT_MOM3)=U_LL(EXT_DENS)*U_LL(EXT_VEL3)
#else
  U_LL(EXT_VEL3)=0.
  U_LL(EXT_MOM3)=0.
#endif
  ! right state: U_R
  U_RR(EXT_DENS)=U_R(DENS,i,j)
  U_RR(EXT_SRHO)=1./U_RR(EXT_DENS)
  U_RR(EXT_ENER)=U_R(ENER,i,j)
  U_RR(EXT_PRES)=UPrim_R(PRES,i,j)
  U_RR(EXT_MUSA)=U_R(MUSA,i,j)
  ! rotate momentum in normal and tangential direction
  U_RR(EXT_VEL1)=DOT_PRODUCT(UPRIM_R(VELV,i,j),nv(:,i,j))
  U_RR(EXT_VEL2)=DOT_PRODUCT(UPRIM_R(VELV,i,j),t1(:,i,j))
  U_RR(EXT_MOM1)=U_RR(EXT_DENS)*U_RR(EXT_VEL1)
  U_RR(EXT_MOM2)=U_RR(EXT_DENS)*U_RR(EXT_VEL2)
#if PP_dim==3
  U_RR(EXT_VEL3)=DOT_PRODUCT(UPRIM_R(VELV,i,j),t2(:,i,j))
  U_RR(EXT_MOM3)=U_RR(EXT_DENS)*U_RR(EXT_VEL3)
#else
  U_RR(EXT_VEL3)=0.
  U_RR(EXT_MOM3)=0.
#endif

  CALL EvalEulerFlux1D_fast(U_LL,F_L)
  CALL EvalEulerFlux1D_fast(U_RR,F_R)

  CALL Riemann_loc(F_L,F_R,U_LL,U_RR,F)

  ! Back Rotate the normal flux into Cartesian direction
  Fout(DENS,i,j)=F(DENS)
  Fout(MOMV,i,j)=nv(:,i,j)*F(MOM1)  &
               + t1(:,i,j)*F(MOM2)  &
#if PP_dim==3
               + t2(:,i,j)*F(MOM3)
#else
               + 0.
#endif
  Fout(ENER,i,j)=F(ENER)
  Fout(MUSA,i,j)=F(MUSA)
END DO; END DO

END SUBROUTINE Riemann_Side


!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
SUBROUTINE Riemann_Point(FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,doBC)
! MODULES
USE MOD_Flux         ,ONLY:EvalEulerFlux1D_fast
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_L        !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_R        !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_L    !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_R    !< primitive solution at right side of the interface
REAL,DIMENSION(3          ),INTENT(IN)  :: nv,t1,t2   !< normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: FOut       !< advective flux
LOGICAL,INTENT(IN)                      :: doBC       !< marker whether side is a BC side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
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
! left state: U_L
U_LL(EXT_DENS)=U_L(DENS)
U_LL(EXT_SRHO)=1./U_LL(EXT_DENS)
U_LL(EXT_ENER)=U_L(ENER)
U_LL(EXT_PRES)=UPrim_L(PRES)
U_LL(EXT_MUSA)=U_L(MUSA)
! rotate velocity in normal and tangential direction
U_LL(EXT_VEL1)=DOT_PRODUCT(UPrim_L(VELV),nv(:))
U_LL(EXT_VEL2)=DOT_PRODUCT(UPrim_L(VELV),t1(:))
U_LL(EXT_MOM1)=U_LL(EXT_DENS)*U_LL(EXT_VEL1)
U_LL(EXT_MOM2)=U_LL(EXT_DENS)*U_LL(EXT_VEL2)
#if PP_dim==3
U_LL(EXT_VEL3)=DOT_PRODUCT(UPrim_L(VELV),t2(:))
U_LL(EXT_MOM3)=U_LL(EXT_DENS)*U_LL(EXT_VEL3)
#else
U_LL(EXT_VEL3)=0.
U_LL(EXT_MOM3)=0.
#endif
! right state: U_R
U_RR(EXT_DENS)=U_R(DENS)
U_RR(EXT_SRHO)=1./U_RR(EXT_DENS)
U_RR(EXT_ENER)=U_R(ENER)
U_RR(EXT_PRES)=UPrim_R(PRES)
U_RR(EXT_MUSA)=U_R(MUSA)
! rotate momentum in normal and tangential direction
U_RR(EXT_VEL1)=DOT_PRODUCT(UPRIM_R(VELV),nv(:))
U_RR(EXT_VEL2)=DOT_PRODUCT(UPRIM_R(VELV),t1(:))
U_RR(EXT_MOM1)=U_RR(EXT_DENS)*U_RR(EXT_VEL1)
U_RR(EXT_MOM2)=U_RR(EXT_DENS)*U_RR(EXT_VEL2)
#if PP_dim==3
U_RR(EXT_VEL3)=DOT_PRODUCT(UPRIM_R(VELV),t2(:))
U_RR(EXT_MOM3)=U_RR(EXT_DENS)*U_RR(EXT_VEL3)
#else
U_RR(EXT_VEL3)=0.
U_RR(EXT_MOM3)=0.
#endif

CALL EvalEulerFlux1D_fast(U_LL,F_L)
CALL EvalEulerFlux1D_fast(U_RR,F_R)

CALL Riemann_loc(F_L,F_R,U_LL,U_RR,F)

! Back Rotate the normal flux into Cartesian direction
Fout(DENS)=F(DENS)
Fout(MOMV)=nv(:)*F(MOM1)     &
                + t1(:)*F(MOM2)  &
#if PP_dim==3
                + t2(:)*F(MOM3)
#else
                + 0.
#endif
Fout(ENER)=F(ENER)
Fout(MUSA)=F(MUSA)

END SUBROUTINE Riemann_Point


#if PARABOLIC
!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux_Side(Nloc,F,UPrim_L,UPrim_R, &
                            gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv)
! MODULES
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
USE MOD_Lifting_Vars ,ONLY: diffFluxX_L,diffFluxY_L,diffFluxZ_L
USE MOD_Lifting_Vars ,ONLY: diffFluxX_R,diffFluxY_R,diffFluxZ_R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                             :: Nloc     !< local polynomial degree
                                                               !> solution in primitive variables at left/right side of interface
REAL,DIMENSION(PP_nVarPrim   ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L,UPrim_R
                                                               !> solution gradients in x/y/z-direction left/right of interface
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DIMENSION(3             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv  !< normal vector
REAL,DIMENSION(PP_nVar       ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: F   !< viscous flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                                       :: p,q
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
CALL EvalDiffFlux3D(Nloc,UPrim_L,   gradUx_L,   gradUy_L,   gradUz_L  &
                                ,diffFluxX_L,diffFluxY_L,diffFluxZ_L  )
CALL EvalDiffFlux3D(Nloc,UPrim_R,   gradUx_R,   gradUy_R,   gradUz_R  &
                                ,diffFluxX_R,diffFluxY_R,diffFluxZ_R  )
! Arithmetic mean of the fluxes
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  F(:,p,q)=0.5*(nv(1,p,q)*(diffFluxX_L(:,p,q)+diffFluxX_R(:,p,q)) &
               +nv(2,p,q)*(diffFluxY_L(:,p,q)+diffFluxY_R(:,p,q)) &
               +nv(3,p,q)*(diffFluxZ_L(:,p,q)+diffFluxZ_R(:,p,q)))
END DO; END DO
END SUBROUTINE ViscousFlux_Side

!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux_Point(F,UPrim_L,UPrim_R, &
                             gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv)
! MODULES
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                           !> solution in primitive variables at left/right side of the interface
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim_L,UPrim_R
                                           !> solution gradients in x/y/z-direction left/right of the interface
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DIMENSION(3             ),INTENT(IN)  :: nv  !< normal vector
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: F   !< viscous flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar)  :: diffFluxX_L,diffFluxY_L,diffFluxZ_L
REAL,DIMENSION(PP_nVar)  :: diffFluxX_R,diffFluxY_R,diffFluxZ_R
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
CALL EvalDiffFlux3D(UPrim_L,   gradUx_L,   gradUy_L,   gradUz_L  &
                           ,diffFluxX_L,diffFluxY_L,diffFluxZ_L  )
CALL EvalDiffFlux3D(UPrim_R,   gradUx_R,   gradUy_R,   gradUz_R  &
                           ,diffFluxX_R,diffFluxY_R,diffFluxZ_R  )
! Arithmetic mean of the fluxes
F(:)=0.5*(nv(1)*(diffFluxX_L(:)+diffFluxX_R(:)) &
         +nv(2)*(diffFluxY_L(:)+diffFluxY_R(:)) &
         +nv(3)*(diffFluxZ_L(:)+diffFluxZ_R(:)))
END SUBROUTINE ViscousFlux_Point
#endif /* PARABOLIC */


!==================================================================================================================================
!> Local Lax-Friedrichs (Rusanov) Riemann solver
!==================================================================================================================================
SUBROUTINE Riemann_LF(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F         !< resulting Riemann flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: LambdaMax
!==================================================================================================================================
! Lax-Friedrichs
LambdaMax = MAX( ABS(U_RR(EXT_VEL1)),ABS(U_LL(EXT_VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )
F = 0.5*((F_L+F_R) - LambdaMax*(U_RR(EXT_CONS) - U_LL(EXT_CONS)))

END SUBROUTINE Riemann_LF


!=================================================================================================================================
!> Roe's approximate Riemann solver using the Harten and Hymen II entropy fix, see
!> Pelanti, Marica & Quartapelle, Luigi & Vigevano, L & Vigevano, Luigi. (2018):
!>  A review of entropy fixes as applied to Roe's linearization.
!=================================================================================================================================
SUBROUTINE Riemann_RoeEntropyFix(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iVar
REAL                    :: c_L,c_R
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                    :: RoeVel(3),RoeH,Roec,RoeDens
REAL,DIMENSION(5)       :: r1,r2,r3,r4,r5,a,al,ar,Delta_U,Alpha  ! Roe eigenvectors
REAL                    :: tmp,da
REAL                    :: LambdaMax
!=================================================================================================================================
c_L       = SPEEDOFSOUND_HE(U_LL)
c_R       = SPEEDOFSOUND_HE(U_RR)
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
RoeDens   = SQRT(U_LL(EXT_DENS)*U_RR(EXT_DENS))
! Roe+Pike version of Roe Riemann solver

! calculate jump
Delta_U(DENS)   = U_RR(EXT_DENS) - U_LL(EXT_DENS)
Delta_U(VELV)   = U_RR(EXT_VELV) - U_LL(EXT_VELV)
Delta_U(PRES)   = U_RR(EXT_PRES) - U_LL(EXT_PRES)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate wave strenghts
tmp      = 0.5/(Roec*Roec)
Alpha(1) = tmp*(Delta_U(5)-RoeDens*Roec*Delta_U(2))
Alpha(2) = Delta_U(1) - Delta_U(5)*2.*tmp
Alpha(3) = RoeDens*Delta_U(3)
Alpha(4) = RoeDens*Delta_U(4)
Alpha(5) = tmp*(Delta_U(5)+RoeDens*Roec*Delta_U(2))

! Harten+Hyman entropy fix (apply only for acoustic waves, don't fix r)

al(1) = U_LL(EXT_VEL1) - c_L
al(2) = U_LL(EXT_VEL1)
al(3) = U_LL(EXT_VEL1)
al(4) = U_LL(EXT_VEL1)
al(5) = U_LL(EXT_VEL1) + c_L
ar(1) = U_RR(EXT_VEL1) - c_R
ar(2) = U_RR(EXT_VEL1)
ar(3) = U_RR(EXT_VEL1)
ar(4) = U_RR(EXT_VEL1)
ar(5) = U_RR(EXT_VEL1) + c_R
! HH1
!IF(ABS(a(1)).LT.da1) a(1)=da1
!IF(ABS(a(5)).LT.da5) a(5)=da5
! HH2
DO iVar=1,5
  da = MAX(0.,a(iVar)-al(iVar),ar(iVar)-a(iVar))

  IF(ABS(a(iVar)).LT.da) THEN
    a(iVar)=0.5*(a(iVar)*a(iVar)/da+da)
  ELSE
    a(iVar) = ABS(a(iVar))
  END IF
END DO

! assemble Roe flux for NS part of the equation system
F(1:5)=0.5*((F_L(1:5)+F_R(1:5)) - &
               Alpha(1)*a(1)*r1 - &
               Alpha(2)*a(2)*r2 - &
               Alpha(3)*a(3)*r3 - &
               Alpha(4)*a(4)*r4 - &
               Alpha(5)*a(5)*r5)

! Revert to LF for the RANS SA equations
LambdaMax = MAX( ABS(U_RR(EXT_VEL1)),ABS(U_LL(EXT_VEL1)) ) + MAX(c_L,c_R)
F(MUSA) = 0.5*((F_L(MUSA)+F_R(MUSA)) - LambdaMax*(U_RR(EXT_MUSA) - U_LL(EXT_MUSA)))
END SUBROUTINE Riemann_RoeEntropyFix


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
