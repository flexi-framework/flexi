!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! Copyright (c) 2016-2017 Gregor Gassner (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Florian Hindenlang (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Andrew Winters (github.com/project-fluxo/fluxo)
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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

ABSTRACT INTERFACE
  PPURE SUBROUTINE RiemannInt(F_L,F_R,U_LL,U_RR,F)
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
    REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
    REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
  END SUBROUTINE RiemannInt
END INTERFACE

PROCEDURE(RiemannInt),POINTER :: Riemann_pointer    !< pointer defining the standard inner Riemann solver
PROCEDURE(RiemannInt),POINTER :: RiemannBC_pointer  !< pointer defining the standard BC    Riemann solver
REAL                   :: RiemannPhi                !< factor for Riemann solver with fox for grid aligned shocks
INTEGER,PARAMETER      :: PRM_RIEMANN_SAME                = -1
INTEGER,PARAMETER      :: PRM_RIEMANN_LF                  = 1
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLC                = 2
INTEGER,PARAMETER      :: PRM_RIEMANN_ROE                 = 3
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEL2               = 32
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEENTROPYFIX       = 33
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEENTROPYGSHOCKFIX = 34
INTEGER,PARAMETER      :: PRM_RIEMANN_HLL                 = 4
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLE                = 5
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLEM               = 6
#ifdef SPLIT_DG
INTEGER,PARAMETER      :: PRM_RIEMANN_CH                  = 7
INTEGER,PARAMETER      :: PRM_RIEMANN_IR                  = 8
INTEGER,PARAMETER      :: PRM_RIEMANN_Average             = 0
#endif

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

PUBLIC:: DefineParametersRiemann
PUBLIC:: InitRiemann
PUBLIC:: Riemann
PUBLIC:: FinalizeRiemann
#if PARABOLIC
PUBLIC:: ViscousFlux
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
CALL prms%CreateIntFromStringOption('Riemann',   "Riemann solver to be used: LF, HLLC, Roe, RoeEntropyFix, HLL, HLLE, HLLEM", &
                                                 "RoeEntropyFix")
CALL addStrListEntry('Riemann','lf',                         PRM_RIEMANN_LF)
CALL addStrListEntry('Riemann','hllc',                       PRM_RIEMANN_HLLC)
CALL addStrListEntry('Riemann','roe',                        PRM_RIEMANN_ROE)
CALL addStrListEntry('Riemann','roeentropyfix',              PRM_RIEMANN_ROEENTROPYFIX)
CALL addStrListEntry('Riemann','roel2',                      PRM_RIEMANN_ROEL2)
CALL addStrListEntry('Riemann','hll',                        PRM_RIEMANN_HLL)
CALL addStrListEntry('Riemann','hlle',                       PRM_RIEMANN_HLLE)
CALL addStrListEntry('Riemann','hllem',                      PRM_RIEMANN_HLLEM)
CALL addStrListEntry('Riemann','roeentropygshockfix',        PRM_RIEMANN_ROEENTROPYGSHOCKFIX)
#ifdef SPLIT_DG
CALL addStrListEntry('Riemann','ch',                         PRM_RIEMANN_CH)
CALL addStrListEntry('Riemann','ir',                         PRM_RIEMANN_IR)
CALL addStrListEntry('Riemann','avg',                        PRM_RIEMANN_Average)
#endif
CALL prms%CreateIntFromStringOption('RiemannBC', "Riemann solver used for boundary conditions: Same, LF, Roe, RoeEntropyFix, "//&
                                                 "HLL, HLLE, HLLEM",&
                                                 "Same")
CALL addStrListEntry('RiemannBC','lf',                       PRM_RIEMANN_LF)
CALL addStrListEntry('RiemannBC','hllc',                     PRM_RIEMANN_HLLC)
CALL addStrListEntry('RiemannBC','roe',                      PRM_RIEMANN_ROE)
CALL addStrListEntry('RiemannBC','roeentropyfix',            PRM_RIEMANN_ROEENTROPYFIX)
CALL addStrListEntry('RiemannBC','roel2',                    PRM_RIEMANN_ROEL2)
CALL addStrListEntry('RiemannBC','hll',                      PRM_RIEMANN_HLL)
CALL addStrListEntry('RiemannBC','hlle',                     PRM_RIEMANN_HLLE)
CALL addStrListEntry('RiemannBC','hllem',                    PRM_RIEMANN_HLLEM)
CALL addStrListEntry('RiemannBC','roeentropygshockfix',      PRM_RIEMANN_ROEENTROPYGSHOCKFIX)
#ifdef SPLIT_DG
CALL addStrListEntry('RiemannBC','ch',                       PRM_RIEMANN_CH)
CALL addStrListEntry('RiemannBC','ir',                       PRM_RIEMANN_IR)
CALL addStrListEntry('RiemannBC','avg',                      PRM_RIEMANN_Average)
#endif
CALL addStrListEntry('RiemannBC','same',                     PRM_RIEMANN_SAME)
! Riemann solver specific parameters
CALL prms%CreateRealOption('RiemannPhiFactor', "Factor to adjust imbalance of advective and acoustic dissipation, fixes grid alligned shock instabilities", "5")
END SUBROUTINE DefineParametersRiemann


!==================================================================================================================================!
!> Initialize Riemann solver routines, read inner and BC Riemann solver parameters and set pointers
!==================================================================================================================================!
SUBROUTINE InitRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR, GETREAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: Riemann
!==================================================================================================================================
#ifndef SPLIT_DG
Riemann = GETINTFROMSTR('Riemann')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_LF)
  Riemann_pointer => Riemann_LF
CASE(PRM_RIEMANN_HLLC)
  Riemann_pointer => Riemann_HLLC
CASE(PRM_RIEMANN_ROE)
  Riemann_pointer => Riemann_Roe
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  Riemann_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_ROEL2)
  Riemann_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_HLL)
  Riemann_pointer => Riemann_HLL
CASE(PRM_RIEMANN_HLLE)
  Riemann_pointer => Riemann_HLLE
CASE(PRM_RIEMANN_HLLEM)
  Riemann_pointer => Riemann_HLLEM
CASE(PRM_RIEMANN_ROEENTROPYGSHOCKFIX)
  Riemann_pointer => Riemann_RoeEntropyGShockFix
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
CASE(PRM_RIEMANN_HLLC)
  RiemannBC_pointer => Riemann_HLLC
CASE(PRM_RIEMANN_ROE)
  RiemannBC_pointer => Riemann_Roe
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  RiemannBC_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_ROEL2)
  RiemannBC_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_HLL)
  RiemannBC_pointer => Riemann_HLL
CASE(PRM_RIEMANN_HLLE)
  RiemannBC_pointer => Riemann_HLLE
CASE(PRM_RIEMANN_HLLEM)
  RiemannBC_pointer => Riemann_HLLEM
CASE(PRM_RIEMANN_ROEENTROPYGSHOCKFIX)
  RiemannBC_pointer => Riemann_RoeEntropyGShockFix
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'RiemannBC solver not defined!')
END SELECT

#else
Riemann = GETINTFROMSTR('Riemann')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_LF)
  Riemann_pointer => Riemann_LF
CASE(PRM_RIEMANN_ROE)
  Riemann_pointer => Riemann_Roe
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  Riemann_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_ROEL2)
  Riemann_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_IR)
  Riemann_pointer => Riemann_IR
CASE(PRM_RIEMANN_CH)
  Riemann_pointer => Riemann_CH
CASE(PRM_RIEMANN_Average)
  Riemann_pointer => Riemann_FluxAverage
CASE(PRM_RIEMANN_ROEENTROPYGSHOCKFIX)
  Riemann_pointer => Riemann_RoeEntropyGShockFix
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
CASE(PRM_RIEMANN_ROE)
  RiemannBC_pointer => Riemann_Roe
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  RiemannBC_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_IR)
  RiemannBC_pointer => Riemann_IR
CASE(PRM_RIEMANN_ROEL2)
  RiemannBC_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_CH)
  RiemannBC_pointer => Riemann_CH
CASE(PRM_RIEMANN_Average)
  RiemannBC_pointer => Riemann_FluxAverage
CASE(PRM_RIEMANN_ROEENTROPYGSHOCKFIX)
  RiemannBC_pointer => Riemann_RoeEntropyGShockFix
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'RiemannBC solver not defined!')
END SELECT
#endif /*SPLIT_DG*/

! Read in Phi for Riemann with grid aligned shock fix
IF (Riemann.EQ.PRM_RIEMANN_ROEENTROPYGSHOCKFIX) &
  RiemannPhi = GETREAL("RiemannRiemannPhi")

END SUBROUTINE InitRiemann


!==================================================================================================================================
!> Computes the numerical flux for a side calling the flux calculation pointwise.
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
SUBROUTINE Riemann_Side(Nloc,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,doBC)
! MODULES
USE MOD_Flux         ,ONLY:EvalEulerFlux1D_fast
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                          :: Nloc      !< local polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L       !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R       !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L   !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_R   !< primitive solution at right side of the interface
REAL,DIMENSION(3          ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv,t1,t2  !> normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: FOut      !< advective flux
LOGICAL,INTENT(IN)                                          :: doBC      !< marker whether side is a BC side
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

DO j=0,ZDIM(Nloc); DO i=0,Nloc
  ! Momentum has to be rotated using the normal system individual for each
  ! left state: U_L
  U_LL(EXT_DENS)=U_L(DENS,i,j)
  U_LL(EXT_SRHO)=1./U_LL(EXT_DENS)
  U_LL(EXT_ENER)=U_L(ENER,i,j)
  U_LL(EXT_PRES)=UPrim_L(PRES,i,j)

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

#ifndef SPLIT_DG
  CALL EvalEulerFlux1D_fast(U_LL,F_L)
  CALL EvalEulerFlux1D_fast(U_RR,F_R)
#endif /*SPLIT_DG*/

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

! Momentum has to be rotated using the normal system individual for each
! left state: U_L
U_LL(EXT_DENS)=U_L(DENS)
U_LL(EXT_SRHO)=1./U_LL(EXT_DENS)
U_LL(EXT_ENER)=U_L(ENER)
U_LL(EXT_PRES)=UPrim_L(PRES)

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

# ifndef SPLIT_DG
CALL EvalEulerFlux1D_fast(U_LL,F_L)
CALL EvalEulerFlux1D_fast(U_RR,F_R)
#endif /*SPLIT_DG*/

CALL Riemann_loc(F_L,F_R,U_LL,U_RR,F)

! Back rotate the normal flux into Cartesian direction
Fout(DENS)=F(DENS)
Fout(MOMV)=nv(:)*F(MOM1)  &
          +t1(:)*F(MOM2)  &
#if PP_dim==3
          +t2(:)*F(MOM3)
#else
          + 0.
#endif
Fout(ENER)=F(ENER)
END SUBROUTINE Riemann_Point


#if PARABOLIC
!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux_Side(Nloc,F,UPrim_L,UPrim_R, &
                            gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
#if EDDYVISCOSITY
                           ,muSGS_L,muSGS_R &
#endif
                           )
! MODULES
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
USE MOD_Lifting_Vars ,ONLY: diffFluxX_L,diffFluxY_L,diffFluxZ_L
USE MOD_Lifting_Vars ,ONLY: diffFluxX_R,diffFluxY_R,diffFluxZ_R
! IMPLICIT VARIABLE HANDLING
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
#if EDDYVISCOSITY
REAL,DIMENSION(1             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: muSGS_L,muSGS_R   !> eddy viscosity left/right of the interface
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                                       :: p,q
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
CALL EvalDiffFlux3D(Nloc,UPrim_L,   gradUx_L,   gradUy_L,   gradUz_L  &
                                ,diffFluxX_L,diffFluxY_L,diffFluxZ_L  &
#if EDDYVISCOSITY
                   ,muSGS_L &
#endif
      )
CALL EvalDiffFlux3D(Nloc,UPrim_R,   gradUx_R,   gradUy_R,   gradUz_R  &
                                ,diffFluxX_R,diffFluxY_R,diffFluxZ_R  &
#if EDDYVISCOSITY
                   ,muSGS_R&
#endif
      )
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
                             gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
#if EDDYVISCOSITY
                            ,muSGS_L,muSGS_R &
#endif
                            )
! MODULES
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                           !> solution in primitive variables at left/right side of the interface
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim_L,UPrim_R
                                           !> solution gradients in x/y/z-direction left/right of the interface
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DIMENSION(3             ),INTENT(IN)  :: nv  !< normal vector
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: F   !< viscous flux
#if EDDYVISCOSITY
REAL,INTENT(IN)                            :: muSGS_L,muSGS_R    !> eddy viscosity left/right of the interface
#endif
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
                           ,diffFluxX_L,diffFluxY_L,diffFluxZ_L  &
#if EDDYVISCOSITY
                   ,muSGS_L &
#endif
      )
CALL EvalDiffFlux3D(UPrim_R,   gradUx_R,   gradUy_R,   gradUz_R  &
                           ,diffFluxX_R,diffFluxY_R,diffFluxZ_R  &
#if EDDYVISCOSITY
                   ,muSGS_R&
#endif
      )
! Arithmetic mean of the fluxes
F(:)=0.5*(nv(1)*(diffFluxX_L(:)+diffFluxX_R(:)) &
         +nv(2)*(diffFluxY_L(:)+diffFluxY_R(:)) &
         +nv(3)*(diffFluxZ_L(:)+diffFluxZ_R(:)))
END SUBROUTINE ViscousFlux_Point
#endif /* PARABOLIC */


!==================================================================================================================================
!> Local Lax-Friedrichs (Rusanov) Riemann solver
!==================================================================================================================================
PPURE SUBROUTINE Riemann_LF(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
! IMPLICIT VARIABLE HANDLING
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
LambdaMax = MAX( ABS(U_RR(EXT_VEL1)),ABS(U_LL(EXT_VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )
#ifndef SPLIT_DG
F = 0.5*((F_L+F_R) - LambdaMax*(U_RR(CONS) - U_LL(CONS)))
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! compute surface flux
F = F - 0.5*LambdaMax*(U_RR(CONS) - U_LL(CONS))
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_LF

!=================================================================================================================================
!> Harten-Lax-Van-Leer Riemann solver resolving contact discontinuity
!=================================================================================================================================
PPURE SUBROUTINE Riemann_HLLC(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: KappaM1!,kappa
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                           !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                           !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F    !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL    :: RoeVel(3),RoeH,Roec,absVel
REAL    :: Ssl,Ssr,SStar
REAL    :: U_Star(PP_nVar),EStar
REAL    :: sMu_L,sMu_R
!REAL    :: c_L,c_R
!=================================================================================================================================
! HLLC flux

! Version A: Basic Davis estimate for wave speed
!Ssl = U_LL(EXT_VEL1) - SPEEDOFSOUND_HE(U_LL)
!Ssr = U_RR(EXT_VEL1) + SPEEDOFSOUND_HE(U_RR)

! Version B: Basic Davis estimate for wave speed
!c_L = SPEEDOFSOUND_HE(U_LL)
!c_R = SPEEDOFSOUND_HE(U_RR)
!Ssl = MIN(U_LL(EXT_VEL1) - c_L,U_RR(EXT_VEL1) - c_R)
!Ssr = MAX(U_LL(EXT_VEL1) + c_L,U_RR(EXT_VEL1) + c_R)

! Version C: Better Roe estimate for wave speeds Davis, Einfeldt
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L       )     * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = SQRT(KappaM1*(RoeH-0.5*absVel))
Ssl       = RoeVel(1) - Roec
Ssr       = RoeVel(1) + Roec

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  sMu_L = Ssl - U_LL(EXT_VEL1)
  sMu_R = Ssr - U_RR(EXT_VEL1)
  SStar = (U_RR(EXT_PRES) - U_LL(EXT_PRES) + U_LL(EXT_MOM1)*sMu_L - U_RR(EXT_MOM1)*sMu_R) / (U_LL(EXT_DENS)*sMu_L - U_RR(EXT_DENS)*sMu_R)
  IF ((Ssl .LE. 0.).AND.(SStar .GE. 0.)) THEN
    EStar  = TOTALENERGY_HE(U_LL) + (SStar-U_LL(EXT_VEL1))*(SStar + U_LL(EXT_PRES)*U_LL(EXT_SRHO)/sMu_L)
    U_Star = U_LL(EXT_DENS) * sMu_L/(Ssl-SStar) * (/ 1., SStar, U_LL(EXT_VEL2:EXT_VEL3), EStar /)
    F=F_L+Ssl*(U_Star-U_LL(CONS))
  ELSE
    EStar  = TOTALENERGY_HE(U_RR) + (SStar-U_RR(EXT_VEL1))*(SStar + U_RR(EXT_PRES)*U_RR(EXT_SRHO)/sMu_R)
    U_Star = U_RR(EXT_DENS) * sMu_R/(Ssr-SStar) * (/ 1., SStar, U_RR(EXT_VEL2:EXT_VEL3), EStar /)
    F=F_R+Ssr*(U_Star-U_RR(CONS))
  END IF
END IF ! subsonic case
END SUBROUTINE Riemann_HLLC

!=================================================================================================================================
!> Roe's approximate Riemann solver
!=================================================================================================================================
PPURE SUBROUTINE Riemann_Roe(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars  ,ONLY: kappaM1
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL                    :: RoeVel(3),RoeH,Roec,absVel
REAL,DIMENSION(PP_nVar) :: a,r1,r2,r3,r4,r5  ! Roe eigenvectors
REAL                    :: Alpha1,Alpha2,Alpha3,Alpha4,Alpha5,Delta_U(PP_nVar+1)
!=================================================================================================================================
! Roe flux
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate differences
Delta_U(CONS)     = U_RR(CONS) - U_LL(CONS)
Delta_U(DELTA_U6) = Delta_U(DELTA_U5)-(Delta_U(DELTA_U3)-RoeVel(DELTA_U2)*Delta_U(DELTA_U1))*RoeVel(2) -&
                    (Delta_U(DELTA_U4)-RoeVel(DELTA_U3)*Delta_U(DELTA_U1))*RoeVel(DELTA_U3)
! calculate factors
Alpha3 = Delta_U(DELTA_U3) - RoeVel(DELTA_U2)*Delta_U(DELTA_U1)
Alpha4 = Delta_U(DELTA_U4) - RoeVel(DELTA_U3)*Delta_U(DELTA_U1)
Alpha2 = ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)
Alpha1 = 0.5/Roec * (Delta_U(DELTA_U1)*(RoeVel(1)+Roec) - Delta_U(DELTA_U2) - Roec*Alpha2)
Alpha5 = Delta_U(DELTA_U1) - Alpha1 - Alpha2
#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R) - &
       Alpha1*ABS(a(1))*r1 - &
       Alpha2*ABS(a(2))*r2 - &
       Alpha3*ABS(a(3))*r3 - &
       Alpha4*ABS(a(4))*r4 - &
       Alpha5*ABS(a(5))*r5)
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! assemble Roe flux
F = F - 0.5*(Alpha1*ABS(a(1))*r1 + &
             Alpha2*ABS(a(2))*r2 + &
             Alpha3*ABS(a(3))*r3 + &
             Alpha4*ABS(a(4))*r4 + &
             Alpha5*ABS(a(5))*r5)
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_Roe


!=================================================================================================================================
!> Roe's approximate Riemann solver using the Harten and Hymen II entropy fix, see
!> Pelanti, Marica & Quartapelle, Luigi & Vigevano, L & Vigevano, Luigi. (2018):
!>  A review of entropy fixes as applied to Roe's linearization.
!=================================================================================================================================
PPURE SUBROUTINE Riemann_RoeEntropyFix(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iVar
REAL                    :: c_L,c_R
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                    :: RoeVel(3),RoeH,Roec,RoeDens
REAL,DIMENSION(PP_nVar) :: r1,r2,r3,r4,r5,a,al,ar,Delta_U,Alpha  ! Roe eigenvectors
REAL                    :: tmp,da
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
Delta_U(DELTA_U1)   = U_RR(EXT_DENS) - U_LL(EXT_DENS)
Delta_U(DELTA_UV)   = U_RR(EXT_VELV) - U_LL(EXT_VELV)
Delta_U(DELTA_U5)   = U_RR(EXT_PRES) - U_LL(EXT_PRES)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate wave strenghts
tmp      = 0.5/(Roec*Roec)
Alpha(1) = tmp*(Delta_U(DELTA_U5)-RoeDens*Roec*Delta_U(DELTA_U2))
Alpha(2) = Delta_U(DELTA_U1) - Delta_U(DELTA_U5)*2.*tmp
Alpha(3) = RoeDens*Delta_U(DELTA_U3)
Alpha(4) = RoeDens*Delta_U(DELTA_U4)
Alpha(5) = tmp*(Delta_U(DELTA_U5)+RoeDens*Roec*Delta_U(DELTA_U2))

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

#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R)        - &
       Alpha(1)*a(1)*r1 - &
       Alpha(2)*a(2)*r2 - &
       Alpha(3)*a(3)*r3 - &
       Alpha(4)*a(4)*r4 - &
       Alpha(5)*a(5)*r5)
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! for KG or PI flux eigenvalues have to be altered to ensure consistent KE dissipation
! assemble Roe flux
F= F - 0.5*(Alpha(1)*a(1)*r1 + &
            Alpha(2)*a(2)*r2 + &
            Alpha(3)*a(3)*r3 + &
            Alpha(4)*a(4)*r4 + &
            Alpha(5)*a(5)*r5)
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_RoeEntropyFix


!=================================================================================================================================
!> Roe's approximate Riemann solver using the Harten and Hymen II entropy fix, see
!> Pelanti, Marica & Quartapelle, Luigi & Vigevano, L & Vigevano, Luigi. (2018):
!>  A review of entropy fixes as applied to Roe's linearization.
!> The GShockFix or Grid Alligned Shock Fix uses the modified Roe flux formulation from the paper:
!> N. Fleischmann, S. Adami, X. Y. Hu, and N. A. Adams, “A low dissipation method to cure the grid-aligned shock instability,”
!> Journal of Computational Physics, vol. 401, p. 109004, Jan. 2020
!=================================================================================================================================
PPURE SUBROUTINE Riemann_RoeEntropyGShockFix(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars  ,ONLY: Kappa,KappaM1
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                                !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR
                                                !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iVar
REAL                                :: c_L,c_R
REAL                                :: H_L,H_R
REAL                                :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                                :: RoeVel(3),RoeH,Roec,RoeDens
REAL,DIMENSION(PP_nVar)             :: r1,r2,r3,r4,r5,a,al,ar,Delta_U,Alpha  ! Roe eigenvectors
REAL                                :: tmp,da
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
Delta_U(DELTA_U1)   = U_RR(EXT_DENS) - U_LL(EXT_DENS)
Delta_U(DELTA_UV)   = U_RR(EXT_VELV) - U_LL(EXT_VELV)
Delta_U(DELTA_U5)   = U_RR(EXT_PRES) - U_LL(EXT_PRES)

! mean eigenvalues and eigenvectors
a  = (/ SIGN(1., RoeVel(1))*MAX(Roec/RiemannPhi, ABS(RoeVel(1))) - Roec, &
        SIGN(1., RoeVel(1))*MAX(Roec/RiemannPhi, ABS(RoeVel(1))), &
        SIGN(1., RoeVel(1))*MAX(Roec/RiemannPhi, ABS(RoeVel(1))), &
        SIGN(1., RoeVel(1))*MAX(Roec/RiemannPhi, ABS(RoeVel(1))), &
        SIGN(1., RoeVel(1))*MAX(Roec/RiemannPhi, ABS(RoeVel(1))) + Roec /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate wave strenghts
tmp      = 0.5/(Roec*Roec)
Alpha(1) = tmp*(Delta_U(DELTA_U5)-RoeDens*Roec*Delta_U(DELTA_U2))
Alpha(2) = Delta_U(DELTA_U1) - Delta_U(DELTA_U5)*2.*tmp
Alpha(3) = RoeDens*Delta_U(DELTA_U3)
Alpha(4) = RoeDens*Delta_U(DELTA_U4)
Alpha(5) = tmp*(Delta_U(DELTA_U5)+RoeDens*Roec*Delta_U(DELTA_U2))

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

#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R)        - &
       Alpha(1)*a(1)*r1 - &
       Alpha(2)*a(2)*r2 - &
       Alpha(3)*a(3)*r3 - &
       Alpha(4)*a(4)*r4 - &
       Alpha(5)*a(5)*r5)
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! for KG or PI flux eigenvalues have to be altered to ensure consistent KE dissipation
! assemble Roe flux
F= F - 0.5*(Alpha(1)*a(1)*r1 + &
            Alpha(2)*a(2)*r2 + &
            Alpha(3)*a(3)*r3 + &
            Alpha(4)*a(4)*r4 + &
            Alpha(5)*a(5)*r5)
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_RoeEntropyGShockFix


!=================================================================================================================================
!> low mach number Roe's approximate Riemann solver according to Oßwald(2015)
!=================================================================================================================================
PPURE SUBROUTINE Riemann_RoeL2(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars  ,ONLY: kappaM1,kappa
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL                    :: RoeVel(3),RoeH,Roec,absVel
REAL                    :: Ma_loc ! local Mach-Number
REAL,DIMENSION(PP_nVar) :: a,r1,r2,r3,r4,r5  ! Roe eigenvectors
REAL                    :: Alpha1,Alpha2,Alpha3,Alpha4,Alpha5,Delta_U(PP_nVar+1)
!=================================================================================================================================
! Roe flux
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate differences
Delta_U(CONS) = U_RR(EXT_CONS) - U_LL(EXT_CONS)
Delta_U(DELTA_U6)   = Delta_U(DELTA_U5)-(Delta_U(DELTA_U3)-RoeVel(2)*Delta_U(DELTA_U1))*RoeVel(2) - &
                      (Delta_U(DELTA_U4)-RoeVel(3)*Delta_U(DELTA_U1))*RoeVel(3)

! low Mach-Number fix
Ma_loc = SQRT(absVel)/(Roec*SQRT(kappa))
Delta_U(DELTA_UV) = Delta_U(DELTA_UV) * Ma_loc

! calculate factors
Alpha3 = Delta_U(DELTA_U3) - RoeVel(2)*Delta_U(DELTA_U1)
Alpha4 = Delta_U(DELTA_U4) - RoeVel(3)*Delta_U(DELTA_U1)
Alpha2 = ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)
Alpha1 = 0.5/Roec * (Delta_U(DELTA_U1)*(RoeVel(1)+Roec) - Delta_U(DELTA_U2) - Roec*Alpha2)
Alpha5 = Delta_U(DELTA_U1) - Alpha1 - Alpha2

#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R) - &
       Alpha1*ABS(a(1))*r1 - &
       Alpha2*ABS(a(2))*r2 - &
       Alpha3*ABS(a(3))*r3 - &
       Alpha4*ABS(a(4))*r4 - &
       Alpha5*ABS(a(5))*r5)
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! assemble Roe flux
F = F - 0.5*(Alpha1*ABS(a(1))*r1 + &
             Alpha2*ABS(a(2))*r2 + &
             Alpha3*ABS(a(3))*r3 + &
             Alpha4*ABS(a(4))*r4 + &
             Alpha5*ABS(a(5))*r5)
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_RoeL2

!=================================================================================================================================
!> Standard Harten-Lax-Van-Leer Riemann solver without contact discontinuity
!=================================================================================================================================
PPURE SUBROUTINE Riemann_HLL(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars, ONLY: KappaM1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roec
REAL    :: Ssl,Ssr
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L)            * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
! HLL flux
! Basic Davis estimate for wave speed
!Ssl = U_LL(EXT_VEL1) - c_L
!Ssr = U_RR(EXT_VEL1) + c_R
! Better Roe estimate for wave speeds Davis, Einfeldt
Ssl = RoeVel(1) - Roec
Ssr = RoeVel(1) + Roec
! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(EXT_CONS)-U_LL(EXT_CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLL


!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt Riemann solver
!=================================================================================================================================
PPURE SUBROUTINE Riemann_HLLE(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roec
REAL    :: Ssl,Ssr,beta
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L)            * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
! HLLE flux (positively conservative)
beta=BETA_RIEMANN_H()
SsL=MIN(RoeVel(1)-Roec,U_LL(EXT_VEL1) - beta*SPEEDOFSOUND_HE(U_LL), 0.)
SsR=MAX(RoeVel(1)+Roec,U_RR(EXT_VEL1) + beta*SPEEDOFSOUND_HE(U_RR), 0.)

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(EXT_CONS)-U_LL(EXT_CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLLE


!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt-Munz Riemann solver
!=================================================================================================================================
PPURE SUBROUTINE Riemann_HLLEM(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                   :: H_L,H_R
REAL                                   :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                                   :: RoeVel(3),RoeH,Roec,RoeDens
REAL                                   :: Ssl,Ssr
REAL                                   :: Alpha(2:4),delta,beta
REAL,DIMENSION(PP_nVar)                :: r2,r3,r4  ! Roe eigenvectors + jump in prims
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L)            * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
RoeDens   = SQRT(U_LL(EXT_DENS)*U_RR(EXT_DENS))
! HLLEM flux (positively conservative)
beta=BETA_RIEMANN_H()
SsL=MIN(RoeVel(1)-Roec,U_LL(EXT_VEL1) - beta*SPEEDOFSOUND_HE(U_LL), 0.)
SsR=MAX(RoeVel(1)+Roec,U_RR(EXT_VEL1) + beta*SPEEDOFSOUND_HE(U_RR), 0.)

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  ! delta
  delta = Roec/(Roec+ABS(0.5*(Ssl+Ssr)))

  ! mean eigenvectors
  Alpha(2)   = (U_RR(EXT_DENS)-U_LL(EXT_DENS))  - (U_RR(EXT_PRES)-U_LL(EXT_PRES))/(Roec*Roec)
  Alpha(3:4) = RoeDens*(U_RR(EXT_VEL2:EXT_VEL3) - U_LL(EXT_VEL2:EXT_VEL3))
  r2 = (/ 1., RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel /)
  r3 = (/ 0., 0.,        1.,        0.,        RoeVel(2)  /)
  r4 = (/ 0., 0.,        0.,        1.,        RoeVel(3)  /)

  F=(Ssr*F_L-Ssl*F_R + Ssl*Ssr* &
     (U_RR(EXT_CONS)-U_LL(EXT_CONS) - delta*(r2*Alpha(2)+r3*Alpha(3)+r4*Alpha(4))))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLLEM


#ifdef SPLIT_DG
!==================================================================================================================================
!> Riemann solver using purely the average fluxes
!==================================================================================================================================
PPURE SUBROUTINE Riemann_FluxAverage(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer
! IMPLICIT VARIABLE HANDLING
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
!==================================================================================================================================
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
END SUBROUTINE Riemann_FluxAverage


!==================================================================================================================================
!> kinetic energy preserving and entropy consistent flux according to Chandrashekar (2012)
!==================================================================================================================================
PPURE SUBROUTINE Riemann_CH(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,sKappaM1
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer,GetLogMean
! IMPLICIT VARIABLE HANDLING
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
REAL                               :: LambdaMax
REAL                               :: beta_LL,beta_RR   ! auxiliary variables for the inverse Temperature
REAL                               :: rhoMean           ! auxiliary variable for the mean density
REAL                               :: uMean,vMean,wMean ! auxiliary variable for the average velocities
REAL                               :: betaLogMean       ! auxiliary variable for the logarithmic mean inverse temperature
!==================================================================================================================================
! Lax-Friedrichs
LambdaMax = MAX( ABS(U_RR(EXT_VEL1)),ABS(U_LL(EXT_VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )

! average quantities
rhoMean = 0.5*(U_LL(EXT_DENS) + U_RR(EXT_DENS))
uMean   = 0.5*(U_LL(EXT_VEL1) + U_RR(EXT_VEL1))
vMean   = 0.5*(U_LL(EXT_VEL2) + U_RR(EXT_VEL2))
wMean   = 0.5*(U_LL(EXT_VEL3) + U_RR(EXT_VEL3))

! inverse temperature
beta_LL = 0.5*U_LL(EXT_DENS)/U_LL(EXT_PRES)
beta_RR = 0.5*U_RR(EXT_DENS)/U_RR(EXT_PRES)

! logarithmic mean
CALL GetLogMean(beta_LL,beta_RR,betaLogMean)

! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)

!compute flux
F(DENS:MOM3) = F(DENS:MOM3) - 0.5*LambdaMax*(U_RR(EXT_DENS:EXT_MOM3)-U_LL(EXT_DENS:EXT_MOM3))
F(ENER)      = F(ENER)      - 0.5*LambdaMax*( &
         (U_RR(EXT_DENS)-U_LL(EXT_DENS))*(0.5*sKappaM1/betaLogMean +0.5*(U_RR(EXT_VEL1)*U_LL(EXT_VEL1)+U_RR(EXT_VEL2)*U_LL(EXT_VEL2)+U_RR(EXT_VEL3)*U_LL(EXT_VEL3))) &
         +rhoMean*uMean*(U_RR(EXT_VEL1)-U_LL(EXT_VEL1)) + rhoMean*vMean*(U_RR(EXT_VEL2)-U_LL(EXT_VEL2)) + rhoMean*wMean*(U_RR(EXT_VEL3)-U_LL(EXT_VEL3)) &
         +0.5*rhoMean*sKappaM1*(1./beta_RR - 1./beta_LL))

END SUBROUTINE Riemann_CH

!=================================================================================================================================
!> Approximate riemann solver for the Ismael and Roe split flux (in entropy variables)
!=================================================================================================================================
PPURE SUBROUTINE Riemann_IR(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1,KappaP1
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer,GetLogMean
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iVar
REAL                    :: sR,sL,rho_pR,rho_pL,vR(5),vL(5)
REAL,DIMENSION(5)       :: z_LL,z_RR
REAL                    :: z1LogMean,z5LogMean,z1Mean
REAL                    :: rhoHat,uHat,vHat,wHat,pHat,p2Hat,hHat,ahat
REAL                    :: Rhat(5,5),dHat(5),vJump(5),diss(5)
!=================================================================================================================================
sR     =  LOG(U_RR(EXT_PRES)) - kappa*LOG(U_RR(EXT_DENS))
sL     =  LOG(U_LL(EXT_PRES)) - kappa*LOG(U_LL(EXT_DENS))
rho_pR =  U_RR(EXT_DENS)/U_RR(EXT_PRES)
rho_pL =  U_LL(EXT_DENS)/U_LL(EXT_PRES)
vR(1)  =  (kappa-sR)/(kappaM1) - 0.5*rho_pR*(SUM(U_RR(EXT_VELV)*U_RR(EXT_VELV)))
vL(1)  =  (kappa-sL)/(kappaM1) - 0.5*rho_pL*(SUM(U_LL(EXT_VELV)*U_LL(EXT_VELV)))
vR(2)  =  rho_pR*U_RR(EXT_VEL1)
vL(2)  =  rho_pL*U_LL(EXT_VEL1)
vR(3)  =  rho_pR*U_RR(EXT_VEL2)
vL(3)  =  rho_pL*U_LL(EXT_VEL2)
vR(4)  =  rho_pR*U_RR(EXT_VEL3)
vL(4)  =  rho_pL*U_LL(EXT_VEL3)
vR(5)  = -rho_pR
vL(5)  = -rho_pL

! Compute parameter vector for left and right state
z_LL(1) = SQRT(U_LL(EXT_DENS)/U_LL(EXT_PRES))
z_LL(2) = z_LL(1)*U_LL(EXT_VEL1)
z_LL(3) = z_LL(1)*U_LL(EXT_VEL2)
z_LL(4) = z_LL(1)*U_LL(EXT_VEL3)
z_LL(5) = SQRT(U_LL(EXT_DENS)*U_LL(EXT_PRES))
z_RR(1) = SQRT(U_RR(EXT_DENS)/U_RR(EXT_PRES))
z_RR(2) = z_RR(1)*U_RR(EXT_VEL1)
z_RR(3) = z_RR(1)*U_RR(EXT_VEL2)
z_RR(4) = z_RR(1)*U_RR(EXT_VEL3)
z_RR(5) = SQRT(U_RR(EXT_DENS)*U_RR(EXT_PRES))

! Compute averaged auxilliary variables
CALL getLogMean(z_LL(1),z_RR(1),z1LogMean)
CALL getLogMean(z_LL(5),z_RR(5),z5LogMean)
z1Mean  = 0.5*(z_LL(1)+z_RR(1))
rhoHat  = z1Mean*z5LogMean
uHat    = 0.5*(z_LL(2)+z_RR(2))/z1Mean
vHat    = 0.5*(z_LL(3)+z_RR(3))/z1Mean
wHat    = 0.5*(z_LL(4)+z_RR(4))/z1Mean
pHat    = 0.5*(z_LL(5)+z_RR(5))/z1Mean
p2Hat   = 1./(2.*Kappa)*(KappaP1 * z5LogMean/z1LogMean + KappaM1 * pHat)
aHat    = SQRT(kappa*p2Hat/rhoHat)
hHat    = Kappa*p2Hat/(rhoHat*KappaM1)+0.5*(uHat**2.+vHat**2.+wHat**2.)

! Matrix of right eigenvectors
RHat(1,:) = (/ 1.               , 1.                                      , 0.   ,  0.  , 1.               /)
RHat(2,:) = (/ uHat - aHat      , uHat                                    , 0.   ,  0.  , uHat + aHat      /)
RHat(3,:) = (/ vHat             , vHat                                    , 1.   ,  0.  , vHat             /)
RHat(4,:) = (/ wHat             , wHat                                    , 0.   ,  1.  , wHat             /)
RHat(5,:) = (/ HHat - uHat*aHat , 0.5*(uHat*uHat + vHat*vHat + wHat*wHat) , vHat , wHat , HHat + uHat*aHat /)
! Diagonal scaling matrix where DHat = ABS(\Lambda)S
DHat(1) = 0.5*ABS(uHat - aHat)*rhoHat/kappa
DHat(2) = ABS(uHat)*(kappaM1/kappa)*rhoHat!*rhoHat*rhoHat
DHat(3) = ABS(uHat)*pHat!*rhoHat
DHat(4) = DHat(3)
DHat(5) = 0.5*ABS(uHat + aHat)*rhoHat/kappa
! Compute the dissipation term RHat*DHat*RHat^T*[v]
vJump = vR - vL
diss  = RHat(1,:)*vJump(1) + RHat(2,:)*vJump(2) + RHat(3,:)*vJump(3) + RHat(4,:)*vJump(4) + RHat(5,:)*vJump(5)
DO iVar = 1,5
  diss(iVar) = DHat(iVar)*diss(iVar)
END DO
diss = RHat(:,1)*diss(1) + RHat(:,2)*diss(2) + RHat(:,3)*diss(3) + RHat(:,4)*diss(4) + RHat(:,5)*diss(5)

! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)

! assemble Roe flux
F= F - 0.5*diss
END SUBROUTINE Riemann_IR
#endif /*SPLIT_DG*/


!==================================================================================================================================
!> Finalize Riemann solver routines
!==================================================================================================================================
SUBROUTINE FinalizeRiemann()
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeRiemann

END MODULE MOD_Riemann
