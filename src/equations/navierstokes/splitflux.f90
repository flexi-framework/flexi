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

!==================================================================================================================================
!> Contains the routines for computing the fluxes for the Split-DG algorithm
!==================================================================================================================================
#include "flexi.h"
#include "eos.h"

MODULE MOD_SplitFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  PURE SUBROUTINE VolumeFlux(URef,UPrimRef,U,UPrim,MRef,M,Flux)
    REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U
    REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim
    REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,M
    REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  PURE SUBROUTINE SurfaceFlux(U_LL,U_RR,F)
    REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F
  END SUBROUTINE
END INTERFACE

PROCEDURE(VolumeFlux),POINTER    :: SplitDGVolume_pointer    !< pointer defining the SpliDG formulation beeing used
PROCEDURE(SurfaceFlux),POINTER   :: SplitDGSurface_pointer   !< pointer defining the SpliDG formulation beeing used
INTEGER                :: SplitIndicator           !< specifying which flux to be used

INTEGER,PARAMETER      :: PRM_SPLITDG_SD          = 0
INTEGER,PARAMETER      :: PRM_SPLITDG_MO          = 1
INTEGER,PARAMETER      :: PRM_SPLITDG_DU          = 2
INTEGER,PARAMETER      :: PRM_SPLITDG_KG          = 3
INTEGER,PARAMETER      :: PRM_SPLITDG_PI          = 4

INTERFACE InitSplitDG
  MODULE PROCEDURE InitSplitDG
END INTERFACE

PUBLIC::InitSplitDG,DefineParametersSplitDG
PUBLIC::SplitDGSurface_pointer,SplitDGVolume_pointer,SplitIndicator
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersSplitDG()
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
CALL prms%SetSection("SplitDG")
CALL prms%CreateIntFromStringOption('SplitDG',"SplitDG formulation to be used: SD, MO, DU, KG, PI")
CALL addStrListEntry('SplitDG','sd',           PRM_SPLITDG_SD)
CALL addStrListEntry('SplitDG','mo',           PRM_SPLITDG_MO)
CALL addStrListEntry('SplitDG','du',           PRM_SPLITDG_DU)
CALL addStrListEntry('SplitDG','kg',           PRM_SPLITDG_KG)
CALL addStrListEntry('SplitDG','pi',           PRM_SPLITDG_PI)

END SUBROUTINE DefineParametersSplitDG

!==================================================================================================================================!
!> Initialize Riemann solver routines, read inner and BC Riemann solver parameters and set pointers
!==================================================================================================================================!
SUBROUTINE InitSplitDG()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: SplitDG
#ifdef DEBUG
REAL,DIMENSION(PP_nVar    ) :: U         ! dummy variables, only to suppress compiler warnings
REAL,DIMENSION(PP_nVarPrim) :: UPrim     ! dummy variables, only to suppress compiler warnings
REAL,DIMENSION(PP_nVar    ) :: f,g,h!,F   ! dummy variables, only to suppress compiler warnings
REAL,DIMENSION(PP_2Var    ) :: U_LL,U_RR ! dummy variables, only to suppress compiler warnings
#endif
!==================================================================================================================================
! check if Gauss-Lobatto-Pointset is beeing used
#if (PP_NodeType==1)
CALL CollectiveStop(__STAMP__,&
  'Wrong Pointset: Gauss-Lobatto-Points are mandatory for using SplitDG !')
#endif
! set pointers
SplitDG = GETINTFROMSTR('SplitDG')
SELECT CASE(SplitDG)
CASE(PRM_SPLITDG_SD)
  SplitDGVolume_pointer  => SplitVolumeFluxSD
  SplitDGSurface_pointer => SplitSurfaceFluxSD
  SplitIndicator = 0
CASE(PRM_SPLITDG_MO)
  SplitDGVolume_pointer  => SplitVolumeFluxMO
  SplitDGSurface_pointer => SplitSurfaceFluxMO
  SplitIndicator = 1
CASE(PRM_SPLITDG_DU)
  SplitDGVolume_pointer  => SplitVolumeFluxDU
  SplitDGSurface_pointer => SplitSurfaceFluxDU
  SplitIndicator = 2
CASE(PRM_SPLITDG_KG)
  SplitDGVolume_pointer  => SplitVolumeFluxKG
  SplitDGSurface_pointer => SplitSurfaceFluxKG
  SplitIndicator = 3
CASE(PRM_SPLITDG_PI)
  SplitDGVolume_pointer  => SplitVolumeFluxPI
  SplitDGSurface_pointer => SplitSurfaceFluxPI
  SplitIndicator = 4
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'SplitDG formulation not defined!')
END SELECT

#ifdef DEBUG
! ===============================================================================
! Following dummy calls do suppress compiler warnings of unused Riemann-functions
! ===============================================================================
IF (0.EQ.1) THEN
  U=1. ;  UPrim=1. ;   U_LL=1. ;   U_RR=1.
  CALL SplitDGVolume_pointer  (U,UPrim,U,UPrim,f,g,h)
  CALL SplitDGSurface_pointer (U_LL,U_RR,F)
  CALL SplitVolumeFluxSD     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxSD    (U_LL,U_RR,F)
  CALL SplitVolumeFluxMO     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxMO    (U_LL,U_RR,F)
  CALL SplitVolumeFluxDU     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxDU    (U_LL,U_RR,F)
  CALL SplitVolumeFluxKG     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxKG    (U_LL,U_RR,F)
  CALL SplitVolumeFluxPI     (U,UPrim,U,UPrim,f,g,h)
  CALL SplitSurfaceFluxPI    (U_LL,U_RR,F)
END IF
#endif
END SUBROUTINE InitSplitDG

!==================================================================================================================================
!> Computes the Split-Flux retaining the standart NS-Equations
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
PURE SUBROUTINE SplitVolumeFluxSD(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,M ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: EpRef,Ep ! auxilery variable for (rho*e+p)
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde,hTilde ! fluxes in physical space
!==================================================================================================================================
! compute auxilery variables
EpRef = URef(5) + UPrimRef(5)
Ep    = U(5) + UPrim(5)

! local Euler fluxes x-direction
  fTilde(1) = (URef(2) + U(2))                                           ! {rho*u}
  fTilde(2) = (URef(2)*UPrimRef(2)+UPrimRef(5) + U(2)*UPrim(2)+UPrim(5)) ! {rho*u²}+{p}
  fTilde(3) = (URef(2)*UPrimRef(3) + U(2)*UPrim(3))                      ! {rho*u*v}
  fTilde(4) = (URef(2)*UPrimRef(4) + U(2)*UPrim(4))                      ! {rho*u*w}
  fTilde(5) = (EpRef*UPrimRef(2) + Ep*UPrim(2))                          ! {(rho*e+p)*u}
! local Euler fluxes y-direction
  gTilde(1) = (URef(3) + U(3))                                           ! {rho*v}
  gTilde(2) = (URef(2)*UPrimRef(3) + U(2)*UPrim(3))                      ! {rho*u*v}
  gTilde(3) = (URef(3)*UPrimRef(3)+UPrimRef(5) + U(3)*UPrim(3)+UPrim(5)) ! {rho*v²}+{p}
  gTilde(4) = (URef(3)*UPrimRef(4) + U(3)*UPrim(4))                      ! {rho*v*w}
  gTilde(5) = (EpRef*UPrimRef(3) + Ep*UPrim(3))                          ! {(rho*e+p)*v}
! local Euler fluxes z-direction
  hTilde(1) = (URef(4) + U(4))                                           ! {rho*w}
  hTilde(2) = (URef(2)*UPrimRef(4) + U(2)*UPrim(4))                      ! {rho*u*w}
  hTilde(3) = (URef(3)*UPrimRef(4) + U(3)*UPrim(4))                      ! {rho*v*w}
  hTilde(4) = (URef(4)*UPrimRef(4)+UPrimRef(5) + U(4)*UPrim(4)+UPrim(5)) ! {rho*v²+p}
  hTilde(5) = (EpRef*UPrimRef(4) + Ep*UPrim(4))                          ! {(rho*e+p)*w}

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
          0.5*(MRef(2)+M(2))*gTilde(:) + &
          0.5*(MRef(3)+M(3))*hTilde(:)

END SUBROUTINE SplitVolumeFluxSD

!==================================================================================================================================
!> Computes the surface flux for the split formulation retaining the standart NS-Equations
!==================================================================================================================================
PURE SUBROUTINE SplitSurfaceFluxSD(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

F(1)= 0.5*(U_LL(MOM1)+U_RR(MOM1))                                                ! {rho*u}
F(2)= 0.5*(U_LL(MOM1)*U_LL(VEL1)+U_LL(PRES)+U_RR(MOM1)*U_RR(VEL1)+U_RR(PRES))    ! {rho*u²}+{p}
F(3)= 0.5*(U_LL(MOM1)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2))                          ! {rho*u*v}
F(4)= 0.5*(U_LL(MOM1)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3))                          ! {rho*u*w}
F(5)= 0.5*((U_LL(ENER)+U_LL(PRES))*U_LL(VEL1)+(U_RR(ENER)+U_RR(PRES))*U_RR(VEL1))! {(rho*e+p)*u}

END SUBROUTINE SplitSurfaceFluxSD

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Ducros
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
PURE SUBROUTINE SplitVolumeFluxDU(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,M ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde,hTilde ! flux in physical space
!==================================================================================================================================

! local Euler fluxes x-direction
  fTilde(1) = 0.5*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))                          ! {rho}*{u}
  fTilde(2) = 0.5*(URef(2)+U(2))*(UPrimRef(2)+UPrim(2)) + (UPrimRef(5)+UPrim(5)) ! {rho*u}*{u}+{p}
  fTilde(3) = 0.5*(URef(3)+U(3))*(UPrimRef(2)+UPrim(2))                          ! {rho*v}*{u}
  fTilde(4) = 0.5*(URef(4)+U(4))*(UPrimRef(2)+UPrim(2))                          ! {rho*w}*{u}
  fTilde(5) = 0.5*(URef(5)+U(5)+UPrimRef(5)+UPrim(5))*(UPrimRef(2)+UPrim(2))     ! ({rho*e}+{p})*{u}
! local Euler fluxes y-direction
  gTilde(1) = 0.5*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))                          ! {rho}*{v}
  gTilde(2) = 0.5*(URef(2)+U(2))*(UPrimRef(3)+UPrim(3))                          ! {rho*u}*{v}
  gTilde(3) = 0.5*(URef(3)+U(3))*(UPrimRef(3)+UPrim(3)) + (UPrimRef(5)+UPrim(5)) ! {rho*v}*{v}+{p}
  gTilde(4) = 0.5*(URef(4)+U(4))*(UPrimRef(3)+UPrim(3))                          ! {rho*w}*{v}
  gTilde(5) = 0.5*(URef(5)+U(5)+UPrimRef(5)+UPrim(5))*(UPrimRef(3)+UPrim(3))     ! ({rho*e}+{p})*{v}
! local Euler fluxes z-direction
  hTilde(1) = 0.5*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))                          ! {rho}*{w}
  hTilde(2) = 0.5*(URef(2)+U(2))*(UPrimRef(4)+UPrim(4))                          ! {rho*u}*{w}
  hTilde(3) = 0.5*(URef(3)+U(3))*(UPrimRef(4)+UPrim(4))                          ! {rho*v}*{w}
  hTilde(4) = 0.5*(URef(4)+U(4))*(UPrimRef(4)+UPrim(4)) + (UPrimRef(5)+UPrim(5)) ! {rho*w}*{w}+{p}
  hTilde(5) = 0.5*(URef(5)+U(5)+UPrimRef(5)+UPrim(5))*(UPrimRef(4)+UPrim(4))     ! ({rho*e}+{p})*{w}

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
          0.5*(MRef(2)+M(2))*gTilde(:) + &
          0.5*(MRef(3)+M(3))*hTilde(:)

END SUBROUTINE SplitVolumeFluxDU

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Ducros
!==================================================================================================================================
PURE SUBROUTINE SplitSurfaceFluxDU(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
F(1)= 0.25*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho}*{u}
F(2)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL1)+U_RR(VEL1)) + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho*u}*{u}+{p}
F(3)= 0.25*(U_LL(MOM2)+U_RR(MOM2))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho*v}*{u}
F(4)= 0.25*(U_LL(MOM3)+U_RR(MOM3))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho*w}*{u}
F(5)= 0.25*(U_LL(ENER)+U_RR(ENER)+U_LL(PRES)+U_RR(PRES))*(U_LL(VEL1)+U_RR(VEL1))         ! ({rho*e}+{p})*{u}

END SUBROUTINE SplitSurfaceFluxDU

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Kennedy and Gruber
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
PURE SUBROUTINE SplitVolumeFluxKG(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,M ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: e,eRef ! auxilery variables for the specific energy
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde,hTilde ! flux in reference space
!==================================================================================================================================

! specific energy
eRef = URef(5)/URef(1)
e    = U(5)/U(1)

! local Euler fluxes x-direction
  fTilde(1) = 0.5* (URef(1)+U(1))*(UPrimRef(2)+UPrim(2))                             ! {rho}*{u}
  fTilde(2) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{u}²+{p}
  fTilde(3) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(UPrimRef(3)+UPrim(3))      ! {rho}*{u}*{v}
  fTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(UPrimRef(4)+UPrim(4))      ! {rho}*{u}*{w}
  fTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(eRef+e) + &
              0.5* (UPrimRef(5)+UPrim(5))*(UPrimRef(2)+UPrim(2))                     ! {rho}*{e}*{u}+{p}*{u}
! local Euler fluxes y-direction
  gTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))                             ! {rho}*{v}
  gTilde(2) = fTilde(3)                                                              ! {rho}*{v}*{u}
  gTilde(3) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{v}²+{p}
  gTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))*(UPrimRef(4)+UPrim(4))      ! {rho}*{v}*{w}
  gTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))*(eRef+e) + &
              0.5* (UPrimRef(5)+UPrim(5))*(UPrimRef(3)+UPrim(3))                     ! {rho}*{e}*{v}+{p}*{v}
! local Euler fluxes z-direction
  hTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))                             ! {rho}*{w}
  hTilde(2) = fTilde(4)                                                              ! {rho}*{w}*{u}
  hTilde(3) = gTilde(4)                                                              ! {rho}*{w}*{v}
  hTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{w}²+{p}
  hTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))*(eRef+e) + &
              0.5 *(UPrimRef(5)+UPrim(5))*(UPrimRef(4)+UPrim(4))                     ! {rho}*{e}*{w}+{p}*{w}

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
          0.5*(MRef(2)+M(2))*gTilde(:) + &
          0.5*(MRef(3)+M(3))*hTilde(:)

END SUBROUTINE SplitVolumeFluxKG

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Kennedy and Gruber
!==================================================================================================================================
PURE SUBROUTINE SplitSurfaceFluxKG(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: e_LL,e_RR ! auxilery variables for the specific energy
!==================================================================================================================================
! specific energy
e_LL = U_LL(ENER)/U_LL(DENS)
e_RR = U_RR(ENER)/U_RR(DENS)
!compute flux
F(1)= 0.25* (U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))                                  ! {rho}*{u}
F(2)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))**2 + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho}*{u}²+{p}
F(3)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL2)+U_RR(VEL2))          ! {rho}*{u}*{v}
F(4)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL3)+U_RR(VEL3))          ! {rho}*{u}*{w}
F(5)= 0.125*(U_LL(DENS)+U_RR(DENS))*(e_LL+e_RR)*(U_LL(VEL1)+U_RR(VEL1)) + &
      0.25 *(U_LL(PRES)+U_RR(PRES))*(U_LL(VEL1)+U_RR(VEL1))                                  ! {rho}*{e}*{u}+{p}*{u}

END SUBROUTINE SplitSurfaceFluxKG

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Morinishi
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
PURE SUBROUTINE SplitVolumeFluxMO(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,M ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: EpRef,Ep ! auxilery variable for inner energy + pressure
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde,hTilde ! flux in physical space
!==================================================================================================================================

! internal energy + pressure
EpRef = URef(5)-0.5*URef(1)*(UPrimRef(2)**2+UPrimRef(3)**2+UPrimRef(4)**2)+UPrimRef(5)
Ep    = U(5)-0.5*U(1)*(UPrim(2)**2+UPrim(3)**2+UPrim(4)**2)+UPrim(5)

! local Euler fluxes x-direction
  fTilde(1) =     (URef(2)+U(2))                                                 ! {rho*u}
  fTilde(2) = 0.5*(URef(2)+U(2))*(UPrimRef(2)+UPrim(2)) + (UPrimRef(5)+UPrim(5)) ! {rho*u}*{u}+{p}
  fTilde(3) = 0.5*(Uref(2)+U(2))*(UPrimRef(3)+UPrim(3))                          ! {rho*u}*{v}
  fTilde(4) = 0.5*(URef(2)+U(2))*(UPrimRef(4)+UPrim(4))                          ! {rho*u}*{w}
  fTilde(5) = (EpRef*UPrimRef(2)+Ep*UPrim(2)) + &                                !{(rho*e_int+p)*u} +
              0.5*(URef(2)*UPrimRef(2)+U(2)*UPrim(2))*(UPrimRef(2)+UPrim(2)) + & !{rho*u²}*{u} +
              0.5*(URef(2)*UPrimRef(3)+U(2)*UPrim(3))*(UPrimRef(3)+UPrim(3)) + & !{rho*u*v}*{v} +
              0.5*(URef(2)*UPrimRef(4)+U(2)*UPrim(4))*(UPrimRef(4)+UPrim(4)) - & !{rho*u*w}*{w} -
              0.5*(URef(2)*UPrimRef(2)*UPrimRef(2)+U(2)*UPrim(2)*UPrim(2))   - & !1/2*({rho*u³} +
              0.5*(URef(2)*UPrimRef(3)*UPrimRef(3)+U(2)*UPrim(3)*UPrim(3))   - & !{rho*u*v²} +
              0.5*(URef(2)*UPrimRef(4)*UPrimRef(4)+U(2)*UPrim(4)*UPrim(4))       !{rho*u*w²})
! local Euler fluxes y-direction
  gTilde(1) =     (URef(3)+U(3))                                                 ! {rho*v}
  gTilde(2) = 0.5*(Uref(3)+U(3))*(UPrimRef(2)+UPrim(2))                          ! {rho*v}*{u}
  gTilde(3) = 0.5*(URef(3)+U(3))*(UPrimRef(3)+UPrim(3)) + (UPrimRef(5)+UPrim(5)) ! {rho*v}*{v}+{p}
  gTilde(4) = 0.5*(URef(3)+U(3))*(UPrimRef(4)+UPrim(4))                          ! {rho*v}*{w}
  gTilde(5) = (EpRef*UPrimRef(3)+Ep*UPrim(3)) + &                                !{(rho*e_int+p)*v} +
              0.5*(URef(3)*UPrimRef(2)+U(3)*UPrim(2))*(UPrimRef(2)+UPrim(2)) + & !{rho*v*u}*{u} +
              0.5*(URef(3)*UPrimRef(3)+U(3)*UPrim(3))*(UPrimRef(3)+UPrim(3)) + & !{rho*v²}*{v} +
              0.5*(URef(3)*UPrimRef(4)+U(3)*UPrim(4))*(UPrimRef(4)+UPrim(4)) - & !{rho*v*w}*{w} -
              0.5*(URef(3)*UPrimRef(2)*UPrimRef(2)+U(3)*UPrim(2)*UPrim(2))   - & !1/2*({rho*v*u²} +
              0.5*(URef(3)*UPrimRef(3)*UPrimRef(3)+U(3)*UPrim(3)*UPrim(3))   - & !{rho*v³} +
              0.5*(URef(3)*UPrimRef(4)*UPrimRef(4)+U(3)*UPrim(4)*UPrim(4))       !{rho*v*w²})
! local Euler fluxes z-direction
  hTilde(1) =     (URef(4)+U(4))                                                 ! {rho*w}
  hTilde(2) = 0.5*(Uref(4)+U(4))*(UPrimRef(2)+UPrim(2))                          ! {rho*w}*{u}
  hTilde(3) = 0.5*(URef(4)+U(4))*(UPrimRef(3)+UPrim(3))                          ! {rho*w}*{v}
  hTilde(4) = 0.5*(URef(4)+U(4))*(UPrimRef(4)+UPrim(4)) + (UPrimRef(5)+UPrim(5)) ! {rho*w}*{w}+{p}
  hTilde(5) = (EpRef*UPrimRef(4)+Ep*UPrim(4)) + &                                !{(rho*e_int+p)*w} +
              0.5*(URef(4)*UPrimRef(2)+U(4)*UPrim(2))*(UPrimRef(2)+UPrim(2)) + & !{rho*w*u}*{u} +
              0.5*(URef(4)*UPrimRef(3)+U(4)*UPrim(3))*(UPrimRef(3)+UPrim(3)) + & !{rho*w*v}*{v} +
              0.5*(URef(4)*UPrimRef(4)+U(4)*UPrim(4))*(UPrimRef(4)+UPrim(4)) - & !{rho*w²}*{w} -
              0.5*(URef(4)*UPrimRef(2)*UPrimRef(2)+U(4)*UPrim(2)*UPrim(2))   - & !1/2*({rho*w*u²} +
              0.5*(URef(4)*UPrimRef(3)*UPrimRef(3)+U(4)*UPrim(3)*UPrim(3))   - & !{rho*w*v²} +
              0.5*(URef(4)*UPrimRef(4)*UPrimRef(4)+U(4)*UPrim(4)*UPrim(4))       !{rho*w³})

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
          0.5*(MRef(2)+M(2))*gTilde(:) + &
          0.5*(MRef(3)+M(3))*hTilde(:)

END SUBROUTINE SplitVolumeFluxMO

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Morinishi
!==================================================================================================================================
PURE SUBROUTINE SplitSurfaceFluxMO(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: Ep_LL,Ep_RR
!==================================================================================================================================
! internal energy + pressure
Ep_LL = U_LL(ENER)-0.5*U_LL(DENS)*(U_LL(VEL1)**2+U_LL(VEL2)**2+U_LL(VEL3)**2)+U_LL(PRES)
Ep_RR = U_RR(ENER)-0.5*U_RR(DENS)*(U_RR(VEL1)**2+U_RR(VEL2)**2+U_RR(VEL3)**2)+U_RR(PRES)

! compute flux
F(1)= 0.5 *(U_LL(MOM1)+U_RR(MOM1))                                                       ! {rho*u}
F(2)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL1)+U_RR(VEL1)) + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho*u}*{u}+{p}
F(3)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL2)+U_RR(VEL2))                               ! {rho*u}*{v}
F(4)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL3)+U_RR(VEL3))                               ! {rho*u}*{w}
F(5)= 0.5 *(Ep_LL*U_LL(VEL1)+Ep_RR*U_RR(VEL1)) +  &                                      !{(rho*e_int+p)*u} +
      0.25*(U_LL(MOM1)*U_LL(VEL1)+U_RR(MOM1)*U_RR(VEL1))*(U_LL(VEL1)+U_RR(VEL1)) + &     !{rho*u²}*{u} +
      0.25*(U_LL(MOM1)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2))*(U_LL(VEL2)+U_RR(VEL2)) + &     !{rho*u*v}*{v} +
      0.25*(U_LL(MOM1)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3))*(U_LL(VEL3)+U_RR(VEL3)) - &     !{rho*u*w}*{w} -
      0.25*(U_LL(MOM1)*U_LL(VEL1)*U_LL(VEL1)+U_RR(MOM1)*U_RR(VEL1)*U_RR(VEL1)) - &       !1/2*({rho*u³} +
      0.25*(U_LL(MOM1)*U_LL(VEL2)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2)*U_RR(VEL2)) - &       !{rho*u*v²} +
      0.25*(U_LL(MOM1)*U_LL(VEL3)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3)*U_RR(VEL3))           !{rho*u*w²})

END SUBROUTINE SplitSurfaceFluxMO

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Pirozzoli
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
PURE SUBROUTINE SplitVolumeFluxPI(URef,UPrimRef,U,UPrim,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U ! conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim ! primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,M ! metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux ! flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: e,eRef ! auxilery variables for the specific enthalpy
REAL,DIMENSION(PP_nVar    )             :: fTilde,gTilde,hTilde ! flux in physical space
!==================================================================================================================================

! specific enthalpy
eRef = (URef(5)+UPrimRef(5))/URef(1)
e    = (U(5)+UPrim(5))/U(1)

! local Euler fluxes x-direction
  fTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))                             ! {rho}*{u}
  fTilde(2) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{u}²+{p}
  fTilde(3) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(UPrimRef(3)+UPrim(3))      ! {rho}*{u}*{v}
  fTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(UPrimRef(4)+UPrim(4))      ! {rho}*{u}*{w}
  fTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(eRef+e)                    ! {rho}*{h}*{u}
! local Euler fluxes y-direction
  gTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))                             ! {rho}*{v}
  gTilde(2) = fTilde(3)                                                              ! {rho}*{v}*{u}
  gTilde(3) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{v}²+{p}
  gTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))*(UPrimRef(4)+UPrim(4))      ! {rho}*{v}*{w}
  gTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))*(eRef+e)                    ! {rho}*{h}*{v}
! local Euler fluxes z-direction
  hTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))                             ! {rho}*{w}
  hTilde(2) = fTilde(4)                                                              ! {rho}*{w}*{u}
  hTilde(3) = gTilde(4)                                                              ! {rho}*{w}*{v}
  hTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{w}²+{p}
  hTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))*(eRef+e)                    ! {rho}*{h}*{w}

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
          0.5*(MRef(2)+M(2))*gTilde(:) + &
          0.5*(MRef(3)+M(3))*hTilde(:)

END SUBROUTINE SplitVolumeFluxPI

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Pirozzoli
!==================================================================================================================================
PURE SUBROUTINE SplitSurfaceFluxPI(U_LL,U_RR,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR ! variables at the left-/right-Surfaces
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F ! resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: e_LL,e_RR ! auxilery variables for the specific energy
!==================================================================================================================================
! specific energy
e_LL = (U_LL(ENER)+U_LL(PRES))/U_LL(DENS)
e_RR = (U_RR(ENER)+U_RR(PRES))/U_RR(DENS)
!compute flux
F(1)= 0.25* (U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))                                  ! {rho}*{u}
F(2)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))**2 + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho}*{u}²+{p}
F(3)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL2)+U_RR(VEL2))          ! {rho}*{u}*{v}
F(4)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL3)+U_RR(VEL3))          ! {rho}*{u}*{w}
F(5)= 0.125*(U_LL(DENS)+U_RR(DENS))*(e_LL+e_RR)*(U_LL(VEL1)+U_RR(VEL1))                      ! {rho}*{h}*{u}

END SUBROUTINE SplitSurfaceFluxPI

END MODULE MOD_SplitFlux
