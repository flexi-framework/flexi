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
!> Contains the (physical) parameters needed for the RANS SA calculation
!==================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL           :: doCalcSource      !< automatically set by calcsource itself
INTEGER           :: IniExactFunc      !< number identifying the used exact function
INTEGER           :: IniRefState       !< RefState for initialization (case IniExactFunc=1 only)
INTEGER           :: nRefState         !< number of refstates defined in parameter file
REAL,ALLOCATABLE  :: RefStatePrim(:,:) !< refstates in primitive variables (as read from ini file)
REAL,ALLOCATABLE  :: RefStateCons(:,:) !< refstates in conservative variables
CHARACTER(LEN=255):: BCStateFile       !< file containing the reference solution on the boundary to be used as BC

! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:) !< array with precomputed BC values (conservative)
REAL,ALLOCATABLE     :: BCDataPrim(:,:,:,:) !< array with precomputed BC values (primitive)
INTEGER,ALLOCATABLE  :: nBCByType(:)   !< number of sides with specific BC type
INTEGER,ALLOCATABLE  :: BCSideID(:,:)  !< array storing side IDs of sides with different BCs

REAL                 :: s43            !< precomputed 4./3.
REAL                 :: s23            !< precomputed 2./3.

! SA-specific variables and parameters
REAL              :: PrTurb            !< Turbulent Prandtl number
REAL, PARAMETER   :: sigma = 2./3.     !< Constant in SA equation
REAL, PARAMETER   :: SAKappa = 0.41    !< Constant in SA equation
REAL, PARAMETER   :: cv1 = 7.1         !< Constant in SA equation
REAL, PARAMETER   :: cv2 = 0.7         !< Constant in SA equation
REAL, PARAMETER   :: cv3 = 0.9         !< Constant in SA equation
REAL, PARAMETER   :: cb1 = 0.1355      !< Constant in SA equation
REAL, PARAMETER   :: cb2 = 0.622       !< Constant in SA equation
REAL, PARAMETER   :: cw1 = (cb1/(SAKappa**2)) + ((1+cb2)/sigma) !< Constant in SA equation
REAL, PARAMETER   :: cw2 = 0.3         !< Constant in SA equation
REAL, PARAMETER   :: cw3 = 2.          !< Constant in SA equation
REAL, PARAMETER   :: cn1 = 16.         !< Constant in SA equation
REAL, PARAMETER   :: rLim = 10.        !< Constant in SA equation
REAL, ALLOCATABLE :: SAd(:,:,:,:,:)    !< Distance from closest wall

LOGICAL           :: includeTrip       !< Swith to include trip terms in SA equation
REAL, PARAMETER   :: ct1 = 1.0         !< Constant in SA equation, specifically for trip
REAL, PARAMETER   :: ct2 = 2.0         !< Constant in SA equation, specifically for trip
REAL, PARAMETER   :: ct3 = 1.2         !< Constant in SA equation, specifically for trip
REAL, PARAMETER   :: ct4 = 0.5         !< Constant in SA equation, specifically for trip
REAL              :: TripX(2)          !< 2D (x,y) location of trip
REAL              :: dxT               !< Mesh width at trip location
REAL              :: omegaT            !< Vorticity at trip
REAL, ALLOCATABLE :: SAdt(:,:,:,:,:)   !< Distance from trip point
INTEGER           :: tripPQ(2)         !< (p,q) index of trip boundary point
INTEGER           :: tripSideId        !< SideID of trip boundary point
LOGICAL           :: tripOnProc        !< Is the trip located on this proc?
#if USE_MPI
INTEGER           :: tripRoot          !< process ID of proc that holds the trip
#endif
REAL,ALLOCATABLE  :: SADebug(:,:,:,:,:) !< Array used for debug output, include Prduction, Destruction, Trip and Diffusion terms
LOGICAL           :: doSADebug          !< Turn debug output on or off (default: off)




CHARACTER(LEN=255),DIMENSION(6),PARAMETER :: StrVarNames =&
  (/ CHARACTER(LEN=255) :: 'Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity','muTilde'/) !< conservative variable names
CHARACTER(LEN=255),DIMENSION(7),PARAMETER :: StrVarNamesPrim=&
  (/ CHARACTER(LEN=255) :: 'Density','VelocityX','VelocityY','VelocityZ','Pressure','Temperature','nuTilde'/) !< primitive variable names

LOGICAL           :: EquationInitIsDone=.FALSE.
!==================================================================================================================================

#ifdef PARABOLIC
INTERFACE fv1
  MODULE PROCEDURE fv1
END INTERFACE
INTERFACE fv2
  MODULE PROCEDURE fv2
END INTERFACE
INTERFACE fw
  MODULE PROCEDURE fw
END INTERFACE
INTERFACE STilde
  MODULE PROCEDURE STilde
END INTERFACE

CONTAINS

PPURE FUNCTION fv1(chi)
!===================================================================================================================================
!> Function fv1 of the Spalart-Allmaras Turbulence model
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: chi !< muTilde/mu
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                           :: fv1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Additional Limiter introduced
IF(CHI.LT.0.)THEN
  fv1=0.
ELSE
  fv1 = chi**3/(chi**3+cv1**3)
END IF

END FUNCTION fv1

PPURE FUNCTION fv2(chi)
!===================================================================================================================================
!> Function fv2 of the Spalart-Allmaras Turbulence model
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: chi !< muTilde/mu
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                           :: fv2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: fv1_loc
!===================================================================================================================================

fv1_loc = chi**3/(chi**3+cv1**3)
fv2 = 1. - (chi/(1.+chi*fv1_loc))

END FUNCTION fv2

PPURE FUNCTION fw(nuTilde, STilde, d)
!===================================================================================================================================
!> Function fw of the negative Spalart-Allmaras Turbulence model
!> See "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Tubulence Model"
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: nuTilde   !< SA kinematic viscosity
REAL,INTENT(IN)                :: STilde    !< modified vorticity
REAL,INTENT(IN)                :: d         !< distance from closest wall
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                           :: fw
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: r         ! auxiliary function
REAL                           :: g         ! auxiliary function
!===================================================================================================================================

IF(nuTilde.GE.0.)THEN
  r = MIN(nuTilde/(STilde*(SAKappa**2)*(d**2)),rLim)
  g = r + cw2*((r**6)-r)

  fw = g*(((1+(cw3**6))/((g**6)+(cw3**6)))**(1./6.))
ELSE
  fw = -1.
END IF

END FUNCTION fw

PPURE FUNCTION fn(chi)
!===================================================================================================================================
!> Function fn of the negative Spalart-Allmaras Turbulence model
!> See "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Tubulence Model"
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: chi       !< muTilde/mu
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                           :: fn
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

fn = (cn1+(chi**3))/(cn1-(chi**3))

END FUNCTION fn

PPURE FUNCTION STilde(nuTilde, d, chi, S)
!===================================================================================================================================
! Modified vorticity of the modifed Spalart-Allmaras Turbulence model
! See "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Tubulence Model"
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: nuTilde   !< SA kinematic viscosity
REAL,INTENT(IN)                :: S         !< vorticity
REAL,INTENT(IN)                :: d         !< distance from closest wall
REAL,INTENT(IN)                :: chi       !< muTilde/mu
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                           :: STilde
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: SBar      ! auxiliary modified vorticity
!===================================================================================================================================

SBar = nuTilde/((SAKappa**2)*(d**2))*fv2(chi)

IF(SBar.GE.(-cv2*S)) THEN
  STilde = S + SBar
ELSE
  STilde = S + (S*((cv2**2)*S+cv3*SBar))/((cv3-2*cv2)*S-SBar)
END IF

END FUNCTION
#endif /*PARABOLIC*/

END MODULE MOD_Equation_Vars
