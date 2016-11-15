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

!==================================================================================================================================
!> Contains the (physical) parameters needed for the Navier Stokes calculation
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
INTEGER           :: IniExactFunc
INTEGER           :: IniRefState       !< RefState for initialization (case IniExactFunc=1 only)
INTEGER           :: nRefState         !< number of refstates defined in parameter file
REAL,ALLOCATABLE  :: RefStatePrim(:,:) !< refstates in primitive variables (as read from ini file)
REAL,ALLOCATABLE  :: RefStateCons(:,:) !< refstates in conservative variables
CHARACTER(LEN=255):: BCStateFile       !< file containing the reference solution on the boundary to be used as BC

! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:) !< array with precomputed BC values (conservative)
REAL,ALLOCATABLE     :: BCDataPrim(:,:,:,:) !< array with precomputed BC values (primitive)
INTEGER,ALLOCATABLE  :: nBCByType(:)   !< number of sides with specific BC type
INTEGER,ALLOCATABLE  :: BCSideID(:,:)

REAL                 :: s43,s23

CHARACTER(LEN=255),DIMENSION(5),PARAMETER :: StrVarNames =&
  (/ CHARACTER(LEN=255) :: 'Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'/) !< conservative variable names
CHARACTER(LEN=255),DIMENSION(6),PARAMETER :: StrVarNamesPrim=&
  (/ CHARACTER(LEN=255) :: 'Density','VelocityX','VelocityY','VelocityZ','Pressure','Temperature'/) !< primitive variable names

LOGICAL           :: EquationInitIsDone=.FALSE.
!==================================================================================================================================

END MODULE MOD_Equation_Vars
