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
!> Contains the constant Advection Velocity Vector used for the linear scalar advection equation
!==================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL           :: doCalcSource             !< Swith to calculate a source term or not, automatically set by calcsource itself
REAL              :: AdvVel(3)                !< Advection velocity
REAL              :: DiffC                    !< Diffusion constant
INTEGER           :: IniExactFunc             !< Number of exact function used for initialization
REAL,ALLOCATABLE  :: RefStatePrim(:,:)        !< Primite reference state
REAL,ALLOCATABLE  :: RefStateCons(:,:)        !< Conservative reference state

! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:)       !< Buffer array for BC data
INTEGER,ALLOCATABLE  :: nBCByType(:)          !< Number of sides for each boundary
INTEGER,ALLOCATABLE  :: BCSideID(:,:)         !< SideIDs for BC types

CHARACTER(LEN=255),PARAMETER :: StrVarNames(1)=(/'Solution'/) !< Variable names for output

LOGICAL           :: EquationInitIsDone=.FALSE.!< Init switch  
!==================================================================================================================================
END MODULE MOD_Equation_Vars
