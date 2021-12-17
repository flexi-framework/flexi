!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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

!===================================================================================================================================
!> Contains global variables provided by the HIT_Analyze routines
!===================================================================================================================================
MODULE MOD_HIT_Analyze_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                  :: mu0               !> Dynamic viscosity
INTEGER               :: N_Filter          !> Cutoff wavenumber for computation of spectrum
INTEGER               :: nVar_HDF5         !> Number of solution variables in statefile
INTEGER               :: N_HDF5            !> Polynomial degree in statefile
INTEGER               :: nElems_HDF5       !> Number of elements in statefile
REAL                  :: Time_HDF5         !> Simulation time in statefile
CHARACTER(LEN=255)    :: NodeType_HDF5     !> Nodetype in statefile
CHARACTER(LEN=255)    :: ProjectName_HDF5  !> Projectname statefile
!===================================================================================================================================
END MODULE MOD_HIT_Analyze_Vars
