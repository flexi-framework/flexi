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
!> Global variables needed by the unittest routines
!==================================================================================================================================
MODULE MOD_Unittest_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER     :: NRef=9        !< Polynomial degree used for reference calculations 
INTEGER,PARAMETER     :: nElemsRef=1   !< Number of elements in reference element
#if PP_dim == 3
INTEGER,PARAMETER     :: nSidesRef=6   !< Number of sides in reference element
INTEGER,PARAMETER     :: NRefZ=NRef    !< Polynomial degree used for reference calculations in third dimension
#else
INTEGER,PARAMETER     :: nSidesRef=4   !< Number of sides in reference element
INTEGER,PARAMETER     :: NRefZ=0       !< Polynomial degree used for reference calculations in third dimension
#endif
!==================================================================================================================================

END MODULE MOD_Unittest_Vars
