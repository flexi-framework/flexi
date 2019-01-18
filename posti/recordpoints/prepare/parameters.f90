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
!===================================================================================================================================
!> Contains the global parameters for the prepare recordpoints tool
!===================================================================================================================================
MODULE MOD_Parameters
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: NSuper                   !< Supersampling polynomial degree, used for initial guess in Newton
REAL                          :: maxTol                   !< Max overshoot in param coords (1+maxTol)
LOGICAL                       :: doVisuRP                 !< Set to True to visualize the record point set for paraview

INTEGER                       :: nCoords                  !< Number of coordinates (physical plus local coordinates)
CHARACTER(LEN=255),ALLOCATABLE:: CoordNames(:)            !< Names of those coordinates
LOGICAL                       :: OutputInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_Parameters
