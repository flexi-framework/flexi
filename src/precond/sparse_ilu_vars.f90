!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
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

!===================================================================================================================================
! Contains global variables used for the ILU(0) factorization and application
!===================================================================================================================================
MODULE MOD_SparseILU_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                     :: nMTriangle
REAL                                        :: epsZero                        !< machine precision
INTEGER,ALLOCATABLE,DIMENSION(:)            :: nUNonZeros,nLNonZeros          !< non-zero entries of L and U matrixes
REAL,ALLOCATABLE,DIMENSION(:,:)             :: Dinv

TYPE tILU                                                                     !< ILU for each element
 REAL,ALLOCATABLE,DIMENSION(:)              :: Entry
 INTEGER,ALLOCATABLE,DIMENSION(:)           :: IEntry,JEntry
END TYPE
TYPE(tILU), ALLOCATABLE                     :: IL(:)                          !< Incomplete Lower matrix
TYPE(tILU), ALLOCATABLE                     :: IU(:)                          !< Incomplete Upper matrix
!===================================================================================================================================
END MODULE MOD_SparseILU_Vars
