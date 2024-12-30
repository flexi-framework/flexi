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
!> Contains global variables used by the preconditioner modules.
!===================================================================================================================================
MODULE MOD_Precond_Vars
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE        :: Ploc(:,:)         !< Non-inverted preconditioner matrix
REAL,ALLOCATABLE        :: Ploc1(:,:)        !< Additional non-inverted matrix for comparision between AD and FD
REAL,ALLOCATABLE,TARGET :: invP(:,:,:)       !< inverse of block Jacobian for each element (1:nDOF_elem,1:nDOFelem,1:nElems)
INTEGER                 :: PrecondType       !< Type of preconditioner, 0: non, 1: analytical, 2: finite difference, 3: compare
INTEGER(DP)             :: PrecondIter       !< Defines how often preconditioner is built
INTEGER                 :: DebugMatrix       !< Set to >0 to get output of preconditioner matrix
LOGICAL                 :: HyperbolicPrecond !< Set TRUE to use the preconditioner of the hyperbolic system, even for
                                             !< parabolic computations
INTEGER                 :: SolveSystem       !< Chose the inversion method (0: exact LU, 1: inexact ILU(0))
#if PARABOLIC
LOGICAL                 :: NoFillIn          !< When building the parabolic precond, use the hyperbolic sparsity pattern
#endif
LOGICAL                 :: DoDisplayPrecond  !< Display building time of preconditioner

LOGICAL                 :: PrecondInitIsDone
!===================================================================================================================================
END MODULE MOD_Precond_Vars
