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

!==================================================================================================================================
!> Variables used for mortars: mortar interpolation and projection matrices
!==================================================================================================================================
MODULE MOD_Mortar_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE,TARGET :: M_0_1(:,:)    !< 1D-Mortar Operator: interpolation full interval 0: [-1,1] to left interval 1: [-1,0]
REAL,ALLOCATABLE,TARGET :: M_0_2(:,:)    !< 1D-Mortar Operator: interpolation full interval 0: [-1,1] to right intervall 2: [0,1]
REAL,ALLOCATABLE,TARGET :: M_1_0(:,:)    !< 1D-Mortar Operator: projection left interval 1: [-1,0]  to full intervall 0: [-1,1]
REAL,ALLOCATABLE,TARGET :: M_2_0(:,:)    !< 1D-Mortar Operator: projection right interval 2: [0,1]  to full intervall 0: [-1,1]
#if FV_ENABLED
REAL,ALLOCATABLE,TARGET :: FV_M_0_1(:,:) !< 1D-Mortar Operator FV: interpolation full interval 0: [-1,1] to left interval 1: [-1,0]
REAL,ALLOCATABLE,TARGET :: FV_M_0_2(:,:) !< 1D-Mortar Operator FV: interpolation full interval 0: [-1,1] to right intervall 2:[0,1]
REAL,ALLOCATABLE,TARGET :: FV_M_1_0(:,:) !< 1D-Mortar Operator FV: projection left interval 1: [-1,0] to full intervall 0: [-1,1]
REAL,ALLOCATABLE,TARGET :: FV_M_2_0(:,:) !< 1D-Mortar Operator FV: projection right interval 2: [0,1] to full intervall 0: [-1,1]
#endif
LOGICAL                 :: MortarInitIsDone=.FALSE. !< marks whether mortar init routines are complete
!==================================================================================================================================
END MODULE MOD_Mortar_Vars
