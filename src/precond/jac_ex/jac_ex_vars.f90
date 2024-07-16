!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz 
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
MODULE MOD_Jac_Ex_Vars
!===================================================================================================================================
! Contains global variables used by the Jac_Ex module.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                               :: Jac_Ex_InitIsDone=.FALSE.
REAL,ALLOCATABLE                      :: LL_plus(:,:)              !< LL_plus(i,j) = L^hat_plus(i)*L_plus(j)
REAL,ALLOCATABLE                      :: LL_minus(:,:)             !< LL_minus(i,j) = L^hat_minus(i)*L_minus(j)
#if PARABOLIC
REAL,ALLOCATABLE                      :: L_mp(:,:)                 !< L_mp(i,iLocSide)=either L_minus(i) or L_plus(i) 
REAL,ALLOCATABLE                      :: R_minus(:,:,:)            !< BR2 lifting surface term
REAL,ALLOCATABLE                      :: R_plus( :,:,:)            !< BR2 lifting surface term
REAL,ALLOCATABLE                      :: JacLiftingFlux(:,:,:,:,:,:)
#endif /*PARABOLIC*/
#if FV_ENABLED && FV_RECONSTRUCT
REAL,ALLOCATABLE                      :: UPrim_extended(:,:,:,:,:) !< extended primitive solution array containing additional
                                                                   !< first layers of neighbouring elements
REAL,ALLOCATABLE                      :: FV_sdx_XI_extended(:,:,:,:)   !< extended inverse of distance between neighboring dofs in x
REAL,ALLOCATABLE                      :: FV_sdx_ETA_extended(:,:,:,:)  !< extended inverse of distance between neighboring dofs in y
REAL,ALLOCATABLE                      :: FV_sdx_ZETA_extended(:,:,:,:) !< extended inverse of distance between neighboring dofs in z
#endif
!===================================================================================================================================
END MODULE MOD_Jac_Ex_Vars
