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
#if PARABOLIC
!==================================================================================================================================
!> Contains global variables used by the lifting procedure.
!==================================================================================================================================
MODULE MOD_Lifting_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! interior face gradients for all elements
REAL,ALLOCATABLE :: gradUx_slave(:,:,:,:)         !< slave side gradients in x-dir
REAL,ALLOCATABLE :: gradUy_slave(:,:,:,:)         !< slave side gradients in x-dir
REAL,ALLOCATABLE :: gradUz_slave(:,:,:,:)         !< slave side gradients in x-dir
REAL,ALLOCATABLE :: gradUx_master(:,:,:,:)        !< master side gradients in x-dir
REAL,ALLOCATABLE :: gradUy_master(:,:,:,:)        !< master side gradients in x-dir
REAL,ALLOCATABLE :: gradUz_master(:,:,:,:)        !< master side gradients in x-dir
REAL,ALLOCATABLE :: FluxX(:,:,:,:)                !< gradient flux in x-dir
REAL,ALLOCATABLE :: FluxY(:,:,:,:)                !< gradient flux in y-dir
REAL,ALLOCATABLE :: FluxZ(:,:,:,:)                !< gradient flux in z-dir
! the lifted gradients
REAL,ALLOCATABLE :: gradUx(:,:,:,:,:)             !< gradients in x-dir at degree N
REAL,ALLOCATABLE :: gradUy(:,:,:,:,:)             !< gradients in y-dir at degree N
REAL,ALLOCATABLE :: gradUz(:,:,:,:,:)             !< gradients in z-dir at degree N
! Diffusive fluxes
REAL,ALLOCATABLE :: diffFluxX_L(:,:,:)
REAL,ALLOCATABLE :: diffFluxX_R(:,:,:)
REAL,ALLOCATABLE :: diffFluxY_L(:,:,:)
REAL,ALLOCATABLE :: diffFluxY_R(:,:,:)
REAL,ALLOCATABLE :: diffFluxZ_L(:,:,:)
REAL,ALLOCATABLE :: diffFluxZ_R(:,:,:)

REAL             :: etaBR2                        !< Lifting penalty for BR2. Increase improves stability
                                                  !< at the cost of performance and reduces jumps between two cells
REAL             :: etaBR2_wall                   !< Lifting penalty for BR2 at the wall. Can be choosen different from etaBR2
                                                  !< to decrease wall velocities
LOGICAL          :: doWeakLifting=.FALSE.         !< marks whether lifting is peformed in weak form, BR2 is always stong
LOGICAL          :: doConservativeLifting=.FALSE. !< marks whether volume contribution to the gradients is computed in
                                                  !< conservative form, i.e. deriving the solution divided by the metric terms,
                                                  !< instead of deriving the solution and multiplying by the metrics
LOGICAL          :: LiftingInitIsDone=.FALSE.     !< marks whether the lifting init routines are complete
!==================================================================================================================================
END MODULE MOD_Lifting_Vars
#endif /*PARABOLIC*/
