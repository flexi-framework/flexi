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
!> Contains global variables used by the filter module
!==================================================================================================================================
MODULE MOD_Filter_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                :: NFilter
INTEGER                :: FilterType            !< filter to be applied 0: no filter, 1: cut-off filter, 2 :Hesthaven filter
REAL                   :: HestFilterParam(3)    !< filter parameters for modal Hesthaven style filter
REAL,ALLOCATABLE       :: FilterMat(:,:)        !< 1D nodal filter matrix
LOGICAL                :: FilterInitIsDone = .FALSE. !< filter routines have been initialized
#if EQNSYSNR==2
REAL,ALLOCATABLE       :: lim(:),eRatio(:),r(:) !< Analysis data for LAF model
REAL,ALLOCATABLE       :: ekin_avg_old(:)       !< cell integral value for ekin avg (LAF)
REAL,ALLOCATABLE       :: ekin_fluc_avg_old(:)  !< cell integral value for ekin fluc avg (LAF))
REAL,ALLOCATABLE       :: Vol(:)                !< cell volume for averaging (LAF) 
REAL,ALLOCATABLE       :: Integrationweight(:,:,:,:)  !< integration weights (LAF)
REAL                   :: normMod               !< spectral normalization (LAF) 
REAL,ALLOCATABLE       :: J_N(:,:,:)            !< Jacobi for volume integral (LAF) 
REAL                   :: LAF_alpha             !< Relaxation factor (LAF)
#endif
!==================================================================================================================================
END MODULE MOD_Filter_Vars
