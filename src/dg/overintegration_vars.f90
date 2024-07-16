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
!> Contains global variables used by the overintegration module
!> Different types / settings for "OverIntegrationType"
!> - 0: no overintegration: collocation DGSEM
!> - 1: Filtering of JU_t by modal projection
!> - 2: Like 1, but collocative division by J on N, not NO. Modal content between N and NO = 0.
!==================================================================================================================================
MODULE MOD_Overintegration_Vars
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER            :: OverintegrationType         !< 0: no overintegration, 1: cutoff, 2: conservative cutoff,

!----------------------------------------------------------------------------------------------------------------------------------
! used for type 1 and 2
INTEGER            :: NUnder                      !< Filter degree for cutoff, i.e. Nunder equals degree of filtered polynomial;
                                                  !< Overintegrationtype 1 and 2 only!

!----------------------------------------------------------------------------------------------------------------------------------
! used for type 1:
REAL,ALLOCATABLE   :: OverintegrationMat(:,:)     !< Overintegration filter matrix, size [0..N,0..N];
                                                  !< Overintegrationtype 1 only!

!----------------------------------------------------------------------------------------------------------------------------------
! used for type 2 only:
REAL,ALLOCATABLE   :: sJNUnder(:,:,:,:)           !< 1/Jacobi on NUnder, size [0..Nunder,0..Nunder,0..Nunder,nElems];
                                                  !< Overintegrationtype 2 only!

REAL,ALLOCATABLE   :: Vdm_NUnder_N(:,:)           !< 1D Vandermonde NUnder->N, size [0..N,0..Nunder];
                                                  !< Overintegrationtype 2 only!

REAL,ALLOCATABLE   :: Vdm_N_NUnder(:,:)           !< 1D Vandermonde N->NUnder, size [0..Nunder,0..N];
                                                  !< Overintegrationtype 2 only!
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL            :: OverintegrationInitIsDone = .FALSE. !< Switch to check if overintegration is initialized
!==================================================================================================================================
END MODULE MOD_Overintegration_Vars
