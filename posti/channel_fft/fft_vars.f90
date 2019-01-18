!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
!> Global variables for the channel fft tool
!===================================================================================================================================
MODULE MOD_FFT_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

! User-defined parameters
INTEGER                 :: OutputFormat      !< Choose the main format for output. 0: Tecplot, 2: HDF5
INTEGER                 :: NCalc             !< Polynomial degree to perform DFFT on
REAL                    :: Re_tau            !< Reynolds number based on friction velocity and channel half height

REAL,ALLOCATABLE        :: VdmGaussEqui(:,:) !< Vandermonde from state to FFT grid

CHARACTER(LEN=255)      :: ProjectName       !< Name of project, read from state
REAL                    :: Time              !< Timestamp of project, read from first state file

! FFT variables
INTEGER(KIND=8)         :: planI
INTEGER(KIND=8)         :: planK
INTEGER                 :: N_FFT(3)
INTEGER                 :: nSamplesI
INTEGER                 :: nSamples_specI
INTEGER                 :: nSamplesK
INTEGER                 :: nSamples_specK
REAL,ALLOCATABLE        :: U_FFT(:,:,:,:)
REAL,ALLOCATABLE        :: X_FFT(:,:,:,:)
REAL,ALLOCATABLE        :: Ex_uu(:,:)
REAL,ALLOCATABLE        :: Ex_vv(:,:)
REAL,ALLOCATABLE        :: Ex_ww(:,:)
REAL,ALLOCATABLE        :: Ex_pp(:,:)
REAL,ALLOCATABLE        :: Ez_uu(:,:)
REAL,ALLOCATABLE        :: Ez_vv(:,:)
REAL,ALLOCATABLE        :: Ez_ww(:,:)
REAL,ALLOCATABLE        :: Ez_pp(:,:)
COMPLEX,ALLOCATABLE     :: inI(:)
COMPLEX,ALLOCATABLE     :: outI(:)
COMPLEX,ALLOCATABLE     :: inK(:)
COMPLEX,ALLOCATABLE     :: outK(:)
REAL,ALLOCATABLE        :: MS_t(:,:)
REAL,ALLOCATABLE        :: MS_PSD(:,:)
REAL,ALLOCATABLE        :: M_t(:,:)
END MODULE MOD_FFT_Vars
