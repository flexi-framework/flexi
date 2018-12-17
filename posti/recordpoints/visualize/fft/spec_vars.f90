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
! Contains global variables provided by spectral analysis routines
!===================================================================================================================================
MODULE MOD_Spec_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: nSamples_spec               !> Number of RP samples used for spectral analysis
REAL,ALLOCATABLE              :: RPData_spec(:,:,:)          !> Spectral data
REAL,ALLOCATABLE              :: RPData_freq(:)              !> List of frequencies
INTEGER                       :: nSamples_oct                !> Number of RP samples used for 1/3 octave analysis
REAL,ALLOCATABLE              :: RPData_Oct(:,:,:)           !> Spectral data for 1/3 octave analysis
REAL,ALLOCATABLE              :: RPData_freqOct(:)           !> List of frequencies for 1/3 octave analysis
!===================================================================================================================================
END MODULE MOD_Spec_Vars
