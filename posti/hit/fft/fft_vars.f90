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

!===================================================================================================================================
!> Contains global variables provided by the HIT_FFT routines
!===================================================================================================================================
MODULE MOD_HIT_FFT_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER               :: N_FFT                !> Number of (global) interpolation points 1D for DFFT
INTEGER               :: NCalc                !> Polynomial degree of equidistant basis in each DG element to build global FFT basis
INTEGER               :: Nc                   !> Nyquist wavenumber Nc=N_FFT/2
INTEGER               :: kmax                 !> Maximum effective wavenumber in 3D
INTEGER               :: endw(3)              !> Max. number of wavenumbers is each of the 3D directions Endw=(Nc,N_FFT,N_FFT)
INTEGER,ALLOCATABLE   :: LocalK(:,:,:,:)      !> 1D Wave number equivalent for each Fourier Mode in 3D
COMPLEX,PARAMETER     :: II = CMPLX(0.,1.0)   !> Complex unit, sqrt(-1)
!===================================================================================================================================
END MODULE MOD_HIT_FFT_Vars
