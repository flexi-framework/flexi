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
!> Defines frequently used variables, that can be either set by the parameter file or precompiled
!> or directly depend on parameters set by the preprocessor
!==================================================================================================================================
MODULE MOD_PreProc
! MODULES
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: DP => REAL64
IMPLICIT NONE
PUBLIC
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,PARAMETER        :: PP_RealTolerance = EPSILON(1.0D0) !< machine precision
REAL,PARAMETER        :: PP_Pi = ACOS(REAL(-1.0,KIND=DP))   !< Pi up to machine accuracy
#if PP_N == N
INTEGER               :: PP_N                              !< polynomial degree
#endif

!==================================================================================================================================
END MODULE MOD_PreProc
