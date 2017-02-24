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
!==================================================================================================================================
!> Contains global variables used by the Testcase Module.
!==================================================================================================================================
MODULE MOD_TestCase_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! REQUIRED VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER :: nAnalyzeTestCase=9999999           !< call AnalyzeTestCase every *th time step. May be adjusted in parameter file
LOGICAL :: doTCSource=.FALSE.                 !< compute source terms for testcase
CHARACTER(LEN=255) :: testcase = "riemann2d"  !< name of testcase
!----------------------------------------------------------------------------------------------------------------------------------
! TESTCASE VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER :: SHOCK         = 1
INTEGER,PARAMETER :: DISCONTINUITY = 2
INTEGER,PARAMETER :: RAREFACTION   = 3
CHARACTER(LEN=13),PARAMETER :: WAVENAMES(3) = (/ 'SHOCK        ', 'DISCONTINUITY', 'RAREFACTION  ' /)
INTEGER           :: RiemannBC_WaveType(4)        ! shock(1), discontinuity(2) or rarefaction(3) wave
REAL              :: RiemannBC_speeds(1:4,1:2)    ! shock, discontinuity, expansion tail, expansion head (for all 4 sides)

INTEGER,PARAMETER :: SideToQuads(2,4) = RESHAPE( (/ 2,1, 3,2, 3,4, 4,1 /), (/2,4/) )
!==================================================================================================================================
END MODULE MOD_TestCase_Vars
