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
MODULE MOD_TestCase_Vars
!==================================================================================================================================
! Contains the variables of your testcase that should be globally accessible!
!==================================================================================================================================
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! REQUIRED VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER :: nAnalyzeTestCase=9999999     !< call AnalyzeTestCase every *th time step. May be adjusted in parameter file
LOGICAL :: doTCSource=.FALSE.           !< compute source terms for testcase
CHARACTER(LEN=255) :: testcase = "taylorgreenvortex"  !< name of testcase
!----------------------------------------------------------------------------------------------------------------------------------
! TESTCASE VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE   :: Time(:)           !< times of log data (nWriteStats)
REAL,ALLOCATABLE   :: writeBuf(:,:)     !< log data (nTGVVars+1,nWriteStats)
REAL               :: Ma0               !< reference Mach number of TGV
REAL,PARAMETER     :: rho0=1.           !< reference density
REAL,PARAMETER     ::   U0=1.           !< magnitude of speed
REAL               ::   T0              !< initial, constant reference temperature
REAL               ::   p0              !< background pressure
LOGICAL            :: IniConstDens      !< flag to indicate whether initial density or rather temperature field should be constant
#if PARABOLIC
INTEGER,PARAMETER  :: nTGVvars=15       !< Number of variables to be evaluated for TGV, time not included
#else
INTEGER,PARAMETER  :: nTGVvars=5        !< Number of variables to be evaluated for TGV, time not included
#endif
INTEGER            :: ioCounter   =0    !< current number of buffer items
INTEGER            :: nWriteStats =-999 !< Write testcase statistics to file at every n-th AnalyzeTestcase step
CHARACTER(LEN=255) :: Filename          !< filename to store testcase log data
!==================================================================================================================================

END MODULE MOD_TestCase_Vars
