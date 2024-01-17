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
!==================================================================================================================================
!> Contains the variables of your testcase that should be globally accessible!
!==================================================================================================================================
MODULE MOD_TestCase_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! REQUIRED VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)        :: testcase = "hit"         !< name of testcase
LOGICAL                   :: doTCSource=.TRUE.        !< compute source terms for testcase
INTEGER                   :: nAnalyzeTestCase=9999999 !< call AnalyzeTestCase every *th time step. dummy variable
!----------------------------------------------------------------------------------------------------------------------------------
! TESTCASE VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                   :: HIT_Forcing              !< Flag to perform HIT forcing
LOGICAL                   :: HIT_Avg                  !< Flag to perform spatial averaging
LOGICAL                   :: HIT_1st                  !< Flag to update forcing only in the 1st RK stage
LOGICAL                   :: HIT_RMS_InitDone=.FALSE. !< Flag to indicate whether HIT_RMS already initialized
REAL                      :: HIT_k                    !< target turbulent kinetic energy
REAL                      :: HIT_rho                  !< density for initialization
REAL                      :: HIT_tauRMS               !< relaxation parameter for linear forcing
REAL,ALLOCATABLE          :: HIT_RMS(:,:,:,:,:)       !< u_rms of current solution
REAL                      :: HIT_tFilter              !< filter width of temporal filter
REAL                      :: A_ILF                    !< forcing coefficient
!----------------------------------------------------------------------------------------------------------------------------------
! ANAYLZE VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE          :: Time(:)                  !< times of log data (nWriteStats)
REAL,ALLOCATABLE          :: writeBuf(:,:)            !< log data (nHITVars+1,nWriteStats)
INTEGER,PARAMETER         :: nHITvars    = 7          !< Number of variables to be evaluated for HIT, time not included
INTEGER                   :: ioCounter   = 0          !< current number of buffer items
INTEGER                   :: nWriteStats =-999        !< Write testcase statistics to file at every n-th AnalyzeTestcase step
CHARACTER(LEN=255)        :: FileName                 !< filename to store testcase log data
!==================================================================================================================================

END MODULE MOD_TestCase_Vars
