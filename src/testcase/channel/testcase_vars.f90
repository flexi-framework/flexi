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
INTEGER :: nAnalyzeTestCase=9999999 !< call AnalyzeTestCase every *th time step. May be adjusted in parameter file
LOGICAL :: doTCSource=.TRUE.        !< compute source terms for testcase
CHARACTER(LEN=255) :: testcase = "channel"  !< name of testcase
!----------------------------------------------------------------------------------------------------------------------------------
! TESTCASE VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! precomputed variables
REAL    :: bulkVel     =0.          !< bulk velocity in domain
REAL    :: dpdx        =0.          !< imposed pressure gradient
REAL    :: Re_tau                   !< target friction Reynolds number
REAL    :: uBulk                    !< bulk velocity in domain
REAL    :: uBulkScale               !< scaling factor for bulk velocity
REAL,ALLOCATABLE   :: writeBuf(:,:) !< buffer to store log testcase data
INTEGER :: ioCounter   =0           !< current number of buffer items
INTEGER :: nWriteStats =-999        !< Write testcase statistics to file at every n-th AnalyzeTestcase step
CHARACTER(LEN=255) :: Filename      !< filename to store testcase log data
!==================================================================================================================================
END MODULE MOD_TestCase_Vars
