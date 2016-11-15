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
!> Contains the variables of phill testcase that should be globally accessible!
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
LOGICAL :: doTCSource=.FALSE.       !< compute source terms for testcase
CHARACTER(LEN=255) :: testcase = "phill"  !< name of testcase
!----------------------------------------------------------------------------------------------------------------------------------
! TESTCASE VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL :: firstTimestep=.TRUE.     !< marks whether timestep is first timestep
REAL    :: tPrev       =0.          !< time of previous forcing computation
REAL    :: dtPrev      =0.          !< dt since previous forcing computation 
REAL    :: massFlowRef =0.          !< reference mass flow for forcing
REAL    :: massFlowPrev=0.          !< massflow at previous forcing computation
REAL    :: bulkVel     =0.          !< bulk velocity in domain
REAL    :: bulkVelPrev =0.          !< bulk velocity in domain prev. forcing computation
REAL    :: dpdx        =0.          !< imposed pressure gradient
REAL    :: dpdxPrev    =0.          !< previous imposed pressure gradient
REAL,ALLOCATABLE   :: writeBuf(:,:) !< buffer to store log testcase data
INTEGER :: massFlowBC  =0           !< index of BC at which massflow is computed
INTEGER :: ioCounter   =0           !< current number of buffer items
INTEGER :: maxBuffer=0              !< max number of time log data
INTEGER :: Buffer=0                 !< expected size of time log buffer entries before IO
CHARACTER(LEN=255) :: Filename      !< filename to store testcase log data
!==================================================================================================================================

END MODULE MOD_TestCase_Vars
