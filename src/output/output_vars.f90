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
!> Contains global variables provided by the output routines
!==================================================================================================================================
MODULE MOD_Output_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                      :: NVisu                    !< Number of points at which solution is sampled for visualization
INTEGER                      :: NOut                     !< Polynomial degree at which solution is written. -1: NOut=N, >0: NOut
REAL,ALLOCATABLE             :: Vdm_GaussN_NVisu(:,:)    !< Vandermonde for direct interpolation from computation grid to visu grid
REAL,ALLOCATABLE             :: Vdm_N_NOut(:,:)          !< output vandermonde
REAL,PARAMETER               :: FileVersion=0.1          !< version written into output file
CHARACTER(LEN=255),PARAMETER :: ProgramName='Flexi'      !< name of program written into output file
CHARACTER(LEN=255)           :: ProjectName              !< Name of the current simulation (mandatory).
INTEGER                      :: outputFormat=0           !< File format for visualization. <=0: no visualization, 1: Tecplot binary,
                                                         !< 2: Tecplot ASCII, 3: Paraview binary. Note: Tecplot output is currently
                                                         !< unavailable due to licensing issues.
INTEGER                      :: ASCIIOutputFormat=0      !< File format for ASCII output. 0: CSV, 1: Tecplot
LOGICAL                      :: OutputInitIsDone=.FALSE. !< marks whether output routines have been initialized
LOGICAL                      :: doPrintStatusLine        !< flag indicating if status line should be printed
LOGICAL                      :: WriteStateFiles=.TRUE.   !< flag indicating if state files should be written
LOGICAL                      :: WriteTimeAvgFiles=.TRUE. !< flag indicating if time average files should be written
INTEGER                      :: userblock_total_len      !< length of userblock file + length of ini-file (with header) in bytes
CHARACTER(LEN=255)           :: UserBlockTmpFile='userblock.tmp' !< name of user block temp file
!==================================================================================================================================
END MODULE MOD_Output_Vars
