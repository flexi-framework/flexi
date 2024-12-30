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
!> Contains global variables used by the Analyze modules.
!==================================================================================================================================
MODULE MOD_Analyze_Vars
! MODULES
USE MOD_Globals, ONLY: DP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER              :: nWriteData=1                      !< data output (writing/visualizeing the solution, timeaverages etc.)
                                                          !< is performed every multiple of analyze_dt
REAL                 :: analyze_dt                        !< time intervall at which analysis routines are called
REAL                 :: WriteData_dt                      !< time intervall at which solution data is written
REAL                 :: tWriteData                        !< actual time at which next solution IO will be performed
! precomputed variables
#if FV_ENABLED == 1
INTEGER(KIND=DP)     :: totalFV_nElems=0                  !< total number of FV elements in domain
#elif FV_ENABLED == 2 || FV_ENABLED == 3
REAL                 :: FV_totalAlpha=0.                  !< average of blending coefficient in domain
#endif
#if PP_LIMITER
INTEGER(KIND=DP)     :: totalPP_nElems=0                  !< number of elements limited with positiy-preserving limiter
#endif
REAL,ALLOCATABLE     :: wGPSurf(:,:)                      !< wGPSurf(i,j)=wGP(i)*wGP(j)
REAL,ALLOCATABLE     :: wGPVol(:,:,:)                     !< wGPVol(i,j,k)=wGP(i)*wGP(j)*wGP(k)
REAL,ALLOCATABLE     :: Surf(:)                           !< surface of each analyze set (e.g. of each boundary condition)
REAL,ALLOCATABLE     :: ElemVol(:)                        !< volume of each element
REAL                 :: Vol                               !< volume of the domain
! Analyze features
LOGICAL              :: doCalcErrorNorms  =.FALSE.        !< marks whether error norms should be computed
LOGICAL              :: doAnalyzeToFile   =.FALSE.        !< marks whether error norms should be written to a file

! Performance features
REAL                 :: PIDTimeStart                      !< start system time for PID calculation
REAL                 :: PIDTimeEnd                        !< end   system time for PID calculation
REAL                 :: PID_kill                          !< kill PID for FLEXI
INTEGER              :: nCalcPID                          !< counter for iterations since last PID calculation
INTEGER              :: nCalcPIDMax                       !< compute PID at least after every Nth iteration
REAL                 :: PID                               !< current PID in FLEXI

! Analyze to file
REAL                 :: iterRestart=0                     !< contains iteration count of previous computation in case a restart is
                                                          !< performed. No restart: 0
REAL                 :: calcTimeRestart=0.                !< contains simulation time at which a restart has been performed
                                                          !< No restart: 0

! Variables for the specific analyze routines
! ErrorNorms
INTEGER              :: NAnalyze                          !< polynomial degree analysis is performed at (e.g. computation of L2
                                                          !< norms). Number of points: NAnalyze+1
INTEGER              :: NAnalyzeZ
INTEGER              :: AnalyzeExactFunc                  !< Exact function used for analyze routines
INTEGER              :: AnalyzeRefState                   !< State used for analyze routines
REAL,ALLOCATABLE     :: wGPVolAnalyze(:,:,:)              !< product of GL integration weights used for analyze routines
REAL,ALLOCATABLE     :: Vdm_GaussN_NAnalyze(:,:)          !< Vandermonde for interpolating the solution to analyze points

#if FV_ENABLED
REAL,ALLOCATABLE     :: FV_Vdm_NAnalyze(:,:)              !< Vandermonde for interpolating the solution to analyze points
#endif

CHARACTER(LEN=255)   :: Filename_ErrNorm                  !< filename into which error norms are written

LOGICAL              :: AnalyzeInitIsDone = .FALSE.       !< marks whether analyze routines have been inittialized
!==================================================================================================================================
END MODULE MOD_Analyze_Vars
