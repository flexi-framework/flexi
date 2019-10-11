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
#include "flexi.h"
MODULE MOD_Implicit_Vars
!===================================================================================================================================
!> Contains global variables used by the implicit modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! DG solution type
REAL,ALLOCATABLE                      :: LinSolverRHS(:,:,:,:,:)    !< right hand side of linear solver
REAL,ALLOCATABLE                      :: Xk(:,:,:,:,:)              !< current value in Newton iteration
REAL,ALLOCATABLE                      :: R_Xk(:,:,:,:,:)            !< DG operator evaluated using current Newton iterated value
! Predictor
INTEGER                               :: PredictorType              !< Choose between four different predictor types
REAL,ALLOCATABLE                      :: U_predictor(:,:,:,:,:)     !< Predicted value for next stage
INTEGER                               :: PredictorOrder             !< Order of the Lagrange polynomial for predictor type 2
REAL,ALLOCATABLE                      :: Upast(:,:,:,:,:,:,:)       !< Values at each stage of last time steps, needed for predictor
                                                                    !< type 2
REAL,ALLOCATABLE                      :: t_old(:,:)                 !< times at the old stages for type 2
REAL,ALLOCATABLE                      :: Un_old(:,:,:,:,:)          !< Solution at old time instance for predictor type 3
REAL,ALLOCATABLE                      :: Ut_old(:,:,:,:,:,:)        !< Time derivative at old stages from last time step for type 3
! number of DOF
INTEGER                               :: nDOFVar1D,nDOFVarElem      !< nVar*(N+1), nVar*(N+1)^3
INTEGER                               :: nDOFVarProc                !< nDOFVarElem*nElems
INTEGER                               :: nDOFVarGlobal              !< nDOFVarElem*nElems summed for all Procs
! Newton variables
INTEGER                               :: nNewtonIter                !< max Newton iterations
INTEGER                               :: nNewtonIterGlobal          !< counter of Newton iterations global for the calculation
INTEGER                               :: nInnerNewton               !< counter of Newton iterations for current stage
REAL                                  :: epsNewton                  !< newton relative epsilon
LOGICAL                               :: adaptepsNewton             !< adaptive eps Newton calculated by the embedded RK
! GMRES variables
INTEGER                               :: nRestarts                  !< max numbers of GMRES restarts
INTEGER                               :: nKDim                      !< max dimension of Krylov space before restart
INTEGER                               :: nGMRESIterGlobal           !< counter of global GMRES inner iterations, nulled at dt_Analyze
INTEGER                               :: nGMRESIterdt               !< counter of GMRES iterations for one time step
INTEGER                               :: nGMRESRestartGlobal        !< counter of GMRES restarts
INTEGER                               :: nInnerGMRES                !< counter of GMRES iterations for single linear solve
LOGICAL                               :: EisenstatWalker=.FALSE.    !< Use Eisenstat Walker abort criterion for GMRES
REAL                                  :: gammaEW                    !< gamma parameter for Eisenstat Walker
REAL                                  :: EpsGMRES                   !< Fixed abort criterion for GMRES (if no Eisenstat Walker)
! FD approximation of A*v
INTEGER                               :: FD_Order                   !< Chose oder of finite difference
REAL                                  :: rEps0,srEps0               !< SQRT(EPSILON(0.0)), for FD step size
REAL                                  :: rEps0_O1,srEps0_O1         !< SQRT(EPSILON(0.0)), for FD step size if FD O1
REAL                                  :: scaleps                    !< scaleps*SQRT(EPSILON(0.0)), used for EpsFD
INTEGER                               :: Eps_Method                 !< Chose method how rEps0 is formed
REAL                                  :: Norm_Xk                    !< Norm required for eps_method=2

LOGICAL                               :: ImplicitInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_Implicit_Vars
