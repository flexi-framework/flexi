#include "flexi.h"
MODULE MOD_Implicit_Vars
!===================================================================================================================================
! Contains global variables used by the DG modules.
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
REAL,ALLOCATABLE                      :: LinSolverRHS(:,:,:,:,:)    ! right hand side of linear solver
REAL,ALLOCATABLE                      :: R_xk(:,:,:,:,:)
REAL,ALLOCATABLE                      :: Xk(:,:,:,:,:)
REAL,ALLOCATABLE                      :: mass(:,:,:,:,:)
! Predictor
INTEGER                               :: PredictorType              ! linear interpolation or Lagrange Interpolation
INTEGER                               :: PredictorCounter=1         ! auxilliary counter for lagrange predictor
INTEGER                               :: PredictorOrder             ! Order of the Lagrange polynomial
REAL,ALLOCATABLE                      :: Upast(:,:,:,:,:,:,:)       ! history of upast, required for predictor
REAL,ALLOCATABLE                      :: U_Predictor(:,:,:,:,:)     ! buffer needed for predictor
REAL,ALLOCATABLE                      :: Un_old(:,:,:,:,:)          ! solution at old time instance
REAL,ALLOCATABLE                      :: Ut_old(:,:,:,:,:,:)        ! time increment at old stages
REAL,ALLOCATABLE                      :: t_old(:,:)                 ! times at the old stages 
! number of DOF
INTEGER                               :: nDOFVar1D,nDOFVarElem      ! nVar*(N+1), nVar*(N+1)^3
INTEGER                               :: nDOFGlobal                 ! nDOFVarElem*nElems
INTEGER                               :: nDOFVarGlobal              ! nDOFVarElem*nElems summed for all Procs
! Newton variables
INTEGER                               :: nNewtonIter                ! max Newton iterations
INTEGER                               :: nNewtonIterGlobal          ! max Newton iterations global for the calculation
INTEGER                               :: nInnerNewton               ! max Newton iterations for actual Newton stage
! GMRES variables
INTEGER                               :: nRestarts                  ! max numbers of GMRES iterations
INTEGER                               :: nKDim                      ! max Krylov space
INTEGER                               :: nGMRESIterGlobal           ! global GMRES inner iterations, nulled at dt_Analyze
INTEGER                               :: nGMRESIterdt               ! GMRES iterations for one time step
INTEGER                               :: nGMRESGlobal               ! max GMRES solves (>1 if restart>1) (number of restarts)
INTEGER                               :: nInnerGMRES                ! max GMRES iterations for actual GMRES restart 
! epsilons
REAL                                  :: rEps0,srEps0               ! SQRT(EPSILON(0.0))
REAL                                  :: scaleps                    ! scaleps*SQRT(EPSILON(0.0)), used for EpsFD
REAL                                  :: EpsNewton                  ! newton relative epsilon
LOGICAL                               :: adaptepsNewton             ! adaptive eps Newton calculated by the embedded RK
REAL                                  :: Eps2Newton                 !  square of newton relative epsilon
REAL                                  :: EpsGMRES
REAL                                  :: gammaEW                    ! gamma parameter for Eisenstat Walker
LOGICAL                               :: EisenstatWalker=.FALSE.
LOGICAL                               :: ImplicitInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_Implicit_Vars

