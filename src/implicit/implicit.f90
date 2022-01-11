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

!===================================================================================================================================
!> Contains the modules needed for implicit time integration. Besides the initaliation, this includes the methods used to solve the
!> non-linear system (Newton's method) and the inner linear system (GMRES). GMRES requires the evaluation of A*v (matrix-vector
!> product), which we approximate using a finite difference approach, which can be found in MatrixVector.
!===================================================================================================================================
MODULE MOD_Implicit
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersImplicit
  MODULE PROCEDURE DefineParametersImplicit
END INTERFACE

INTERFACE InitImplicit
  MODULE PROCEDURE InitImplicit
END INTERFACE

INTERFACE Newton
  MODULE PROCEDURE Newton
END INTERFACE

INTERFACE FinalizeImplicit
  MODULE PROCEDURE FinalizeImplicit
END INTERFACE

PUBLIC::InitImplicit,Newton,FinalizeImplicit,DefineParametersImplicit
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters of implicit scheme
!==================================================================================================================================
SUBROUTINE DefineParametersImplicit()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Implicit")
CALL prms%CreateLogicalOption('adaptepsNewton', "Adaptive Newton eps by Runge-Kutta error estimation", value='.FALSE.')
CALL prms%CreateRealOption(   'EpsNewton',      "Newton tolerance, only used if adaptepsNewton=F", value='1.E-3')
CALL prms%CreateIntOption(    'nNewtonIter',    "Maximum amount of Newton iterations", value='50')
CALL prms%CreateLogicalOption('EisenstatWalker',"Adaptive abort criterion for GMRES", value='.FALSE.')
CALL prms%CreateRealOption(   'gammaEW',        "Parameter for Eisenstat Walker adaptation", value='0.9')
CALL prms%CreateRealOption(   'EpsGMRES',       "GMRES Tolerance, only used of EisenstatWalker=F", value='1.E-3')
CALL prms%CreateIntOption(    'nRestarts',      "Maximum number of GMRES Restarts", value='10')
CALL prms%CreateIntOption(    'nKDim',          "Maximum number of Krylov subspaces for GMRES, after that a restart is performed", &
                                                 value='30')
CALL prms%CreateIntOption(    'Eps_Method',     "Method of determining the step size of FD approximation of A*v in GMRES, &
                                                &1: sqrt(machineAccuracy)*scaleps, 2: take norm of solution into account", value='2')
CALL prms%CreateRealOption(   'scaleps',        "Scaling factor for step size in FD, mainly used in Eps_Method=1", value='1.')
CALL prms%CreateIntOption(    'FD_Order',       "Order of FD approximation (1/2)", value='1')
CALL prms%CreateIntOption(    'PredictorType',  "Type of predictor to be used, 0: use current U, 1: use right hand side, 2: &
                                                 &polynomial extrapolation, 3: dense output formula of RK scheme",value='0')
CALL prms%CreateIntOption(    'PredictorOrder', "Order of predictor to be used (PredictorType=2)",value='0')

END SUBROUTINE DefineParametersImplicit

!===================================================================================================================================
!> Initialize implicit time integration. Mainly read in parameters for the Newton and GMRES solvers and call preconditioner init.
!===================================================================================================================================
SUBROUTINE InitImplicit()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars
USE MOD_Mesh_Vars,         ONLY:MeshInitIsDone,nElems,nGlobalElems
#if PP_dim==3
USE MOD_Mesh_Vars,         ONLY:firstInnerSide,lastInnerSide,SideToElem
#endif
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone
USE MOD_ReadInTools,       ONLY:GETINT,GETREAL,GETLOGICAL
USE MOD_DG_Vars,           ONLY:nDOFElem
USE MOD_TimeDisc_Vars,     ONLY:TimeDiscType,RKb_embedded
USE MOD_Precond,           ONLY:InitPrecond
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if PP_dim==3
INTEGER                   :: i
#endif
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.ImplicitInitIsDone)THEN
   CALL abort(__STAMP__,'InitImplicit not ready to be called or already called.')
END IF
IF(TimeDiscType.EQ.'ESDIRK') THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' INIT Implicit...'

  nDOFVar1D     = PP_nVar*(PP_N+1)
  nDOFVarElem   = PP_nVar*nDOFElem
  nDOFVarProc   = nDOFVarElem*nElems
  nDOFVarGlobal = nDOFVarElem*nGlobalElems

  !Abort for 2D periodic meshes, when compiling in 3D. Preconditioner is not working in that case
#if PP_dim==3
  DO i=firstInnerSide,lastInnerSide
    IF(SideToElem(S2E_ELEM_ID,i).EQ.SideToElem(S2E_NB_ELEM_ID,i)) THEN
      CALL CollectiveStop(__STAMP__,'ERROR - This is a 2D mesh.')
    ENDIF
  END DO
#endif

  !========================================================================================
  ! Variables used for Newton
  ! Root function: F_Xk = Xk - Q - alpha * dt* R_Xk(t+beta*dt,Xk) = 0!
  !========================================================================================
  ! the constant vector part of non-linear system (Q)
  ALLOCATE(LinSolverRHS(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  ! Xk is the root of the Newton method
  ALLOCATE(Xk(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  ! R_Xk is the DG Operator depending on the solution of the current time (implicit)
  ! In Newton: DG Operator depending on the actual Newton iteration value "Xk"
  ALLOCATE(R_Xk(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  LinSolverRHS = 0.
  Xk           = 0.
  R_Xk         = 0.

  ! Abort criterion for Newton
  EpsNewton      = GETREAL('EpsNewton','0.001')
  ! Adaptive abort criterion, we need specific coefficients from the embedded RK scheme for this to work. Not available for all of
  ! them!
  adaptepsNewton = GETLOGICAL('adaptepsNewton','.FALSE.')
  IF((.NOT.ALLOCATED(RKb_embedded)).AND.(adaptepsNewton.EQV..TRUE.))THEN
    SWRITE(UNIT_stdOut,'(A)') ' Attention, Newton abort criterium is not adaptive!'
    adaptepsNewton = .FALSE.
  END IF
  nNewtonIter    = GETINT('nNewtonIter','50')
  ! initialize identifier if Newton's method has failed
  NewtonConverged = .TRUE.

  ! Parameters for the finite difference approximation of A*v (matrix free GMRES)
  scaleps        = GETREAL('scaleps','1.') ! (A*v = (R(Xk+eps)-R(Xk))/(scaleps*eps))
  ! Choose method how eps is calculated:
  ! 1: According to Qin,Ludlow,Shaw: A matrix-free preconditioned Newton/GMRES method for unsteady Navier-Stokes solutions (Eq. 13),
  !                                  Int.J.Numer.Meth.Fluids 33 (2000) 223-248
  ! 2: According to Knoll,Keyes: Jacobian-free Newton-Krylov methods: A survey of approaches and applications (Eq. (14)),
  !                              JCP 193 (2004) 357-397
  Eps_Method     = GETINT('Eps_Method','2')
  ! Choose order of finite difference
  FD_Order       = GETINT('FD_Order','1')

  ! Adapt machine epsilon to order of finite difference according to: An, Weng, Feng: On finite difference approximation of a
  ! matrix-vector product in the Jacobian-free Newton-Krylov method (Eqs. (11)-(13)), J.Comp.Appl.Math. 263 (2011) 1399-1409
  SELECT CASE(FD_Order)
  CASE(1)
    rEps0        = scaleps*SQRT(EPSILON(0.0))
    rEps0_O1     = rEps0
  CASE(2)
    rEps0        = (scaleps**2*0.5*EPSILON(0.0))**(1./3.)
    rEps0_O1     = scaleps*SQRT(EPSILON(0.0))
  END SELECT
  srEps0           =1./rEps0
  srEps0_O1        =1./rEps0_O1

  nNewtonIterGlobal=0
  nInnerNewton     =0

  !========================================================================================
  ! Variables using for GMRES
  ! GMRES solving the LES: A * y = b
  ! Matrix          A = I - alpha*dt * dR/dx
  ! Solution        y = dx
  ! Right Hand Side b = -F_Xk
  !========================================================================================
  nRestarts       =GETINT    ('nRestarts','10')
  nKDim           =GETINT    ('nKDim','30')
  EisenstatWalker =GETLOGICAL('EisenstatWalker','.FALSE.')
  gammaEW         =GETREAL   ('gammaEW','0.9')
  IF (.NOT.EisenstatWalker) EpsGMRES = GETREAL('EpsGMRES','0.001')

  nGMRESIterGlobal    = 0
  nGMRESIterdt        = 0
  nGMRESRestartGlobal = 0

  ! Preconditioner
  CALL InitPrecond()

  ImplicitInitIsDone=.TRUE.
  SWRITE(UNIT_stdOut,'(A)')' INIT Implicit DONE!'
  SWRITE(UNIT_stdOut,'(132("-"))')
END IF
END SUBROUTINE InitImplicit

!===================================================================================================================================
!> Solves the non-linear system with Newton
!> Root function: F_Xk = Xk - Q- alpha * dt* R_Xk(t+beta*dt,Xk) = 0!
!> Newton algorithm:
!> dF/dX|Xk * DeltaX = -F_Xk
!> X_K+1 = X_k + DeltaX
!> Attention: we use actual U as X0
!===================================================================================================================================
SUBROUTINE Newton(t,Alpha)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mathtools     ,ONLY: GlobalVectorDotProduct
USE MOD_DG            ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars       ,ONLY: U,Ut
USE MOD_Mesh_Vars     ,ONLY: nElems
USE MOD_Implicit_Vars ,ONLY: EpsNewton,Xk,R_Xk,nNewtonIter,nNewtonIterGlobal,nInnerNewton,nDOFVarGlobal,nDOFVarProc
USE MOD_Implicit_Vars ,ONLY: gammaEW,EpsGMRES,EisenstatWalker,LinSolverRHS,U_predictor
USE MOD_Implicit_Vars ,ONLY: PredictorType,Eps_Method,Norm_Xk,NewtonConverged
USE MOD_TimeDisc_Vars ,ONLY: dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)      :: t       !< current simulation time
REAL,INTENT(IN)      :: Alpha   !< coefficient of RK stage
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: AbortCritNewton
REAL             :: F_X0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL             :: F_Xk(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL             :: Norm_F_X0,Norm2_F_X0,Norm2_F_Xk,Norm2_F_Xk_old,Norm_F_Xk
REAL             :: DeltaX(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL             :: eta_k,etaA,etaB,etaMax
INTEGER          :: i,j,k,iElem,iVar
!===================================================================================================================================
NewtonConverged = .TRUE.
!Initialization for the Newton iterations

CALL DGTimeDerivative_weakForm(t)
DO iElem=1,nElems
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          F_X0(iVar,i,j,k,iElem) = U(iVar,i,j,k,iElem)-LinSolverRHS(iVar,i,j,k,iElem)-alpha*dt*Ut(iVar,i,j,k,iElem)
        END DO
      END DO
    END DO
  END DO
END DO
IF(PredictorType.NE.0)THEN ! U is not already U_predictor
  U = U_predictor
  CALL DGTimeDerivative_weakForm(t)
END IF
DO iElem=1,nElems
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          F_Xk(iVar,i,j,k,iElem) = U(iVar,i,j,k,iElem)-LinSolverRHS(iVar,i,j,k,iElem)-alpha*dt*Ut(iVar,i,j,k,iElem)
          Xk(iVar,i,j,k,iElem)   = U(iVar,i,j,k,iElem)
          R_Xk(iVar,i,j,k,iElem) = Ut(iVar,i,j,k,iElem)
        END DO
      END DO
    END DO
  END DO
END DO

! Eps_method 2 takes the norm of Xk into account for the FD step size
IF(Eps_Method.EQ.2)THEN
  CALL GlobalVectorDotProduct(Xk,Xk,nDOFVarProc,Norm_Xk)
  Norm_Xk = SQRT(Norm_Xk)
END IF

! Preparation for the Abort Criteria of Newton
! |F_Xk| < epsNewton * |F_Xk|
CALL GlobalVectorDotProduct(F_X0,F_X0,nDOFVarProc,Norm2_F_X0)
Norm_F_X0=SQRT(Norm2_F_X0)
IF (Norm_F_X0.LE.1.E-9*nDOFVarGlobal) THEN ! do not iterate, as U is already the implicit solution
  Norm_F_Xk=0.
ELSE ! we need iterations
  IF(PredictorType.NE.0)THEN
    CALL GlobalVectorDotProduct(F_Xk,F_Xk,nDOFVarProc,Norm2_F_Xk)
  ELSE
    Norm2_F_Xk=Norm2_F_X0
  END IF
  Norm_F_Xk=SQRT(Norm2_F_Xk)
END IF

!Eisenstat and Walker controls the Abort Criteria for GMRES
IF (EisenstatWalker.EQV..TRUE.) THEN
  etaMax = 0.9999
END IF

AbortCritNewton=EpsNewton*Norm_F_X0

!===================================================================================================================================
! Newton loop
!===================================================================================================================================
nInnerNewton=0 ! counts the newton steps in one implicit solve

DO WHILE((Norm_F_Xk.GT.AbortCritNewton).AND.(nInnerNewton.LT.nNewtonIter))

  ! Computation of the forcing terms eta_k for the Abort Criteria of GMRES via Eisenstat and Walker
  ! S. C. Eisenstat and H. F. Walker. “Choosing the forcing terms in an inexact Newton method”.
  ! In: SIAM Journal on Scientific Computing 17.1 (1996), pp. 16–32.
  IF(EisenstatWalker.EQV..TRUE.) THEN
    IF (nInnerNewton.EQ.0) THEN
      eta_k=etaMax
    ELSE
      etaA = gammaEW*(Norm2_F_Xk)/(Norm2_F_Xk_old)
      IF (gammaEW*eta_k*eta_k .LE. 0.1) THEN
        etaB = min(etaMax,etaA)
      ELSE
        etaB = min(etaMax, max(etaA,gammaEW*eta_k*eta_k))
      ENDIF
      eta_k = min(etaMax,max(etaB,0.5*AbortCritNewton/SQRT(Norm2_F_Xk)))
    END IF
    Norm2_F_Xk_old=Norm2_F_Xk
  ELSE
    eta_k=EpsGMRES
  END IF !EisenstatWalker
  nInnerNewton=nInnerNewton+1

  !======================================================================
  ! CALL of the linear solver GMRES
  CALL GMRES_M(t,Alpha,-F_Xk,Norm_F_Xk,eta_k,DeltaX)
  !======================================================================

  ! Update to next iteration
  Xk  =Xk+DeltaX
  U   =Xk
  ! Calculate new norm
  CALL DGTimeDerivative_weakForm(t)
  R_Xk=Ut
  F_Xk=(U-LinSolverRHS-alpha*dt*R_Xk)
  CALL GlobalVectorDotProduct(F_Xk,F_Xk,nDOFVarProc,Norm2_F_Xk)
  Norm_F_Xk=SQRT(Norm2_F_Xk)
END DO

! Global Newton iteration counter
nNewtonIterGlobal=nNewtonIterGlobal+nInnerNewton
! Check if we left the loop because the maximum number of newton iterations have been reached
IF((Norm_F_Xk.GT.EpsNewton*Norm_F_X0).AND.(nInnerNewton.EQ.nNewtonIter)) THEN
  NewtonConverged = .FALSE.
END IF

END SUBROUTINE Newton

!===================================================================================================================================
!> Matrix free Linear Krylov subspace solver
!> A * y = b
!> Matrix          A = I - Alpha*dt * dR/dx
!> Solution        y = delta x
!> Right Hand Side b = -F_Xk
!> Inital guess: Delta x = 0
!> Y. Saad and M. H. Schultz. “GMRES: A generalized minimal residual  algorithm for solving nonsymmetric linear systems”
!> In: SIAM Journalon scientific and statistical computing 7.3 (1986), pp. 856–869.
!===================================================================================================================================
SUBROUTINE GMRES_M(t,Alpha,B,Norm_B,AbortCrit,DeltaX)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mathtools     ,ONLY:GlobalVectorDotProduct
USE MOD_Implicit_Vars ,ONLY:nKDim,nRestarts,nGMRESIterGlobal,nGMRESRestartGlobal,nInnerGMRES,nGMRESIterdt
USE MOD_Implicit_Vars ,ONLY:nDOFVarProc
USE MOD_Precond       ,ONLY:ApplyPrecond
USE MOD_Precond_Vars  ,ONLY:PrecondType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t
REAL,INTENT(IN)   :: Alpha,Norm_B
REAL,INTENT(IN)   :: B(nDOFVarProc)
REAL,INTENT(IN)   :: AbortCrit
REAL,INTENT(OUT)  :: DeltaX(nDOFVarProc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: V(nDOFVarProc,1:nKDim)
REAL              :: W(nDOFVarProc)
REAL              :: Z(nDOFVarProc,nKDim)
REAL              :: R0(nDOFVarProc)
REAL              :: Gam(1:nKDim+1),C(1:nKDim),S(1:nKDim),H(1:nKDim+1,1:nKDim+1),Alp(1:nKDim)
REAL              :: Norm_R0,Resu,Temp,Bet
INTEGER           :: Restart
INTEGER           :: m,nn,o
REAL              :: tol
!===================================================================================================================================
!Initializations
R0         =B
Norm_R0    =Norm_B
DeltaX     =0.
Restart    =0
nInnerGMRES=0

! |dF/dU|U_k *delta U + F(U_K) | < eta_k * |F(U_K)|
! eta_k = AbortCrit forcing terms
tol=Norm_B*AbortCrit

! GMRES Loop
DO WHILE (Restart<nRestarts)
! GMRES(m)
  V(:,1)=R0/Norm_R0
  Gam(1)=Norm_R0
  DO m=1,nKDim
    nInnerGMRES=nInnerGMRES+1
    ! Preconditioner
    IF(PrecondType.NE.0) THEN
      CALL ApplyPrecond(V(:,m),Z(:,m))
    ELSE
      Z(:,m)=V(:,m)
    END IF
    ! matrix vector
    CALL MatrixVector(t,Alpha,Z(:,m),W)
    ! modified Gram-Schmidt
    DO nn=1,m
      CALL GlobalVectorDotProduct(V(:,nn),W,nDOFVarProc,H(nn,m))
      W=W-H(nn,m)*V(:,nn)
    END DO !nn
    CALL GlobalVectorDotProduct(W,W,nDOFVarProc,Resu)
    H(m+1,m)=SQRT(Resu)
    ! Givens Rotation
    DO nn=1,m-1
      Temp     =   C(nn)*H(nn,m) + S(nn)*H(nn+1,m)
      H(nn+1,m) = - S(nn)*H(nn,m) + C(nn)*H(nn+1,m)
      H(nn,m)   =   Temp
    END DO !nn
    Bet=SQRT(H(m,m)*H(m,m)+H(m+1,m)*H(m+1,m))
    S(m)=H(m+1,m)/Bet
    C(m)=H(m,m)/Bet
    H(m,m)=Bet
    Gam(m+1)=-S(m)*Gam(m)
    Gam(m)=C(m)*Gam(m)
    IF ((ABS(Gam(m+1)).LE.tol) .OR. (m.EQ.nKDim)) THEN !converge or max Krylov reached
      DO nn=m,1,-1
         Alp(nn)=Gam(nn)
         DO o=nn+1,m
           Alp(nn)=Alp(nn) - H(nn,o)*Alp(o)
         END DO !o
         Alp(nn)=Alp(nn)/H(nn,nn)
      END DO !nn
      ! Preconditioner back is not needed, as we used the basis vectors z
      DO nn=1,m
        DeltaX=DeltaX+Alp(nn)*Z(:,nn)
      END DO !nn
      IF (ABS(Gam(m+1)).LE.tol) THEN !converged
        nGMRESRestartGlobal = nGMRESRestartGlobal+Restart+1
        nGMRESIterGlobal    = nGMRESIterGlobal+nInnerGMRES
        nGMRESIterdt        = nGMRESIterdt + nInnerGMRES
        RETURN
      END IF  ! converged
    ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.tol) .OR. (m.EQ.nKDim))
      V(:,m+1)=W/H(m+1,m)
    END IF ! ((ABS(Gam(m+1)).LE.tol) .OR. (m.EQ.nKDim))
  END DO ! m
  ! Restart needed
  Restart=Restart+1
  CALL MatrixVector(t,Alpha,DeltaX,R0)
  R0=B-R0
  CALL GlobalVectorDotProduct(R0,R0,nDOFVarProc,Norm_R0)
  Norm_R0=SQRT(Norm_R0)
END DO ! While Restart

! convergence criterion not reached, nevertheless continue with Newton iterations
nGMRESRestartGlobal = nGMRESRestartGlobal+Restart+1
nGMRESIterGlobal    = nGMRESIterGlobal+nInnerGMRES
nGMRESIterdt        = nGMRESIterdt + nInnerGMRES
END SUBROUTINE GMRES_M


!===================================================================================================================================
!> Computes Matrix Vector Product using the spatial operator and finite difference approach, see Dissertation Serena Vangelatos.
!> Computes resu=A*v
!> A is operator at linearization state xk (Newton iteration)
!> Important: needs definition of xk before calling subroutine
!>            needs computation of R_xk before calling subroutine
!===================================================================================================================================
SUBROUTINE MatrixVector(t,Alpha,V,Resu)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mathtools     ,ONLY:GlobalVectorDotProduct
USE MOD_DG_Vars       ,ONLY:U,Ut
USE MOD_Mesh_Vars     ,ONLY:nElems
USE MOD_DG            ,ONLY:DGTimeDerivative_weakForm
USE MOD_Implicit_Vars ,ONLY:Xk,R_Xk,rEps0,Eps_Method,Norm_Xk,FD_Order
USE MOD_Implicit_Vars ,ONLY:nDOFVarProc
USE MOD_TimeDisc_Vars ,ONLY:dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t
REAL,INTENT(IN)   :: Alpha
REAL,INTENT(IN)   :: V(   1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL,INTENT(OUT)  :: Resu(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: V_abs,EpsFD
REAL             :: Ut_plus(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
!===================================================================================================================================
! needed for FD matrix vector approximation
CALL GlobalVectorDotProduct(V,V,nDOFVarProc,V_abs)

! Different approaches to compute an appropriate step size for the FD
SELECT CASE(Eps_Method)
CASE(1) ! Qin, Ludlow, Shaw
  EpsFD = rEps0/SQRT(V_abs) !EpsFD= sqrt(Machine accuracy) / |Delta x|
CASE(2) ! Knoll, Keyes (Eq. (14))
  EpsFD = rEps0/SQRT(V_abs)*SQRT(1.+Norm_Xk)
END SELECT
! Compute pertubated state and DG operator
U = Xk + EpsFD*V
CALL DGTimeDerivative_weakForm(t)

! Actual computation of FD
SELECT CASE(FD_Order)
CASE(1) ! first order FD for approximation of Jacobian
  Resu = V - (Alpha*dt/EpsFD)*(Ut - R_Xk)
CASE(2) ! second order FD for approximation of Jacobian, needs a second pertubation
  Ut_plus = Ut
  U = Xk - EpsFD*V
  CALL DGTimeDerivative_weakForm(t)
  Resu = V - (Alpha*dt/(2.*EpsFD))*(Ut_plus - Ut)
END SELECT
END SUBROUTINE MatrixVector


!===================================================================================================================================
!> Deallocate global variables
!===================================================================================================================================
SUBROUTINE FinalizeImplicit()
! MODULES
USE MOD_Implicit_Vars
USE MOD_Predictor,     ONLY:FinalizePredictor
USE MOD_Precond,       ONLY:FinalizePrecond
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(LinSolverRHS)
SDEALLOCATE(R_Xk)
SDEALLOCATE(Xk)
CALL FinalizePredictor()
CALL FinalizePrecond()
ImplicitInitIsDone = .FALSE.
END SUBROUTINE FinalizeImplicit

END MODULE MOD_Implicit
