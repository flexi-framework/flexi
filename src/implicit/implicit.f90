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
! Contains the initialization of the Implicit global variables
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

PUBLIC::InitImplicit,Newton,FinalizeImplicit,DefineParametersImplicit,VectorDotProduct
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersImplicit()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Implicit")
CALL prms%CreateRealOption(  'EpsNewton',      "Newton tolerance", value='1.E-6')
CALL prms%CreateIntOption(   'nNewtonIter',    "Amounts of Newton iterations", value='50')
CALL prms%CreateLogicalOption('adaptepsNewton',      "adaptive Newton eps", value='.FALSE.')
CALL prms%CreateLogicalOption('EisenstatWalker',    "", value='.FALSE.')
CALL prms%CreateRealOption(  'gammaEW',        "", value='0.9')
CALL prms%CreateRealOption(  'EpsGMRES',       "GMRES Tolerance", value='1.E-4')
CALL prms%CreateIntOption(   'nRestarts',      "GMRES Restarts", value='1')
CALL prms%CreateIntOption(   'nKDim',          "Number of Krylov subspaces K^m for GMRES", value='10')
CALL prms%CreateRealOption(  'scaleps',        "scale of eps", value='1.')
CALL prms%CreateIntOption(   'Eps_Method',     "how eps of finite difference is calulated", value='2')
CALL prms%CreateLogicalOption('withmass',    "", value='.FALSE.')
CALL prms%CreateIntOption( 'PredictorType',  "Type of predictor to be used",value='0')
CALL prms%CreateIntOption( 'PredictorOrder', "Order of predictor to be used",value='0')

END SUBROUTINE DefineParametersImplicit

!===================================================================================================================================
!> Allocate global variable 
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
USE MOD_Analyze_Vars,      ONLY:wGPVol
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
INTEGER                   :: i,j,k
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.ImplicitInitIsDone)THEN
   CALL abort(__STAMP__,'InitImplicit not ready to be called or already called.')
END IF
IF(TimeDiscType.EQ.'ESDIRK') THEN
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' INIT Implicit...'

  nDOFVar1D=PP_nVar*(PP_N+1)
  nDOFVarElem=PP_nVar*nDOFElem
  nDOFGlobal=nDOFVarElem*nElems
  nDOFVarGlobal=nDOFElem*nGlobalElems

  !Abort for 2D periodic meshes, when compiling in 3D. Preconditioner is not working in that case
#if PP_dim==3
  DO i=firstInnerSide,lastInnerSide
    IF(SideToElem(S2E_ELEM_ID,i).EQ.SideToElem(S2E_NB_ELEM_ID,i)) THEN
      CALL CollectiveStop(__STAMP__,'ERROR - This is a 2D mesh.')
    ENDIF
  END DO
#endif

  !========================================================================================
  ! Variables using for Newton
  ! Root function: F_Xk = Xk - Q- alpha * dt* R_Xk(t+beta*dt,Xk) = 0!
  !========================================================================================
  ! the constant vector part of non-linear system
  ALLOCATE(LinSolverRHS(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  ! Xk is the root of the Newton method
  ALLOCATE(Xk(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  ! R_Xk is the DG Operator depending on the solution of the current time (implicit)
  ! In Newton: DG Operator depending on the actual Newton iteration value "xk"
  ALLOCATE(R_Xk(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  LinSolverRHS = 0.
  Xk           = 0.
  R_Xk         = 0.

  EpsNewton      = GETREAL('EpsNewton','0.001')
  adaptepsNewton = GETLOGICAL('adaptepsNewton','.FALSE.')
  IF((.NOT.ALLOCATED(RKb_embedded)).AND.(adaptepsNewton.EQV..TRUE.))THEN
    SWRITE(UNIT_stdOut,'(A)') ' Attention, Newton abort criterium is not adaptive!'
    adaptepsNewton = .FALSE.
  END IF
  scaleps        = GETREAL('scaleps','1.')
  ! Choose method how eps is calculated:
  ! 1: According to Qin,Ludlow,Shaw: A matrix-free preconditioned Newton/GMRES method for unsteady Navier-Stokes solutions (Eq. 13),
  !                                  Int.J.Numer.Meth.Fluids 33 (2000) 223-248
  ! 2: According to Knoll,Keyes: Jacobian-free Newton-Krylov methods: A survey of approaches and applications (Eq. (14)),
  !                              JCP 193 (2004) 357-397
  Eps_Method     = GETINT('Eps_Method','2')
  nNewtonIter    = GETINT('nNewtonIter','50')

  ! Newton takes the quadratic norm into account
  Eps2Newton=EpsNewton**2
  rEps0            =scaleps*SQRT(EPSILON(0.0))
  srEps0           =1./rEps0
  nInnerNewton     =0
  nNewtonIterGlobal=0
  nGMRESGlobal     =0

  !========================================================================================
  ! Variables using for GMRES 
  ! GMRES solving the LES: A * y = b
  ! Matrix          A = I - alpha*dt * dR/dx
  ! Solution        y = dx
  ! Right Hand Side b = -F_Xk
  !========================================================================================
  nRestarts       =GETINT    ('nRestarts','1')
  nKDim           =GETINT    ('nKDim','10')
  EpsGMRES        =GETREAL   ('EpsGMRES','0.0001')
  EisenstatWalker =GETLOGICAL('EisenstatWalker','.FALSE.')
  gammaEW         =GETREAL   ('gammaEW','0.9')


  nGMRESIterGlobal=0
  nGMRESIterdt=0

  !If CG solver is used, we have to take Mass into account
  ALLOCATE(Mass(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  IF(.NOT.GETLOGICAL('withmass','F'))THEN
    mass=1.
  ELSE
    DO k=0,PP_NZ
      DO j=0,PP_N
        DO i=0,PP_N
          Mass(1:PP_nVar,i,j,k,:)=wGPVol(i,j,k)
        END DO ! i
      END DO ! j
    END DO !k
  END IF

  CALL InitPrecond()

  ImplicitInitIsDone=.TRUE.
  SWRITE(UNIT_stdOut,'(A)')' INIT Implicit DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')
END IF
END SUBROUTINE InitImplicit

!===================================================================================================================================
!> Solves the non-linear system with Newton
!> Root function: F_Xk = Xk - Q- alpha * dt* R_Xk(t+beta*dt,Xk) = 0!
!> Newton algorithm: 
!> dF/dX|Xk * DeltaX = F_Xk
!> X_K+1 = X_k + DeltaX
!> Attention: we use actual U as X0
!===================================================================================================================================
SUBROUTINE Newton(t,Alpha)
! MODULES
USE MOD_DG_Vars       ,ONLY: U,Ut
USE MOD_PreProc
USE MOD_Mesh_Vars     ,ONLY: nElems
USE MOD_DG            ,ONLY: DGTimeDerivative_weakForm
USE MOD_Implicit_Vars ,ONLY: EpsNewton,Xk,R_Xk,nNewtonIter,nNewtonIterGlobal,nInnerNewton,nDOFVarGlobal!,Eps2Newton
USE MOD_Implicit_Vars ,ONLY: gammaEW,EpsGMRES,EisenstatWalker,LinSolverRHS,U_predictor
USE MOD_TimeDisc_Vars ,ONLY: dt
USE MOD_Implicit_Vars ,ONLY: Mass,PredictorType,Eps_Method,Norm_Xk
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)      :: t
REAL,INTENT(IN)      :: Alpha
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

IF(Eps_Method.EQ.2)THEN
  CALL VectorDotProduct(Xk,Xk,Norm_Xk)
  Norm_Xk = SQRT(Norm_Xk)
END IF
! Preparation for the Abort Criteria of Newton
! |F_Xk| < epsNewton * |F_Xk|
CALL VectorDotProduct(F_X0,F_X0,Norm2_F_X0)
Norm_F_X0=SQRT(Norm2_F_X0)
IF (Norm_F_X0.LE.1.E-13*nDOFVarGlobal) THEN ! do not iterate, as U is already the implicit solution
  Norm_F_Xk=0.
ELSE ! we need iterations
  IF(PredictorType.NE.0)THEN
    CALL VectorDotProduct(F_Xk,F_Xk,Norm2_F_Xk)
  ELSE
    Norm2_F_Xk=Norm2_F_X0
  END IF 
  Norm_F_Xk=SQRT(Norm2_F_Xk)
END IF

!Eisenstat and Walker controls the Abort Criteria for GMRES
IF (EisenstatWalker.EQV..TRUE.) THEN
  etaMax   =0.9999
END IF

!AbortCritNewton=Eps2Newton*Norm2_F_X0
AbortCritNewton=EpsNewton*Norm_F_X0

!===================================================================================================================================
! Newton loop
!===================================================================================================================================
nInnerNewton=0
!DO WHILE((Norm2_F_Xk.GT.AbortCritNewton).AND. (nInnerNewton.LT.nNewtonIter))
DO WHILE((Norm_F_Xk.GT.AbortCritNewton).AND. (nInnerNewton.LT.nNewtonIter))

  ! Computation of the forcing terms eta_k for the Abort Criteria of GMRES via Eisenstat and Walker
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
  !CALL of the linear solver GMRES
  CALL GMRES_M(t,Alpha,-F_Xk,Norm_F_Xk,eta_k,DeltaX)
  !======================================================================

  Xk  =Xk+DeltaX
  U   =Xk
  CALL DGTimeDerivative_weakForm(t)
  R_Xk=Ut
  F_Xk=mass*(U-LinSolverRHS-alpha*dt*R_Xk)
  CALL VectorDotProduct(F_Xk,F_Xk,Norm2_F_Xk)
  Norm_F_Xk=SQRT(Norm2_F_Xk)
END DO

nNewtonIterGlobal=nNewtonIterGlobal+nInnerNewton
!IF((Norm2_F_Xk.GT.Eps2Newton*Norm2_F_X0).AND. (nInnerNewton.EQ.nNewtonIter)) THEN
IF((Norm_F_Xk.GT.EpsNewton*Norm_F_X0).AND. (nInnerNewton.EQ.nNewtonIter)) THEN
  CALL abort(__STAMP__, &
  !'NEWTON NOT CONVERGED WITH NEWTON ITERATIONS AND RESIDUAL REDUCTION F_Xk/F_X0:',nInnerNewton,SQRT(Norm2_F_Xk/Norm2_F_X0))
  'NEWTON NOT CONVERGED WITH NEWTON ITERATIONS AND RESIDUAL REDUCTION F_Xk/F_X0:',nInnerNewton,SQRT(Norm_F_Xk/Norm_F_X0))
END IF !nInnerNewton

END SUBROUTINE Newton

!===================================================================================================================================
!> Matrix free Linear Krylov subspace solver
!> A * y = b
!> Matrix          A = I - Alpha*dt * dR/dx
!> Solution        y = delta x
!> Right Hand Side b = -F_Xk
!> Inital guess: Delta x = 0
!===================================================================================================================================
SUBROUTINE GMRES_M(t,Alpha,B,Norm_B,AbortCrit,DeltaX)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Implicit_Vars ,ONLY:nKDim,nRestarts,nGMRESIterGlobal,nGMRESGlobal,nInnerGMRES,nGMRESIterdt
USE MOD_Implicit_Vars ,ONLY:nDOFGlobal
USE MOD_Precond       ,ONLY:ApplyPrecond
USE MOD_Precond_Vars  ,ONLY:PrecondType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t
REAL,INTENT(IN)   :: Alpha,Norm_B
REAL,INTENT(IN)   :: B(nDOFGlobal)
REAL,INTENT(IN)   :: AbortCrit
REAL,INTENT(OUT)  :: DeltaX(nDOFGlobal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: V(nDOFGlobal,1:nKDim)
REAL              :: W(nDOFGlobal)
REAL              :: Z(nDOFGlobal,nKDim)
REAL              :: R0(nDOFGlobal)
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
      CALL VectorDotProduct(V(:,nn),W,H(nn,m))
      W=W-H(nn,m)*V(:,nn)
    END DO !nn
    CALL VectorDotProduct(W,W,Resu)
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
        nGMRESGlobal=nGMRESGlobal+Restart+1 
        nGMRESIterGlobal=nGMRESIterGlobal+nInnerGMRES 
        nGMRESIterdt=nGMRESIterdt + nInnerGMRES
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
  CALL VectorDotProduct(R0,R0,Norm_R0)
  Norm_R0=SQRT(Norm_R0)
END DO ! While Restart
CALL abort(__STAMP__, &
     'GMRES_M NOT CONVERGED WITH RESTARTS AND GMRES ITERATIONS:',Restart,REAL(nInnerGMRES))
END SUBROUTINE GMRES_M


!===================================================================================================================================
!> Computes Matrix Vector Product using the spatiall operator and finite difference approach, see Imperator.pdf
!> Computes resu=A*v 
!> A is operator at linearization state xk (Newton iteration)
!> Important: needs definition of xk before calling subroutine
!>            needs computation of R_xk before calling subroutine
!===================================================================================================================================
SUBROUTINE MatrixVector(t,Alpha,V,Resu)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars       ,ONLY:U,Ut
USE MOD_Mesh_Vars     ,ONLY:nElems
USE MOD_DG            ,ONLY:DGTimeDerivative_weakForm
USE MOD_Implicit_Vars ,ONLY:Xk,R_Xk,rEps0,Eps_Method,Norm_Xk
USE MOD_Implicit_Vars ,ONLY:mass
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
!===================================================================================================================================
! needed for FD matrix vector approximation
CALL VectorDotProduct(V,V,V_abs)

SELECT CASE(Eps_Method)
CASE(1) ! Qin, Ludlow, Shaw
  EpsFD = rEps0/SQRT(V_abs) !EpsFD= sqrt(Machine accuracy) / |Delta x|
CASE(2) ! Knoll, Keyes (Eq. (14))
  EpsFD = rEps0/SQRT(V_abs)*SQRT(1.+Norm_Xk)
END SELECT
U = Xk + EpsFD*V
CALL DGTimeDerivative_weakForm(t)
Resu = mass*(V - (Alpha*dt/EpsFD)*(Ut - R_Xk))
END SUBROUTINE MatrixVector

!===================================================================================================================================
!> Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
SUBROUTINE VectorDotProduct(A,B,Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars ,ONLY:nDOFGlobal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: A(nDOFGlobal)
REAL,INTENT(IN)   :: B(nDOFGlobal)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
!===================================================================================================================================
Resu=0.
DO i=1,nDOFGlobal
  Resu=Resu + A(i)*B(i)
END DO

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
#endif

END SUBROUTINE VectorDotProduct

!===================================================================================================================================
!> Deallocate global variable Q (explicit part of the stages before) and R_Xk (dg time derivative), Xk (newton iteration)
!> mass ( wGP, needed for CG solver)
!===================================================================================================================================
SUBROUTINE FinalizeImplicit()
! MODULES
USE MOD_Implicit_Vars, ONLY:LinSolverRHS,Xk,mass,R_Xk
USE MOD_Implicit_Vars, ONLY:ImplicitInitIsDone
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
SDEALLOCATE(mass)
CALL FinalizePredictor()
CALL FinalizePrecond()
ImplicitInitIsDone = .FALSE.
END SUBROUTINE FinalizeImplicit

END MODULE MOD_Implicit
