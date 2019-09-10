#include "flexi.h"

MODULE MOD_Implicit
!===================================================================================================================================
! Contains the initialization of the Implicit global variables
!===================================================================================================================================
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

#if EQNSYSNR==1
INTERFACE ConjugateGradients
  MODULE PROCEDURE ConjugateGradients
END INTERFACE
#endif

INTERFACE FinalizeImplicit
  MODULE PROCEDURE FinalizeImplicit
END INTERFACE

PUBLIC::InitImplicit,Newton,FinalizeImplicit,DefineParametersImplicit,VectorDotProduct
#if EQNSYSNR==1
PUBLIC::ConjugateGradients
#endif
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
CALL prms%CreateLogicalOption('adapteps',      "adaptive Newton eps", value='.FALSE.')
CALL prms%CreateLogicalOption('adaptRef',      "epsNewton from parameter file with first timesptep CFL=1", value='.FALSE.')
CALL prms%CreateRealOption(  'EpsGMRES',       "GMRES Tolerance", value='1.E-4')
CALL prms%CreateRealOption(  'gammaEW',        "", value='0.9')
CALL prms%CreateRealOption(  'scaleps',        "scale of eps", value='1.')
CALL prms%CreateIntOption(   'nRestarts',      "GMRES Restarts", value='1')
CALL prms%CreateIntOption(   'nKDim',          "Number of Krylov subspaces K^m for GMRES", value='10')
CALL prms%CreateIntOption(   'nNewtonIter',    "Amounts of Newton iterations", value='50')
CALL prms%CreateLogicalOption('EisenstatWalker',    "", value='.FALSE.')
CALL prms%CreateLogicalOption('withmass',    "", value='.FALSE.')
CALL prms%CreateIntOption( 'PredictorType',  "Type of predictor to be used",value='0')
CALL prms%CreateIntOption( 'PredictorOrder', "Order of predictor to be used",value='0')


END SUBROUTINE DefineParametersImplicit



SUBROUTINE InitImplicit()
!===================================================================================================================================
! Allocate global variable 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars
USE MOD_Mesh_Vars,         ONLY:MeshInitIsDone,nElems,nGlobalElems,firstInnerSide,lastInnerSide,SideToElem
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone
USE MOD_ReadInTools,       ONLY:GETINT,GETREAL,GETLOGICAL
USE MOD_DG_Vars,           ONLY:nDOFElem
USE MOD_Analyze_Vars,      ONLY:wGPVol
USE MOD_Predictor,         ONLY:InitPredictor
USE MOD_TimeDisc_Vars,     ONLY:TimeDiscMode
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
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Implicit...'

nDOFVar1D=PP_nVar*(PP_N+1)
nDOFVarElem=PP_nVar*nDOFElem
nDOFGlobal=nDOFVarElem*nElems
nDOFVarGlobal=nDOFElem*nGlobalElems

!Abort for 2D periodic meshes, when compiling in 3D. Preconditioner is not working in that case
#if PP_dim==3
IF(TimeDiscMode.EQ.'Implicit') THEN
  DO i=firstInnerSide,lastInnerSide
    IF(SideToElem(S2E_ELEM_ID,i).EQ.SideToElem(S2E_NB_ELEM_ID,i)) THEN
      CALL CollectiveStop(__STAMP__,'ERROR - This is a 2D mesh.')
    ENDIF
  END DO
END IF !TimeDiscMode
#endif

!========================================================================================
! Variables using for Newton
! Root function: F_Xk = Xk - Q- alpha * dt* R_Xk(t+beta*dt,Xk) = 0!
!========================================================================================
! the constant vector part of non-linear system
! Q contains only solution at the stages before (explicit)
ALLOCATE(Q(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
! Xk is the root of the Newton method
ALLOCATE(Xk(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
! R_Xk is the DG Operator depending on the solution of the current time (implicit)
! In Newton: DG Operator depending on the actual Newton iteration value "xk"
ALLOCATE(R_Xk(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

Q   =0.
Xk  =0.
R_Xk=0.

EpsNewton  =GETREAL('EpsNewton','0.000001')
adapteps   =GETLOGICAL('adapteps','.FALSE.')
adaptRef   =GETLOGICAL('adaptRef','.FALSE.')
scaleps    =GETREAL('scaleps','1.')
nNewtonIter=GETINT('nNewtonIter','50')

IF(adapteps) THEN
  adaptRef=.TRUE.
END IF
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
ApplyPrecondTime=0.
MatrixVectorTime=0.
GMREStime=0.

!If CG solver is used, we have to take Mass into account
ALLOCATE(Mass(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
DO k=0,PP_NZ
  DO j=0,PP_N
    DO i=0,PP_N
      Mass(1:PP_nVar,i,j,k,:)=wGPVol(i,j,k)
    END DO ! i
  END DO ! j
END DO !k
IF(.NOT.GETLOGICAL('withmass','F')) mass=1.

!Predictor for Newton
CALL InitPredictor


ImplicitInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Implicit DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitImplicit


SUBROUTINE Newton(Alpha,imex)
!===================================================================================================================================
! Solves the non-linear system with Newton
! Root function: F_Xk = Xk - Q- alpha * dt* R_Xk(t+beta*dt,Xk) = 0!
! Newton algorithm: 
! dF/dX|Xk * DeltaX = F_Xk
! X_K+1 = X_k + DeltaX
! Attention: we use actual U as X0
!            Q (explicit) is determined outside of Newton by time discretization method
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars       ,ONLY: U,Ut
USE MOD_PreProc
USE MOD_Mesh_Vars     ,ONLY: nElems
USE MOD_DG            ,ONLY: DGTimeDerivative_weakForm
USE MOD_Implicit_Vars ,ONLY: timeStage,Q,EpsNewton,Eps2Newton,Xk,R_Xk,nNewtonIter,nNewtonIterGlobal,nInnerNewton,nDOFVarGlobal
USE MOD_Implicit_Vars ,ONLY: gammaEW,GMRESTime,EpsGMRES,EisenstatWalker,U_stage_old
USE MOD_TimeDisc_Vars ,ONLY: dt
USE MOD_Implicit_Vars ,ONLY: Mass,PredictorType,nGMRESIterGlobal
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)      :: Alpha
INTEGER,INTENT(IN)   :: imex
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: AbortCritNewton,EpsNewtonloc
LOGICAL          :: UsePredictor_loc
REAL             :: U_store(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL             :: F_X0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL             :: F_Xk(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL             :: Norm_F_X0,Norm2_F_X0,Norm2_F_Xk,Norm2_F_Xk_old,Norm_F_Xk
REAL             :: DeltaX(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL             :: eta_k,etaA,etaB,etaMax,eta_k_old
REAL             :: timestart,timeend,Time
!REAL             :: timestartGMRES,timeendGMRES
!===================================================================================================================================
!Initialization for the Newton iterations

IF(PredictorType.NE.0)THEN
  UsePredictor_loc = .TRUE.
ELSE
  UsePredictor_loc = .FALSE.
END IF
IF(UsePredictor_loc)THEN
  !U of old stage (unpredicted)
  Xk        =U_stage_old(:,:,:,:,:,0)
  U_store   =U
  U         =U_stage_old(:,:,:,:,:,0)
  CALL DGTimeDerivative_weakForm(timeStage,imex)
  R_Xk      =Ut
  F_X0      =Mass*(U-Q-Alpha*dt*R_Xk)
  U         =U_store
END IF

!U predicted
Xk        =U
CALL DGTimeDerivative_weakForm(timeStage,imex)
R_Xk      =Ut
F_Xk      =Mass*(U-Q-Alpha*dt*R_Xk)
IF(.NOT.UsePredictor_loc) F_X0=F_Xk

! Preparation for the Abort Criteria of Newton
! |F_Xk| < epsNewton * |F_Xk|
CALL VectorDotProduct(F_X0,F_X0,Norm2_F_X0)
Norm_F_X0=sqrt(Norm2_F_X0)
IF (Norm_F_X0.LE.1.E-13*nDOFVarGlobal) THEN ! do not iterate, as U is already the implicit solution
  !Norm_F_Xk=TINY(1.)
  Norm_F_Xk=0.
ELSE ! we need iterations
  IF(UsePredictor_loc) THEN
    CALL VectorDotProduct(F_Xk,F_Xk,Norm2_F_Xk)
  ELSE
    Norm2_F_Xk=Norm2_F_X0
  END IF 
  Norm_F_Xk=sqrt(Norm2_F_Xk)
END IF
!Norm_F_Xk=sqrt(Norm2_F_Xk)

!Eisenstat and Walker controls the Abort Criteria for GMRES
IF (EisenstatWalker.EQV..TRUE.) THEN
  etaMax   =0.9999
END IF
!AbortCritNewton=Eps2Newton*Norm2_F_X0
AbortCritNewton=EpsNewton*Norm_F_X0
!EpsNewtonloc=SQRT(AbortCritNewton) + 2.7E-16

!===================================================================================================================================
! Newton loop
!===================================================================================================================================
nInnerNewton=0
!DO WHILE((Norm2_F_Xk.GT.AbortCritNewton).AND. (nInnerNewton.LT.nNewtonIter))
!Time = 0.
!SWRITE(*,*) 'Norm F_X0 ', Norm_F_X0
!IF(MPIROOT) THEN
!OPEN(UNIT=42,  &
       !FILE='NewtonResiduum.dat',      &
       !STATUS='UNKNOWN',  &
       !ACTION='WRITE',    &
       !POSITION='APPEND', &
       !IOSTAT=i)
 !END IF
  !SWRITE(42,*) 'Outer Residual,  ', 'CPU time', 'GMRESIter', 'epsGMRES'
  !SWRITE(42,*) Norm_F_Xk, 0, 0, 0

DO WHILE((Norm_F_Xk.GT.AbortCritNewton).AND. (nInnerNewton.LT.nNewtonIter))

  !Timestart=FLEXITIME()
  !DO i=1,20
  ! Computation of the forcing terms eta_k for the Abort Criteria of GMRES via Eisenstat and Walker
  IF(EisenstatWalker.EQV..TRUE.) THEN
    IF (nInnerNewton.EQ.0) THEN
      eta_k=etaMax
    ELSE
      etaA = gammaEW*(Norm2_F_Xk)/(Norm2_F_Xk_old)
      !SWRITE(*,*) 'old ', Norm2_F_Xk_old
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
  !Norm2_F_Xk_old=Norm2_F_Xk
  nInnerNewton=nInnerNewton+1
  !SWRITE(*,*) eta_k

  !======================================================================
  !CALL of the linear solver GMRES
  !TimeStartGMRES=FLEXITIME()
  CALL GMRES_M(Alpha,-F_Xk,Norm_F_Xk,eta_k,DeltaX,imex)
  !TimeEndGMRES=FLEXITIME()
  !GMRESTime = GMRESTime+(TimeEndGMRES-TimeStartGMRES)
  !======================================================================

  Xk  =Xk+DeltaX
  U   =Xk
  CALL DGTimeDerivative_weakForm(timeStage,imex)
  R_Xk=Ut
  F_Xk=mass*(U-Q-alpha*dt*R_Xk)
  CALL VectorDotProduct(F_Xk,F_Xk,Norm2_F_Xk)
  Norm_F_Xk=sqrt(Norm2_F_Xk)

  !if(i<20) THEN
    !R_Xk=R_Xk_old
    !U=U_old
    !Xk=Xk_old
    !F_Xk=F_Xk_old
    !Norm_F_Xk=Norm_F_Xk_old
    !eta_k=eta_k_old
    !Norm2_F_Xk=Norm2_F_Xkoldold
  !END IF
  !SWRITE(*,*) eta_k
  !SWRITE(*,*) Norm2_F_Xk
  !END DO
  !TimeEnd=FLEXITIME()
  !nInnerNewton=nInnerNewton+1
  !Norm2_F_Xk_old=Norm2_F_Xk

  !Time = Time + (TimeEnd-Timestart)/20
  !Time = Time + (TimeEnd-Timestart)
  !Time = TimeEnd-Timestart
  !SWRITE(42,*) Norm_F_Xk, Time, nGMRESIterGlobal, eta_k
  !SWRITE(*,*) Norm_F_Xk, Time, nGMRESIterGlobal, eta_k
  !stop
END DO
!CLOSE(42)

!SWRITE(*,*) 'Time for one Newton Step without GMRES Iterations:', TimeEnd-TimeStart-(TimeEndGMRES-TimeStartGMRES)
nNewtonIterGlobal=nNewtonIterGlobal+nInnerNewton
!IF((Norm2_F_Xk.GT.Eps2Newton*Norm2_F_X0).AND. (nInnerNewton.EQ.nNewtonIter)) THEN
IF((Norm_F_Xk.GT.EpsNewton*Norm_F_X0).AND. (nInnerNewton.EQ.nNewtonIter)) THEN
  CALL abort(__STAMP__, &
  !'NEWTON NOT CONVERGED WITH NEWTON ITERATIONS AND RESIDUAL REDUCTION F_Xk/F_X0:',nInnerNewton,SQRT(Norm2_F_Xk/Norm2_F_X0))
  'NEWTON NOT CONVERGED WITH NEWTON ITERATIONS AND RESIDUAL REDUCTION F_Xk/F_X0:',nInnerNewton,SQRT(Norm_F_Xk/Norm_F_X0))
END IF !nInnerNewton

END SUBROUTINE Newton

#if EQNSYSNR==1
SUBROUTINE ConjugateGradients(alpha,beta)
!===================================================================================================================================
! Matrix free Linear Krylov subspace solver
! A * y = b
! Matrix          A = I - Alpha*dt * dR/dx
! Solution        y = delta x
! Right Hand Side b = -F_Xk
! Inital guess: Delta x = 0
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars       ,ONLY: U,Ut
USE MOD_DG            ,ONLY: DGTimeDerivative_weakForm
USE MOD_Mesh_Vars     ,ONLY: nElems
USE MOD_Implicit_Vars ,ONLY: Q,nNewtonIter,nNewtonIterGlobal
USE MOD_Implicit_Vars ,ONLY: timeStage,Eps2Newton,R_Xk,nDOFVarGlobal
USE MOD_TimeDisc_Vars ,ONLY: dt,t
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: alpha
REAL,INTENT(IN)   :: beta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: x(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL              :: d(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL              :: r(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL              :: z(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL              :: a,b,dz
REAL              :: Norm_r,Norm2_r,Norm2_r_old
REAL              :: timestart,timeend
REAL              :: AbortCrit
INTEGER           :: nInnerCG
!===================================================================================================================================
nInnerCG   =0

!Usage of the Newton abort criterium for the cg method
AbortCrit=Eps2Newton
timeStage =t+beta*dt

x=U
!R_Xkz=d*Laplace*x
CALL DGTimeDerivative_weakForm(timeStage,-1)
R_Xk      =Ut
!z=Ax=x-alpha*dt*R_Xk
z=x-alpha*dt*R_Xk
!Residuum b=Q
r=Q-z
CALL VectorDotProduct(r,r,Norm2_r)

IF (Norm2_r.LE.(1.E-8)**2*nDOFVarGlobal) THEN ! do not iterate, as U is already the implicit solution
  Norm2_r=TINY(1.)
END IF
Norm_r=sqrt(Norm2_r)

!initialization search direction
d=r
!===================================================================================================================================
! CG loop
!===================================================================================================================================
DO WHILE((Norm2_r.GT.AbortCrit).AND. (nInnerCG.LT.nNewtonIter))
  !Timestart=FLEXITIME()

  Norm2_r_old=Norm2_r
  nInnerCG=nInnerCG+1

  U=d
  CALL DGTimeDerivative_weakForm(timeStage,-1)
  R_Xk=Ut
  !z=A*d
  z   =d-alpha*dt*R_Xk
  CALL VectorDotProduct(d,z,dz)
  a=Norm2_r/dz
  x=x+a*d
  r=r-a*z

  CALL VectorDotProduct(r,r,Norm2_r)
  IF ((Norm2_r.LE.AbortCrit)) THEN !converged
    U=x
    CALL DGTimeDerivative_weakForm(timeStage,-1)
    R_Xk=Ut
    nNewtonIterGlobal=nNewtonIterGlobal+nInnerCG
    RETURN
  ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) 
    b=Norm2_r/Norm2_r_old
    d=r+b*d
  END IF ! ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim))
  !TimeEnd=FLEXITIME()
END DO

!SWRITE(*,*) 'Time for one Newton Step without GMRES Iterations:', TimeEnd-TimeStart-(TimeEndGMRES-TimeStartGMRES)
IF((Norm2_r.GT.Eps2Newton*Norm2_r).AND. (nInnerCG.EQ.nNewtonIter)) THEN
  CALL abort(__STAMP__, &
  'CG NOT CONVERGED WITH CG ITERATIONS AND RESIDUAL F_r:',nInnerCG,Norm_r)
END IF !nInnerNewton

END SUBROUTINE ConjugateGradients 
#endif
  

SUBROUTINE GMRES_M(Alpha,B,Norm_B,AbortCrit,DeltaX,imex)
!===================================================================================================================================
! Matrix free Linear Krylov subspace solver
! A * y = b
! Matrix          A = I - Alpha*dt * dR/dx
! Solution        y = delta x
! Right Hand Side b = -F_Xk
! Inital guess: Delta x = 0
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Implicit_Vars ,ONLY:nKDim,nRestarts,nGMRESIterGlobal,nGMRESGlobal,nInnerGMRES,nGMRESIterdt
USE MOD_Implicit_Vars ,ONLY:MatrixVectorTime
USE MOD_Implicit_Vars ,ONLY:ApplyPrecondTime
USE MOD_Implicit_Vars ,ONLY:nDOFGlobal
USE MOD_Precond       ,ONLY:ApplyPrecond
USE MOD_Precond_Vars  ,ONLY:PrecondType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: Alpha,Norm_B
REAL,INTENT(IN)   :: B(nDOFGlobal)
INTEGER,INTENT(IN):: imex
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
REAL              :: timestart,timeend,TotalTime1,TotalTime2
INTEGER           :: m,nn,o,i
REAL              :: tol
REAL              :: deltaz(nDOFGlobal),deltau(nDOFGlobal)
!===================================================================================================================================
!Initializations
R0         =B
Norm_R0    =Norm_B
DeltaX     =0.
Restart    =0
nInnerGMRES=0
TotalTime1 =0.
!TotalTime2 =0.

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
      !TimeStart=FLEXITIME()
      CALL ApplyPrecond(V(:,m),Z(:,m))  
      !TimeEnd=FLEXITIME()
      !TotalTime1=TotalTime1+(TimeEnd-TimeStart)
    ELSE
      Z(:,m)=V(:,m)
    END IF
    ! matrix vector
    !TimeStart=FLEXITIME() 
    CALL MatrixVector(Alpha,Z(:,m),W,imex)
    !TimeEnd=FLEXITIME()
    !TotalTime2=TotalTime2+(TimeEnd-TimeStart)
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
        !ApplyPrecondTime = ApplyPrecondTime+TotalTime1
        !MatrixVectorTime = MatrixVectorTime+TotalTime2
        RETURN
      END IF  ! converged
    ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.tol) .OR. (m.EQ.nKDim)) 
      V(:,m+1)=W/H(m+1,m)
    END IF ! ((ABS(Gam(m+1)).LE.tol) .OR. (m.EQ.nKDim))
  END DO ! m 
  ! Restart needed
  Restart=Restart+1
  CALL MatrixVector(Alpha,DeltaX,R0,imex)
  R0=B-R0
  CALL VectorDotProduct(R0,R0,Norm_R0)
  Norm_R0=SQRT(Norm_R0)
END DO ! While Restart
CALL abort(__STAMP__, &
     'GMRES_M NOT CONVERGED WITH RESTARTS AND GMRES ITERATIONS:',Restart,REAL(nInnerGMRES))
 WRITE(*,*) ' You are using single precision for the preconditioner application!! Try the version with &
              double precision'
END SUBROUTINE GMRES_M


SUBROUTINE MatrixVector(Alpha,V,Resu,imex)
!===================================================================================================================================
! Computes Matrix Vector Product using the spatiall operator and finite difference approach, see Imperator.pdf
! Computes resu=A*v 
! A is operator at linearization state xk (Newton iteration)
! Important: needs definition of xk before calling subroutine
!            needs computation of R_xk before calling subroutine
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars       ,ONLY:U,Ut
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY:nElems
USE MOD_DG            ,ONLY:DGTimeDerivative_weakForm
#if EQNSYSNR==2
USE MOD_Implicit_Vars ,ONLY:Xk,R_Xk,rEps0
#endif
USE MOD_Implicit_Vars ,ONLY:mass,timeStage 
USE MOD_TimeDisc_Vars ,ONLY:dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: Alpha
REAL,INTENT(IN)   :: V(   1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
INTEGER,INTENT(IN):: imex
REAL,INTENT(OUT)  :: Resu(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if EQNSYSNR==2
REAL             :: V_abs,EpsFD
#endif
!===================================================================================================================================
#if EQNSYSNR==1
U=V
CALL DGTimeDerivative_weakForm(timeStage,imex)
Resu = Mass*(V - alpha*dt*Ut)
#endif

#if EQNSYSNR==2
! needed for FD matrix vector approximation
CALL VectorDotProduct(V,V,V_abs)

!EpsFD= sqrt(Machine accuracy) / |Delta x|
EpsFD= rEps0/SQRT(V_abs)
U = Xk+EpsFD*V
CALL DGTimeDerivative_weakForm(timeStage,imex)
Resu = mass*(V - (Alpha*dt/EpsFD)*(Ut - R_Xk))
#endif
END SUBROUTINE MatrixVector

SUBROUTINE VectorDotProduct(A,B,Resu)
!===================================================================================================================================
! Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
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
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif

END SUBROUTINE VectorDotProduct



SUBROUTINE FinalizeImplicit()
!===================================================================================================================================
! Deallocate global variable Q (explicit part of the stages before) and R_Xk (dg time derivative), Xk (newton iteration)
!    mass ( wGP, needed for CG solver)
!===================================================================================================================================
! MODULES
USE MOD_Implicit_Vars,ONLY:Q,Xk,mass,R_Xk
USE MOD_Implicit_Vars,ONLY:ImplicitInitIsDone
USE MOD_Predictor,    ONLY:FinalizePredictor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

SDEALLOCATE(Q)
SDEALLOCATE(R_Xk)
SDEALLOCATE(Xk)
SDEALLOCATE(mass)
CALL FinalizePredictor
ImplicitInitIsDone = .FALSE.
END SUBROUTINE FinalizeImplicit


END MODULE MOD_Implicit
