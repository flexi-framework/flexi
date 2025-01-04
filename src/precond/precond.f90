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
#if USE_PRECOND
#include "flexi.h"

!===================================================================================================================================
!> Module for the Block-Jacobi Preconditioner
!===================================================================================================================================
MODULE MOD_Precond
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: DefineParametersPrecond
PUBLIC:: InitPrecond
PUBLIC:: BuildPrecond
PUBLIC:: ApplyPrecond
PUBLIC:: FinalizePrecond
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for preconditioner
!==================================================================================================================================
SUBROUTINE DefineParametersPrecond()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Preconditioner")
CALL prms%CreateIntOption(    'PrecondType',     "Preconditioner Type (0: no Preconditioner, 1: analytic, 2: finite difference, 3: &
                                                 & compute both and compare)"                                     , value='1')
CALL prms%CreateIntOption(    'PrecondIter',     "Defines how often preconditioner is built"                      , value='1')
CALL prms%CreateIntOption(    'SolveSystem',     "Solver of the preconditioned system (0: exact LU inversion, 1: inexact ILU(0)    &
                                                 &inversion, always with NoFillIn=T)"                             , value='1')
CALL prms%CreateIntOption(    'DebugMatrix',     "Write Jacobians to file for debug purposes (0: no output, 1: non-inverted matrix,&
                                                 &2: additionally inverted matrix, 3: additionally check inversion accuracy)"      &
                                                                                                                  , value='0')
CALL prms%CreateLogicalOption('HyperbolicPrecond',"Preconditioner only for the hyperbolic flux"                   , value='.FALSE.')
#if PARABOLIC
CALL prms%CreateLogicalOption('NoFillIn'   ,     "Precond for parabolic system forced to have the same sparsity as the &
                                                 &Euler Precond"                                                  , value='.FALSE.')
#endif
CALL prms%CreateLogicalOption('DoDisplayPrecond',"Display building time of preconditioner"                        , value='.FALSE.')

END SUBROUTINE DefineParametersPrecond


!===================================================================================================================================
!> Initialize preconditioner and call initialize of type of preconditioner
!===================================================================================================================================
SUBROUTINE InitPrecond()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Precond_Vars
USE MOD_SparseILU     ,ONLY:InitSparseILU
USE MOD_Mesh_Vars     ,ONLY:nElems
USE MOD_Implicit_Vars ,ONLY:nDOFVarElem
USE MOD_Jac_ex        ,ONLY:InitJac_ex
USE MOD_ReadInTools   ,ONLY:GETINT,GETLOGICAL
#if PP_dim==3
USE MOD_Mesh_Vars     ,ONLY:firstInnerSide,lastInnerSide,SideToElem
#endif
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
IF(PrecondInitIsDone)THEN
   SWRITE(*,*) "InitPrecond already called."
   RETURN
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PRECONDITIONER...'

PrecondType       = GETINT(    'PrecondType')
PrecondIter       = GETINT(    'PrecondIter')
SolveSystem       = GETINT(    'SolveSystem')
DebugMatrix       = GETINT(    'DebugMatrix')
HyperbolicPrecond = GETLOGICAL('HyperbolicPrecond')
#if PARABOLIC
NoFillIn          = GETLOGICAL('NoFillIn')
#endif
DoDisplayPrecond  = GETLOGICAL('DoDisplayPrecond')

! Allocate the preconditioner matrix
ALLOCATE(Ploc(1:nDOFVarElem,1:nDOFVarElem))
IF(PrecondType.EQ.3) ALLOCATE(Ploc1(1:nDOFVarElem,1:nDOFVarElem))

! Initialization of the analytical Preconditioner
IF ((PrecondType.EQ.2).OR.(PrecondType.EQ.3)) THEN
#if PP_dim==3
  ! Abort for 2D periodic meshes, when compiling in 3D. Preconditioner is not working in that case
  DO i=firstInnerSide,lastInnerSide
    IF(SideToElem(S2E_ELEM_ID,i).EQ.SideToElem(S2E_NB_ELEM_ID,i)) THEN
      CALL CollectiveStop(__STAMP__,'ERROR - This is a 2D mesh.')
    ENDIF
  END DO
#endif
END IF

! Initialization of the analytical Preconditioner
IF ((PrecondType.EQ.1).OR.(PrecondType.EQ.3)) CALL InitJac_Ex()

! Method how to solve the preconditioned linear system
IF(PrecondType.NE.0)THEN
  SELECT CASE(SolveSystem)
  CASE(0)
    ! Exact inversion by LU
    ALLOCATE(invP(1:nDOFVarElem,1:nDOFVarElem,nElems))
  CASE(1)
    CALL InitSparseILU()
#if PARABOLIC
    ! Always use the Euler-type sparsity for ILU!
    NoFillIn = .TRUE.
#endif
  END SELECT
END IF

PrecondInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PRECONDITIONER DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitPrecond


!===================================================================================================================================
!> Build preconditioner for each element, calls a type of preconditioner. The block Jacobi preconditioner only takes into account
!> the dependencies of the fluxes in a single element w.r.t. the DOFs of that element, and not the dependencies on DOFs from
!> neighbouring elements! The main advantage is that building and applying the preconditioner is a cell-local operation.
!> For FV, this means we treat subcells on DG cell boundaries differently from the inner cells, although the dependencies would be
!> the same, to keep the block structure in all cases.
!===================================================================================================================================
SUBROUTINE BuildPrecond(t,alpha,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG            ,ONLY:DGTimeDerivative_WeakForm
USE MOD_Implicit_Vars ,ONLY:nDOFVarElem
USE MOD_Jac_ex        ,ONLY:Jac_ex
USE MOD_Jac_FD        ,ONLY:Jac_FD
USE MOD_Mathtools     ,ONLY:Inverse
USE MOD_Mesh_Vars     ,ONLY:nElems
USE MOD_SparseILU     ,ONLY:BuildILU0
USE MOD_Precond_Vars  ,ONLY:invP,DebugMatrix,PrecondType,SolveSystem,DoDisplayPrecond,Ploc,Ploc1
#if USE_MPI
USE MOD_DG_Vars       ,ONLY:U_master,UPrim_master
USE MOD_Mesh_Vars     ,ONLY:nSides
USE MOD_MPI_Vars
USE MOD_MPI           ,ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#if PARABOLIC
USE MOD_Lifting_Vars  ,ONLY:gradUx_master,gradUy_master
#if PP_dim==3
USE MOD_Lifting_Vars  ,ONLY:gradUz_master
#endif /*PP_dim*/
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars ,ONLY:muSGS_master
#endif /*EDDYVISCOSITY*/
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_FV_Vars       ,ONLY:FV_Elems_Sum,FV_Elems_master,FV_Elems_slave
USE MOD_MPI           ,ONLY:StartExchange_FV_Elems
#if FV_RECONSTRUCT
USE MOD_DG_Vars       ,ONLY:U_slave,UPrim_slave
USE MOD_FV_Vars       ,ONLY:FV_dx_master
#endif /*FV_ENABLED*/
#endif /*FV_RECONSTRUCT*/
#endif /*USE_MPI*/
#if FV_ENABLED && FV_RECONSTRUCT
USE MOD_DG_Vars       ,ONLY:UPrim
USE MOD_Jac_Ex_Reconstruction,ONLY:Fill_ExtendedState
USE MOD_Jac_Ex_Vars   ,ONLY:UPrim_extended,FV_sdx_XI_extended,FV_sdx_ETA_extended
#if PP_dim == 3
USE MOD_Jac_Ex_Vars   ,ONLY:FV_sdx_ZETA_extended
#endif /*PP_dim*/
#endif /*FV_ENABLED && FV_RECONSTRUCT*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: alpha   !< RK coefficient for current stage
REAL,INTENT(IN)    :: dt      !< current time step
REAL,INTENT(IN)    :: t       !< current simulation (stage) time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,r,s
REAL               :: Time
INTEGER            :: ind(2)
#if USE_MPI
REAL               :: TimeMPI
#endif /*USE_MPI*/
!===================================================================================================================================
IF(PrecondType.EQ.0) RETURN !NO PRECONDITIONER

! Output of building time for preconditioner
IF(DoDisplayPrecond)THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' BUILD PRECONDITIONER...'
END IF

CALL DGTimeDerivative_WeakForm(t)
#if USE_MPI
! For all side U_master,UPrim_master,U_slave,Uprim_slave,gradU_master,gradU_slave are needed
! U_slave and Uprim_slave were already sent in DGTimeDerivative_WeakForm
! Send MINE - receive YOUR
CALL StartReceiveMPIData(U_master    ,DataSizeSide    ,1,nSides,MPIRequest_U(:,RECV),SendID=1) ! Receive YOUR / U_master
CALL StartSendMPIData(   U_master    ,DataSizeSide    ,1,nSides,MPIRequest_U(:,SEND),SendID=1) ! SEND MINE    / U_master
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)                                            ! U_master: master -> slave
CALL StartReceiveMPIData(UPrim_master,DataSizeSidePrim,1,nSides,MPIRequest_U(:,SEND),SendID=1) ! Receive YOUR / UPrim_master
CALL StartSendMPIData(   UPrim_master,DataSizeSidePrim,1,nSides,MPIRequest_U(:,RECV),SendID=1) ! Send MINE    / UPrim_master
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)                                            ! UPrim_master: master -> slave

#if PARABOLIC
! Sending gradU_master
! gradU_slave was already sent in DGTimeDerivative_WeakForm
! Send MINE - receive YOUR
CALL StartReceiveMPIData(gradUx_master,DataSizeSideGrad,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=1) ! Receive YOUR / gradUx_master
CALL StartReceiveMPIData(gradUy_master,DataSizeSideGrad,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=1) ! Receive YOUR / gradUy_master
#if PP_dim==3
CALL StartReceiveMPIData(gradUz_master,DataSizeSideGrad,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=1) ! Receive YOUR / gradUz_master
#endif
CALL StartSendMPIData(   gradUx_master,DataSizeSideGrad,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=1) ! SEND MINE    / gradUx_master
CALL StartSendMPIData(   gradUy_master,DataSizeSideGrad,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=1) ! SEND MINE    / gradUy_master
#if PP_dim==3
CALL StartSendMPIData(   gradUz_master,DataSizeSideGrad,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=1) ! SEND MINE    / gradUz_master
#endif
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU) ! gradU(x,y,z)_master: master -> slave
#if EDDYVISCOSITY
CALL StartReceiveMPIData(muSGS_master,DataSizeSideSGS,1,nSides,MPIRequest_SGS(:,RECV),SendID=1)
CALL StartSendMPIData   (muSGS_master,DataSizeSideSGS,1,nSides,MPIRequest_SGS(:,SEND),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_SGS)   ! muSGS_master: master -> slave
#endif
#endif /*PARABOLIC*/

#if FV_ENABLED
CALL StartExchange_FV_Elems(FV_Elems_master,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems)                                    ! FV_Elems_master: master -> slave
! Build FV_Elems_Sum on sides (FV_Elems_slave has already been sent to master)
FV_Elems_Sum = FV_Elems_master + 2 * FV_Elems_slave
#if FV_RECONSTRUCT
CALL StartReceiveMPIData(FV_dx_master, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest_Rec_MS(:,SEND),SendID=1)
CALL StartSendMPIData(   FV_dx_master, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest_Rec_MS(:,RECV),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Rec_MS)                                      ! FV_dx_master: master -> slave
! slave states have to be communicated again as they are overwritten (MY Sides) by FV reconstruction with reconstructed values
CALL StartReceiveMPIData(U_slave,    DataSizeSide,    1,nSides,MPIRequest_U(:,SEND),SendID=1) ! Send MINE    / U_slave
CALL StartSendMPIData(   U_slave,    DataSizeSide,    1,nSides,MPIRequest_U(:,RECV),SendID=1) ! Receive YOUR / U_slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)                                           ! U_slave: master -> slave
CALL StartReceiveMPIData(UPrim_slave,DataSizeSidePrim,1,nSides,MPIRequest_U(:,SEND),SendID=1) ! Send MINE    / UPrim_slave
CALL StartSendMPIData(   UPrim_slave,DataSizeSidePrim,1,nSides,MPIRequest_U(:,RECV),SendID=1) ! Receive YOUR / UPrim_slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)                                           ! UPrim_slave: master -> slave
#endif /*FV_RECONSTRUCT*/
#endif /*FV_ENABLED*/
#endif /*USE_MPI*/

#if FV_ENABLED && FV_RECONSTRUCT
IF((PrecondType.EQ.1).OR.(PrecondType.EQ.3))THEN
  ! Build a state vector that contains the primitive values for each element and additionally the first layer of primitive values
  ! from the neighbouring elements. The same is done for the sDx values (lengths for reconstruction).
  CALL Fill_ExtendedState(t,UPrim,UPrim_extended,FV_sdx_XI_extended,FV_sdx_ETA_extended &
#if PP_dim == 3
                          ,FV_sdx_ZETA_extended &
#endif
                           )
END IF
#endif

IF(DoDisplayPrecond) Time=FLEXITIME(MPI_COMM_FLEXI)

DO iElem=1,nElems
  Ploc=0.
  ! compute derivative of DG residual, dRdU
  SELECT CASE(PrecondType)
  CASE(1) ! analytic block Jacobian
    CALL Jac_ex(t,iElem,Ploc,.TRUE.,.TRUE.)
  CASE(2) ! finite difference block Jacobian: use only for debugging
    CALL Jac_FD(t,iElem,Ploc)
  CASE(3) ! debug: difference between analytic and FD preconditioner
    ! Only usable for for LU inversion, since non-exact inversion methods can never reproduce the FD result
    SolveSystem=0
    CALL Jac_FD(t,iElem,Ploc)
    Ploc1=0.
#if FV_ENABLED && FV_RECONSTRUCT && USE_MPI
    ! slave states have to be communicated again as they are overwritten (MY Sides) by fv reconstruction
    CALL StartReceiveMPIData(U_slave,    DataSizeSide,    1,nSides,MPIRequest_U(:,SEND),SendID=1) !Send MINE    / U_slave
    CALL StartSendMPIData(   U_slave,    DataSizeSide,    1,nSides,MPIRequest_U(:,RECV),SendID=1) !Receive YOUR / U_slave
    CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)                   ! U_slave: master -> slave
    CALL StartReceiveMPIData(UPrim_slave,DataSizeSidePrim,1,nSides,MPIRequest_U(:,SEND),SendID=1) !Send MINE    / UPrim_slave
    CALL StartSendMPIData(   UPrim_slave,DataSizeSidePrim,1,nSides,MPIRequest_U(:,RECV),SendID=1) !Receive YOUR / UPrim_slave
    CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)                   ! UPrim_slave: master -> slave
#endif
    CALL Jac_ex(t,iElem,Ploc1,.TRUE.,.TRUE.)
    WRITE(*,*) 'FD= ', SUM(ABS(Ploc))
    WRITE(*,*) 'AD= ', SUM(ABS(Ploc1))
    Ploc1=Ploc1-Ploc
    ind = MAXLOC(ABS(Ploc1))
    IPWRITE(*,*) 'difference: ',MAXVAL(ABS(Ploc1)), &
                ' absolute value: ',ABS(Ploc(ind(1),ind(2))), &
                ' relative difference: ', MAXVAL(ABS(Ploc1))/MAX(ABS(Ploc(ind(1),ind(2))),1.E-16)
    IF(MAXVAL(ABS(Ploc1)).GT.1.0E-4) STOP
  CASE DEFAULT
    CALL Abort(__STAMP__,'No valid preconditioner chosen!')
  END SELECT

  ! add contibution I-alpha*dt*dRdU
  DO s=1,nDOFVarElem
    DO r=1,nDOFVarElem
      Ploc(r,s)=-alpha*dt*Ploc(r,s)
    END DO !r
    Ploc(s,s)=Ploc(s,s)+1.
  END DO !s

  ! Invert the matrix
  SELECT CASE(SolveSystem)
  CASE(0)
    ! Block inverse with Lapack LU factorization
    invP(:,:,iElem)=Inverse(Ploc)
  CASE(1)
    CALL BuildILU0(Ploc,iElem)
  CASE DEFAULT
    CALL Abort(__STAMP__,'No valid linear solver for inverting preconditioner chosen!')
  END SELECT
  IF(DebugMatrix.NE.0) CALL CheckBJPrecond(Ploc,invP(:,:,iElem),iElem)
END DO !iElem

IF(DoDisplayPrecond)THEN
  Time=FLEXITIME(MPI_COMM_FLEXI)-Time
#if USE_MPI
  CALL MPI_REDUCE(Time,TimeMPI,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
  IF(MPIRoot) THEN
    Time=TimeMPI
  END IF
#endif /*USE_MPI*/
  SWRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL DERIVATING & INVERTING TIME =[',Time,' ]'
  SWRITE(UNIT_stdOut,'(A)')' BUILD PRECONDITIONER DONE!'
  SWRITE(UNIT_stdOut,'(132("-"))')
END IF

END SUBROUTINE  BuildPrecond

!===================================================================================================================================
!> Apply the preconditioner to the state vector from the linear solver
!===================================================================================================================================
SUBROUTINE ApplyPrecond(v,z)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY:nElems
USE MOD_Precond_Vars  ,ONLY:SolveSystem,InvP
USE MOD_Implicit_Vars ,ONLY:nDOFVarElem
USE MOD_SparseILU     ,ONLY:ApplyILU
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: v(nDOFVarElem,nElems)    !< non-preconditioned state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: z(nDOFVarElem,nElems)    !< state vector after application of precondtioner
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iElem
!===================================================================================================================================
SELECT CASE(SolveSystem)
CASE(0)
  IF (PP_N.LT.4) THEN
    ! For small N<4 DGEMV is faster than MATMUL
    DO iElem=1,nElems
      CALL DGEMV('N',nDOFVarElem,nDOFVarElem,1.,invP(:,:,iElem),nDOFVarElem,v(:,iElem),1,0.,z(:,iElem), 1)
    END DO !iElem
  ELSE
    DO iElem=1,nElems
      z(:,iElem)=MATMUL(invP(:,:,iElem),v(:,iElem))
    END DO !iElem
  END IF
CASE(1)
  CALL ApplyILU(v,z)
END SELECT
END SUBROUTINE ApplyPrecond


!===================================================================================================================================
!> Debug routine for checking block Jacobian preconditioners. Output options include the non-inverted and the inverted
!> preconditioner matrix. Can also check the inversion.
!===================================================================================================================================
SUBROUTINE CheckBJPrecond(Ploc,invPloc,iElem)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars,      ONLY:nDOFVarElem
USE MOD_Precond_Vars,       ONLY:DebugMatrix,PrecondType
#if USE_MPI
USE MOD_Mesh_Vars,          ONLY:offsetElem
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:nDOFVarElem,1:nDOFVarElem),INTENT(IN) :: Ploc      !< local preconditioner matrix (before inversion)
REAL,DIMENSION(1:nDOFVarElem,1:nDOFVarElem),INTENT(IN) :: invPloc   !< inverted preconditioner
INTEGER,INTENT(IN)                                     :: iElem     !< current element counter
!-----------------------------------------------------------------------------------------------------------------------------------
!LOCAL VARIABLES
REAL               :: dummy
CHARACTER(LEN=255) :: Filename
CHARACTER(LEN=17)  :: strfmt
INTEGER            :: r,s
!===================================================================================================================================
! output of BlockPrecond, no inverse
IF(DebugMatrix.GE.1)THEN
#if USE_MPI
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_Mat_', offsetElem + iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_Mat_', iElem,'.dat'
#endif
  WRITE(strfmt,'(A1,I4,A12)')'(',nDOFVarElem,'(1X,E23.16))'
  WRITE(UNIT_stdOut,*)'Debug Precond (no Inverse) to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
  DO r=1,nDOFVarElem
    WRITE(103,strfmt)Ploc(r,:)
  END DO
  CLOSE(103)
END IF !DebugMatrix >=1

! output of Inverse
IF(DebugMatrix.GE.2)THEN
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_InvMat_', iElem,'.dat'
  WRITE(strfmt,'(A1,I4,A12)')'(',nDOFVarElem,'(1X,E23.16))'
  WRITE(UNIT_stdOut,*)'Debug Precond to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
  DO r=1,nDOFVarElem
    WRITE(103,strfmt)invPloc(r,:)
  END DO
  CLOSE(103)
END IF !DebugMatrix >= 2

! sanity check of inverse
IF((DebugMatrix.GE.3))THEN
  !CHECK INVERSE
  WRITE(UNIT_stdOut,'(A)')'    =>Check invert... '
  dummy=0.
  DO s=1,nDOFVarElem
    DO r=1,nDOFVarElem
      dummy=dummy+ABS(SUM(invPloc(r,:)*Ploc(:,s)))
    END DO
  END DO
  dummy=(dummy-REAL(nDOFVarElem))/REAL(nDOFVarElem*nDOFVarElem) !relative error per matrix entry

  IF(dummy.GT. 1.0E-08) THEN
    IPWRITE(UNIT_stdOut,*)'WARNING!!! accuracy problems in with preconditioner inverse..',dummy
  END IF
  IF(dummy.NE. dummy) THEN  !NAN
    IPWRITE(UNIT_stdOut,*)'WARNING!!! NAN problem in with preconditioner inverse..',dummy
    STOP
  END IF
END IF !DebugMatrix >= 3

END SUBROUTINE CheckBJPrecond


!===================================================================================================================================
!> Finalizes variables
!===================================================================================================================================
SUBROUTINE FinalizePrecond()
! MODULES
USE MOD_Precond_Vars
USE MOD_Jac_Ex,         ONLY:FinalizeJac_Ex
USE MOD_SparseILU,      ONLY:FinalizeSparseILU
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(invP)
SDEALLOCATE(Ploc)
SDEALLOCATE(Ploc1)
IF (SolveSystem.EQ.1) CALL FinalizeSparseILU()
IF ((PreCondType.EQ.1).OR.(PreCondType.EQ.3)) CALL FinalizeJac_Ex()
PrecondInitIsDone = .FALSE.
END SUBROUTINE FinalizePrecond

END MODULE MOD_Precond
#endif /*USE_PRECOND*/
