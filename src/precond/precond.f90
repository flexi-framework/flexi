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
!> Module for the Block-Jacobi Preconditioner  
!===================================================================================================================================
MODULE MOD_Precond
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersPrecond
  MODULE PROCEDURE DefineParametersPrecond 
END INTERFACE

INTERFACE InitPrecond
  MODULE PROCEDURE InitPrecond
END INTERFACE

INTERFACE BuildPrecond
  MODULE PROCEDURE BuildPrecond 
END INTERFACE

INTERFACE FinalizePrecond
  MODULE PROCEDURE FinalizePrecond
END INTERFACE

PUBLIC :: InitPrecond,BuildPrecond,DefineParametersPrecond
PUBLIC :: ApplyPrecond
PUBLIC :: FinalizePrecond
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for preconditioner
!==================================================================================================================================
SUBROUTINE DefineParametersPrecond()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Preconditioner")
CALL prms%CreateIntOption(    'PrecondType',    "Preconditioner Type", value='1')
CALL prms%CreateIntOption(    'PrecondIter',    "Defines how often preconditioner is built", value='1')
CALL prms%CreateIntOption(    'SolveSystem',    "Solver of the preconditioned system", value='0')
CALL prms%CreateIntOption(    'DebugMatrix',    "Debug Matrix", value='0')
CALL prms%CreateLogicalOption('EulerPrecond',   "Preconditioner only for the advective flux", value='.FALSE.')
CALL prms%CreateLogicalOption('NoFillIN'   ,    "Precond has the same Sparsity as the Euler Precond", value='.FALSE.')
CALL prms%CreateLogicalOption('DoDisplayPrecond',"Display building time of preconditioner",'.FALSE.')

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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(PrecondInitIsDone)THEN
   SWRITE(*,*) "InitPrecond already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PRECONDITIONER...'

PrecondType  = GETINT(    'PrecondType','1' )
PrecondIter  = GETINT(    'PrecondIter','1' )
SolveSystem  = GETINT(    'SolveSystem','0' )
DebugMatrix  = GETINT(    'DebugMatrix','0' )
EulerPrecond = GETLOGICAL('EulerPrecond','F')
NoFillIN     = GETLOGICAL('NoFillIN','F'    )
DoDisplayPrecond = GETLOGICAL('DoDisplayPrecond','.FALSE.')

!Initialization of the different Preconditioner Types
SELECT CASE(PrecondType) 
CASE(1,3)
  CALL InitJac_Ex()
END SELECT

! Method how to solve the precontitioned linear system
IF(PrecondType.NE.0)THEN
  SELECT CASE(SolveSystem)
  CASE(0)
    ALLOCATE(invP(1:nDOFVarElem,1:nDOFVarElem,nElems))
  CASE(1)
    CALL InitSparseILU()
    NoFillIN = .TRUE.
  END SELECT
END IF


PrecondInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PRECONDITIONER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPrecond

!===================================================================================================================================
!> Build preconditioner for each element, calls a type of preconditioner 
!===================================================================================================================================
SUBROUTINE BuildPrecond(t,alpha,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mathtools     ,ONLY:Inverse
USE MOD_Mesh_Vars     ,ONLY:nElems,nSides
USE MOD_Implicit_Vars ,ONLY:nDOFVarElem
USE MOD_SparseILU     ,ONLY:BuildILU0
USE MOD_Precond_Vars  ,ONLY:invP,DebugMatrix,PrecondType,SolveSystem,DoDisplayPrecond
USE MOD_Jac_ex        ,ONLY:Jac_ex
USE MOD_Jac_FD        ,ONLY:Jac_FD
USE MOD_DG            ,ONLY:DGTimeDerivative_WeakForm
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI           ,ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_DG_Vars       ,ONLY:U_master,UPrim_master
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
USE MOD_MPI           ,ONLY:StartExchange_FV_Elems
USE MOD_FV_Vars       ,ONLY:FV_Elems_Sum,FV_Elems_master,FV_Elems_slave
#if FV_RECONSTRUCT
USE MOD_DG_Vars       ,ONLY:U_slave,UPrim_slave
USE MOD_FV_Vars       ,ONLY:FV_dx_master
#endif /*FV_ENABLED*/
#endif /*FV_RECONSTRUCT*/
#endif /*MPI*/
#if FV_ENABLED && FV_RECONSTRUCT
USE MOD_DG_Vars       ,ONLY:UPrim
USE MOD_Jac_Reconstruction,ONLY:Fill_ExtendedState
USE MOD_Jac_Ex_Vars   ,ONLY:UPrim_extended,FV_sdx_XI_extended,FV_sdx_ETA_extended
#if PP_dim == 3
USE MOD_Jac_Ex_Vars   ,ONLY:FV_sdx_ZETA_extended
#endif /*PP_dim*/
#endif /*FV_ENABLED && FV_RECONSTRUCT*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: alpha,dt,t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,r,s
REAL,ALLOCATABLE   :: Ploc(:,:)
REAL,ALLOCATABLE   :: Ploc1(:,:)
REAL               :: TimeStart,TimeEnd,TotalTime
INTEGER            :: ind(2)
#if USE_MPI
REAL               :: TotalTimeMPI
#endif /*MPI*/
!===================================================================================================================================
IF(PrecondType.EQ.0) RETURN !NO PRECONDITIONER
IF(DoDisplayPrecond)THEN
  TotalTime=0.
  TimeStart=0.
  TimeEnd=0.
  IF(MPIroot)THEN
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' BUILD PRECONDITIONER...'
  END IF
END IF

CALL DGTimeDerivative_WeakForm(t)
#if USE_MPI
! For all side U_master,UPrim_master,U_slave,Uprim_slave,gradU_master,gradU_slave are needed
! Sending U_master and computing then UPrim_master
! U_slave and Uprim_slave were already sent in DG_TimeDerivative_WeakForm
! Send MINE - receive YOUR
CALL StartReceiveMPIData(U_master    ,DataSizeSide    ,1,nSides,MPIRequest_U(:,RECV),SendID=1) ! Receive YOUR / U_master
CALL StartSendMPIData(   U_master    ,DataSizeSide    ,1,nSides,MPIRequest_U(:,SEND),SendID=1) ! SEND MINE    / U_master
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)        ! U_master: master -> slave 
CALL StartReceiveMPIData(UPrim_master,DataSizeSidePrim,1,nSides,MPIRequest_U(:,SEND),SendID=1) ! Receive YOUR / UPrim_master
CALL StartSendMPIData(   UPrim_master,DataSizeSidePrim,1,nSides,MPIRequest_U(:,RECV),SendID=1) ! Send MINE    / UPrim_master 
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)        ! UPrim_master: master -> slave

#if PARABOLIC
! Sending gradU_master 
! gradU_slave was already sent in DG_TimeDerivative_WeakForm
! Send MINE - receive YOUR
MPIRequest_gradU=MPI_REQUEST_NULL
CALL StartReceiveMPIData(gradUx_master,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,1,RECV),SendID=1) ! Receive YOUR / gradUx_master
CALL StartReceiveMPIData(gradUy_master,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,2,RECV),SendID=1) ! Receive YOUR / gradUy_master
#if PP_dim==3
CALL StartReceiveMPIData(gradUz_master,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,3,RECV),SendID=1) ! Receive YOUR / gradUz_master
#endif

CALL StartSendMPIData(   gradUx_master,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,1,SEND),SendID=1) ! SEND MINE    / gradUx_master
CALL StartSendMPIData(   gradUy_master,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,2,SEND),SendID=1) ! SEND MINE    / gradUy_master
#if PP_dim==3
CALL StartSendMPIData(   gradUz_master,DataSizeSidePrim,1,nSides,MPIRequest_gradU(:,3,SEND),SendID=1) ! SEND MINE    / gradUz_master
#endif
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU)    ! gradU(x,y,z)_master: master -> slave
#if EDDYVISCOSITY
CALL StartReceiveMPIData(muSGS_master,DataSizeSideSGS,1,nSides,MPIRequest_SGS(:,RECV),SendID=1)
CALL StartSendMPIData   (muSGS_master,DataSizeSideSGS,1,nSides,MPIRequest_SGS(:,SEND),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_SGS)  ! muSGS_slave: slave -> master
#endif
#endif /*PARABOLIC*/

#if FV_ENABLED
CALL StartExchange_FV_Elems(FV_Elems_master,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=1) 
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_FV_Elems) ! FV_Elems_master: master -> slave
! Build FV_Elems_Sum on sides (FV_Elems_slave has already been sent to master)
FV_Elems_Sum = FV_Elems_master + 2 * FV_Elems_slave
#if FV_RECONSTRUCT
CALL StartReceiveMPIData(FV_dx_master, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest_Rec_MS(:,SEND),SendID=1)
CALL StartSendMPIData(   FV_dx_master, (PP_N+1)*(PP_NZ+1), 1,nSides,MPIRequest_Rec_MS(:,RECV),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Rec_MS) ! FV_dx_master: master -> slave
! slave states have to be communicated again as they are overwritten (MY Sides) by fv reconstruction with reconstructed values
CALL StartReceiveMPIData(U_slave,    DataSizeSide,    1,nSides,MPIRequest_U(:,SEND),SendID=1) !Send MINE    / U_slave 
CALL StartSendMPIData(   U_slave,    DataSizeSide,    1,nSides,MPIRequest_U(:,RECV),SendID=1) !Receive YOUR / U_slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)                   ! U_slave: master -> slave
CALL StartReceiveMPIData(UPrim_slave,DataSizeSidePrim,1,nSides,MPIRequest_U(:,SEND),SendID=1) !Send MINE    / UPrim_slave 
CALL StartSendMPIData(   UPrim_slave,DataSizeSidePrim,1,nSides,MPIRequest_U(:,RECV),SendID=1) !Receive YOUR / UPrim_slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)                   ! UPrim_slave: master -> slave
#endif
#endif
#endif /*USE_MPI*/

#if FV_ENABLED && FV_RECONSTRUCT
IF((PrecondType.EQ.1).OR.(PrecondType.EQ.3))THEN
  CALL Fill_ExtendedState(t,UPrim,UPrim_extended,FV_sdx_XI_extended,FV_sdx_ETA_extended &
#if PP_dim == 3
                          ,FV_sdx_ZETA_extended &
#endif
                           )
END IF
#endif

ALLOCATE(Ploc(1:nDOFVarElem,1:nDOFVarElem))
IF(PrecondType.EQ.3) ALLOCATE(Ploc1(1:nDOFVarElem,1:nDOFVarElem))

IF(DoDisplayPrecond) TimeStart=FLEXITIME(MPI_COMM_FLEXI)
DO iElem=1,nElems
  Ploc=0.
  !compute derivative of DG residual, dRdU  
  SELECT CASE(PrecondType)
  CASE(1) ! analytic block Jacobian
    CALL Jac_ex(t,iElem,Ploc,.TRUE.,.TRUE.)
  CASE(2) ! finite difference block Jacobian: use only for debugging
    CALL Jac_FD(t,iElem,Ploc)
  CASE(3) ! debug: difference between analytic and fd preconditioner
    !test for LU preconditioner, since the other preconditioners use NoFillIn
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
    IPWRITE(*,*)'Differenz:',MAXVAL(ABS(Ploc1)), &
                'Absolutwert:',ABS(Ploc(ind(1),ind(2))),'Abweichung rel.:', MAXVAL(ABS(Ploc1))/MAX(ABS(Ploc(ind(1),ind(2))),1.E-16)
    IF(MAXVAL(ABS(Ploc1)).GT.1.0E-3) STOP
  CASE DEFAULT
    CALL abort(__STAMP__,'No valid preconditioner chosen!')
  END SELECT

  ! add contibution I-alpha*dt*dRdU
  DO s=1,nDOFVarElem
    DO r=1,nDOFVarElem
      Ploc(r,s)=-alpha*dt*Ploc(r,s)
    END DO !r
    Ploc(s,s)=Ploc(s,s)+1.
  END DO !s
  SELECT CASE(SolveSystem)
  CASE(0)
    ! Block inverse with Lapack LU factorization
    invP(:,:,iElem)=Inverse(Ploc)
  CASE(1)
    CALL BuildILU0(Ploc,iElem)
  CASE DEFAULT
    CALL abort(__STAMP__,'No valid linear solver for inverting preconditioner chosen!')
  END SELECT
  IF(DebugMatrix.NE.0)THEN
     CALL CheckBJPrecond(Ploc,Ploc,invP(:,:,iElem),iElem)
  END IF ! DebugMatrix
END DO !iElem

IF(DoDisplayPrecond)THEN
  TimeEnd=FLEXITIME(MPI_COMM_FLEXI)
  TotalTime=TotalTime+(TimeEnd-TimeStart)
#if USE_MPI
  CALL MPI_REDUCE(TotalTime,TotalTimeMPI ,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,IERROR)
  IF(MPIRoot) THEN
    TotalTime=TotalTimeMPI
  END IF
#endif /*MPI*/
  IF(MPIRoot)THEN
    WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL DERIVATING & INVERTING TIME =[',TotalTime,' ]'
    WRITE(UNIT_stdOut,'(A)')' BUILD PRECONDITIONER DONE!'
    WRITE(UNIT_StdOut,'(132("-"))')
  END IF
END IF

DEALLOCATE(Ploc)
IF(PrecondType.EQ.3) DEALLOCATE(Ploc1)
END SUBROUTINE  BuildPrecond

!===================================================================================================================================
!> Build preconditioner for each element, calls a type of preconditioner
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
REAL,INTENT(IN)  :: v(nDOFVarElem,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: z(nDOFVarElem,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iElem
!===================================================================================================================================
SELECT CASE(SolveSystem)
CASE(0)
  ! For small N<4 DGEMV is faster than MATMUL
  !DO iElem=1,nElems
  !CALL DGEMV('N',nDOFVarElem, nDOFVarElem,1.,invP(:,:,iElem),nDOFVarElem,v(:,iElem),1,0.,z(:,iElem), 1)
  !END DO !iElem
  DO iElem=1,nElems
    z(:,iElem)=MATMUL(invP(:,:,iElem),v(:,iElem))
  END DO !iElem
CASE(1)
  CALL ApplyILU(v,z)
END SELECT
END SUBROUTINE  ApplyPrecond

!===================================================================================================================================
!> Debug routine for checking block jacobian preconditioners
!===================================================================================================================================
SUBROUTINE CheckBJPrecond(Ploc1,Ploc,invPloc,iElem)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Implicit_Vars,      ONLY:nDOFVarElem
USE MOD_Precond_Vars,       ONLY:DebugMatrix,PrecondType
USE MOD_Mesh_Vars,          ONLY:offsetElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:nDOFVarElem,1:nDOFVarElem),INTENT(IN) :: Ploc,invPloc,Ploc1
INTEGER,INTENT(IN)                                     :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
!LOCAL VARIABLES
!REAL               :: diff(1:nDOFVarElem,1:nDOFVarElem)
REAL               :: dummy
CHARACTER(LEN=255) :: Filename
CHARACTER(LEN=17)  :: strfmt
INTEGER            :: r,s,counter
REAL               :: sparsity
REAL,DIMENSION(1:nDOFVarElem,1:nDOFVarElem)  :: Ploc_loc
!===================================================================================================================================

! output of BlockJac  | only derivativ
IF(DebugMatrix.GT.0)THEN
#if USE_MPI
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond',PreCondType,'_Mat_', offsetElem + iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond',PreCondType,'_Mat_', iElem,'.dat'
#endif
  WRITE(strfmt,'(A1,I4,A12)')'(',nDOFVarElem,'(1X,E23.16))'

  Ploc_loc=Ploc
  IF(iElem.EQ.4) THEN
  !IF(iElem.EQ.15) THEN
    counter=0
    DO r=1,nDOFVarElem
      DO s=1,nDOFVarElem
        IF(abs(Ploc(r,s)).GT.1.E-15) THEN
        !IF(abs(Ploc(r,s)).LE.1.E-15) THEN
          !Ploc_loc(r,s) = 0.
          counter = counter +1
        END IF
      END DO
    END DO
    WRITE(Filename,'(A,I2.2,A)')'Storage_Precond',PreCondType,'.dat'
    OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND')
      WRITE(103,*) PP_N, counter*0.0078125
    CLOSE(103)
    sparsity= REAL(counter)/(nDOFVarElem**2)
    !WRITE(*,*) 'sparsity= ', sparsity 
    !WRITE(*,*) 'density= ', 1.-sparsity
    ! WRITE(UNIT_stdOut,*)'Debug Block Jacobian to:',TRIM(Filename)
    !OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
    !DO r=1,nDOFVarElem
      !WRITE(103,strfmt)Ploc_loc(r,:)
    !END DO
    !CLOSE(103)
  END IF
END IF !DebugMatrix >0

! output of BlockPrecond, no inverse
IF(DebugMatrix.GT.1)THEN
#if USE_MPI
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_Mat_', offsetElem + iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_Mat_', iElem,'.dat'
#endif
  WRITE(strfmt,'(A1,I4,A12)')'(',nDOFVarElem,'(1X,E23.16))'
  WRITE(UNIT_stdOut,*)'Debug Precond (no Inverse) to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
  DO r=1,nDOFVarElem
    WRITE(103,strfmt)Ploc1(r,:)
  END DO
  CLOSE(103)
END IF !DebugMatrix >1

! output of Inverse
IF(DebugMatrix.GT.2)THEN
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_InvMat_', iElem,'.dat'
  WRITE(strfmt,'(A1,I4,A12)')'(',nDOFVarElem,'(1X,E23.16))'
  WRITE(UNIT_stdOut,*)'Debug Precond to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
  DO r=1,nDOFVarElem
    WRITE(103,strfmt)invPloc(r,:)
  END DO
  CLOSE(103)
END IF !DebugMatrix>2

! sanity check of inverse
IF((DebugMatrix.GT.3).OR.(DebugMatrix.LT.0))THEN
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
END IF !DebugMatrix>3

END SUBROUTINE CheckBJPrecond

SUBROUTINE FinalizePrecond()
!===================================================================================================================================
! Finalizes variables 
!===================================================================================================================================
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
CALL FinalizeJac_Ex()
SDEALLOCATE(invP)
IF(SolveSystem.EQ.1) CALL FinalizeSparseILU()
PrecondInitIsDone = .FALSE.
END SUBROUTINE FinalizePrecond

END MODULE MOD_Precond
