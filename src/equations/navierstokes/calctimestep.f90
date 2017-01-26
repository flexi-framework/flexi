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
#include "eos.h"

!==================================================================================================================================
!> This module contains the routines to calculate the equation system specific allowable timestep.
!==================================================================================================================================
MODULE MOD_CalcTimeStep
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitCalctimestep
  MODULE PROCEDURE InitCalctimestep
END INTERFACE

INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
END INTERFACE

INTERFACE FinalizeCalctimestep
  MODULE PROCEDURE FinalizeCalctimestep
END INTERFACE


PUBLIC :: InitCalctimestep,CALCTIMESTEP,FinalizeCalctimestep
!==================================================================================================================================

REAL,ALLOCATABLE :: MetricsAdv(:,:,:,:,:,:)  !< support variable: NORM2(Metricsfgh)/J
#if PARABOLIC
REAL,ALLOCATABLE :: MetricsVisc(:,:,:,:,:,:) !< support variable: kappa/Pr*(SUM((Metricsfgh/J)**2))
#endif

CONTAINS

!==================================================================================================================================
!> Precompute some metric support variables
!==================================================================================================================================
SUBROUTINE InitCalctimestep()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
#if PARABOLIC
USE MOD_EOS_Vars ,ONLY:KappasPr
#endif /*PARABOLIC*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem,FVE
#if PARABOLIC
REAL                         :: KappasPr_max
#endif /*PARABOLIC*/
!==================================================================================================================================

ALLOCATE(MetricsAdv(3,0:PP_N,0:PP_N,0:PP_N,nElems,0:FV_ENABLED))
DO FVE=0,FV_ENABLED
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      MetricsAdv(1,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_fTilde(:,i,j,k,iElem,FVE))
      MetricsAdv(2,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_gTilde(:,i,j,k,iElem,FVE))
      MetricsAdv(3,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_hTilde(:,i,j,k,iElem,FVE))
    END DO; END DO; END DO
  END DO
END DO
#if PARABOLIC
ALLOCATE(MetricsVisc(3,0:PP_N,0:PP_N,0:PP_N,nElems,0:FV_ENABLED))
KappasPr_max=KAPPASPR_MAX_TIMESTEP_H()
DO FVE=0,FV_ENABLED
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      MetricsVisc(1,i,j,k,iElem,FVE)=KappasPR_max*(SUM((Metrics_fTilde(:,i,j,k,iElem,FVE)*sJ(i,j,k,iElem,FVE))**2))
      MetricsVisc(2,i,j,k,iElem,FVE)=KappasPR_max*(SUM((Metrics_gTilde(:,i,j,k,iElem,FVE)*sJ(i,j,k,iElem,FVE))**2))
      MetricsVisc(3,i,j,k,iElem,FVE)=KappasPR_max*(SUM((Metrics_hTilde(:,i,j,k,iElem,FVE)*sJ(i,j,k,iElem,FVE))**2))
    END DO; END DO; END DO
  END DO
END DO
#endif /*PARABOLIC*/
END SUBROUTINE


!==================================================================================================================================
!> Compute the time step for the current update of U for the Navier-Stokes-Equations
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
USE MOD_DG_Vars      ,ONLY:U
USE MOD_Mesh_Vars    ,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP,nElems
USE MOD_EOS_Vars
#if PARABOLIC
USE MOD_TimeDisc_Vars,ONLY:DFLScale
USE MOD_Viscosity
#endif /*PARABOLIC*/
USE MOD_TimeDisc_Vars,ONLY:CFLScale,ViscousTimeStep,dtElem
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
#endif
#ifdef EDDYVISCOSITY
USE MOD_EddyVisc_Vars, ONLY: muSGSMax
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: CalcTimeStep
INTEGER,INTENT(OUT)          :: errType
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL,DIMENSION(PP_2Var)      :: UE
REAL                         :: TimeStepConv, TimeStepVisc, TimeStep(3)
REAL                         :: Max_Lambda(3),c,vsJ(3)
#if PARABOLIC
REAL                         :: Max_Lambda_v(3),mu,prim(PP_nVarPrim)
#endif /*PARABOLIC*/
INTEGER                      :: FVE
!==================================================================================================================================
errType=0

TimeStepConv=HUGE(1.)
TimeStepVisc=HUGE(1.)
DO iElem=1,nElems
  FVE = FV_Elems(iElem)
  Max_Lambda=0.
#if PARABOLIC
  Max_Lambda_v=0.
#endif /*PARABOLIC*/
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    ! TODO: ATTENTION: Temperature of UE not filled!!!
    UE(CONS)=U(:,i,j,k,iElem)
    UE(SRHO)=1./UE(DENS)
    UE(VELV)=VELOCITY_HE(UE)
    UE(PRES)=PRESSURE_HE(UE)
    UE(TEMP)=TEMPERATURE_HE(UE)
    ! Convective Eigenvalues
    IF(IEEE_IS_NAN(UE(DENS)))THEN
      ERRWRITE(*,'(A,3ES16.7)')'Density NaN, Position= ',Elem_xGP(:,i,j,k,iElem)
      errType=1
    END IF
    c=SPEEDOFSOUND_HE(UE)
    vsJ=UE(VELV)*sJ(i,j,k,iElem,FVE)
    Max_Lambda(1)=MAX(Max_Lambda(1),ABS(SUM(Metrics_fTilde(:,i,j,k,iElem,FVE)*vsJ)) + &
                                              c*MetricsAdv(1,i,j,k,iElem,FVE))
    Max_Lambda(2)=MAX(Max_Lambda(2),ABS(SUM(Metrics_gTilde(:,i,j,k,iElem,FVE)*vsJ)) + &
                                              c*MetricsAdv(2,i,j,k,iElem,FVE))
    Max_Lambda(3)=MAX(Max_Lambda(3),ABS(SUM(Metrics_hTilde(:,i,j,k,iElem,FVE)*vsJ)) + &
                                              c*MetricsAdv(3,i,j,k,iElem,FVE))
#if PARABOLIC
    ! Viscous Eigenvalues
    prim = UE(PRIM)
    mu=(VISCOSITY_PRIM(prim) &
#ifdef EDDYVISCOSITY
      +50.*muSGSMax(iElem) &
#endif
    )
    Max_Lambda_v=MAX(Max_Lambda_v,mu*UE(SRHO)*MetricsVisc(:,i,j,k,iElem,FVE))
#endif /* PARABOLIC*/
  END DO; END DO; END DO ! i,j,k

  dtElem(iElem)=CFLScale(FVE)*2./SUM(Max_Lambda)
  TimeStepConv=MIN(TimeStepConv,dtElem(iElem))
  IF(IEEE_IS_NAN(TimeStepConv))THEN
    ERRWRITE(*,'(A,I0,A,I0)')'Convective timestep NaN on proc',myRank,' for element: ',iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_conv=',TimeStepConv,' dt_visc=',TimeStepVisc
    errType=2
  END IF

#if PARABOLIC
  IF(SUM(Max_Lambda_v).GT.0.)THEN
    dtElem(iElem)=MIN(dtElem(iElem),DFLScale(FVE)*4./SUM(Max_Lambda_v))
    TimeStepVisc= MIN(TimeStepVisc, DFLScale(FVE)*4./SUM(Max_Lambda_v))
  END IF
  IF(IEEE_IS_NAN(TimeStepVisc))THEN
    ERRWRITE(*,'(A,I0,A,I0)')'Viscous timestep NaN on proc ',myRank,' for element: ', iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_visc=',TimeStepVisc,' dt_conv=',TimeStepConv
    errType=3
  END IF
#endif /* PARABOLIC*/

END DO ! iElem=1,nElems

TimeStep(1)=TimeStepConv
TimeStep(2)=TimeStepVisc
#if MPI
TimeStep(3)=-errType ! reduce with timestep, minus due to MPI_MIN
CALL MPI_ALLREDUCE(MPI_IN_PLACE,TimeStep,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
errType=INT(-TimeStep(3))
#endif /*MPI*/
ViscousTimeStep=(TimeStep(2) .LT. TimeStep(1))
CalcTimeStep=MINVAL(TimeStep(1:2))
END FUNCTION CALCTIMESTEP


!==================================================================================================================================
!> Deallocate CalcTimeStep arrays
!==================================================================================================================================
SUBROUTINE FinalizeCalctimestep()
! MODULES
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(MetricsAdv)
#if PARABOLIC
SDEALLOCATE(MetricsVisc)
#endif
END SUBROUTINE

END MODULE MOD_CalcTimeStep
