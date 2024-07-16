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
#include "flexi.h"

!==================================================================================================================================
!> Routines providing initialization and initial solutions for the linear advection-diffusion equation
!==================================================================================================================================
MODULE MOD_Exactfunc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersExactFunc
  MODULE PROCEDURE DefineParametersExactFunc
END INTERFACE

INTERFACE InitExactFunc
  MODULE PROCEDURE InitExactFunc
END INTERFACE

INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc
END INTERFACE

INTERFACE CalcSource
  MODULE PROCEDURE CalcSource
END INTERFACE


PUBLIC::DefineParametersExactFunc
PUBLIC::InitExactFunc
PUBLIC::ExactFunc
PUBLIC::CalcSource
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersExactFunc()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Exactfunc")
CALL prms%CreateIntFromStringOption(      'IniExactFunc', "Number of exactfunction to be used, to initialize the solution. "//&
                                                          "1: linear, 2: inviscid convtest sine")
CALL addStrListEntry('IniExactFunc','linear0'  ,1)
CALL addStrListEntry('IniExactFunc','linear10' ,21)
CALL addStrListEntry('IniExactFunc','constant' ,101)
CALL addStrListEntry('IniExactFunc','sine'     ,2)
CALL addStrListEntry('IniExactFunc','sinex'    ,41)
CALL addStrListEntry('IniExactFunc','sinetime' ,31)
CALL addStrListEntry('IniExactFunc','quadratic',4)
CALL addStrListEntry('IniExactFunc','advdiff'  ,5)
CALL addStrListEntry('IniExactFunc','dissdisp' ,6)
CALL prms%CreateRealOption('OmegaRef',     "Angular frequency for one-dimensional cosine wave (IniExactFunc 6).")
END SUBROUTINE DefineParametersExactFunc

!==================================================================================================================================
!> Get some parameters needed for exact function
!==================================================================================================================================
SUBROUTINE InitExactFunc()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools,   ONLY: GETINTFROMSTR,GETREAL,GETINT
USE MOD_ExactFunc_Vars
USE MOD_Equation_Vars, ONLY: AdvVel,IniExactFunc,IniRefState
#if PARABOLIC
USE MOD_Equation_Vars, ONLY: DiffC
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT EXACT FUNCTION...'

! Read in boundary parameters
IniExactFunc = GETINTFROMSTR('IniExactFunc')
IniRefState  = -1 ! only dummy for linadv

! Read in parameters specific to certain init functions
SELECT CASE (IniExactFunc)
CASE(6) ! One-Dimensional cosine wave with angular frequency set by user
  omegaRef = GETREAL('OmegaRef')
  ! Set advection velocity to (1,0,0) and diffusion coefficient to 0
  AdvVel(1)   = 1.
  AdvVel(2:3) = 0.
#if PARABOLIC
  DiffC       = 0.
#endif
END SELECT

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT FUNCTION DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitExactFunc

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu,RefStateOpt)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars, ONLY: AdvVel
USE MOD_Exactfunc_Vars,ONLY: OmegaRef
USE MOD_Timedisc_Vars, ONLY: fullBoundaryOrder,CurrentStage,dt,RKb,RKc,t
#if PARABOLIC
USE MOD_Equation_Vars, ONLY: DiffC
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: tIn                    !< input time (either time at RK stage or time at the beginning of
                                                          !< timestep if full boundary order is used (only with RK3)
REAL,INTENT(IN)                 :: x(3)                   !< coordinates to evaluate exact function
INTEGER,INTENT(IN)              :: ExactFunction          !< specifies the exact function to be used
REAL,INTENT(OUT)                :: Resu(PP_nVar)          !< output state in conservative variables
INTEGER,INTENT(IN),OPTIONAL     :: RefStateOpt            !< refstate to be used for exact func (dummy)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: tEval
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar) ! temporal state deriv in conservative variables
REAL                            :: Frequency,Amplitude,Omega
REAL                            :: Cent(3),x0(3)
REAL                            :: Pi
REAL                            :: k
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration

Pi = ACOS(-1.)

Resu   =0.
Resu_t =0.
Resu_tt=0.

Cent=x-AdvVel*tEval
SELECT CASE (ExactFunction)
CASE(1) !linear
  Resu=0.+SUM(Cent)
  IF(fullBoundaryOrder)THEN
    Resu_t=-SUM(Advvel)
  END IF
CASE(21) !linear
  Resu=10.+SUM(Cent)
  IF(fullBoundaryOrder)THEN
    Resu_t=-SUM(Advvel)
  END IF
CASE(101) !constant
  Resu=7.7
CASE(2) !sinus
  Frequency=0.5
  Amplitude=1.
  Omega=2.*PP_PI*Frequency
  Resu=Amplitude*SIN(Omega*SUM(Cent))
  IF(fullBoundaryOrder)THEN
    Resu_t =-Amplitude*COS(Omega*SUM(Cent))*Omega*SUM(AdvVel)
    Resu_tt=-Amplitude*SIN(Omega*SUM(Cent))*Omega*SUM(AdvVel)*Omega*SUM(AdvVel)
  END IF
CASE(22) ! two-dimensional sinus with source for convergence test
  Frequency=1.
  Amplitude=1.
  Omega=2.*PP_PI*Frequency
  Resu=Amplitude*SIN(Omega*SUM(Cent))
  IF(fullBoundaryOrder)THEN
    Resu_t =-Amplitude*COS(Omega*SUM(Cent))*Omega*SUM(AdvVel)
    Resu_tt=-Amplitude*SIN(Omega*SUM(Cent))*Omega*SUM(AdvVel)*Omega*SUM(AdvVel)
  END IF
CASE(41) ! sinus in x direction
  Frequency=0.5
  Amplitude=1.
  Omega=2.*PP_PI*Frequency
  Resu=Amplitude*SIN(Omega*Cent(1))
  IF(fullBoundaryOrder)THEN
    Resu_t =-Amplitude*COS(Omega*Cent(1))*Omega*AdvVel(1)
    Resu_tt=-Amplitude*SIN(Omega*Cent(1))*Omega*AdvVel(1)*Omega*AdvVel(1)
  END IF
CASE(31)
  !Resu=Cent(1)
  Resu=SIN(Pi*(x(1)-AdvVel(1)/2.*t))
CASE(4) ! quadratic
  Resu=5.*(x(1)-AdvVel(1)*tEval)**2
  IF(fullBoundaryOrder)THEN
    Resu_t=-10.*AdvVel(1)*(x(1)-AdvVel(1)*tEval)
    Resu_tt=10.*AdvVel(1)*AdvVel(1)
  END IF
#if PARABOLIC
CASE(5) ! Kopriva page 200, advection-diffusion, but for 3D with 1/( (4t+1)^(3/2) )
  x0 = (/-0.5,-0.5,-0.5/)
  Resu   =1./((4.*t+1.)**(1.5))*EXP(-(SUM((x(:)-AdvVel(:)*t-x0(:))**2))/(DiffC*(4.*t+1.)))
  Resu_t =Resu  *(-6./(4.*t+1.) &
                  +2./(DiffC*(4.*t+1.)   )*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))) &
                  +4./(DiffC*(4.*t+1.)**2)*SUM((x(:)-AdvVel(:)*t-x0(:))**2) )
  Resu_tt=Resu_t*(-6./(4.*t+1.) &
                  +2./(DiffC*(4.*t+1.))*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))) &
                  +4./(DiffC*(4.*t+1.)**2)*SUM((x(:)-AdvVel(:)*t-x0(:))**2)) &
         + Resu*( 24./(4.*t+1.)**2 &
                  -8./(DiffC*(4.*t+1.)**2)*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))) &
                  -2./(DiffC*(4.*t+1.))*SUM(AdvVel(:)*AdvVel(:))    &
                 -32./(DiffC*(4.*t+1.)**3)*SUM((x(:)-AdvVel(:)*t-x0(:))**2)   &
                  -8./(DiffC*(4.*t+1.)**2)*SUM(AdvVel(:)*(x(:)-AdvVel(:)*t-x0(:))))
#endif
CASE(6) ! One-dimensional cosine wave with amplitude 1 and angular frequency specified by user.
  ! Used for the dissipation and dispersion analysis following Gassner & Kopriva 2011.
  k = OmegaRef/AdvVel(1)
  Resu = COS(k*x(1)-OmegaRef*tEval)
  IF(fullBoundaryOrder)THEN
    Resu_t  = OmegaRef*SIN(k*x(1)-OmegaRef*tEval)
    Resu_tt = -1.*(OmegaRef**2)*COS(k*x(1)-OmegaRef*tEval)
  END IF
CASE DEFAULT
  CALL Abort(__STAMP__,&
             'Specified exactfuntion not implemented!')
END SELECT ! ExactFunction

! For O3 LS 3-stage RK, we have to define proper time dependent BC
IF(fullBoundaryOrder)THEN ! add resu_t, resu_tt if time dependant
  SELECT CASE(CurrentStage)
  CASE(1)
    ! resu = g(t)
  CASE(2)
    ! resu = g(t) + dt/3*g'(t)
    Resu=Resu + dt*RKc(2)*Resu_t
  CASE(3)
    ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
    Resu=Resu + RKc(3)*dt*Resu_t + RKc(2)*RKb(2)*dt*dt*Resu_tt
  CASE DEFAULT
    ! Stop, works only for 3 Stage O3 LS RK
    CALL Abort(__STAMP__,&
               'Time-dependant exactfuntion works only for 3 Stage O3 LS RK!')
  END SELECT
END IF
END SUBROUTINE ExactFunc



!==================================================================================================================================
!> Compute source terms for some specific testcases and adds it to DG time derivative
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,t)
! MODULES
USE MOD_Globals,       ONLY:Abort
USE MOD_Equation_Vars, ONLY:IniExactFunc,AdvVel,doCalcSource
USE MOD_Mesh_Vars,     ONLY:nElems,Elem_xGP,sJ
USE MOD_PreProc
#if PARABOLIC
USE MOD_Equation_Vars, ONLY:DiffC
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k
REAL                :: Omega
REAL                :: Ut_src(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
!==================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(22)
#if PARABOLIC
  Omega=2.*PP_PI
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut_src(:,i,j,k)=-2.*DiffC*Omega**2.*SIN(Omega*(AdvVel(1)*t-Elem_xGP(1,i,j,k,iElem) + &
                                                     AdvVel(2)*t-Elem_xGP(2,i,j,k,iElem)))
    END DO; END DO; END DO ! i,j,k
    ! Add contribution and take jacobian into account
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/(sJ(i,j,k,iElem,0))
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
#endif
CASE(31)
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut_src(:,i,j,k) = AdvVel(1)/2.*COS(PP_PI*(Elem_xGP(1,i,j,k,iElem)-AdvVel(1)/2.*t))*PP_PI
#if PARABOLIC
      Ut_src(:,i,j,k) = Ut_src(:,i,j,k) + PP_PI*PP_Pi*SIN(PP_PI*(Elem_xGP(1,i,j,k,iElem)-AdvVel(1)/2.*t))*diffC
#endif
    END DO; END DO; END DO ! i,j,k
    ! Add contribution and take jacobian into account
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/(sJ(i,j,k,iElem,0))
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
CASE DEFAULT
  doCalcSource=.FALSE.
END SELECT ! ExactFunction
END SUBROUTINE CalcSource

END MODULE MOD_Exactfunc
