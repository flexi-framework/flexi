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

!==================================================================================================================================
!> Soubroutines providing exactly evaluated functions used in initialization or boundary conditions.
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
!> Define parameters of exact functions
!==================================================================================================================================
SUBROUTINE DefineParametersExactFunc()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Exactfunc")
CALL prms%CreateIntFromStringOption('IniExactFunc', "Exact function to be used for computing initial solution.")
CALL addStrListEntry('IniExactFunc','testcase' ,-1)
CALL addStrListEntry('IniExactFunc','testcase' ,0)
CALL addStrListEntry('IniExactFunc','refstate' ,1)
CALL addStrListEntry('IniExactFunc','sinedens' ,2)
CALL addStrListEntry('IniExactFunc','sinedensx',21)
CALL addStrListEntry('IniExactFunc','lindens'  ,3)
CALL addStrListEntry('IniExactFunc','sinevel'  ,4)
CALL addStrListEntry('IniExactFunc','sinevelx' ,41)
CALL addStrListEntry('IniExactFunc','sinevely' ,42)
CALL addStrListEntry('IniExactFunc','sinevelz' ,43)
CALL addStrListEntry('IniExactFunc','roundjet' ,5)
CALL addStrListEntry('IniExactFunc','cylinder' ,6)
CALL addStrListEntry('IniExactFunc','shuvortex',7)
CALL addStrListEntry('IniExactFunc','couette'  ,8)
CALL addStrListEntry('IniExactFunc','cavity'   ,9)
CALL addStrListEntry('IniExactFunc','shock'    ,10)
CALL addStrListEntry('IniExactFunc','sod'      ,11)
CALL addStrListEntry('IniExactFunc','dmr'      ,13)
#if PARABOLIC
CALL addStrListEntry('IniExactFunc','blasius'  ,1338)
#endif
CALL prms%CreateRealArrayOption(    'AdvVel',       "Advection velocity (v1,v2,v3) required for exactfunction CASE(2,21,4,8)")
CALL prms%CreateRealOption(         'MachShock',    "Parameter required for CASE(10)", '1.5')
CALL prms%CreateRealOption(         'PreShockDens', "Parameter required for CASE(10)", '1.0')
CALL prms%CreateRealArrayOption(    'IniCenter',    "Shu Vortex CASE(7) (x,y,z)")
CALL prms%CreateRealArrayOption(    'IniAxis',      "Shu Vortex CASE(7) (x,y,z)")
CALL prms%CreateRealOption(         'IniAmplitude', "Shu Vortex CASE(7)", '0.2')
CALL prms%CreateRealOption(         'IniHalfwidth', "Shu Vortex CASE(7)", '0.2')
#if PARABOLIC
CALL prms%CreateRealOption(         'delta99_in',   "Blasius boundary layer CASE(1338)")
CALL prms%CreateRealArrayOption(    'x_in',         "Blasius boundary layer CASE(1338)")
#endif

END SUBROUTINE DefineParametersExactFunc

!==================================================================================================================================
!> Get some parameters needed for exact function
!==================================================================================================================================
SUBROUTINE InitExactFunc()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_ExactFunc_Vars
USE MOD_Equation_Vars      ,ONLY: IniExactFunc,IniRefState

! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT EXACT FUNCTION...'

IniExactFunc = GETINTFROMSTR('IniExactFunc')
IniRefState  = GETINT('IniRefState', "-1")
! Read in boundary parameters
SELECT CASE (IniExactFunc)
CASE(2,21,3,4,41,42,43) ! synthetic test cases
  AdvVel       = GETREALARRAY('AdvVel',3)
CASE(7) ! Shu Vortex
  IniCenter    = GETREALARRAY('IniCenter',3,'(/0.,0.,0./)')
  IniAxis      = GETREALARRAY('IniAxis',3,'(/0.,0.,1./)')
  IniAmplitude = GETREAL('IniAmplitude','0.2')
  IniHalfwidth = GETREAL('IniHalfwidth','0.2')
CASE(8) ! couette-poiseuille flow
  P_Parameter  = GETREAL('P_Parameter','0.0')
  U_Parameter  = GETREAL('U_Parameter','0.01')
CASE(10) ! shock
  MachShock    = GETREAL('MachShock','1.5')
  PreShockDens = GETREAL('PreShockDens','1.0')
#if PARABOLIC
CASE(1338) ! Blasius boundary layer solution
  delta99_in      = GETREAL('delta99_in')
  x_in            = GETREALARRAY('x_in',2,'(/0.,0./)')
  BlasiusInitDone = .TRUE. ! Mark Blasius init as done so we don't read the parameters again in BC init
#endif 
CASE DEFAULT
END SELECT ! IniExactFunc

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT FUNCTION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitExactFunc

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!> t is the actual time
!> dt is only needed to compute the time dependent boundary values for the RK scheme
!> for each function resu and the first and second time derivative resu_t and resu_tt have to be defined (is trivial for constants)
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu)
! MODULES
USE MOD_Preproc        ,ONLY: PP_PI
USE MOD_Globals        ,ONLY: Abort
USE MOD_Mathtools      ,ONLY: CROSS
USE MOD_Eos_Vars       ,ONLY: Kappa,sKappaM1,KappaM1,KappaP1,R
USE MOD_Exactfunc_Vars ,ONLY: IniCenter,IniHalfwidth,IniAmplitude,IniAxis,AdvVel
USE MOD_Exactfunc_Vars ,ONLY: MachShock,PreShockDens
USE MOD_Exactfunc_Vars ,ONLY: P_Parameter,U_Parameter
USE MOD_Equation_Vars  ,ONLY: IniRefState,RefStateCons,RefStatePrim
USE MOD_Timedisc_Vars  ,ONLY: fullBoundaryOrder,CurrentStage,dt,RKb,RKc,t
USE MOD_TestCase       ,ONLY: ExactFuncTestcase
USE MOD_EOS            ,ONLY: PrimToCons,ConsToPrim
#if PARABOLIC
USE MOD_Eos_Vars       ,ONLY: mu0
USE MOD_Exactfunc_Vars ,ONLY: delta99_in,x_in
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: ExactFunction          !< determines the exact function
REAL,INTENT(IN)                 :: x(3)                   !< physical coordinates
REAL,INTENT(IN)                 :: tIn                    !< solution time (Runge-Kutta stage)
REAL,INTENT(OUT)                :: Resu(5)                !< state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: tEval
REAL                            :: Resu_t(5),Resu_tt(5),ov ! state in conservative variables
REAL                            :: Frequency,Amplitude
REAL                            :: Omega
REAL                            :: Vel(3),Cent(3),a
REAL                            :: Prim(PP_nVarPrim)
REAL                            :: r_len
REAL                            :: Ms,xs
REAL                            :: Resul(5),Resur(5)
REAL                            :: random
REAL                            :: du, dTemp, RT, r2       ! aux var for SHU VORTEX,isentropic vortex case 12
REAL                            :: pi_loc,phi,radius       ! needed for cylinder potential flow
REAL                            :: h,sRT,pexit,pentry   ! needed for Couette-Poiseuille
#if PARABOLIC
! needed for blasius BL
INTEGER                         :: nSteps,i
REAL                            :: eta,deta,deta2,f,fp,fpp,fppp,fbar,fpbar,fppbar,fpppbar
REAL                            :: x_eff(3),x_offset(3)
#endif
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration

Resu   =0.
Resu_t =0.
Resu_tt=0.

! Determine the value, the first and the second time derivative
SELECT CASE (ExactFunction)
CASE DEFAULT
  CALL ExactFuncTestcase(tEval,x,Resu,Resu_t,Resu_tt)
CASE(0)
  CALL ExactFuncTestcase(tEval,x,Resu,Resu_t,Resu_tt)
CASE(1) ! constant
  Resu = RefStateCons(:,IniRefState)
CASE(2) ! sinus
  Frequency=0.5
  Amplitude=0.3
  Omega=2.*PP_Pi*Frequency
  ! base flow
  prim(1)   = 1.
  prim(2:4) = AdvVel
  prim(5)   = 1.
  Vel=prim(2:4)
  cent=x-Vel*tEval
  prim(1)=prim(1)*(1.+Amplitude*SIN(Omega*SUM(cent(1:3))))
  ! g(t)
  Resu(1)=prim(1) ! rho
  Resu(2:4)=prim(1)*prim(2:4) ! rho*vel
  Resu(5)=prim(5)*sKappaM1+0.5*prim(1)*SUM(prim(2:4)*prim(2:4)) ! rho*e

  IF(fullBoundaryOrder)THEN
    ov=Omega*SUM(vel)
    ! g'(t)
    Resu_t(1)=-Amplitude*cos(Omega*SUM(cent(1:3)))*ov
    Resu_t(2:4)=Resu_t(1)*prim(2:4) ! rho*vel
    Resu_t(5)=0.5*Resu_t(1)*SUM(prim(2:4)*prim(2:4))
    ! g''(t)
    Resu_tt(1)=-Amplitude*sin(Omega*SUM(cent(1:3)))*ov**2.
    Resu_tt(2:4)=Resu_tt(1)*prim(2:4)
    Resu_tt(5)=0.5*Resu_tt(1)*SUM(prim(2:4)*prim(2:4))
  END IF
CASE(21) ! sinus x
  Frequency=0.5
  Amplitude=0.3
  Omega=2.*PP_Pi*Frequency
  ! base flow
  prim(1)   = 1.
  prim(2:4) = AdvVel
  prim(5)   = 1.
  Vel=prim(2:4)
  cent=x-Vel*tEval
  prim(1)=prim(1)*(1.+Amplitude*SIN(Omega*cent(1)))
  ! g(t)
  Resu(1)=prim(1) ! rho
  Resu(2:4)=prim(1)*prim(2:4) ! rho*vel
  Resu(5)=prim(5)*sKappaM1+0.5*prim(1)*SUM(prim(2:4)*prim(2:4)) ! rho*e

  IF(fullBoundaryOrder)THEN
    ov=Omega*SUM(vel)
    ! g'(t)
    Resu_t(1)=-Amplitude*cos(Omega*cent(1))*ov
    Resu_t(2:4)=Resu_t(1)*prim(2:4) ! rho*vel
    Resu_t(5)=0.5*Resu_t(1)*SUM(prim(2:4)*prim(2:4))
    ! g''(t)
    Resu_tt(1)=-Amplitude*sin(Omega*cent(1))*ov**2.
    Resu_tt(2:4)=Resu_tt(1)*prim(2:4)
    Resu_tt(5)=0.5*Resu_tt(1)*SUM(prim(2:4)*prim(2:4))
  END IF
CASE(3) ! linear in rho
  ! base flow
  prim(1)   = 100.
  prim(2:4) = AdvVel
  prim(5)   = 1.
  Vel=prim(2:4)
  cent=x-Vel*tEval
  prim(1)=prim(1)+SUM(AdvVel*cent)
  ! g(t)
  Resu(1)=prim(1) ! rho
  Resu(2:4)=prim(1)*prim(2:4) ! rho*vel
  Resu(5)=prim(5)*sKappaM1+0.5*prim(1)*SUM(prim(2:4)*prim(2:4)) ! rho*e
  IF(fullBoundaryOrder)THEN
    ! g'(t)
    Resu_t(1)=-SUM(Vel)
    Resu_t(2:4)=Resu_t(1)*prim(2:4) ! rho*vel
    Resu_t(5)=0.5*Resu_t(1)*SUM(prim(2:4)*prim(2:4))
  END IF
CASE(4) ! exact function
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi

  ! g(t)
  Resu(1:4)=2.+ Amplitude*sin(Omega*SUM(x) - a*tEval)
  Resu(5)=Resu(1)*Resu(1)
  IF(fullBoundaryOrder)THEN
    ! g'(t)
    Resu_t(1:4)=(-a)*Amplitude*cos(Omega*SUM(x) - a*tEval)
    Resu_t(5)=2.*Resu(1)*Resu_t(1)
    ! g''(t)
    Resu_tt(1:4)=-a*a*Amplitude*sin(Omega*SUM(x) - a*tEval)
    Resu_tt(5)=2.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))
  END IF
CASE(41) ! SINUS in x
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(1:2)=2.+ Amplitude*sin(Omega*x(1) - a*tEval)
  Resu(5)=Resu(1)*Resu(1)
  IF(fullBoundaryOrder)THEN
    Resu_t = 0.
    Resu_tt = 0.
    ! g'(t)
    Resu_t(1:2)=(-a)*Amplitude*cos(Omega*x(1) - a*tEval)
    Resu_t(5)=2.*Resu(1)*Resu_t(1)
    ! g''(t)
    Resu_tt(1:2)=-a*a*Amplitude*sin(Omega*x(1) - a*tEval)
    Resu_tt(5)=2.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))
  END IF
CASE(42) ! SINUS in y
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(2)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(1)=2.+ Amplitude*sin(Omega*x(2) - a*tEval)
  Resu(3)=Resu(1)
  Resu(5)=Resu(1)*Resu(1)
  IF(fullBoundaryOrder)THEN
    Resu_tt = 0.
    ! g'(t)
    Resu_t(1)=(-a)*Amplitude*cos(Omega*x(2) - a*tEval)
    Resu_t(3) = Resu_t(1)
    Resu_t(5)=2.*Resu(1)*Resu_t(1)
    ! g''(t)
    Resu_tt(1)=-a*a*Amplitude*sin(Omega*x(2) - a*tEval)
    Resu_t(3) = Resu_t(1)
    Resu_tt(5)=2.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))
  END IF
CASE(43) ! SINUS in z
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(3)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(1)=2.+ Amplitude*sin(Omega*x(3) - a*tEval)
  Resu(4)=Resu(1)
  Resu(5)=Resu(1)*Resu(1)
  IF(fullBoundaryOrder)THEN
    Resu_tt = 0.
    ! g'(t)
    Resu_t(1)=(-a)*Amplitude*cos(Omega*x(3) - a*tEval)
    Resu_t(4) = Resu_t(1)
    Resu_t(5)=2.*Resu(1)*Resu_t(1)
    ! g''(t)
    Resu_tt(1)=-a*a*Amplitude*sin(Omega*x(3) - a*tEval)
    Resu_t(4) = Resu_t(1)
    Resu_tt(5)=2.*(Resu_t(1)*Resu_t(1) + Resu(1)*Resu_tt(1))
  END IF
CASE(5) !Roundjet Bogey Bailly 2002, Re=65000, x-axis is jet axis
  prim(1)  =1.
  prim(2:4)=0.
  prim(5)  =1./Kappa
  prim(6) = prim(5)/(prim(1)*R)
  ! Jet inflow (from x=0, diameter 2.0)
  ! Initial jet radius: rj=1.
  ! Momentum thickness: delta_theta0=0.05=1/20
  ! Re=65000
  ! Uco=0.
  ! Uj=0.9
  r_len=SQRT((x(2)*x(2)+x(3)*x(3)))
  prim(2)=0.9*0.5*(1.+TANH((1.-r_len)*10.))
  CALL RANDOM_NUMBER(random)
  ! Random disturbance +-5%
  random=0.05*2.*(random-0.5)
  prim(2)=prim(2)+random*prim(2)
  prim(3)=x(2)/r_len*0.5*random*prim(2)
  prim(4)=x(3)/r_len*0.5*random*prim(2)
  CALL PrimToCons(prim,ResuL)
  prim(1)  =1.
  prim(2:4)=0.
  prim(5)  =1./Kappa
  prim(6) = prim(5)/(prim(1)*R)
  CALL PrimToCons(prim,ResuR)
!   after x=10 blend to ResuR
  Resu=ResuL+(ResuR-ResuL)*0.5*(1.+tanh(x(1)-10.))
CASE(6)  ! Cylinder flow
  IF(tEval .EQ. 0.)THEN   ! Initialize potential flow
    prim(1)=RefStatePrim(1,IniRefState)  ! Density
    prim(4)=0.           ! VelocityZ=0. (2D flow)
    ! Calculate cylinder coordinates (0<phi<Pi/2)
    pi_loc=ASIN(1.)*2.
    IF(x(1) .LT. 0.)THEN
      phi=ATAN(ABS(x(2))/ABS(x(1)))
      IF(x(2) .LT. 0.)THEN
        phi=pi_loc+phi
      ELSE
        phi=pi_loc-phi
      END IF
    ELSEIF(x(1) .GT. 0.)THEN
      phi=ATAN(ABS(x(2))/ABS(x(1)))
      IF(x(2) .LT. 0.) phi=2.*pi_loc-phi
    ELSE
      IF(x(2) .LT. 0.)THEN
        phi=pi_loc*1.5
      ELSE
        phi=pi_loc*0.5
      END IF
    END IF
    ! Calculate radius**2
    radius=x(1)*x(1)+x(2)*x(2)
    ! Calculate velocities, radius of cylinder=0.5
    prim(2)=RefStatePrim(2,IniRefState)*(COS(phi)**2*(1.-0.25/radius)+SIN(phi)**2*(1.+0.25/radius))
    prim(3)=RefStatePrim(2,IniRefState)*(-2.)*SIN(phi)*COS(phi)*0.25/radius
    ! Calculate pressure, RefState(2)=u_infinity
    prim(5)=RefStatePrim(5,IniRefState) + &
            0.5*prim(1)*(RefStatePrim(2,IniRefState)*RefStatePrim(2,IniRefState)-prim(2)*prim(2)-prim(3)*prim(3))
    prim(6) = prim(5)/(prim(1)*R)
  ELSE  ! Use RefState as BC
    prim=RefStatePrim(:,IniRefState)
  END IF  ! t=0
  CALL PrimToCons(prim,resu)
CASE(7) ! SHU VORTEX,isentropic vortex
  ! base flow
  prim=RefStatePrim(:,IniRefState)  ! Density
  ! ini-Parameter of the Example
  vel=prim(2:4)
  RT=prim(PP_nVar)/prim(1) !ideal gas
  cent=(iniCenter+vel*tEval)!centerpoint time dependant
  cent=x-cent ! distance to centerpoint
  cent=CROSS(iniAxis,cent) !distance to axis, tangent vector, length r
  cent=cent/iniHalfWidth !Halfwidth is dimension 1
  r2=SUM(cent*cent) !
  du = iniAmplitude/(2.*PP_Pi)*exp(0.5*(1.-r2)) ! vel. perturbation
  dTemp = -kappaM1/(2.*kappa*RT)*du**2 ! adiabatic
  prim(1)=prim(1)*(1.+dTemp)**(1.*skappaM1) !rho
  prim(2:4)=prim(2:4)+du*cent(:) !v
  prim(PP_nVar)=prim(PP_nVar)*(1.+dTemp)**(kappa/kappaM1) !p
  prim(6) = prim(5)/(prim(1)*R)
  CALL PrimToCons(prim,resu)
CASE(8) !Couette-Poiseuille flow between plates: exact steady lamiar solution with height=1 !
        !(Gao, hesthaven, Warburton)
  RT=1. ! Hesthaven: Absorbing layers for weakly compressible flows
  sRT=1/RT
  ! size of domain must be [-0.5, 0.5]^2 -> (x*y)
  h=0.5
  prim(2)= U_parameter*(0.5*(1.+x(2)/h) + P_parameter*(1-(x(2)/h)**2))
  prim(3:4)= 0.
  pexit=0.9996
  pentry=1.0004
  prim(5)= ( ( x(1) - (-0.5) )*( pexit - pentry) / ( 0.5 - (-0.5)) ) + pentry
  prim(1)=prim(5)*sRT
  prim(6)=prim(5)/(prim(1)*R)
  CALL PrimToCons(prim,Resu)
CASE(9) !lid driven cavity flow from Gao, Hesthaven, Warburton
        !"Absorbing layers for weakly compressible flows", to appear, JSC, 2016
        ! Special "regularized" driven cavity BC to prevent singularities at corners
        ! top BC assumed to be in x-direction from 0..1
  Prim = RefStatePrim(:,IniRefState)
  IF (x(1).LT.0.2) THEN 
    prim(2)=1000*4.9333*x(1)**4-1.4267*1000*x(1)**3+0.1297*1000*x(1)**2-0.0033*1000*x(1)
  ELSEIF (x(1).LE.0.8) THEN
    prim(2)=1.0
  ELSE  
    prim(2)=1000*4.9333*x(1)**4-1.8307*10000*x(1)**3+2.5450*10000*x(1)**2-1.5709*10000*x(1)+10000*0.3633
  ENDIF
  CALL PrimToCons(prim,Resu)
CASE(10) ! shock
  prim=0.

  ! pre-shock
  prim(1) = PreShockDens
  Ms      = MachShock

  prim(5)=prim(1)/Kappa
  prim(6)=prim(5)/(prim(1)*R)
  CALL PrimToCons(prim,Resur)

  ! post-shock
  prim(3)=prim(1) ! temporal storage of pre-shock density
  prim(1)=prim(1)*((KappaP1)*Ms*Ms)/(KappaM1*Ms*Ms+2.)
  prim(5)=prim(5)*(2.*Kappa*Ms*Ms-KappaM1)/(KappaP1)
  prim(6)=prim(5)/(prim(1)*R)
  IF (prim(2) .EQ. 0.0) THEN
    prim(2)=Ms*(1.-prim(3)/prim(1))
  ELSE
    prim(2)=prim(2)*prim(3)/prim(1)
  END IF
  prim(3)=0. ! reset temporal storage
  CALL PrimToCons(prim,Resul)
  xs=5.+Ms*tEval ! 5. bei 10x10x10 Rechengebiet
  ! Tanh boundary
  Resu=-0.5*(Resul-Resur)*TANH(5.0*(x(1)-xs))+Resur+0.5*(Resul-Resur)
CASE(11) ! Sod Shock tube
  xs = 0.5
  IF (X(1).LE.xs) THEN
    Resu = RefStateCons(:,1)
  ELSE
    Resu = RefStateCons(:,2)
  END IF
CASE(12) ! Shu Osher density fluctuations shock wave interaction 
  IF (x(1).LT.-4.0) THEN
    prim(1) = 3.857143
    prim(2) = 2.629369
    prim(3:4) = 0.
    prim(5) = 10.33333
  ELSE
    prim(1) = 1.+0.2*SIN(5.*x(1))
    prim(2:4) = 0.
    prim(5) = 1.
  END IF
  CALL PrimToCons(prim,resu)

CASE(13) ! DoubleMachReflection (see e.g. http://www.astro.princeton.edu/~jstone/Athena/tests/dmr/dmr.html )
  IF (x(1).EQ.0.) THEN
    prim = RefStatePrim(:,1)
  ELSE IF (x(1).EQ.4.0) THEN
    prim = RefStatePrim(:,2)
  ELSE
    IF (x(1).LT.1./6.+(x(2)+20.*t)*1./3.**0.5) THEN
      prim = RefStatePrim(:,1)
    ELSE
      prim = RefStatePrim(:,2)
    END IF
  END IF
  CALL PrimToCons(prim,resu)
#if PARABOLIC
CASE(1338) ! blasius
  prim=RefStatePrim(:,IniRefState)
  ! calculate equivalent x for Blasius flat plate to have delta99_in at x_in
  x_offset(1)=(delta99_in/5)**2*prim(1)*prim(2)/mu0-x_in(1)
  x_offset(2)=-x_in(2)
  x_offset(3)=0.
  x_eff=x+x_offset
  IF(x_eff(2).GE.0 .AND. x_eff(1).GT.0) THEN
    ! scale bl position in physical space to reference space, eta=5 is ~99% bl thickness
    eta=x_eff(2)*(prim(1)*prim(2)/(mu0*x_eff(1)))**0.5

    deta=0.02 ! step size
    nSteps=CEILING(eta/deta)
    deta =eta/nSteps
    deta2=0.5*deta

    f=0.
    fp=0.
    fpp=0.332 ! default literature value, don't change if you don't know what you're doing
    fppp=0.
    !Blasius boundary layer
    DO i=1,nSteps
      ! predictor 
      fbar    = f   + deta * fp
      fpbar   = fp  + deta * fpp
      fppbar  = fpp + deta * fppp
      fpppbar = -0.5*fbar*fppbar
      ! corrector 
      f       = f   + deta2 * (fp   + fpbar)
      fp      = fp  + deta2 * (fpp  + fppbar)
      fpp     = fpp + deta2 * (fppp + fpppbar)
      fppp    = -0.5*f*fpp
    END DO
    prim(3)=0.5*(mu0*prim(2)/prim(1)/x_eff(1))**0.5*(fp*eta-f)
    prim(2)=RefStatePrim(2,IniRefState)*fp
  ELSE
    IF(x_eff(2).LT.0) THEN
      prim(2)=0.
    END IF 
  END IF
  CALL PrimToCons(prim,resu)
#endif
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
    CALL abort(__STAMP__,&
               'Exactfuntion works only for 3 Stage O3 LS RK!')
  END SELECT
END IF
END SUBROUTINE ExactFunc

!==================================================================================================================================
!> Compute source terms for some specific testcases and adds it to DG time derivative
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:IniExactFunc,doCalcSource
USE MOD_Eos_Vars,ONLY:Kappa,KappaM1
USE MOD_Exactfunc_Vars,ONLY:AdvVel
#if PARABOLIC
USE MOD_Eos_Vars,ONLY:mu0,Pr
#endif
USE MOD_Mesh_Vars,    ONLY:Elem_xGP,sJ,nElems
#if FV_ENABLED
USE MOD_ChangeBasis,  ONLY:ChangeBasis3D
USE MOD_FV_Vars,      ONLY:FV_Vdm,FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< current solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
REAL                :: Ut_src(5,0:PP_N,0:PP_N,0:PP_N)
REAL                :: Frequency,Amplitude,Omega,a
REAL                :: sinXGP,sinXGP2,cosXGP,at
REAL                :: tmp(6)
REAL                :: C
#if FV_ENABLED
REAL                :: Ut_src2(5,0:PP_N,0:PP_N,0:PP_N)
#endif
!==================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(4) ! exact function
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi
  tmp(1)=-a+3*Omega
  tmp(2)=-a+0.5*Omega*(1.+kappa*5.)
  tmp(3)=Amplitude*Omega*KappaM1
  tmp(4)=0.5*((9.+Kappa*15.)*Omega-8.*a)
  tmp(5)=Amplitude*(3.*Omega*Kappa-a)
#if PARABOLIC
  tmp(6)=3.*mu0*Kappa*Omega*Omega/Pr
#else
  tmp(6)=0.
#endif
  tmp=tmp*Amplitude
  at=a*t
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      cosXGP=COS(omega*SUM(Elem_xGP(:,i,j,k,iElem))-at)
      sinXGP=SIN(omega*SUM(Elem_xGP(:,i,j,k,iElem))-at)
      sinXGP2=2.*sinXGP*cosXGP !=SIN(2.*(omega*SUM(Elem_xGP(:,i,j,k,iElem))-a*t))
      Ut_src(1  ,i,j,k) = tmp(1)*cosXGP
      Ut_src(2:4,i,j,k) = tmp(2)*cosXGP + tmp(3)*sinXGP2
      Ut_src(5  ,i,j,k) = tmp(4)*cosXGP + tmp(5)*sinXGP2 + tmp(6)*sinXGP
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED    
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem     
      CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
      END DO; END DO; END DO ! i,j,k
    ELSE
#endif      
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
#if FV_ENABLED    
    END IF
#endif
  END DO ! iElem
CASE(41) ! Sinus in x
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi
  C = 2.0

  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
#if PARABOLIC      
      Ut_src(1,i,j,k) = (-Amplitude*a+Amplitude*omega)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)

      Ut_src(2,i,j,k) = (-Amplitude**2*omega+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(1,i,j,k,iElem)-2.*a*t)+&
                      (-Amplitude*a+2.*Amplitude*omega*kappa*C-1./2.*Amplitude*omega*kappa+&
                      3./2.*Amplitude*omega-2.*Amplitude*omega*C)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)

      Ut_src(3:4,i,j,k) = 0.0

      Ut_src(5,i,j,k) = mu0*kappa*Amplitude*sin(omega*Elem_xGP(1,i,j,k,iElem)-a*t)*omega**2/Pr+&
                      (-Amplitude**2*a+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(1,i,j,k,iElem)-2.*a*t)+&
                      1./2.*(-4.*Amplitude*a*C+4.*Amplitude*omega*kappa*C-Amplitude*omega*kappa+&
                           Amplitude*omega)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)
#else
      Ut_src(1,i,j,k) = (-amplitude*a+amplitude*omega)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t) 
      Ut_src(2,i,j,k) = (-amplitude**2*omega+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(1,i,j,k,iElem)-2.*a*t)+ &
                      (-amplitude*a+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+ &
                      3./2.*amplitude*omega-2.*amplitude*omega*C)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)
                      
      Ut_src(3:4,i,j,k) = 0.0
      Ut_src(5,i,j,k) = (-amplitude**2*a+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(1,i,j,k,iElem)-2.*a*t)+&
                      (-2.*amplitude*a*C+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+&
                      1./2.*amplitude*omega)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)
                      
#endif
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED    
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
      CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
      END DO; END DO; END DO ! i,j,k
    ELSE
#endif
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
#if FV_ENABLED    
    END IF
#endif
  END DO
CASE(42) ! Sinus in y
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(2)*2.*PP_Pi
  C = 2.0

  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
#if PARABOLIC      
      Ut_src(1,i,j,k) = (-Amplitude*a+Amplitude*omega)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)
      Ut_src(2,i,j,k) = 0.0
      Ut_src(3,i,j,k) = (-Amplitude**2*omega+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(2,i,j,k,iElem)-2.*a*t)+&
                      (-Amplitude*a+2.*Amplitude*omega*kappa*C-1./2.*Amplitude*omega*kappa+&
                      3./2.*Amplitude*omega-2.*Amplitude*omega*C)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)

      Ut_src(4,i,j,k) = 0.0

      Ut_src(5,i,j,k) = mu0*kappa*Amplitude*sin(omega*Elem_xGP(2,i,j,k,iElem)-a*t)*omega**2/Pr+&
                      (-Amplitude**2*a+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(2,i,j,k,iElem)-2.*a*t)+&
                      1./2.*(-4.*Amplitude*a*C+4.*Amplitude*omega*kappa*C-Amplitude*omega*kappa+&
                           Amplitude*omega)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)
#else
      Ut_src(1,i,j,k) = (-amplitude*a+amplitude*omega)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t) 
      Ut_src(2,i,j,k) = 0.0
      Ut_src(3,i,j,k) = (-amplitude**2*omega+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(2,i,j,k,iElem)-2.*a*t)+ &
                      (-amplitude*a+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+ &
                      3./2.*amplitude*omega-2.*amplitude*omega*C)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)
                      
      Ut_src(4,i,j,k) = 0.0
      Ut_src(5,i,j,k) = (-amplitude**2*a+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(2,i,j,k,iElem)-2.*a*t)+&
                      (-2.*amplitude*a*C+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+&
                      1./2.*amplitude*omega)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)
                      
#endif
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED    
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
      CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
      END DO; END DO; END DO ! i,j,k
    ELSE
#endif
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
#if FV_ENABLED    
    END IF
#endif
  END DO

CASE(43) ! Sinus in z
  Frequency=1.
  Amplitude=0.1
  Omega=PP_Pi*Frequency
  a=AdvVel(3)*2.*PP_Pi
  C = 2.0

  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
#if PARABOLIC      
      Ut_src(1,i,j,k) = (-Amplitude*a+Amplitude*omega)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)
      Ut_src(2:3,i,j,k) = 0.0
      Ut_src(4,i,j,k) = (-Amplitude**2*omega+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(3,i,j,k,iElem)-2.*a*t)+&
                      (-Amplitude*a+2.*Amplitude*omega*kappa*C-1./2.*Amplitude*omega*kappa+&
                      3./2.*Amplitude*omega-2.*Amplitude*omega*C)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)


      Ut_src(5,i,j,k) = mu0*kappa*Amplitude*sin(omega*Elem_xGP(3,i,j,k,iElem)-a*t)*omega**2/Pr+&
                      (-Amplitude**2*a+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(3,i,j,k,iElem)-2.*a*t)+&
                      1./2.*(-4.*Amplitude*a*C+4.*Amplitude*omega*kappa*C-Amplitude*omega*kappa+&
                           Amplitude*omega)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)
#else
      Ut_src(1,i,j,k) = (-amplitude*a+amplitude*omega)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t) 
      Ut_src(2:3,i,j,k) = 0.0
      Ut_src(4,i,j,k) = (-amplitude**2*omega+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(3,i,j,k,iElem)-2.*a*t)+ &
                      (-amplitude*a+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+ &
                      3./2.*amplitude*omega-2.*amplitude*omega*C)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)
                      
      Ut_src(5,i,j,k) = (-amplitude**2*a+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(3,i,j,k,iElem)-2.*a*t)+&
                      (-2.*amplitude*a*C+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+&
                      1./2.*amplitude*omega)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)
                      
#endif
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED    
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
      CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
      END DO; END DO; END DO ! i,j,k
    ELSE
#endif
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
#if FV_ENABLED    
    END IF
#endif
  END DO
CASE DEFAULT
  ! No source -> do nothing and set marker to not run again
  doCalcSource=.FALSE.
END SELECT ! ExactFunction
END SUBROUTINE CalcSource

END MODULE MOD_Exactfunc
