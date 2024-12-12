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
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Subroutines providing exactly evaluated functions used in initialization or boundary conditions.
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
CALL addStrListEntry('IniExactFunc','testcase'          ,-1)
CALL addStrListEntry('IniExactFunc','testcase'          ,0)
CALL addStrListEntry('IniExactFunc','refstate'          ,1)
CALL addStrListEntry('IniExactFunc','sinedens'          ,2)
CALL addStrListEntry('IniExactFunc','sinedensx'         ,21)
CALL addStrListEntry('IniExactFunc','lindens'           ,3)
CALL addStrListEntry('IniExactFunc','sinevel'           ,4)
CALL addStrListEntry('IniExactFunc','sinevelx'          ,41)
CALL addStrListEntry('IniExactFunc','sinevely'          ,42)
CALL addStrListEntry('IniExactFunc','sinevelz'          ,43)
CALL addStrListEntry('IniExactFunc','roundjet'          ,5)
CALL addStrListEntry('IniExactFunc','parabjet'          ,51)
CALL addStrListEntry('IniExactFunc','cylinder'          ,6)
CALL addStrListEntry('IniExactFunc','shuvortex'         ,7)
CALL addStrListEntry('IniExactFunc','couette'           ,8)
CALL addStrListEntry('IniExactFunc','cavity'            ,9)
CALL addStrListEntry('IniExactFunc','shock'             ,10)
CALL addStrListEntry('IniExactFunc','sod'               ,11)
CALL addStrListEntry('IniExactFunc','dmr'               ,13)
CALL addStrListEntry('IniExactFunc','harmonicgausspulse',14)
CALL addStrListEntry('IniExactFunc','blast_shock'       ,111)
#if PARABOLIC
CALL addStrListEntry('IniExactFunc','blasius'  ,1338)
#endif
CALL addStrListEntry('IniExactFunc','sedov'    ,1342)
CALL addStrListEntry('IniExactFunc','leblanc'  ,1344)
CALL addStrListEntry('IniExactFunc','hui'      ,1346)

CALL prms%CreateRealArrayOption(    'AdvVel',       "Advection velocity (v1,v2,v3) required for exactfunction CASE(2,21,4,8)")
CALL prms%CreateRealOption(         'IniAmplitude', "Amplitude for synthetic test case")
CALL prms%CreateRealOption(         'IniFrequency', "Frequency for synthetic test case")
CALL prms%CreateRealOption(         'MachShock',    "Parameter required for CASE(10)", '1.5')
CALL prms%CreateRealOption(         'PreShockDens', "Parameter required for CASE(10)", '1.0')
CALL prms%CreateRealOption(         'ShockPos',     "Parameter required for CASE(11)", '0.5')
CALL prms%CreateRealArrayOption(    'IniCenter',    "Shu Vortex CASE(7) (x,y,z)")
CALL prms%CreateRealArrayOption(    'IniAxis',      "Shu Vortex CASE(7) (x,y,z)")
CALL prms%CreateRealOption(         'IniHalfwidth', "Shu Vortex CASE(7)", '0.2')
CALL prms%CreateRealOption(         'JetRadius',    "Roundjet CASE(5,51,33)", '1.0')
CALL prms%CreateRealOption(         'JetEnd',       "Roundjet CASE(5,51,33)", '10.0')
CALL prms%CreateRealOption(         'JetAmplitude', "Roundjet CASE(5,51,33)", '1.0')
CALL prms%CreateRealOption(         'Ramping',      "Subsonic mass inflow CASE(28)"  , '1.0')
CALL prms%CreateRealOption(         'P_Parameter',  "Couette-Poiseuille flow CASE(8)", '0.0')
CALL prms%CreateRealOption(         'U_Parameter',  "Couette-Poiseuille flow CASE(8)", '0.01')
CALL prms%CreateRealOption(         'AmplitudeFactor',         "Harmonic Gauss Pulse CASE(14)", '0.1')
CALL prms%CreateRealOption(         'HarmonicFrequency',       "Harmonic Gauss Pulse CASE(14)", '400')
CALL prms%CreateRealOption(         'SigmaSqr',                "Harmonic Gauss Pulse CASE(14)", '0.1')
#if PARABOLIC
CALL prms%CreateRealOption(         'delta99_in',              "Blasius boundary layer CASE(1338)")
CALL prms%CreateRealArrayOption(    'x_in',                    "Blasius boundary layer CASE(1338)")
#endif
CALL prms%CreateRealOption(         'eradius_in',   "Radius of high energy area CASE(1342)")
CALL prms%CreateRealOption(         'leblanc_x0',   "Location of initial discontinity CASE(1344)")
CALL prms%CreateRealArrayOption(    'Hui_Data',    "Orientation and axis intercept data for Hui problem CASE(1346)")


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
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT EXACT FUNCTION...'

IniExactFunc = GETINTFROMSTR('IniExactFunc')
IniRefState  = GETINT('IniRefState', "-1")
! Read in boundary parameters
SELECT CASE (IniExactFunc)
  CASE(2,21) ! sinus
    AdvVel          = GETREALARRAY('AdvVel',3)
    IniFrequency    = GETREAL('IniFrequency','0.5')
    IniAmplitude    = GETREAL('IniAmplitude','0.3')
  CASE(3) ! synthetic test cases
    AdvVel          = GETREALARRAY('AdvVel',3)
  CASE(4,41,42,43) ! synthetic test cases
    AdvVel          = GETREALARRAY('AdvVel',3)
    IniFrequency    = GETREAL('IniFrequency','1.0')
    IniAmplitude    = GETREAL('IniAmplitude','0.1')
  CASE(5,51)
    JetRadius       = GETREAL('JetRadius')
    JetEnd          = GETREAL('JetEnd')
    JetAmplitude    = GETREAL('JetAmplitude')
  CASE(7) ! Shu Vortex
    IniCenter       = GETREALARRAY('IniCenter',3,'(/0.,0.,0./)')
    IniAxis         = GETREALARRAY('IniAxis',3,'(/0.,0.,1./)')
    IniAmplitude    = GETREAL('IniAmplitude','0.2')
    IniHalfwidth    = GETREAL('IniHalfwidth','0.2')
  CASE(8) ! couette-poiseuille flow
    P_Parameter     = GETREAL('P_Parameter')
    U_Parameter     = GETREAL('U_Parameter')
  CASE(10) ! shock
    MachShock       = GETREAL('MachShock')
    PreShockDens    = GETREAL('PreShockDens')
  CASE(11,111) ! 1D shock tube problem
    ShockPos        = GETREAL('ShockPos')
  CASE(14)
    HarmonicFrequency = GETREAL('HarmonicFrequency')
    AmplitudeFactor   = GETREAL('AmplitudeFactor')
    SiqmaSqr          = GETREAL('SigmaSqr')
#if PARABOLIC
  CASE(1338) ! Blasius boundary layer solution
    delta99_in      = GETREAL('delta99_in')
    x_in            = GETREALARRAY('x_in',2,'(/0.,0./)')
    BlasiusInitDone = .TRUE. ! Mark Blasius init as done so we don't read the parameters again in BC init
#endif
  CASE(1342) ! Energy radius in for Sedov
    SWRITE(UNIT_stdOut,'(A)')' READING ERADIUS!'
    eradius_in      = GETREAL('eradius_in')
    SWRITE(UNIT_stdOut,'(A)')' INIT ERADIUS DONE!'
  CASE(1344) ! Discontinuity location for LeBlanc
    SWRITE(UNIT_stdOut,'(A)')' READING LEBLANC_X0!'
    leblanc_x0      = GETREAL('leblanc_x0')
    SWRITE(UNIT_stdOut,'(A)')' INIT LEBLANC_X0 DONE!'
  CASE(1346) ! Axis angle (theta in rads) and Interface location (l) for Hui problem
    SWRITE(UNIT_stdOut,'(A)')' READING HUI_INTERCEPT!'
    Hui_Data    = GETREALARRAY('Hui_Data',2,'(/0.,1./)')
    SWRITE(UNIT_stdOut,'(A)')' INIT HUI_INTERCEPT DONE!'
  CASE DEFAULT
    ! Everything defined, do nothing
END SELECT ! IniExactFunc

#if PP_dim==2
SELECT CASE (IniExactFunc)
CASE(43) ! synthetic test cases
  CALL CollectiveStop(__STAMP__,'The selected exact function is not available in 2D!')
CASE(2,3,4,41,42) ! synthetic test cases
  IF(AdvVel(3).NE.0.) THEN
    CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set AdvVel(3) = 0!')
  END IF
END SELECT
#endif

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT FUNCTION DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitExactFunc

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!> t is the actual time
!> dt is only needed to compute the time dependent boundary values for the RK scheme
!> for each function resu and the first and second time derivative resu_t and resu_tt have to be defined (is trivial for constants)
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu,RefStateOpt)
! MODULES
USE MOD_Preproc        ,ONLY: PP_PI
USE MOD_Globals        ,ONLY: Abort
USE MOD_Mathtools      ,ONLY: CROSS
USE MOD_Eos_Vars       ,ONLY: Kappa,sKappaM1,KappaM1,KappaP1,R
USE MOD_Exactfunc_Vars ,ONLY: IniCenter,IniHalfwidth,IniAmplitude,IniFrequency,IniAxis,AdvVel
USE MOD_Exactfunc_Vars ,ONLY: MachShock,PreShockDens,ShockPos
USE MOD_Exactfunc_Vars ,ONLY: P_Parameter,U_Parameter
USE MOD_Exactfunc_Vars ,ONLY: JetRadius,JetEnd,JetAmplitude
USE MOD_Equation_Vars  ,ONLY: IniRefState,RefStateCons,RefStatePrim
USE MOD_Timedisc_Vars  ,ONLY: fullBoundaryOrder,CurrentStage,dt,RKb,RKc,t
USE MOD_TestCase       ,ONLY: ExactFuncTestcase
USE MOD_EOS            ,ONLY: PrimToCons,ConsToPrim
#if PARABOLIC
USE MOD_Eos_Vars       ,ONLY: mu0
USE MOD_Exactfunc_Vars ,ONLY: delta99_in,x_in
#endif
USE MOD_Exactfunc_Vars ,ONLY: eradius_in
USE MOD_Exactfunc_Vars ,ONLY: leblanc_x0
USE MOD_Exactfunc_Vars ,ONLY: Hui_Data

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: ExactFunction          !< determines the exact function
REAL,INTENT(IN)                 :: x(3)                   !< physical coordinates
REAL,INTENT(IN)                 :: tIn                    !< solution time (Runge-Kutta stage)
REAL,INTENT(OUT)                :: Resu(PP_nVar)          !< state in conservative variables
INTEGER,INTENT(IN),OPTIONAL     :: RefStateOpt            !< refstate to be used for exact func
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: RefState
REAL                            :: tEval
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar),ov ! state in conservative variables
REAL                            :: Frequency,Amplitude
REAL                            :: Omega
REAL                            :: Vel(3),Cent(3),a
REAL                            :: Prim(PP_nVarPrim)
REAL                            :: r_len
REAL                            :: Ms,xs
REAL                            :: Resul(PP_nVar),Resur(PP_nVar)
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
IF (PRESENT(RefStateOpt)) THEN
  RefState = RefStateOpt
ELSE
  RefState = IniRefState
END IF

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
  Resu = RefStateCons(:,RefState)
CASE(2) ! sinus
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=2.*PP_Pi*Frequency
  ! base flow
  prim(DENS)   = 1.
  prim(VELV) = AdvVel
  prim(PRES)   = 1.
  Vel=prim(VELV)
  cent=x-Vel*tEval
  prim(DENS)=prim(DENS)*(1.+Amplitude*SIN(Omega*SUM(cent(1:3))))
  ! g(t)
  Resu(DENS)=prim(DENS) ! rho
  Resu(MOMV)=prim(DENS)*prim(VELV) ! rho*vel
  Resu(ENER)=prim(PRES)*sKappaM1+0.5*SUM(Resu(MOMV)*prim(VELV)) ! rho*e

  IF(fullBoundaryOrder)THEN
    ov=Omega*SUM(vel)
    ! g'(t)
    Resu_t(DENS)=-Amplitude*cos(Omega*SUM(cent(1:3)))*ov
    Resu_t(MOMV)=Resu_t(DENS)*prim(VELV) ! rho*vel
    Resu_t(ENER)=0.5*SUM(Resu_t(MOMV)*prim(VELV))
    ! g''(t)
    Resu_tt(DENS)=-Amplitude*sin(Omega*SUM(cent(1:3)))*ov**2.
    Resu_tt(MOMV)=Resu_tt(DENS)*prim(VELV)
    Resu_tt(ENER)=0.5*SUM(Resu_tt(MOMV)*prim(VELV))
  END IF
CASE(21) ! sinus x
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=2.*PP_Pi*Frequency
  ! base flow
  prim(DENS)   = 1.
  prim(VELV) = AdvVel
  prim(PRES)   = 1.
  Vel=prim(VELV)
  cent=x-Vel*tEval
  prim(DENS)=prim(DENS)*(1.+Amplitude*SIN(Omega*cent(1)))
  ! g(t)
  Resu(DENS)=prim(DENS) ! rho
  Resu(MOMV)=prim(DENS)*prim(VELV) ! rho*vel
  Resu(ENER)=prim(PRES)*sKappaM1+0.5*SUM(Resu(MOMV)*prim(VELV)) ! rho*e

  IF(fullBoundaryOrder)THEN
    ov=Omega*SUM(vel)
    ! g'(t)
    Resu_t(DENS)=-Amplitude*cos(Omega*cent(1))*ov
    Resu_t(MOMV)=Resu_t(DENS)*prim(VELV) ! rho*vel
    Resu_t(ENER)=0.5*SUM(Resu_t(MOMV)*prim(VELV))
    ! g''(t)
    Resu_tt(DENS)=-Amplitude*sin(Omega*cent(1))*ov**2.
    Resu_tt(MOMV)=Resu_tt(DENS)*prim(VELV)
    Resu_tt(ENER)=0.5*SUM(Resu_tt(MOMV)*prim(VELV))
  END IF
CASE(3) ! linear in rho
  ! base flow
  prim(DENS)   = 100.
  prim(VELV)   = AdvVel
  prim(PRES)   = 1.
  Vel=prim(VELV)
  cent=x-Vel*tEval
  prim(DENS)=prim(DENS)+SUM(AdvVel*cent)
  ! g(t)
  Resu(DENS)=prim(DENS) ! rho
  Resu(MOMV)=prim(DENS)*prim(VELV) ! rho*vel
  Resu(ENER)=prim(PRES)*sKappaM1+0.5*SUM(Resu(MOMV)*prim(VELV)) ! rho*e
  IF(fullBoundaryOrder)THEN
    ! g'(t)
    Resu_t(DENS)=-SUM(Vel)
    Resu_t(MOMV)=Resu_t(DENS)*prim(VELV) ! rho*vel
    Resu_t(ENER)=0.5*SUM(Resu_t(MOMV)*prim(VELV))
  END IF
CASE(4) ! oblique sine wave (in x,y,z for 3D calculations, and x,y for 2D)
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi

  ! g(t)
#if (PP_dim == 3)
  Resu(DENS:MOM3) = 2.+ Amplitude*sin(Omega*SUM(x(1:PP_dim)) - a*tEval)
#else
  Resu(DENS:MOM2) = 2.+ Amplitude*sin(Omega*SUM(x(1:PP_dim)) - a*tEval)
  Resu(MOM3)   = 0.
#endif
  Resu(ENER)=Resu(DENS)*Resu(DENS)
  IF(fullBoundaryOrder)THEN
    ! g'(t)
#if (PP_dim == 3)
    Resu_t(DENS:MOM3)=(-a)*Amplitude*cos(Omega*SUM(x(1:PP_dim)) - a*tEval)
#else
    Resu_t(DENS:MOM2)=(-a)*Amplitude*cos(Omega*SUM(x(1:PP_dim)) - a*tEval)
    Resu_t(MOM3)=0.
#endif
    Resu_t(ENER)=2.*Resu(DENS)*Resu_t(DENS)
    ! g''(t)
#if (PP_dim == 3)
    Resu_tt(DENS:MOM3)=-a*a*Amplitude*sin(Omega*SUM(x(1:PP_dim)) - a*tEval)
#else
    Resu_tt(DENS:MOM2)=-a*a*Amplitude*sin(Omega*SUM(x(1:PP_dim)) - a*tEval)
    Resu_tt(MOM3)=0.
#endif
    Resu_tt(ENER)=2.*(Resu_t(DENS)*Resu_t(DENS) + Resu(DENS)*Resu_tt(DENS))
  END IF
CASE(41) ! SINUS in x
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(DENS:MOM1)=2.+ Amplitude*sin(Omega*x(1) - a*tEval)
  Resu(ENER)=Resu(DENS)*Resu(DENS)
  IF(fullBoundaryOrder)THEN
    Resu_t = 0.
    Resu_tt = 0.
    ! g'(t)
    Resu_t(DENS:MOM1)=(-a)*Amplitude*cos(Omega*x(1) - a*tEval)
    Resu_t(ENER)=2.*Resu(1)*Resu_t(1)
    ! g''(t)
    Resu_tt(DENS:MOM1)=-a*a*Amplitude*sin(Omega*x(1) - a*tEval)
    Resu_tt(ENER)=2.*(Resu_t(DENS)*Resu_t(DENS) + Resu(DENS)*Resu_tt(DENS))
  END IF
CASE(42) ! SINUS in y
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(2)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(DENS)=2.+ Amplitude*sin(Omega*x(2) - a*tEval)
  Resu(MOM2)=Resu(DENS)
  Resu(ENER)=Resu(DENS)*Resu(DENS)
  IF(fullBoundaryOrder)THEN
    Resu_tt = 0.
    ! g'(t)
    Resu_t(DENS)=(-a)*Amplitude*cos(Omega*x(2) - a*tEval)
    Resu_t(MOM2)=Resu_t(DENS)
    Resu_t(ENER)=2.*Resu(DENS)*Resu_t(DENS)
    ! g''(t)
    Resu_tt(DENS)=-a*a*Amplitude*sin(Omega*x(2) - a*tEval)
    Resu_tt(MOM2)=Resu_tt(DENS)
    Resu_tt(ENER)=2.*(Resu_t(DENS)*Resu_t(DENS) + Resu(DENS)*Resu_tt(DENS))
  END IF
#if PP_dim==3
CASE(43) ! SINUS in z
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(3)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(DENS)=2.+ Amplitude*sin(Omega*x(3) - a*tEval)
  Resu(MOM3)=Resu(DENS)
  Resu(ENER)=Resu(DENS)*Resu(DENS)
  IF(fullBoundaryOrder)THEN
    Resu_tt = 0.
    ! g'(t)
    Resu_t(DENS)=(-a)*Amplitude*cos(Omega*x(3) - a*tEval)
    Resu_t(MOM3)=Resu_t(DENS)
    Resu_t(ENER)=2.*Resu(DENS)*Resu_t(DENS)
    ! g''(t)
    Resu_tt(DENS)=-a*a*Amplitude*sin(Omega*x(3) - a*tEval)
    Resu_tt(MOM3)=Resu_tt(DENS)
    Resu_tt(ENER)=2.*(Resu_t(DENS)*Resu_t(DENS) + Resu(DENS)*Resu_tt(DENS))
  END IF
#endif
CASE(5) !Roundjet Bogey Bailly 2002, Re=65000, x-axis is jet axis
  prim(DENS)  =1.
  prim(VELV)  =0.
  prim(PRES)  =1./Kappa
  prim(TEMP)  = prim(PRES)/(prim(DENS)*R)
  ! Jet inflow (from x=0, diameter 2.0)
  ! Initial jet radius: rj=1.
  ! Momentum thickness: delta_theta0=0.05=1/20
  ! Re=65000
  ! Uco=0.
  ! Uj=0.9
  r_len=SQRT((x(2)*x(2)+x(3)*x(3)))
  prim(VEL1)=0.9*0.5*(1.+TANH((JetRadius-r_len)/JetRadius*10.))
  CALL RANDOM_NUMBER(random)
  ! Random disturbance +-5%; uniform distribution between -1,1
  random=0.05*2.*(random-0.5)
  prim(VEL1)=prim(VEL1)+random*prim(VEL1)
  prim(VEL2)=x(2)/r_len*0.5*random*prim(VEL1)
  prim(VEL3)=x(3)/r_len*0.5*random*prim(VEL1)
  CALL PrimToCons(prim,ResuL)
  prim(VELV)  =0.
  CALL PrimToCons(prim,ResuR)
  ! after x/r0=10 blend to ResuR
  Resu=ResuL+(ResuR-ResuL)*0.5*(1.+tanh(x(1)/JetRadius-JetEnd))
CASE(51)
  ! Parabolic velocity distribution following
  ! "Benchmark Computations of Laminar Flow Around a Cylinder", Schäfer and Turek, 1996.
  ! https://doi.org/10.1007/978-3-322-89849-4_39
  ! ATTENTION: In contrast to paper, velocity profile is defined around y,z=0. and NOT y,z=H/2.
  prim = RefStatePrim(:,RefState)
  prim(VELV) = 0.
  IF (NORM2(x(2:PP_dim)).LE.JetRadius) THEN
    prim(VEL1) = JetAmplitude*4.*(JetRadius-(x(2)))*(JetRadius+(x(2)))/(2.*JetRadius)**2 ! 2D-1, 2D-2
#if PP_dim == 3
    ! Multiply contribution of z-dimension
    prim(VEL1) = prim(VEL1)  *4.*(JetRadius-(x(3)))*(JetRadius+(x(3)))/(2.*JetRadius)**2 ! 3D-1, 3D-2
#endif
  END IF
  CALL PrimToCons(prim,resu)
CASE(6)  ! Cylinder flow
  IF(tEval .EQ. 0.)THEN   ! Initialize potential flow
    prim(DENS)=RefStatePrim(DENS,RefState)  ! Density
    prim(VEL3)=0.                           ! VelocityZ=0. (2D flow)
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
    prim(VEL1)=RefStatePrim(VEL1,RefState)*(COS(phi)**2*(1.-0.25/radius)+SIN(phi)**2*(1.+0.25/radius))
    prim(VEL2)=RefStatePrim(VEL1,RefState)*(-2.)*SIN(phi)*COS(phi)*0.25/radius
    ! Calculate pressure, RefState(2)=u_infinity
    prim(PRES)=RefStatePrim(PRES,RefState) + &
            0.5*prim(DENS)*(RefStatePrim(VEL1,RefState)*RefStatePrim(VEL1,RefState)-prim(VEL1)*prim(VEL1)-prim(VEL2)*prim(VEL2))
    prim(TEMP) = prim(PRES)/(prim(DENS)*R)
  ELSE  ! Use RefState as BC
    prim=RefStatePrim(:,RefState)
  END IF  ! t=0
  CALL PrimToCons(prim,resu)
CASE(7) ! SHU VORTEX,isentropic vortex
  ! base flow
  prim=RefStatePrim(:,RefState)  ! Density
  ! ini-Parameter of the Example
  vel=prim(VELV)
  RT=prim(PRES)/prim(DENS) !ideal gas
  cent=(iniCenter+vel*tEval)    !centerpoint time dependant
  cent=x-cent                   ! distance to centerpoint
  cent=CROSS(iniAxis,cent)      !distance to axis, tangent vector, length r
  cent=cent/iniHalfWidth        !Halfwidth is dimension 1
  r2=SUM(cent*cent) !
  du = IniAmplitude/(2.*PP_Pi)*exp(0.5*(1.-r2))   ! vel. perturbation
  dTemp = -kappaM1/(2.*kappa*RT)*du**2            ! adiabatic
  prim(DENS)=prim(DENS)*(1.+dTemp)**(1.*skappaM1) !rho
  prim(VELV)=prim(VELV)+du*cent(:)                !v
#if PP_dim == 2
  prim(VEL3)=0.
#endif
  prim(PRES) = prim(PRES)*(1.+dTemp)**(kappa/kappaM1) !p
  prim(TEMP) = prim(PRES)/(prim(DENS)*R)
  CALL PrimToCons(prim,resu)
CASE(8) !Couette-Poiseuille flow between plates: exact steady lamiar solution with height=1 !
        !(Gao, hesthaven, Warburton)
  RT=1. ! Hesthaven: Absorbing layers for weakly compressible flows
  sRT=1./RT
  ! size of domain must be [-0.5, 0.5]^2 -> (x*y)
  h=0.5
  prim(VEL1)     = U_parameter*(0.5*(1.+x(2)/h) + P_parameter*(1-(x(2)/h)**2))
  prim(VEL2:VEL3)= 0.
  pexit=0.9996
  pentry=1.0004
  prim(PRES)= ( ( x(1) - (-0.5) )*( pexit - pentry) / ( 0.5 - (-0.5)) ) + pentry
  prim(DENS)=prim(PRES)*sRT
  prim(TEMP)=prim(PRES)/(prim(DENS)*R)
  CALL PrimToCons(prim,Resu)
CASE(9) !lid driven cavity flow from Gao, Hesthaven, Warburton
        !"Absorbing layers for weakly compressible flows", to appear, JSC, 2016
        ! Special "regularized" driven cavity BC to prevent singularities at corners
        ! top BC assumed to be in x-direction from 0..1
  Prim = RefStatePrim(:,RefState)
  IF (x(1).LT.0.2) THEN
    prim(VEL1)=1000*4.9333*x(1)**4-1.4267*1000*x(1)**3+0.1297*1000*x(1)**2-0.0033*1000*x(1)
  ELSEIF (x(1).LE.0.8) THEN
    prim(VEL1)=1.0
  ELSE
    prim(VEL1)=1000*4.9333*x(1)**4-1.8307*10000*x(1)**3+2.5450*10000*x(1)**2-1.5709*10000*x(1)+10000*0.3633
  ENDIF
  CALL PrimToCons(prim,Resu)
CASE(10) ! shock
  prim=0.

  ! pre-shock
  prim(DENS) = PreShockDens
  Ms         = MachShock

  prim(PRES)=prim(DENS)/Kappa
  prim(TEMP)=prim(PRES)/(prim(DENS)*R)
  CALL PrimToCons(prim,Resur)

  ! post-shock
  prim(VEL2)=prim(DENS) ! temporal storage of pre-shock density
  prim(DENS)=prim(DENS)*((KappaP1)*Ms*Ms)/(KappaM1*Ms*Ms+2.)
  prim(PRES)=prim(PRES)*(2.*Kappa*Ms*Ms-KappaM1)/(KappaP1)
  prim(TEMP)=prim(PRES)/(prim(DENS)*R)
  IF (prim(VEL1) .EQ. 0.0) THEN
    prim(VEL1)=Ms*(1.-prim(VEL2)/prim(DENS))
  ELSE
    prim(VEL1)=prim(VEL1)*prim(VEL2)/prim(DENS)
  END IF
  prim(3)=0. ! reset temporal storage
  CALL PrimToCons(prim,Resul)
  xs=5.+Ms*tEval ! 5. bei 10x10x10 Rechengebiet
  ! Tanh boundary
  Resu=-0.5*(Resul-Resur)*TANH(5.0*(x(1)-xs))+Resur+0.5*(Resul-Resur)
CASE(11) ! Sod Shock tube
  IF (X(1).LE.ShockPos) THEN
    Resu = RefStateCons(:,1)
  ELSE
    Resu = RefStateCons(:,2)
  END IF
CASE(111) ! Sedov blast wave
  Ms = SQRT(SUM((X(1:3)-1.5)**2))
  IF ((Ms.LE.ShockPos).AND.(Ms.NE.0)) THEN
    prim(DENS)      = 1.3416
    prim(VEL1:VEL3) = 0.3615*(X(1:3)-1.5)/Ms
    prim(PRES)      = 1.5133
  ELSE
    prim(DENS)      = 1.
    prim(VELV)      = 0.
    prim(PRES)      = 1.
  END IF
  prim(TEMP)=prim(PRES)/(prim(DENS)*R)
  CALL PrimToCons(prim,resu)
CASE(12) ! Shu Osher density fluctuations shock wave interaction
  IF (x(1).LT.-4.0) THEN
    prim(DENS)      = 3.857143
    prim(VEL1)      = 2.629369
    prim(VEL2:VEL3) = 0.
    prim(PRES)      = 10.33333
  ELSE
    prim(DENS)      = 1.+0.2*SIN(5.*x(1))
    prim(VELV)      = 0.
    prim(PRES)      = 1.
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
CASE(14) ! harmonic gauss pulse
  Resu = RefStateCons(:,RefState)
#if PARABOLIC
CASE(1338) ! blasius
  prim=RefStatePrim(:,RefState)
  ! calculate equivalent x for Blasius flat plate to have delta99_in at x_in
  x_offset(1)=(delta99_in/5)**2*prim(DENS)*prim(VEL1)/mu0-x_in(1)
  x_offset(2)=-x_in(2)
  x_offset(3)=0.
  x_eff=x+x_offset
  IF(x_eff(2).GT.0 .AND. x_eff(1).GT.0) THEN
    ! scale bl position in physical space to reference space, eta=5 is ~99% bl thickness
    eta=x_eff(2)*(prim(DENS)*prim(VEL1)/(mu0*x_eff(1)))**0.5

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
    prim(VEL2)=0.5*(mu0*prim(VEL1)/prim(DENS)/x_eff(1))**0.5*(fp*eta-f)
    prim(VEL1)=RefStatePrim(VEL1,RefState)*fp
  ELSE
    IF(x_eff(2).LE.0) THEN
      prim(VEL1)=0.
    END IF
  END IF
  CALL PrimToCons(prim,resu)
#endif
CASE(1342) ! sedov
  IF (X(1)**2.0 + X(2)**2.0.GE.eradius_in**2.0) THEN
    Resu = RefStateCons(:,1)
  ELSE
    Resu = RefStateCons(:,2)
  END IF
CASE(1344) ! LeBlanc Shock Tube
  IF (X(1).LE.leblanc_x0) THEN
    Resu = RefStateCons(:,1)
  ELSE
    Resu = RefStateCons(:,2)
  END IF
CASE(1346) ! Hui Shock 2D
  IF (X(2) - ((sin(Hui_Data(1)) + cos(Hui_Data(1)))*X(1) - Hui_Data(2))/(sin(Hui_Data(1)) - cos(Hui_Data(1))) .GT. 0 ) THEN
    Resu = RefStateCons(:,1)
  ELSE
    Resu = RefStateCons(:,2)
  END IF
END SELECT ! ExactFunction
#if PP_dim==2
Resu(MOM3)=0.
#endif

! For O3 LS 3-stage RK, we have to define proper time dependent BC
IF(fullBoundaryOrder)THEN ! add resu_t, resu_tt if time dependant
#if PP_dim==2
  Resu_t=0.
  Resu_tt=0.
#endif
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
USE MOD_EOS_Vars         ,ONLY: Kappa,KappaM1
USE MOD_Equation_Vars    ,ONLY: IniExactFunc,doCalcSource
USE MOD_Exactfunc_Vars   ,ONLY: AdvVel,IniAmplitude,IniFrequency
USE MOD_Exactfunc_Vars   ,ONLY: HarmonicFrequency,AmplitudeFactor,SiqmaSqr
USE MOD_Mesh_Vars        ,ONLY: Elem_xGP,sJ,nElems
#if PARABOLIC
USE MOD_EOS_Vars         ,ONLY: mu0,Pr
#endif
#if FV_ENABLED
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisVolume
USE MOD_FV_Vars          ,ONLY: FV_Vdm,FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< current solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
REAL                :: Ut_src(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL                :: Frequency,Amplitude,Omega,a
REAL                :: sinXGP,sinXGP2,cosXGP,at
REAL                :: tmp(6)
REAL                :: C
#if FV_ENABLED
REAL                :: Ut_src2(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#endif
!==================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(4) ! exact function
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi
  tmp(1)=-a+REAL(PP_dim)*Omega
#if (PP_dim == 3)
  tmp(2)=-a+0.5*Omega*(1.+kappa*5.)
#else
  tmp(2)=-a+Omega*(-1.+kappa*3.)
#endif
  tmp(3)=Amplitude*Omega*KappaM1
#if (PP_dim == 3)
  tmp(4)=0.5*((9.+Kappa*15.)*Omega-8.*a)
  tmp(5)=Amplitude*(3.*Omega*Kappa-a)
#else
  tmp(4)=((2.+Kappa*6.)*Omega-4.*a)
  tmp(5)=Amplitude*(2.*Omega*Kappa-a)
#endif
#if PARABOLIC
  tmp(6)=REAL(PP_dim)*mu0*Kappa*Omega*Omega/Pr
#else
  tmp(6)=0.
#endif
  tmp=tmp*Amplitude
  at=a*t
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      cosXGP=COS(omega*SUM(Elem_xGP(1:PP_dim,i,j,k,iElem))-at)
      sinXGP=SIN(omega*SUM(Elem_xGP(1:PP_dim,i,j,k,iElem))-at)
      sinXGP2=2.*sinXGP*cosXGP !=SIN(2.*(omega*SUM(Elem_xGP(:,i,j,k,iElem))-a*t))
      Ut_src(DENS      ,i,j,k) = tmp(1)*cosXGP
#if (PP_dim == 3)
      Ut_src(MOM1:MOM3,i,j,k)  = tmp(2)*cosXGP + tmp(3)*sinXGP2
#else
      Ut_src(MOM1:MOM2,i,j,k)  = tmp(2)*cosXGP + tmp(3)*sinXGP2
      Ut_src(MOM3      ,i,j,k) = 0.
#endif
      Ut_src(ENER,i,j,k)       = tmp(4)*cosXGP + tmp(5)*sinXGP2 + tmp(6)*sinXGP
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
      CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
      END DO; END DO; END DO ! i,j,k
    ELSE
#endif
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
    END IF
#endif
  END DO ! iElem
CASE(14) ! Harmonic Gausspulse
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut_src(1,i,j,k) = AmplitudeFactor*cos(2.*PP_Pi*HarmonicFrequency*t)*1/sqrt(((2*PP_Pi)**2)*2*SiqmaSqr)*EXP(-0.5*SUM(Elem_xGP(1:2,i,j,k,iElem)**2)/SiqmaSqr)
      Ut_src(2:5,i,j,k) = 0.0
    END DO; END DO; END DO ! i,j,k
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
  END DO
CASE(41) ! Sinus in x
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi
  C = 2.0

  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
#if PARABOLIC
      Ut_src(DENS,i,j,k) = (-Amplitude*a+Amplitude*omega)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)

      Ut_src(MOM1,i,j,k) = (-Amplitude**2*omega+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(1,i,j,k,iElem)-2.*a*t)+&
                           (-Amplitude*a+2.*Amplitude*omega*kappa*C-1./2.*Amplitude*omega*kappa+&
                           3./2.*Amplitude*omega-2.*Amplitude*omega*C)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)

      Ut_src(MOM2:MOM3,i,j,k) = 0.0

      Ut_src(ENER,i,j,k) = mu0*kappa*Amplitude*sin(omega*Elem_xGP(1,i,j,k,iElem)-a*t)*omega**2/Pr+&
                           (-Amplitude**2*a+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(1,i,j,k,iElem)-2.*a*t)+&
                           1./2.*(-4.*Amplitude*a*C+4.*Amplitude*omega*kappa*C-Amplitude*omega*kappa+&
                                Amplitude*omega)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)
#else
      Ut_src(DENS,i,j,k) = (-amplitude*a+amplitude*omega)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)
      Ut_src(MOM1,i,j,k) = (-amplitude**2*omega+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(1,i,j,k,iElem)-2.*a*t)+ &
                           (-amplitude*a+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+ &
                           3./2.*amplitude*omega-2.*amplitude*omega*C)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)

      Ut_src(MOM2:MOM3,i,j,k) = 0.0
      Ut_src(ENER,i,j,k) = (-amplitude**2*a+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(1,i,j,k,iElem)-2.*a*t)+&
                           (-2.*amplitude*a*C+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+&
                           1./2.*amplitude*omega)*cos(omega*Elem_xGP(1,i,j,k,iElem)-a*t)

#endif
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
      CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
      END DO; END DO; END DO ! i,j,k
    ELSE
#endif
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
    END IF
#endif
  END DO
CASE(42) ! Sinus in y
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(2)*2.*PP_Pi
  C = 2.0

  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
#if PARABOLIC
      Ut_src(DENS,i,j,k) = (-Amplitude*a+Amplitude*omega)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)
      Ut_src(MOM1,i,j,k) = 0.0
      Ut_src(MOM2,i,j,k) = (-Amplitude**2*omega+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(2,i,j,k,iElem)-2.*a*t)+&
                           (-Amplitude*a+2.*Amplitude*omega*kappa*C-1./2.*Amplitude*omega*kappa+&
                           3./2.*Amplitude*omega-2.*Amplitude*omega*C)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)

      Ut_src(MOM3,i,j,k) = 0.0

      Ut_src(ENER,i,j,k) = mu0*kappa*Amplitude*sin(omega*Elem_xGP(2,i,j,k,iElem)-a*t)*omega**2/Pr+&
                           (-Amplitude**2*a+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(2,i,j,k,iElem)-2.*a*t)+&
                           1./2.*(-4.*Amplitude*a*C+4.*Amplitude*omega*kappa*C-Amplitude*omega*kappa+&
                                Amplitude*omega)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)
#else
      Ut_src(DENS,i,j,k) = (-amplitude*a+amplitude*omega)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)
      Ut_src(MOM1,i,j,k) = 0.0
      Ut_src(MOM2,i,j,k) = (-amplitude**2*omega+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(2,i,j,k,iElem)-2.*a*t)+ &
                           (-amplitude*a+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+ &
                           3./2.*amplitude*omega-2.*amplitude*omega*C)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)

      Ut_src(MOM3,i,j,k) = 0.0
      Ut_src(ENER,i,j,k) = (-amplitude**2*a+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(2,i,j,k,iElem)-2.*a*t)+&
                           (-2.*amplitude*a*C+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+&
                           1./2.*amplitude*omega)*cos(omega*Elem_xGP(2,i,j,k,iElem)-a*t)

#endif
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
      CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
      END DO; END DO; END DO ! i,j,k
    ELSE
#endif
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
    END IF
#endif
  END DO

#if PP_dim==3
CASE(43) ! Sinus in z
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(3)*2.*PP_Pi
  C = 2.0

  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
#if PARABOLIC
      Ut_src(DENS,i,j,k) = (-Amplitude*a+Amplitude*omega)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)
      Ut_src(MOM1:MOM2,i,j,k) = 0.0
      Ut_src(MOM3,i,j,k) = (-Amplitude**2*omega+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(3,i,j,k,iElem)-2.*a*t)+&
                           (-Amplitude*a+2.*Amplitude*omega*kappa*C-1./2.*Amplitude*omega*kappa+&
                           3./2.*Amplitude*omega-2.*Amplitude*omega*C)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)


      Ut_src(ENER,i,j,k) = mu0*kappa*Amplitude*sin(omega*Elem_xGP(3,i,j,k,iElem)-a*t)*omega**2/Pr+&
                           (-Amplitude**2*a+Amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(3,i,j,k,iElem)-2.*a*t)+&
                           1./2.*(-4.*Amplitude*a*C+4.*Amplitude*omega*kappa*C-Amplitude*omega*kappa+&
                                Amplitude*omega)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)
#else
      Ut_src(DENS,i,j,k) = (-amplitude*a+amplitude*omega)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)
      Ut_src(MOM1:MOM2,i,j,k) = 0.0
      Ut_src(MOM3,i,j,k) = (-amplitude**2*omega+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(3,i,j,k,iElem)-2.*a*t)+ &
                           (-amplitude*a+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+ &
                           3./2.*amplitude*omega-2.*amplitude*omega*C)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)

      Ut_src(ENER,i,j,k) = (-amplitude**2*a+amplitude**2*omega*kappa)*sin(2.*omega*Elem_xGP(3,i,j,k,iElem)-2.*a*t)+&
                           (-2.*amplitude*a*C+2.*amplitude*omega*kappa*C-1./2.*omega*kappa*amplitude+&
                           1./2.*amplitude*omega)*cos(omega*Elem_xGP(3,i,j,k,iElem)-a*t)

#endif
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
      CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
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
#endif
CASE DEFAULT
  ! No source -> do nothing and set marker to not run again
  doCalcSource=.FALSE.
END SELECT ! ExactFunction
END SUBROUTINE CalcSource

END MODULE MOD_Exactfunc
