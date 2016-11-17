!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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
!> Riemann2D testcase 
!==================================================================================================================================
MODULE MOD_Testcase
! MODULES
USE MOD_Testcase_Vars
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE DefineParametersTestcase
  MODULE PROCEDURE DO_NOTHING
End INTERFACE

INTERFACE InitTestcase
  MODULE PROCEDURE InitTestcase
END INTERFACE

INTERFACE FinalizeTestcase
  MODULE PROCEDURE DO_NOTHING
END INTERFACE

INTERFACE ExactFuncTestcase
  MODULE PROCEDURE ExactFuncTestcase
END INTERFACE

INTERFACE CalcForcing
  MODULE PROCEDURE DO_NOTHING
END INTERFACE

INTERFACE AnalyzeTestCase
  MODULE PROCEDURE DO_NOTHING
END INTERFACE

INTERFACE GetBoundaryFluxTestcase
  MODULE PROCEDURE GetBoundaryFluxTestcase
END INTERFACE

INTERFACE GetBoundaryFVgradientTestcase
  MODULE PROCEDURE GetBoundaryFVgradientTestcase
END INTERFACE

INTERFACE Lifting_GetBoundaryFluxTestcase
  MODULE PROCEDURE Lifting_GetBoundaryFluxTestcase
END INTERFACE

PUBLIC:: DefineParametersTestcase
PUBLIC:: InitTestcase
PUBLIC:: FinalizeTestcase
PUBLIC:: ExactFuncTestcase
PUBLIC:: TestcaseSource
PUBLIC:: CalcForcing
PUBLIC:: AnalyzeTestCase
PUBLIC:: GetBoundaryFluxTestcase
PUBLIC:: GetBoundaryFVgradientTestcase
PUBLIC:: Lifting_GetBoundaryFluxTestcase

CONTAINS

!==================================================================================================================================
!> Empty placeholder routine
!==================================================================================================================================
SUBROUTINE DO_NOTHING(optionalREAL,optionalREAL2)
IMPLICIT NONE
REAL,OPTIONAL,INTENT(IN)  :: optionalREAL,optionalREAL2
END SUBROUTINE DO_NOTHING


!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE InitTestcase()
! MODULES
USE MOD_Globals
USE MOD_Testcase_Vars
USE MOD_Equation_Vars ,ONLY: IniExactFunc
USE MOD_Equation_Vars ,ONLY: RefStatePrim
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TESTCASE Riemann2D...'

CALL CalcIniStates()

CALL Riemann_Speeds(1,RiemannBC_WaveType(1),RefStatePrim(2,:),RefStatePrim(1,:),RiemannBC_speeds(1,:))
CALL Riemann_Speeds(2,RiemannBC_WaveType(2),RefStatePrim(3,:),RefStatePrim(2,:),RiemannBC_speeds(2,:))
CALL Riemann_Speeds(1,RiemannBC_WaveType(3),RefStatePrim(3,:),RefStatePrim(4,:),RiemannBC_speeds(3,:))
CALL Riemann_Speeds(2,RiemannBC_WaveType(4),RefStatePrim(4,:),RefStatePrim(1,:),RiemannBC_speeds(4,:))

IF (MPIRoot) THEN
  WRITE(*,*) " Calculate states, shock, discontinuity and expansion speeds for all 4 sides"
  WRITE(*,"(A,I3)") "  Configuration",-IniExactFunc
  WRITE(*,'(31(" "),A13," speed=",2F8.4)') WAVENAMES(RiemannBC_WaveType(1)), RiemannBC_speeds(1,:)
  WRITE(*,'(14(" "),"|",21("-"),"|",21("-")"|")') 
  WRITE(*,'(14(" "),"|",21(" "),"|",21(" ")"|")') 
  WRITE(*,'(14(" "),"| p=",F7.4," r=",F7.4," | p=",F7.4," r=",F7.4," |")')  RefStatePrim(2,5), RefStatePrim(2,1),&
      RefStatePrim(1,5), RefStatePrim(1,1) 
  WRITE(*,'(14(" "),"| u=",F7.4," v=",F7.4," | u=",F7.4," v=",F7.4," |")')  RefStatePrim(2,2), RefStatePrim(2,3),&
      RefStatePrim(1,2), RefStatePrim(1,3) 
  WRITE(*,'(14(" "),"|",21(" "),"|",21(" ")"|")') 
  WRITE(*,'(A13," |",21("-"),"|",21("-")"| ",A13)') ADJUSTR(TRIM(WAVENAMES(RiemannBC_WaveType(2)))), WAVENAMES(RiemannBC_WaveType(4))
  WRITE(*,'("speed=",8(" "),"|",21(" "),"|",21(" ")"| speed=")') 
  WRITE(*,'(F13.4," | p=",F7.4," r=",F7.4," | p=",F7.4," r=",F7.4," | ",F7.4)') RiemannBC_speeds(2,1), &
      RefStatePrim(3,5), RefStatePrim(3,1), RefStatePrim(4,5), RefStatePrim(4,1), RiemannBC_speeds(4,1)
  WRITE(*,'(F13.4," | u=",F7.4," v=",F7.4," | u=",F7.4," v=",F7.4," | ",F7.4)') RiemannBC_speeds(2,2), &
      RefStatePrim(3,2), RefStatePrim(3,3), RefStatePrim(4,2), RefStatePrim(4,3), RiemannBC_speeds(4,2)
  WRITE(*,'(14(" "),"|",21(" "),"|",21(" ")"|")') 
  WRITE(*,'(14(" "),"|",21("-"),"|",21("-")"|")') 
  WRITE(*,'(31(" "),A13," speed=",2F8.4)') WAVENAMES(RiemannBC_WaveType(3)), RiemannBC_speeds(3,:)
ENDIF

SWRITE(UNIT_stdOut,'(A)')' INIT TESTCASE Riemann2D DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitTestcase

FUNCTION GetPHI(rhoL, rhoR, pL, pR) 
USE MOD_EOS_Vars      ,ONLY: Kappa, KappaM1
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN) :: rhoL, rhoR, pL, pR
REAL            :: GetPHI
!===================================================================================================================================
  GetPHI =  2.*SQRT(Kappa)/KappaM1 * (SQRT(pL/rhoL) - SQRT(pR/rhoR))
END FUNCTION GetPHI

FUNCTION GetPSI(rhoL,rhoR,pL,pR) 
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN) :: rhoL, rhoR, pL, pR
REAL            :: GetPSI
!===================================================================================================================================
  GetPSI =  SQRT( ((pL-pR) * (rhoL-rhoR)) / &
      (rhoL*rhoR) )
END FUNCTION GetPSI

FUNCTION GetPI(pL, pR) 
USE MOD_EOS_Vars      ,ONLY: KappaM1, KappaP1
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN) :: pL, pR
REAL            :: GetPI
!===================================================================================================================================
  GetPI =  ( pL/pR + KappaM1/KappaP1 ) / &
      ( 1. + KappaM1/KappaP1 * pL /  pR )
END FUNCTION GetPI


FUNCTION GetRho_RankineHugoniot(rhoL, pL, pR) 
USE MOD_EOS_Vars      ,ONLY: KappaM1, KappaP1
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN) :: rhoL, pL, pR
REAL            :: GetRho_RankineHugoniot
!===================================================================================================================================
  GetRho_RankineHugoniot = rhoL * (pL*KappaM1 + pR*KappaP1) / (pL*KappaP1 + pR*KappaM1)
END FUNCTION GetRho_RankineHugoniot

FUNCTION GetRho_Rarefaction(rhoL, pL, pR) 
USE MOD_EOS_Vars      ,ONLY: Kappa
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN) :: rhoL, pL, pR
REAL            :: GetRho_Rarefaction
!===================================================================================================================================
  GetRho_Rarefaction = (pr/pL)**(1./Kappa) * rhoL
END FUNCTION GetRho_Rarefaction

FUNCTION GetP_RankineHugoniot(rhoL, rhoR, pL) 
USE MOD_EOS_Vars      ,ONLY: KappaM1, KappaP1
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN) :: rhoL, rhoR, pL
REAL            :: GetP_RankineHugoniot
!===================================================================================================================================
  GetP_RankineHugoniot = pL* (KappaP1 - KappaM1*rhoL/rhoR) / (KappaP1*rhoL/rhoR - KappaM1)
END FUNCTION GetP_RankineHugoniot


SUBROUTINE Calc_p_rho_RankineHugoniot(rhoL, vL, vR, pL, rhoR, pR) 
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaP1
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)  :: rhoL
REAL,INTENT(IN)  :: vL
REAL,INTENT(IN)  :: vR
REAL,INTENT(IN)  :: pL
REAL,INTENT(OUT) :: rhoR(2)
REAL,INTENT(OUT) :: pR(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: z1,z2,z3
!===================================================================================================================================
  z1 = KappaP1*(vL-vR)*rhoL**2
  z2 = sqrt(Kappa**2*(vL-vR)**2*rhoL + 2.*Kappa*(vL-vR)**2*rhoL + (vL-vR)**2 *rhoL + 16.*Kappa*pL)
  z3 = (Kappa-3.)*(vL-vR)*rhoL

  rhoR(1)=(z1 - rhoL**1.5 * z2) / (z3 - sqrt(rhoL)*z2)
  rhoR(2)=(z1 + rhoL**1.5 * z2) / (z3 + sqrt(rhoL)*z2)

  z1 = 4.*pL+KappaP1*(vL-vR)**2*rhoL
  z2 = (vL-vR)*sqrt((1.+2.*Kappa+Kappa**2)*(vL-vR)**2*rhoL**2.+16.*Kappa*pL*rhoL)

  pR(1) = (z1 - z2) / 4.
  pR(2) = (z1 + z2) / 4.

END SUBROUTINE Calc_p_rho_RankineHugoniot

SUBROUTINE Calc_v_RankineHugoniot(rhoL, rhoR,  vL, pL, pR, vR) 
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)  :: rhoL
REAL,INTENT(IN)  :: rhoR
REAL,INTENT(IN)  :: vL
REAL,INTENT(IN)  :: pL
REAL,INTENT(IN)  :: pR
REAL,INTENT(OUT) :: vR(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  vR(1) = -(sqrt((pR-pL)*rhoL*rhoR**2+(pL-pR)*rhoL**2*rhoR)-rhoL*rhoR*vL)/(rhoL*rhoR)
  vR(2) =  (sqrt((pR-pL)*rhoL*rhoR**2+(pL-pR)*rhoL**2*rhoR)+rhoL*rhoR*vL)/(rhoL*rhoR)
END SUBROUTINE Calc_v_RankineHugoniot


SUBROUTINE CalcIniStates() 
USE MOD_Globals
USE MOD_PreProc
USE MOD_Testcase_Vars
USE MOD_Equation_Vars ,ONLY: IniExactFunc
USE MOD_Equation_Vars ,ONLY: RefStatePrim,RefStateCons, nRefState
USE MOD_EOS_Vars      ,ONLY: Kappa,R
USE MOD_EOS           ,ONLY: PrimToCons
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: p(4)
REAL    :: rho(4)
REAL    :: u(4)
REAL    :: v(4)
REAL    :: rho_vec(2)
REAL    :: p_vec(2)
REAL    :: v_vec(2)
INTEGER :: i
REAL    :: UE(PP_2Var)
!===================================================================================================================================

SELECT CASE(-IniExactFunc)
CASE(3)  ! Schulz-Rinne: Cfg 3
  RiemannBC_WaveType(1) = SHOCK
  RiemannBC_WaveType(2) = SHOCK
  RiemannBC_WaveType(3) = SHOCK
  RiemannBC_WaveType(4) = SHOCK
  p(1) = 1.5
  p(2) = 0.3
  rho(1) = 1.5
  u(1) = 0.
  v(1) = 0.

  rho(2) = GetRho_RankineHugoniot(rho(1), p(1), p(2))
  u(2) = GetPSI(rho(2), rho(1), p(2), p(1))+u(1)
  v(2) = v(1)


  p(4) = p(2)
  rho(4) = GetRho_RankineHugoniot(rho(1), p(1), p(2))
  u(4) = u(1)
  v(4) = GetPSI(rho(4), rho(1), p(4), p(1))+v(1)

  u(3) = u(2)
  v(3) = v(4)

  CALL Calc_p_rho_RankineHugoniot(rho(4), u(4), u(3), p(4), rho_vec, p_vec)
  rho(3) = MINVAL(rho_vec)
  p(3) = MINVAL(p_vec)

CASE(4) ! Schulz-Rinne: Cfg 4
  RiemannBC_WaveType(1) = SHOCK
  RiemannBC_WaveType(2) = SHOCK
  RiemannBC_WaveType(3) = SHOCK
  RiemannBC_WaveType(4) = SHOCK
  p(1) = 1.1
  p(2) = 0.35
  rho(1) = 1.1
  u(1) = 0.
  v(1) = 0.

  rho(2) = GetRho_RankineHugoniot(rho(1), p(1), p(2))
  CALL Calc_v_RankineHugoniot(rho(1), rho(2), u(1), p(1), p(2), v_vec)
  u(2) = MAXVAL(v_vec)

  p(4) = p(2)
  rho(4) = rho(2)
  CALL Calc_v_RankineHugoniot(rho(1), rho(4), v(1), p(1), p(4), v_vec)
  v(4) = MAXVAL(v_vec)

  p(3) = p(1)
  rho(3) = rho(1)
  CALL Calc_v_RankineHugoniot(rho(4), rho(3), u(4), p(4), p(3), v_vec)
  u(3) = MAXVAL(v_vec)
  CALL Calc_v_RankineHugoniot(rho(2), rho(3), v(2), p(2), p(3), v_vec)
  v(3) = MAXVAL(v_vec)

CASE(11) ! Schulz-Rinne: Cfg E
  RiemannBC_WaveType(1) = SHOCK
  RiemannBC_WaveType(2) = DISCONTINUITY
  RiemannBC_WaveType(3) = DISCONTINUITY
  RiemannBC_WaveType(4) = SHOCK
  p(1) = 1.
  p(2) = 0.4
  rho(1) = 1.
  rho(3) = 0.8
  u(1) = 0.1
  v(1) = 0.

  rho(2) = GetRho_RankineHugoniot(rho(1), p(1), p(2))
  u(2) = GetPSI(rho(2), rho(1), p(2), p(1)) + u(1)
  v(2) = v(1)
  
  p(4) = p(2)
  rho(4) = GetRho_RankineHugoniot(rho(1), p(1), p(4))
  u(4) = u(1)
  v(4) = GetPSI(rho(4), rho(1), p(4), p(1)) + v(1)

  p(3) = p(2)
  u(3) = u(1)
  v(3) = v(1)

CASE(12) ! Schulz-Rinne: Cfg F
  RiemannBC_WaveType(1) = SHOCK
  RiemannBC_WaveType(2) = DISCONTINUITY
  RiemannBC_WaveType(3) = DISCONTINUITY
  RiemannBC_WaveType(4) = SHOCK
  p(1) = 0.4
  p(2) = 1.
  rho(2) = 1.
  rho(3) = 0.8
  u(1) = 0.
  v(1) = 0.

  rho(1) = GetRho_RankineHugoniot(rho(2), p(2), p(1))

  u(2) = GetPSI(rho(2), rho(1), p(2), p(1)) + u(1)
  v(2) = v(1)
  
  p(4) = p(2)
  rho(4) = GetRho_RankineHugoniot(rho(1), p(1), p(4))
  u(4) = u(1)
  v(4) = GetPSI(rho(4), rho(1), p(4), p(1)) + v(1)

  p(3) = p(2)
  u(3) = u(1)
  v(3) = v(1)

CASE(13) ! Schulz-Rinne: Cfg J
  RiemannBC_WaveType(1) = DISCONTINUITY
  RiemannBC_WaveType(2) = SHOCK
  RiemannBC_WaveType(3) = DISCONTINUITY
  RiemannBC_WaveType(4) = SHOCK
  p(1) = 1.
  p(3) = 0.4
  rho(1) = 1.
  rho(2) = 2.
  u(1) =  0.
  v(1) = -0.3
  v(2) =  0.3

  p(2) = p(1)
  u(2) = u(1)

  rho(3) = GetRho_RankineHugoniot(rho(2), p(2), p(3))
  u(3) = u(1)
  v(3) = GetPSI(rho(3), rho(2), p(3), p(2))+v(2)

  p(4) = p(3)
  rho(4) = GetRho_RankineHugoniot(rho(1), p(1), p(4))
  u(4) = u(1)
  v(4) = GetPSI(rho(4), rho(1), p(4), p(1))+v(1)

!CASE(14) ! Schulz-Rinne: like Cfg J, but different initial conditions?



CASE(15) ! Schulz-Rinne: Cfg G
  RiemannBC_WaveType(1) = RAREFACTION
  RiemannBC_WaveType(2) = DISCONTINUITY
  RiemannBC_WaveType(3) = DISCONTINUITY
  RiemannBC_WaveType(4) = SHOCK
  p(1) = 1.
  p(2) = 0.4
  rho(1) = 1.
  rho(3) = 0.8
  u(1) = 0.1
  v(1) = -0.3

  p(4) = p(2)
  rho(4) = GetRho_RankineHugoniot(rho(1), p(1), p(4))
  u(4) = u(1)
  v(4) = GetPSI(rho(4), rho(1), p(4), p(1))+v(1)

  rho(2) = GetRho_Rarefaction(rho(1), p(1), p(2))
  u(2)=GetPHI(rho(2), rho(1), p(2), p(1)) + u(1)
  v(2) = v(1)
 
  p(3) = p(2)
  u(3) = u(1)
  v(3) = v(1)
 
CASE(16) ! Schulz-Rinne: Cfg H
  RiemannBC_WaveType(1) = RAREFACTION
  RiemannBC_WaveType(2) = DISCONTINUITY
  RiemannBC_WaveType(3) = DISCONTINUITY
  RiemannBC_WaveType(4) = SHOCK
  p(1) = 0.4
  p(2) = 1.0
  rho(3) = 0.8
  rho(4) = 1.0
  u(1) = 0.1
  v(1) = 0.1

  p(3) = p(2)
  p(4) = p(2)

  rho(1) = GetRho_RankineHugoniot(rho(4), p(4), p(1))
  rho(2) = (p(2)/p(1))**(1./Kappa) * rho(1)

  u(2) = u(1) - GetPHI(rho(2),rho(1),p(2),p(1))
  u(3) = u(1)
  u(4) = u(1)

  v(2) = v(1)
  v(3) = v(1)
  v(4) = v(1) + GetPSI(rho(4),rho(1),p(4),p(1))

CASE(17) ! Schulz-Rinne: Cfg K (but with different v(1) = 0.3)
  RiemannBC_WaveType(1) = DISCONTINUITY
  RiemannBC_WaveType(2) = SHOCK
  RiemannBC_WaveType(3) = DISCONTINUITY
  RiemannBC_WaveType(4) = RAREFACTION
  p(1) = 1.
  p(3) = 0.4
  rho(1) = 1.
  rho(2) = 2.
  u(1) = 0.
  v(1) = -0.4
  v(2) = -0.3

  p(2) = p(1)
  u(2) = u(1)

  rho(3) = GetRho_RankineHugoniot(rho(2), p(2), p(3))
  u(3) = u(1)
  v(3) = GetPSI(rho(3),rho(2),p(3),p(2)) + v(2)

  p(4) = p(3)
  rho(4) = GetRho_Rarefaction(rho(1), p(1), p(4))
  u(4) = u(1)
  v(4) = GetPHI(rho(4),rho(1),p(4),p(1)) + v(1)

CASE(18) ! Schulz-Rinne: Cfg K (but with different v(1) = 0.3)
  RiemannBC_WaveType(1) = DISCONTINUITY
  RiemannBC_WaveType(2) = SHOCK
  RiemannBC_WaveType(3) = DISCONTINUITY
  RiemannBC_WaveType(4) = RAREFACTION
  p(1) = 1.
  p(3) = 0.4
  rho(1) = 1.
  rho(2) = 2.
  u(1) = 0.
  v(1) = 1.0
  v(2) = -0.3

  p(2) = p(1)
  u(2) = u(1)

  rho(3) = GetRho_RankineHugoniot(rho(2), p(2), p(3))
  u(3) = u(1)
  v(3) = GetPSI(rho(3),rho(2),p(3),p(2)) + v(2)

  p(4) = p(3)
  rho(4) = GetRho_Rarefaction(rho(1), p(1), p(4))
  u(4) = u(1)
  v(4) = GetPHI(rho(4),rho(1),p(4),p(1)) + v(1)

CASE(19) ! Schulz-Rinne: Cfg K 
  RiemannBC_WaveType(1) = DISCONTINUITY
  RiemannBC_WaveType(2) = SHOCK
  RiemannBC_WaveType(3) = DISCONTINUITY
  RiemannBC_WaveType(4) = RAREFACTION
  p(1) = 1.
  p(3) = 0.4
  rho(1) = 1.
  rho(2) = 2.
  u(1) = 0.
  v(1) = 0.3
  v(2) = -0.3

  p(2) = p(1)
  u(2) = u(1)

  rho(3) = GetRho_RankineHugoniot(rho(2), p(2), p(3))
  u(3) = u(1)
  v(3) = GetPSI(rho(3),rho(2),p(3),p(2)) + v(2)

  p(4) = p(3)
  rho(4) = GetRho_Rarefaction(rho(1), p(1), p(4))
  u(4) = u(1)
  v(4) = GetPHI(rho(4),rho(1),p(4),p(1)) + v(1)
CASE  DEFAULT 
  CALL Abort(__STAMP__, &
      "Riemann2D problem for this IniExactFunc not implemented.")
END SELECT

nRefState = 4
SDEALLOCATE(RefStatePrim)
ALLOCATE(RefStatePrim(4,PP_nVarPrim))
RefStatePrim(:,1) = rho
RefStatePrim(:,2) = u
RefStatePrim(:,3) = v
RefStatePrim(:,4) = 0.
RefStatePrim(:,5) = p

SDEALLOCATE(RefStateCons)
ALLOCATE(RefStateCons(4,PP_nVar))
DO i=1,4
  ! TODO: ATTENTION only sRho and Pressure of UE filled!!!
  UE(SRHO) = 1./RefStatePrim(i,1)
  UE(PRES) = RefStatePrim(i,5)
  RefStatePrim(i,6) = TEMPERATURE_HE(UE)
  CALL PrimToCons(RefStatePrim(i,:), RefStateCons(i,:))
END DO 
END SUBROUTINE CalcIniStates

!==================================================================================================================================
!> Specifies all the initial conditions.
!==================================================================================================================================
SUBROUTINE ExactFuncTestcase(tIn,x,Resu,Resu_t,Resu_tt)
! MODULES
USE MOD_Globals       ,ONLY: Abort
USE MOD_Equation_Vars ,ONLY: RefStateCons
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: x(3)        !< position in physical coordinates
REAL,INTENT(IN)                 :: tIn         !< current simulation time
REAL,INTENT(OUT)                :: Resu(5)     !< exact fuction evaluated at tIn, returning state in conservative variables
REAL,INTENT(OUT)                :: Resu_t(5)   !< first time deriv of exact fuction
REAL,INTENT(OUT)                :: Resu_tt(5)  !< second time deriv of exact fuction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
!(see: Solution of Two-Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers
!      Alexander Kurganov, Eitan Tadmor)
IF ((x(1).GE.0.0).AND.(x(2).GE.0.0)) THEN
    Resu = RefStateCons(1,:)
ELSE IF ((x(1).LT.0.0).AND.(x(2).GE.0.0)) THEN
    Resu = RefStateCons(2,:)
ELSE IF ((x(1).LT.0.0).AND.(x(2).LT.0.0)) THEN
    Resu = RefStateCons(3,:)
ELSE IF ((x(1).GE.0.0).AND.(x(2).LT.0.0)) THEN
    Resu = RefStateCons(4,:)
END IF

Resu_t =0.
Resu_tt=0.

END SUBROUTINE ExactFuncTestcase

SUBROUTINE GetBoundaryFluxTestcase(SideID,t,Nloc,Flux,UPrim_master,                   &
#if PARABOLIC
                           gradUx_master,gradUy_master,gradUz_master,&
#endif
                           NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Testcase_Vars
USE MOD_Mesh_Vars     ,ONLY: nSides
USE MOD_Mesh_Vars     ,ONLY: BoundaryType, BC
USE MOD_FV_Vars       ,ONLY: FV_Elems_master
USE MOD_Riemann       ,ONLY: Riemann
USE MOD_EOS           ,ONLY: PrimToCons
USE MOD_Equation_Vars ,ONLY: RefStatePrim,RefStateCons
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: SideID  
REAL,INTENT(IN)      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)   :: Nloc    !< polynomial degree
REAL,INTENT(IN)      :: UPrim_master( PP_nVarPrim,0:Nloc,0:Nloc,1:nSides) !< inner surface solution
#if PARABOLIC
                                                           !> inner surface solution gradients in x/y/z-direction
REAL,INTENT(IN)      :: gradUx_master(PP_nVarPrim,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(IN)      :: gradUy_master(PP_nVarPrim,0:Nloc,0:Nloc,1:nSides)
REAL,INTENT(IN)      :: gradUz_master(PP_nVarPrim,0:Nloc,0:Nloc,1:nSides)
#endif /*PARABOLIC*/
                                                           !> normal and tangential vectors on surfaces
REAL,INTENT(IN)      :: NormVec (3,0:Nloc,0:Nloc,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)      :: TangVec1(3,0:Nloc,0:Nloc,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)      :: TangVec2(3,0:Nloc,0:Nloc,0:FV_ENABLED,1:nSides)
REAL,INTENT(IN)      :: Face_xGP(3,0:Nloc,0:Nloc,0:FV_ENABLED,1:nSides) !< positions of surface flux points
REAL,INTENT(OUT)     :: Flux(PP_nVar,0:Nloc,0:Nloc,1:nSides)  !< resulting boundary fluxes
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: BCState, dir, p, q, iL, iR
REAL              :: pos(2), x, rel
REAL              :: UCons_boundary(PP_nVar    ,0:PP_N,0:PP_N)
REAL              :: UPrim_boundary(PP_nVarPrim,0:PP_N,0:PP_N)
REAL              :: UCons_master  (PP_nVar    ,0:Nloc,0:Nloc)
!==================================================================================================================================
BCState = Boundarytype(BC(SideID),BC_STATE)
pos = RiemannBC_Speeds(BCState,:) * t 
dir = MOD(BCState+1,2)+1
iL = SideToQuads(1,BCState)
iR = SideToQuads(2,BCState)
DO q=0,PP_N; DO p=0,PP_N
  x = Face_xGP(dir,p,q,FV_Elems_master(SideID),SideID)
  SELECT CASE (RiemannBC_WaveType(BCState))
  CASE(SHOCK,DISCONTINUITY)
    IF (x.LT.pos(1)) THEN
      UCons_boundary(:,p,q)=RefStateCons(iL,:)
      UPrim_boundary(:,p,q)=RefStatePrim(iL,:)
    ELSE 
      UCons_boundary(:,p,q)=RefStateCons(iR,:)
      UPrim_boundary(:,p,q)=RefStatePrim(iR,:)
    END IF
  CASE(RAREFACTION)

    IF (x.LE.pos(1)) THEN
      UCons_boundary(:,p,q)=RefStateCons(iL,:)
      UPrim_boundary(:,p,q)=RefStatePrim(iL,:)
    ELSE IF( x.GE.pos(2)) THEN
      UCons_boundary(:,p,q)=RefStateCons(iR,:)
      UPrim_boundary(:,p,q)=RefStatePrim(iR,:)
    ELSE
      rel = (x-pos(1))/(pos(2)-pos(1))
      UCons_boundary(:,p,q)= (1-rel)*RefStateCons(iL,:) + rel*RefStateCons(iR,:)
      UPrim_boundary(:,p,q)= (1-rel)*RefStatePrim(iL,:) + rel*RefStatePrim(iR,:)
    END IF
  END SELECT
END DO; END DO
DO q=0,PP_N; DO p=0,PP_N
  CALL PrimToCons(UPrim_master(:,p,q,SideID), UCons_master(:,p,q)) 
END DO; END DO ! p,q=0,PP_N

!WRITE (*,*) Face_xGP(:,0,0,FV_Elems_master(SideID),SideID)
!WRITE (*,*) x, pos
!WRITE (*,*) UPrim_master (:,0,0,SideID)
!WRITE (*,*) UPrim_boundary(:,0,0)
!read*
CALL Riemann(Nloc,Flux(:,:,:,SideID),UCons_master,UCons_boundary,UPrim_master(:,:,:,SideID),UPrim_boundary,&
    NormVec (:,:,:,FV_Elems_master(SideID),SideID),&
    TangVec1(:,:,:,FV_Elems_master(SideID),SideID),&
    TangVec2(:,:,:,FV_Elems_master(SideID),SideID),&
    doBC=.TRUE.)

END SUBROUTINE GetBoundaryFluxTestcase


SUBROUTINE GetBoundaryFVgradientTestcase(SideID,t,gradU,UPrim_master)
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: nSides
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(OUT)   :: gradU       (PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
!==================================================================================================================================
gradU(:,:,:,SideID) = 0.
END SUBROUTINE GetBoundaryFVgradientTestcase


SUBROUTINE Lifting_GetBoundaryFluxTestcase(SideID,t,UPrim_master,Flux)
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: nSides
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID
REAL,INTENT(IN)    :: t                                    !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< primitive solution from the inside
REAL,INTENT(OUT)   :: Flux(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides) !< lifting boundary flux
!==================================================================================================================================
END SUBROUTINE Lifting_GetBoundaryFluxTestcase

SUBROUTINE Riemann_Speeds(dir, wavetype, prim_L, prim_R, speeds)
!=================================================================================================================================
!=================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_EOS_Vars, ONLY: Kappa,kappaM1,KappaP1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: dir
INTEGER,INTENT(IN)                 :: wavetype
REAL,DIMENSION(PP_nVar),INTENT(IN) :: prim_L, prim_R
REAL,DIMENSION(1:2),INTENT(OUT)    :: speeds
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                                  :: r1,r2,r3,r4
REAL                                  :: p1,p2,p3,p4
REAL                                  :: u1,u2,u3,u4
REAL                                  :: a1,a3,a4
REAL                                  :: t2t1,r2r1,p2p1,wsp
!=================================================================================================================================

IF (wavetype == DISCONTINUITY) THEN 
  speeds = (prim_L(1)*prim_L(dir+1) - prim_R(1)*prim_R(dir+1))/(prim_L(1)-prim_R(1))
  RETURN
END IF

IF (prim_L(5).LT.prim_R(5)) THEN
  r1 =  prim_L(1)
  u1 = -prim_L(dir+1)
  p1 =  prim_L(5) 
  r4 =  prim_R(1)
  u4 = -prim_R(dir+1)
  p4 =  prim_R(5)
ELSE 
  r1 =  prim_R(1)
  u1 =  prim_R(dir+1)
  p1 =  prim_R(5) 
  r4 =  prim_L(1)
  u4 =  prim_L(dir+1)
  p4 =  prim_L(5)
END IF

a1 = sqrt( Kappa * p1 / r1 )
a4 = sqrt( Kappa * p4 / r4 )
call exact_riemann(Kappa, r4,r1,r2, u4,u1,u2, p4,p1,p2, a4,a1)

p2p1 = p2/p1
t2t1 = p2p1 * ( KappaP1/KappaM1 + p2p1 ) / ( 1.0 + KappaP1 * p2p1 / KappaM1 )
r2r1 = ( 1.0 + KappaP1 * p2p1 / KappaM1 ) / ( KappaP1 / KappaM1 + p2p1 )

! shock-wave speed.
wsp  = u1 + a1 * sqrt( KappaP1 * ( p2p1 - 1.0 ) / ( 2.0 * Kappa ) + 1.0 )

IF (wavetype == SHOCK) THEN
  IF (prim_L(5).LT.prim_R(5)) THEN
    speeds = -wsp
  ELSE
    speeds = wsp
  END IF
  RETURN
END IF

! State 3.
p3 = p2

! Isentropic between 3 and 4.
r3 = r4 * ( p3 / p4 )**(1.0/Kappa)
a3 = sqrt( Kappa * p3 / r3 )

! Speed of contact discontinuity.
u3 = u2

IF (wavetype == RAREFACTION) THEN
  IF (prim_L(5).LT.prim_R(5)) THEN
    speeds(1) = -(u3-a3)
    speeds(2) = -(u4-a4)
  ELSE
    speeds(1) = u3-a3
    speeds(2) = u4-a4
  END IF
  speeds = (/ MINVAL(speeds), MAXVAL(speeds) /)
END IF

END SUBROUTINE Riemann_Speeds

SUBROUTINE exact_riemann(gamma,rhol,rhor,rho,ul,ur,u,pl,pr,p,al,ar)
!-------------------------------------------------------
IMPLICIT NONE
!-------------------------------------------------------
REAL                                  :: rhol,rhor,ul,ur,pl,pr,al,ar,&
    du, ducrit, p, p0, fr,             &
    fl, frd, fld, cha, rho, u, pm, um, &
    gamma, s
REAL                                  :: G(1:9),tol
INTEGER                               :: KK,nIter
!-------------------------------------------------------
INTENT(IN)                                          :: rhol,rhor,ul,ur,pl,pr,al,ar
INTENT(OUT)                                         :: rho,u,p
!-------------------------------------------------------
s = 0.0
G     = (/ (gamma-1.0)/(2.0*gamma), &
    (gamma+1.0)/(2.0*gamma), &
    2.0*gamma/(gamma-1.0),   &
    2.0/(gamma-1.0),         &
    2.0/(gamma+1.0),         &
    (gamma-1.0)/(gamma+1.0), &
    0.5*(gamma-1.0),         &
    1.0/gamma,               &
    gamma-1.0                /)

du = ur-ul
ducrit = G(4)*(al+ar)-du

IF (ducrit .LE. 0.0) THEN
  WRITE(*,*)'Vacuum is generated by given data!!'
  WRITE(*,*)'rhol=',rhol
  WRITE(*,*)'rhor=',rhor
  WRITE(*,*)'ul=',ul
  WRITE(*,*)'ur=',ur
  WRITE(*,*)'pl=',pl
  WRITE(*,*)'pr=',pr
  WRITE(*,*)'al=',al
  WRITE(*,*)'ar=',ar
  WRITE(*,*)'du=',du
  WRITE(*,*)'G(4), al+ar',G(4), al+ar
  STOP
END IF

tol = 1.0e-10

CALL STARTE(G,tol, p,rhol,rhor,ul,ur,pl,pr,al,ar)

p0  = p
cha = 2.0*tol
KK  = 0

nIter = 1000
DO WHILE ((cha.GT.tol).AND.(KK.LT.nIter))

  KK = KK+1

  CALL PREFUN(G,fl,fld,p,rhol,pl,al)
  CALL PREFUN(G,fr,frd,p,rhor,pr,ar)

  p   = p - (fl+fr+du)/(fld+frd)
  cha = 2.0 * ABS((p-p0)/(p+p0))

  IF (p .LT. 0.0) THEN
    p = tol
  END IF

  p0=p

END DO

IF (KK.GE. nIter) THEN
  WRITE(*,*)'WARNING: Divergence in Newton-Raphson scheme'
END IF

u = 0.5*(ul+ur+fr-fl)

pm = p
um = u

CALL SAMPLE(G,s,p,u,rho,rhol,rhor,ul,ur,um,pl,pr,pm,al,ar)

END SUBROUTINE exact_riemann


SUBROUTINE STARTE(G,tol,p,rhol,rhor,ul,ur,pl,pr,al,ar)
!-------------------------------------------------------
IMPLICIT NONE
!-------------------------------------------------------
REAL                                  :: G(1:9),tol,rhol,rhor,ul,ur,pl,pr,&
    al, ar, p, pv, pmin,pmax,   &
    pnu, pde, qmax,         &
    qrat, gel, ger
!-------------------------------------------------------
INTENT(IN)                                          :: G,tol,rhol,rhor,ul,ur,pl,pr,al,ar
INTENT(OUT)                                         :: p
!-------------------------------------------------------

qmax = 2.0

pv   = 0.5*(pl+pr)-0.125*(ur-ul)*(rhol+rhor)*(al+ar)
pmin = MIN(pl,pr)
pmax = MAX(pl,pr)
qrat = pmax/pmin

IF ((qrat.LE.qmax).AND.((pmin.LE.pv).AND.(pv.LE.pmax))) THEN
  p = MAX(tol,pv)
ELSE
  IF (pv.LT.pmin) THEN
    pnu = al+ar-G(7)*(ur-ul)
    pde = al/pl**G(1) + ar/pr**G(1)
    p   = (pnu/pde)**G(3)
  ELSE
    gel = SQRT((G(5)/rhol)/(G(6)*pl+MAX(tol,pv)))
    ger = SQRT((G(5)/rhor)/(G(6)*pr+MAX(tol,pv)))
    p   = (gel*pl+ger*pr-(ur-ul))/(gel+ger)
    p   = MAX(tol,p)
  END IF
END IF
END SUBROUTINE STARTE


SUBROUTINE PREFUN(G,f,fd,p,rhok,pk,ck)
!-------------------------------------------------------
IMPLICIT NONE
!-------------------------------------------------------
REAL                                  :: G(1:9),f,fd,p,rhok,pk,ck,prat,ak,qrt,bk
!-------------------------------------------------------
INTENT(IN)                                          :: p,rhok,pk,ck
INTENT(OUT)                                         :: f,fd
!-------------------------------------------------------    
IF (p.LE.pk) THEN
  prat = p/pk
  f    = G(4)*ck*(prat**G(1) - 1.0)
  fd   = (1.0/(rhok*ck))*prat**(-G(2))
ELSE
  ak  = G(5)/rhok
  bk  = G(6)*pk
  qrt = SQRT(ak/(bk+p))
  f   = (p-pk) * qrt
  fd  = (1.0-0.5*(p-pk)/(bk+p))*qrt
END IF
END SUBROUTINE PREFUN

SUBROUTINE SAMPLE(G,s,p,u,rho,rhol,rhor,ul,ur,um,pl,pr,pm,al,ar)
!-------------------------------------------------------
IMPLICIT NONE
!-------------------------------------------------------
REAL                                  :: G(1:9),s,p,u,rho,rhol,rhor,ul,ur,um,pl,pr,&
    pm,al,ar,pml,pmr,str,stl,cmr,      &
    cml,shr,shl,sr,sl,c
!-------------------------------------------------------
INTENT(IN)                                          :: G,s,rhol,rhor,ul,ur,um,pl,pr,pm,al,ar
INTENT(OUT)                                         :: p,u,rho
!-------------------------------------------------------    
IF(s.LE.um) THEN
  IF(pm.LE.pl) THEN
    shl=ul-al
    IF(s.LE.shl) THEN
      rho=rhol
      u=ul
      p=pl
    ELSE
      cml = al*(pm/pl)**G(1)
      stl = um-cml
      IF(s.GT.stl) THEN
        rho=rhol*(pm/pl)**G(8)
        u=um
        p=pm
      ELSE
        u=G(5)*(al+G(7)*ul+s)
        c=G(5)*(al+G(7)*(ul-s))
        rho=rhol*(c/al)**G(4)
        p=pl*(c/al)**G(3)
      END IF
    END IF
  ELSE
    pml = pm/pl
    sl  = ul-al*SQRT(G(2)*pml+G(1))
    IF(s.LE.sl) THEN
      rho=rhol
      u=ul
      p=pl
    ELSE
      rho=rhol*(pml+G(6))/(pml*G(6)+1.0)
      u=um
      p=pm
    END IF
  END IF
ELSE
  IF (pm.GT.pr) THEN
    pmr=pm/pr
    sr =ur+ar*SQRT(G(2)*pmr+G(1))
    IF (s.GE.sr) THEN
      rho=rhor
      u=ur
      p=pr
    ELSE
      rho=rhor*(pmr+G(6))/(pmr*G(6)+1.0)
      u=um
      p=pm
    END IF
  ELSE
    shr=ur+ar
    IF (s.GE.SHR) THEN
      rho=rhor
      u=ur
      p=pr
    ELSE
      cmr=ar*(pm/pr)**G(1)
      str=um+cmr
      IF(s.LE.str) THEN
        rho=rhor*(pm/pr)**G(8)
        u=um
        p=pm
      ELSE
        u=G(5)*(-ar+G(7)*ur+s)
        c=G(5)*(ar-G(7)*(ur-s))
        rho=rhor*(c/ar)**G(4)
        p=pr*(c/ar)**G(3)
      END IF
    END IF
  END IF
END IF
END SUBROUTINE SAMPLE

!==================================================================================================================================
!> Add testcases source term to solution time derivative
!==================================================================================================================================
SUBROUTINE TestcaseSource(Ut)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(*),INTENT(IN) :: Ut                        !< solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE TestcaseSource

END MODULE MOD_Testcase

