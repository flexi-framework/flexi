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
!> Contains global variables used by the Timedisc modules.
!==================================================================================================================================
MODULE MOD_TimeDisc_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

INTERFACE SetTimeDiscCoefs
  MODULE PROCEDURE SetTimeDiscCoefs
END INTERFACE

! > Dummy interface for time step function pointer
ABSTRACT INTERFACE
  SUBROUTINE TimeIntegrator(t)
    REAL,INTENT(IN) :: t
  END SUBROUTINE
END INTERFACE


!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL             :: t=0.                      !< current physical time
REAL             :: dt                        !< current timestep
REAL             :: TEnd                      !< End time of simulation
REAL             :: TAnalyze                  !< Analyze time intervall
REAL             :: CFLScale(0:FV_ENABLED)    !< Convective CFL number
REAL             :: DFLScale(0:FV_ENABLED)    !< Viscous CFL number (only if PARABOLIC)
REAL,ALLOCATABLE :: dtElem(:)                 !< Timestep for each element
INTEGER          :: CurrentStage=1            !< Current Runge-Kutta stage within timestep
INTEGER          :: nCalcTimeStepMax          !< Compute dt at least after every Nth timestep
INTEGER(KIND=8)  :: maxIter                   !< Maximum permitted number of timesteps
LOGICAL          :: fullBoundaryOrder=.FALSE. !< temporal order degradation, occuring for
                                              !< time-dependant BCs, can easily be fixed when
                                              !< using 3 stage 3rd order RK schemes (no others!)
LOGICAL          :: ViscousTimeStep=.FALSE.   !< Info wether we have convection of viscous dominated timestep
LOGICAL          :: TimediscInitIsDone=.FALSE.!< Indicate wheter InitTimeDisc routine has been run

!----------------------------------------------------------------------------------------------------------------------------------
! TIME INTEGRATION: RUNGE_KUTTA COEFFICIENTS AND STABILITY NUMBERS
!----------------------------------------------------------------------------------------------------------------------------------

! NOTE: using simple arrays to store coefs, if using classes or types performance seems to degrade
CHARACTER(LEN=50)   :: TimeDiscName  !< name of specific time discretization scheme
CHARACTER(LEN=50)   :: TimeDiscType  !< general type of time discretization scheme
INTEGER             :: nRKStages     !< number of stages of Runge-Kutta method
PROCEDURE(TimeIntegrator),POINTER :: TimeStep !< pointer to timestepping routine, depends on td

!> Runge-Kutta low storage coefficients for Williamson 2 register and Ketcheson 3 register schemes
REAL,ALLOCATABLE    :: RKA(:),RKb(:),RKc(:),RKg1(:),RKg2(:),RKg3(:),RKdelta(:)

REAL                :: CFLScaleAlpha(1:15)  !< timestep scaling factor for CFL number, depends on scheme
#if PARABOLIC
!> DFL in DG depends on the polynomial degree
!> since DFl is only on real axis, stability numbers are defined for RK3 and then scaled for RK4
REAL,PARAMETER      :: DFLScaleAlpha(1:10) = &
#if (PP_NodeType==1)
(/ 1.12, 0.76, 0.55, 0.41, 0.32, 0.25, 0.20, 0.17, 0.14, 0.12 /)
#elif (PP_NodeType==2)
(/ 4.50, 1.95, 1.18, 0.79, 0.56, 0.42, 0.32, 0.25, 0.20, 0.17 /)
#endif /*PP_NodeType*/
REAL                :: RelativeDFL          !< scaling factor for DFL defined by scheme
#endif /*PARABOLIC*/


CONTAINS

!==================================================================================================================================
!> Function returning the coefficients of the timestepping scheme (e.g. Runge-Kutta)
!==================================================================================================================================
SUBROUTINE SetTimeDiscCoefs(TimeDiscMethod)
! MODULES
USE MOD_Globals             ,ONLY:CollectiveStop
#if (PP_NodeType==2)
USE MOD_Overintegration_Vars,ONLY:OverintegrationType
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: TimeDiscMethod !< name of time discretization to be used
!==================================================================================================================================

! Set the type of the Runge Kutta scheme:
! We currently have Williamson 2 stage and Ketcheson 3 stage schemes
SELECT CASE (TRIM(TimeDiscMethod))
CASE('standardrk3-3','carpenterrk4-5','niegemannrk4-14',&
     'toulorgerk4-8c','toulorgerk3-7c','toulorgerk4-8f')
  TimeDiscType='LSERKW2'
CASE('ketchesonrk4-20','ketchesonrk4-18')
  TimeDiscType='LSERKK3'
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
                      'Unknown method of time discretization: '//TRIM(TimeDiscMethod))
END SELECT


SELECT CASE (TimeDiscMethod)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Low-storage Runge-Kutta 3, 3 stages, Kopriva,Algorithm 42
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CASE('standardrk3-3')

  TimeDiscName = 'Standard RK3-3'
  nRKStages=3
#if PARABOLIC
  RelativeDFL=1.
#endif
  fullBoundaryOrder=.TRUE.

#if (PP_NodeType==1)
  CFLScaleAlpha(1:15) = &
  (/ 1.2285, 1.0485, 0.9101, 0.8066, 0.7268, 0.6626, 0.6109, 0.5670, 0.5299, 0.4973, 0.4703, 0.4455, 0.4230, 0.4039, 0.3859 /)
#elif (PP_NodeType==2)
  IF (OverintegrationType.GT.0) THEN
    ! Overintegration with Gauss-Lobatto nodes results in a projection DG formulation, i.e. we have to use the Gauss nodes timestep
    CFLScaleAlpha(1:15) = &
    (/ 1.2285, 1.0485, 0.9101, 0.8066, 0.7268, 0.6626, 0.6109, 0.5670, 0.5299, 0.4973, 0.4703, 0.4455, 0.4230, 0.4039, 0.3859 /)
  ELSE
    CFLScaleAlpha(1:15) = &
    (/ 3.1871, 2.2444, 1.7797, 1.5075, 1.3230, 1.1857, 1.0800, 0.9945, 0.9247, 0.8651, 0.8134, 0.7695, 0.7301, 0.6952, 0.6649 /)
  ENDIF
#endif /*PP_NodeType*/

  ALLOCATE(RKA(2:nRKStages),RKb(1:nRKStages),RKc(2:nRKStages))
  RKA(2:nRKStages) = (/ 5./9.,153./128. /)
  RKb(1:nRKStages) = (/ 1./3., 15./16., 8./15. /)
  RKc(2:nRKStages) = (/ 1./3., 0.75 /)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Runge-Kutta 4 - Carpenter 1994 NASA Report
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CASE('carpenterrk4-5')

  TimeDiscName = 'Carpenter RK4-5'
  nRKStages=5
#if (PP_NodeType==1)
  CFLScaleAlpha(1:15) = &
  (/ 2.0351, 1.7595, 1.5401, 1.3702, 1.2375, 1.1318, 1.0440, 0.9709, 0.9079, 0.8539, 0.8066, 0.7650, 0.7290, 0.6952, 0.6660 /)
#elif (PP_NodeType==2)
  IF (OverintegrationType.GT.0) THEN
    ! Overintegration with Gauss-Lobatto nodes results in a projection DG formulation, i.e. we have to use the Gauss nodes timestep
    CFLScaleAlpha(1:15) = &
    (/ 2.0351, 1.7595, 1.5401, 1.3702, 1.2375, 1.1318, 1.0440, 0.9709, 0.9079, 0.8539, 0.8066, 0.7650, 0.7290, 0.6952, 0.6660 /)
  ELSE
    CFLScaleAlpha(1:15) = &
    (/ 4.7497, 3.4144, 2.8451, 2.4739, 2.2027, 1.9912, 1.8225, 1.6830, 1.5682, 1.4692, 1.3849, 1.3106, 1.2454, 1.1880, 1.1362 /)
  END IF
#endif /*PP_NodeType*/
#if PARABOLIC
  RelativeDFL=1.853
#endif
  ALLOCATE(RKA(2:nRKStages),RKb(1:nRKStages),RKc(2:nRKStages))
  RKA(2:nRKStages) = (/  567301805773.0/  1357537059087.0,&
                        2404267990393.0/  2016746695238.0,&
                        3550918686646.0/  2091501179385.0,&
                        1275806237668.0/   842570457699.0 /)

  RKb(1:nRKStages) = (/ 1432997174477.0/  9575080441755.0,&
                        5161836677717.0/ 13612068292357.0,&
                        1720146321549.0/  2090206949498.0,&
                        3134564353537.0/  4481467310338.0,&
                        2277821191437.0/ 14882151754819.0 /)

  RKc(2:nRKStages) = (/ 1432997174477.0/  9575080441755.0,&
                        2526269341429.0/  6820363962896.0,&
                        2006345519317.0/  3224310063776.0,&
                        2802321613138.0/  2924317926251.0 /)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Low storage Runge-Kutta 4, 14 stages version - Niegemann et al 2012
! Fastest RK4 scheme implemented, but less accurate then Carpenter RK4
! Very good performance for high N
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CASE('niegemannrk4-14')

  TimeDiscName = 'Niegemann RK4-14'
  nRKStages=14
#if (PP_NodeType==1)
  CFLScaleAlpha(1:15) = &
  (/ 6.9716, 5.7724, 5.1863, 4.7880, 4.4741, 4.2120, 3.9836, 3.7811, 3.5617, 3.3705, 3.1995, 3.0488, 2.9137, 2.7900, 2.6786 /)
#elif (PP_NodeType==2)
  IF (OverintegrationType.GT.0) THEN
    ! Overintegration with Gauss-Lobatto nodes results in a projection DG formulation, i.e. we have to use the Gauss nodes timestep
    CFLScaleAlpha(1:15) = &
    (/ 6.9716, 5.7724, 5.1863, 4.7880, 4.4741, 4.2120, 3.9836, 3.7811, 3.5617, 3.3705, 3.1995, 3.0488, 2.9137, 2.7900, 2.6786 /)
  ELSE
    CFLScaleAlpha(1:15) = &
    (/ 14.7882, 9.5906, 7.9447, 7.0965, 6.5486, 6.1436, 5.8185, 5.5440, 5.3055, 5.0940, 4.9028, 4.7295, 4.5697, 4.4235, 4.2885 /)
  END IF
#endif /*PP_NodeType*/
#if PARABOLIC
  RelativeDFL=7.379
#endif
  ALLOCATE(RKA(2:nRKStages),RKb(1:nRKStages),RKc(2:nRKStages))
  RKA = (/ 0.7188012108672410,&
           0.7785331173421570,&
           0.0053282796654044,&
           0.8552979934029281,&
           3.9564138245774565,&
           1.5780575380587385,&
           2.0837094552574054,&
           0.7483334182761610,&
           0.7032861106563359,&
          -0.0013917096117681,&
           0.0932075369637460,&
           0.9514200470875948,&
           7.1151571693922548 /)

  RKb = (/ 0.0367762454319673,&
           0.3136296607553959,&
           0.1531848691869027,&
           0.0030097086818182,&
           0.3326293790646110,&
           0.2440251405350864,&
           0.3718879239592277,&
           0.6204126221582444,&
           0.1524043173028741,&
           0.0760894927419266,&
           0.0077604214040978,&
           0.0024647284755382,&
           0.0780348340049386,&
           5.5059777270269628 /)

  RKc = (/ 0.0367762454319673,&
           0.1249685262725025,&
           0.2446177702277698,&
           0.2476149531070420,&
           0.2969311120382472,&
           0.3978149645802642,&
           0.5270854589440328,&
           0.6981269994175695,&
           0.8190890835352128,&
           0.8527059887098624,&
           0.8604711817462826,&
           0.8627060376969976,&
           0.8734213127600976 /)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CASE('toulorgerk4-8c')
! Low storage Runge-Kutta 4, 8 stages version - Toulorge et al 2012 - RKC84 version
! Slightly slower than Niegemann 14 scheme (performance gap grows with N), but error is nearly on par with Carpenter RK4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TimeDiscName = 'Toulorge RK4-8C'
  nRKStages=8
#if (PP_NodeType==1)
  CFLScaleAlpha(1:15) = &
    (/ 3.7726, 3.2306, 2.8926, 2.6346, 2.3946, 2.1886, 2.0196, 1.8776, 1.7576, 1.6536, 1.5626, 1.4816, 1.4106, 1.3466, 1.2886 /)
#elif (PP_NodeType==2)
  IF (OverintegrationType.GT.0 ) THEN
    ! Overintegration with Gauss-Lobatto nodes results in a projection DG formulation, i.e. we have to use the Gauss nodes timestep
    CFLScaleAlpha(1:15) = &
    (/ 3.7726, 3.2306, 2.8926, 2.6346, 2.3946, 2.1886, 2.0196, 1.8776, 1.7576, 1.6536, 1.5626, 1.4816, 1.4106, 1.3466, 1.2886 /)
  ELSE
    CFLScaleAlpha(1:15) = &
    (/ 8.2506, 5.7006, 4.7876, 4.2616, 3.8866, 3.5946, 3.3536, 3.1516, 2.9766, 2.8246, 2.6786, 2.5356, 2.4096, 2.2976, 2.1966 /)
  END IF
#endif /*PP_NodeType*/
#if PARABOLIC
  RelativeDFL=3.328
#endif
  ALLOCATE(RKA(2:nRKStages),RKb(1:nRKStages),RKc(2:nRKStages))
  RKA = (/0.7212962482279240,&
          0.0107733657161298,&
          0.5162584698930970,&
          1.7301002866322010,&
          5.2001293044030760,&
         -0.7837058945416420,&
          0.5445836094332190 /)

  RKb = (/0.2165936736758085,&
          0.1773950826411583,&
          0.0180253861162329,&
          0.0847347637254149,&
          0.8129106974622483,&
          1.9034160304227600,&
          0.1314841743399048,&
          0.2082583170674149 /)

  RKc = (/0.2165936736758085,&
          0.2660343487538170,&
          0.2840056122522720,&
          0.3251266843788570,&
          0.4555149599187530,&
          0.7713219317101170,&
          0.9199028964538660 /)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Low storage Runge-Kutta 3, 7 stages version - Toulorge et al 2012 - RKC73 version
! For Gauss: fastest scheme up to N~9. For GL: generally fastest scheme
! POOR ACCURACY (but still on par with standard RK3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CASE('toulorgerk3-7c')

  TimeDiscName = 'Toulorge RK3-7C'
  nRKStages=7
#if (PP_NodeType==1)
  CFLScaleAlpha(1:15) = &
  (/ 3.7416, 3.2036, 2.8686, 2.6136, 2.4066, 2.2116, 2.0416, 1.8986, 1.7756, 1.6676, 1.5746, 1.4916, 1.4186, 1.3536, 1.2946 /)
#elif (PP_NodeType==2)
  IF (OverintegrationType.GT.0) THEN
    ! Overintegration with Gauss-Lobatto nodes results in a projection DG formulation, i.e. we have to use the Gauss nodes timestep
    CFLScaleAlpha(1:15) = &
    (/ 3.7416, 3.2036, 2.8686, 2.6136, 2.4066, 2.2116, 2.0416, 1.8986, 1.7756, 1.6676, 1.5746, 1.4916, 1.4186, 1.3536, 1.2946 /)
  ELSE
    CFLScaleAlpha(1:15) = &
    (/ 8.1776, 5.6366, 4.7326, 4.2156, 3.8486, 3.5606, 3.3246, 3.1246, 2.9526, 2.8026, 2.6696, 2.5516, 2.4336, 2.3216, 2.2196 /)
  END IF
#endif /*PP_NodeType*/
#if PARABOLIC
  RelativeDFL=3.342
#endif
  ALLOCATE(RKA(2:nRKStages),RKb(1:nRKStages),RKc(2:nRKStages))
  RKA = (/0.808316387498383,&
          1.503407858773331,&
          1.053064525050744,&
          1.463149119280508,&
          0.659288128108783,&
          1.667891931891068 /)

  RKb = (/0.0119705267309784,&
          0.8886897793820711,&
          0.4578382089261419,&
          0.5790045253338471,&
          0.3160214638138484,&
          0.2483525368264122,&
          0.0677123095940884 /)

  RKc = (/0.0119705267309784,&
          0.182317794036199 ,&
          0.5082168062551849,&
          0.653203122014859 ,&
          0.853440138567825 ,&
          0.998046608462379  /)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Low storage Runge-Kutta 4, 8 stages version - Toulorge et al 2012 - RKF84 version
! Slightly slower than Niegemann 14 and RKC84 scheme (performance gap grows with N), but error is nearly on par with Carpenter RK4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CASE('toulorgerk4-8f')

  TimeDiscName = 'Toulorge RK4-8F'
  nRKStages=8
#if (PP_NodeType==1)
  CFLScaleAlpha(1:15) = &
  (/ 3.5856, 3.0966, 2.7566, 2.4856, 2.2666, 2.0756, 1.9116, 1.7746, 1.6586, 1.5586, 1.4716, 1.3946, 1.3266, 1.2656, 1.2106 /)
#elif (PP_NodeType==2)
  IF (OverintegrationType.GT.0) THEN
    ! Overintegration with Gauss-Lobatto nodes results in a projection DG formulation, i.e. we have to use the Gauss nodes timestep
    CFLScaleAlpha(1:15) = &
    (/ 3.5856, 3.0966, 2.7566, 2.4856, 2.2666, 2.0756, 1.9116, 1.7746, 1.6586, 1.5586, 1.4716, 1.3946, 1.3266, 1.2656, 1.2106 /)
  ELSE
    CFLScaleAlpha(1:15) = &
    (/ 7.9706, 5.5346, 4.6486, 4.1376, 3.7686, 3.4726, 3.2246, 3.0116, 2.8276, 2.6676, 2.5266, 2.4026, 2.2866, 2.1776, 2.0806 /)
  ENDIF 
#endif /*PP_NodeType*/
#if PARABOLIC
  RelativeDFL=3.142
#endif
  ALLOCATE(RKA(2:nRKStages),RKb(1:nRKStages),RKc(2:nRKStages))
  RKA = (/0.5534431294501569,&
         -0.0106598757020349,&
          0.5515812888932000,&
          1.885790377558741 ,&
          5.701295742793264 ,&
         -2.113903965664793 ,&
          0.5339578826675280 /)

  RKb = (/0.0803793688273695,&
          0.5388497458569843,&
          0.0197497440903196,&
          0.0991184129733997,&
          0.7466920411064123,&
          1.6795842456188940,&
          0.2433728067008188,&
          0.1422730459001373 /)

  RKc = (/0.0803793688273695,&
          0.3210064250338430,&
          0.3408501826604660,&
          0.3850364824285470,&
          0.5040052477534100,&
          0.6578977561168540,&
          0.9484087623348481 /)


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Here come the Ketcheson type 3 register schemes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Low storage Runge-Kutta 5, 20 stages version - Ketcheson et al 2012
! Error slightly better then Carpenter RK5 and by far better then Niegemann RK14
! Per stage slightly faster then Carpenter RK5, much slower then Niegemann RK14
! => Recommended as default RK scheme!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CASE('ketchesonrk4-20')

  TimeDiscName = 'Ketcheson RK4-20'
  nRKStages=20

  !TODO: CFL/DFL coefs unknown, use coef from Niegemann RK14
#if (PP_NodeType==1)
  CFLScaleAlpha(1:15) = &
  (/ 6.9716, 5.7724, 5.1863, 4.7880, 4.4741, 4.2120, 3.9836, 3.7811, 3.5617, 3.3705, 3.1995, 3.0488, 2.9137, 2.7900, 2.6786 /)
#elif (PP_NodeType==2)
  IF (OverintegrationType.GT.0) THEN
    ! Overintegration with Gauss-Lobatto nodes results in a projection DG formulation, i.e. we have to use the Gauss nodes timestep
    CFLScaleAlpha(1:15) = &
    (/ 6.9716, 5.7724, 5.1863, 4.7880, 4.4741, 4.2120, 3.9836, 3.7811, 3.5617, 3.3705, 3.1995, 3.0488, 2.9137, 2.7900, 2.6786 /)
  ELSE
    CFLScaleAlpha(1:15) = &
    (/ 14.7882, 9.5906, 7.9447, 7.0965, 6.5486, 6.1436, 5.8185, 5.5440, 5.3055, 5.0940, 4.9028, 4.7295, 4.5697, 4.4235, 4.2885 /)
  ENDIF
#endif /*PP_NodeType*/
#if PARABOLIC
  RelativeDFL=7.379
#endif
  ALLOCATE(RKg1(nRKStages),RKg2(nRKStages),RKg3(nRKStages),RKb(nRKStages),RKc(nRKStages),RKdelta(nRKStages))
  RKg1 = (/0.0000000000000000e+0,&
          -1.1682479703229380e+0,&
          -2.5112155037089772e+0,&
          -5.5259960154735988e-1,&
           2.9243033509511740e-3,&
          -4.7948973385386493e+0,&
          -5.3095533497183016e+0,&
          -2.3624194456630736e+0,&
           2.0068995756589547e-1,&
          -1.4985808661597710e+0,&
           4.8941228502377687e-1,&
          -1.0387512755259576e-1,&
          -1.3287664273288191e-1,&
           7.5858678822837511e-1,&
          -4.3321586294096939e+0,&
           4.8199700138402146e-1,&
          -7.0924756614960671e-3,&
          -8.8422252029506054e-1,&
          -8.9129367099545231e-1,&
           1.5297157134040762e+0/)

  RKg2 = (/1.0000000000000000e+0,&
           8.8952052154583572e-1,&
           8.8988129100385194e-1,&
           3.5701564494677057e-1,&
           2.4232462479216824e-1,&
           1.2727083024258155e+0,&
           1.1126977210342681e+0,&
           5.1360709645409097e-1,&
           1.1181089682044856e-1,&
           2.7881272382085232e-1,&
           4.9032886260666715e-2,&
           4.1871051065897870e-2,&
           4.4602463796686219e-2,&
           1.4897271251154750e-2,&
           2.6244269699436817e-1,&
          -4.7486056986590294e-3,&
           2.3219312682036197e-2,&
           6.2852588972458059e-2,&
           5.4473719351268962e-2,&
           2.4345446089014514e-2 /)


  RKg3 = (/0.0000000000000000e+0,&
           0.0000000000000000e+0,&
           0.0000000000000000e+0,&
           1.9595487007932735e-1,&
          -6.9871675039100595e-5,&
           1.0592231169810050e-1,&
           1.0730426871909635e+0,&
           8.9257826744389124e-1,&
          -1.4078912484894415e-1,&
          -2.6869890558434262e-1,&
          -6.5175753568318007e-2,&
           4.9177812903108553e-1,&
           4.6017684776493678e-1,&
          -6.4689512947008251e-3,&
           4.4034728024115377e-1,&
           6.1086885767527943e-1,&
           5.0546454457410162e-1,&
           5.4668509293072887e-1,&
           7.1414182420995431e-1,&
          -1.0558095282893749e+0/)

  RKb = (/1.7342385375780556e-1,&
          2.8569004728564801e-1,&
          6.8727044379779589e-1,&
          1.2812121060977319e-1,&
          4.9137180740403122e-4,&
          4.7033584446956857e-2,&
          4.4539998128170821e-1,&
          1.2259824887343720e+0,&
          2.0616463985024421e-2,&
          1.5941162575324802e-1,&
          1.2953803678226099e+0,&
          1.7287352967302603e-3,&
          1.1660483420536467e-1,&
          7.7997036621815521e-2,&
          3.2563250234418012e-1,&
          1.0611520488333197e+0,&
          6.5891625628040993e-4,&
          8.3534647700054046e-2,&
          9.8972579458252483e-2,&
          4.3010116145097040e-2/)

  RKc = (/0.0000000000000000e+0,&
          1.7342385375780556e-1,&
          3.0484982420032158e-1,&
          5.5271395645729193e-1,&
          4.7079204549750037e-2,&
          1.5652540451324129e-1,&
          1.8602224049074517e-1,&
          2.8426620035751449e-1,&
          9.5094727548792268e-1,&
          6.8046501070096010e-1,&
          5.9705366562360063e-1,&
          1.8970821645077285e+0,&
          2.9742664004529606e-1,&
          6.0813463700134940e-1,&
          7.3080004188477765e-1,&
          9.1656999044951792e-1,&
          1.4309687554614530e+0,&
          4.1043824968249148e-1,&
          8.4898255952298962e-1,&
          3.3543896258348421e-1 /)

  RKdelta = (/1.0000000000000000e+0,&
              1.4375468781258596e+0,&
              1.5081653637261594e+0,&
             -1.4575347066062688e-1,&
              3.1495761082838158e-1,&
              3.5505919368536931e-1,&
              2.3616389374566960e-1,&
              1.0267488547302055e-1,&
              3.5991243524519438e+0,&
              1.5172890003890782e+0,&
              1.8171662741779953e+0,&
              2.8762263521436831e+0,&
              4.6350154228218754e-1,&
              1.5573122110727220e+0,&
              2.0001066778080254e+0,&
              9.1690694855534305e-1,&
              2.0474618401365854e+0,&
             -3.2336329115436924e-1,&
              3.2899060754742177e-1,&
              0.0000000000000000e+0/)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 3 register low storage Runge-Kutta 4, 18 stages version - Ketcheson et al 2012
! Error slightly better then Carpenter RK5 and by far better then Niegemann RK14
! Per stage slightly faster then Carpenter RK5, much slower then Niegemann RK14
! => Recommended as default RK scheme!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CASE('ketchesonrk4-18')

  TimeDiscName = 'Ketcheson RK4-18'
  nRKStages=18

  !TODO: CFL/DFL coefs unknown, use coef from Niegemann RK14
#if (PP_NodeType==1)
  CFLScaleAlpha(1:15) = &
  (/ 6.9716, 5.7724, 5.1863, 4.7880, 4.4741, 4.2120, 3.9836, 3.7811, 3.5617, 3.3705, 3.1995, 3.0488, 2.9137, 2.7900, 2.6786 /)
#elif (PP_NodeType==2)
  IF (OverintegrationType.GT.0) THEN
    ! Overintegration with Gauss-Lobatto nodes results in a projection DG formulation, i.e. we have to use the Gauss nodes timestep
    CFLScaleAlpha(1:15) = &
    (/ 6.9716, 5.7724, 5.1863, 4.7880, 4.4741, 4.2120, 3.9836, 3.7811, 3.5617, 3.3705, 3.1995, 3.0488, 2.9137, 2.7900, 2.6786 /)
  ELSE
    CFLScaleAlpha(1:15) = &
    (/ 14.7882, 9.5906, 7.9447, 7.0965, 6.5486, 6.1436, 5.8185, 5.5440, 5.3055, 5.0940, 4.9028, 4.7295, 4.5697, 4.4235, 4.2885 /)
  ENDIF
#endif /*PP_NodeType*/
#if PARABOLIC
  RelativeDFL=7.379
#endif
  ALLOCATE(RKg1(nRKStages),RKg2(nRKStages),RKg3(nRKStages),RKb(nRKStages),RKc(nRKStages),RKdelta(nRKStages))
  RKg1 = (/ 0.0000000000000000e+0,&
            1.1750819811951678e+0,&
            3.0909017892654811e-1,&
            1.4409117788115862e+0,&
           -4.3563049445694069e-1,&
            2.0341503014683893e-1,&
            4.9828356971917692e-1,&
            3.5307737157745489e+0,&
           -7.9318790975894626e-1,&
            8.9120513355345166e-1,&
            5.7091009196320974e-1,&
            1.6912188575015419e-2,&
            1.0077912519329719e+0,&
           -6.8532953752099512e-1,&
            1.0488165551884063e+0,&
            8.3647761371829943e-1,&
            1.3087909830445710e+0,&
            9.0419681700177323e-1 /)

  RKg2 = (/ 1.0000000000000000e+0,&
           -1.2891068509748144e-1,&
            3.5609406666728954e-1,&
           -4.0648075226104241e-1,&
            6.0714786995207426e-1,&
            1.0253501186236846e+0,&
            2.4411240760769423e-1,&
           -1.2813606970134104e+0,&
            8.1625711892373898e-1,&
            1.0171269354643386e-1,&
            1.9379378662711269e-1,&
            7.4408643544851782e-1,&
           -1.2591764563430008e-1,&
            1.1996463179654226e+0,&
            4.5772068865370406e-2,&
            8.3622292077033844e-1,&
           -1.4179124272450148e+0,&
            1.3661459065331649e-1 /)

  RKg3 = (/ 0.0000000000000000e+0,&
            0.0000000000000000e+0,&
            0.0000000000000000e+0,&
            2.5583378537249163e-1,&
            5.2676794366988289e-1,&
           -2.5648375621792202e-1,&
            3.1932438003236391e-1,&
           -3.1106815010852862e-1,&
            4.7631196164025996e-1,&
           -9.8853727938895783e-2,&
            1.9274726276883622e-1,&
            3.2389860855971508e-2,&
            7.5923980038397509e-2,&
            2.0635456088664017e-1,&
           -8.9741032556032857e-2,&
            2.6899932505676190e-2,&
            4.1882069379552307e-2,&
            6.2016148912381761e-2 /)

  RKb = (/ 1.2384169480626298e-1,&
           1.0176262534280349e+0,&
          -6.9732026387527429e-2,&
           3.4239356067806476e-1,&
           1.8177707207807942e-2,&
          -6.1188746289480445e-3,&
           7.8242308902580354e-2,&
          -3.7642864750532951e-1,&
          -4.5078383666690258e-2,&
          -7.5734228201432585e-1,&
          -2.7149222760935121e-1,&
           1.1833684341657344e-3,&
           2.8858319979308041e-2,&
           4.6005267586974657e-1,&
           1.8014887068775631e-2,&
          -1.5508175395461857e-2,&
          -4.0095737929274988e-1,&
           1.4949678367038011e-1 /)

  RKc = (/ 0.0000000000000000e+0,&
           1.2384169480626298e-1,&
           1.1574324659554065e+0,&
           5.4372099141546926e-1,&
           8.8394666834280744e-1,&
          -1.2212042176605774e-1,&
           4.4125685133082082e-1,&
           3.8039092095473748e-1,&
           5.4591107347528367e-2,&
           4.8731855535356028e-1,&
          -2.3007964303896034e-1,&
          -1.8907656662915873e-1,&
           8.1059805668623763e-1,&
           7.7080875997868803e-1,&
           1.1712158507200179e+0,&
           1.2755351018003545e+0,&
           8.0422507946168564e-1,&
           9.7508680250761848e-1 /)

  RKdelta = (/ 1.0000000000000000e+0,&
               3.5816500441970289e-1,&
               5.8208024465093577e-1,&
              -2.2615285894283538e-1,&
              -2.1715466578266213e-1,&
              -4.6990441450888265e-1,&
              -2.7986911594744995e-1,&
               9.8513926355272197e-1,&
              -1.1899324232814899e-1,&
               4.2821073124370562e-1,&
              -8.2196355299900403e-1,&
               5.8113997057675074e-2,&
              -6.1283024325436919e-1,&
               5.6800136190634054e-1,&
              -3.3874970570335106e-1,&
              -7.3071238125137772e-1,&
               8.3936016960374532e-2,&
               0.0000000000000000e+0 /)

CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
                      'Unknown method of time discretization: '//TRIM(TimeDiscMethod))
END SELECT

END SUBROUTINE SetTimeDiscCoefs

!==================================================================================================================================
END MODULE MOD_TimeDisc_Vars
