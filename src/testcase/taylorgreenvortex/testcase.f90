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
!> Subroutines defining the Taylor-Green isentropic vortex testcase
!==================================================================================================================================
MODULE MOD_TestCase
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE DefineParametersTestcase
  MODULE PROCEDURE DefineParametersTestcase
End INTERFACE

INTERFACE InitTestcase
  MODULE PROCEDURE InitTestcase
END INTERFACE

INTERFACE FinalizeTestcase
  MODULE PROCEDURE FinalizeTestcase
END INTERFACE

INTERFACE ExactFuncTestcase
  MODULE PROCEDURE ExactFuncTestcase
END INTERFACE

INTERFACE CalcForcing
  MODULE PROCEDURE DO_NOTHING
END INTERFACE

!INTERFACE TestcaseSource
!  MODULE PROCEDURE TestcaseSource
!END INTERFACE

INTERFACE AnalyzeTestcase
  MODULE PROCEDURE AnalyzeTestcase
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
PUBLIC:: AnalyzeTestcase
PUBLIC:: GetBoundaryFluxTestcase
PUBLIC:: GetBoundaryFVgradientTestcase
PUBLIC:: Lifting_GetBoundaryFluxTestcase

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersTestcase()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Testcase")
CALL prms%CreateIntOption('nWriteStats',      "Write testcase statistics to file at every n-th AnalyzeTestcase step.", '100')
CALL prms%CreateIntOption('nAnalyzeTestCase', "Call testcase specific analysis routines every n-th timestep. "//&
                                              "(Note: always called at global analyze level)"                         , '10')
CALL prms%CreateRealOption('MachNumber',      "Reference Mach number for TGV testcase", '0.1')
CALL prms%CreateLogicalOption('IniConstDens', "True:  Initial Density     field is constant, Temperature entails fluctuations.\n&
                                              &False: Initial Temperature field is constant, Density     entails fluctuations." &
                                              ,'T')
END SUBROUTINE DefineParametersTestcase


!==================================================================================================================================
!> Initialize testcase variables
!==================================================================================================================================
SUBROUTINE InitTestcase()
! MODULES
USE MOD_Globals
USE MOD_TestCase_Vars
USE MOD_ReadInTools,    ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_Output_Vars,    ONLY: ProjectName
USE MOD_Output,         ONLY: InitOutputToFile
USE MOD_EOS_Vars,       ONLY: Kappa,R
#if (PP_VISC==1)
USE MOD_EOS_Vars,       ONLY: Tref,Ts
#elif (PP_VISC==2)
USE MOD_EOS_Vars,       ONLY: Tref,ExpoSuth,mu0
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=31)        :: varnames(nTGVVars)
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TESTCASE TAYLOR-GREEN VORTEX...'

! Length of Buffer for TGV output
nWriteStats      = GETINT( 'nWriteStats')
nAnalyzeTestCase = GETINT( 'nAnalyzeTestCase')

! Check whether initial density or rather temperature field should be constant
IniConstDens = GETLOGICAL('IniConstDens')

! Get reference Mach number of TGV
Ma0 = GETREAL('MachNumber')

! Compute reference temperature and pressure of TGV via ideal gas relation and Mach number
p0 = (U0/Ma0)**2/Kappa*rho0
T0 = p0/(rho0*R)

! Viscosity computation has to be adapted to yield correct results with dimensionless temperature as input
#if PP_VISC == 1
! Reformulate Sutherland's law (ratio Ts/Tref is kept constant, hence, Ts doesnt change)
Tref = 1./T0  ! ATTENTION: Tref = 1./Tref | Ts = Ts/Tref
! Provide user warning
CALL PrintWarning("  Viscosity via Sutherland's Law is computed in dimensionless form for the TGV testcase.\n&
                  &  Hence, Tref is set to dimensionless reference temperature T0 = 1/(Kappa*R*MachNumber^2)")
SWRITE(UNIT_stdOut,'(A,ES13.7)') ' |   T0   = ',T0
SWRITE(UNIT_stdOut,'(A,ES13.7)') ' |   Tref = ',1./Tref
SWRITE(UNIT_stdOut,'(A,ES13.7)') ' |   Ts   = ',Ts*T0
#elif PP_VISC == 2
! Adapt precomputed viscosity with new reference temperature
mu0 = mu0 * (Tref/T0)**ExpoSuth
Tref= T0
! Provide user warning
CALL PrintWarning("  Viscosity via Power Law is computed in dimensionless form for the TGV testcase.\n&
                  &  Hence, Tref is set to dimensionless reference temperature T0 = 1/(Kappa*R*MachNumber^2)")
SWRITE(UNIT_stdOut,'(A,ES13.7)') ' |   T0   = ',T0
SWRITE(UNIT_stdOut,'(A,ES13.7)') ' |   Tref = ',Tref
#endif


IF(MPIRoot)THEN
  ALLOCATE(Time(nWriteStats))
  ALLOCATE(writeBuf(nTGVvars,nWriteStats))
  Filename = TRIM(ProjectName)//'_TGVAnalysis'

#if PARABOLIC
  varnames(1) ="Dissipation Rate Incompressible"
  varnames(2) ="Dissipation Rate Compressible"
  varnames(3) ="Ekin incomp"
  varnames(4) ="Ekin comp"
  varnames(5) ="Enstrophy comp"
  varnames(6) ="DR_u"
  varnames(7) ="DR_S"
  varnames(8) ="DR_Sd"
  varnames(9) ="DR_p"
  varnames(10)="Maximum Vorticity"
  varnames(11)="Mean Temperature"
  varnames(12)="uprime"
  varnames(13)="Mean Entropy"
  varnames(14)="ED_S"
  varnames(15)="ED_D"
#else
  varnames(1) ="Ekin incomp"
  varnames(2) ="Ekin comp"
  varnames(3)= "Mean Temperature"
  varnames(4)= "uprime"
  varnames(5)= "Mean Entropy"
#endif
  CALL InitOutputToFile(FileName,'Taylor-Green Vortex Analysis Data',nTGVVars,varnames)
END IF

SWRITE(UNIT_stdOut,'(A)')' INIT TESTCASE TAYLOR-GREEN VORTEX DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitTestcase


!==================================================================================================================================
!> Specifies the initial conditions for the TGV testcase in a (weakly) compressible formulation with size [0,2*PI]^3 following:
!>    "Comparison of high-order numerical methodologies for the simulation of the supersonic Taylor-Green Vortex flow",
!>    Chapelier et al., Physics of Fluids, 2024.
!>
!> Two versions of initial conditions are implemented for the (weakly) compressible form. Either set initial density field as
!> constant and compute the temperature field to be thermodynamically consistent to the intial pressure (containing fluctuations)
!> or vice versa.
!==================================================================================================================================
SUBROUTINE ExactFuncTestcase(tIn,x,Resu,Resu_t,Resu_tt)
! MODULES
USE MOD_Globals
USE MOD_EOS,           ONLY: PrimToCons
USE MOD_EOS_Vars,      ONLY: R
USE MOD_Testcase_Vars, ONLY: rho0,U0,p0,T0,IniConstDens
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
REAL            :: prim(PP_nVarPrim)
!==================================================================================================================================
! Initial velocity and pressure with reference quantities:
!   -U0,rho0 parameters set to unity in Testcase_Vars
!   -T0,  p0 precomputed in InitTestcase based on U0, rho0 and user-specified Ma0
prim(VEL1) = U0*SIN(x(1))*COS(x(2))*COS(x(3))                                        ! (6)
prim(VEL2) =-U0*COS(x(1))*SIN(x(2))*COS(x(3))                                        ! (6)
prim(VEL3) = 0.                                                                      ! (6)
! Expand product in (7) analytically to reduce influence of limited floating point precision
prim(PRES) = p0 + (rho0*U0**2)/16.*( COS(2*x(1))*COS(2.*x(3)) + 2.*COS(2.*x(2)) + 2.*COS(2.*x(1)) + COS(2*x(2))*COS(2.*x(3)) ) ! (7)


! Two different variations of initialization are possible:
! Either set intial density field as constant and compute local temperature thermodynamically consistent or vice versa.
IF(IniConstDens) THEN
  prim(DENS) = rho0 ! Constant initial density
  prim(TEMP) = prim(PRES)/(prim(DENS)*R)
ELSE
  prim(TEMP) = T0   ! Constant initial temperature
  prim(DENS) = prim(PRES)/(prim(TEMP)*R)
END IF

CALL PrimToCons(prim,Resu)

Resu_t =0.
Resu_tt=0.
END SUBROUTINE ExactFuncTestcase


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


!==================================================================================================================================
!> Perform TGV-specific analysis: compute dissipation rates, kinetic energy and enstrophy
!==================================================================================================================================
SUBROUTINE AnalyzeTestcase(t,doFlush)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TestCase_Vars
USE MOD_DG_Vars,        ONLY: U
USE MOD_EOS,            ONLY: ConsToPrim
#if PARABOLIC
USE MOD_Lifting_Vars,   ONLY: GradUx,GradUy,GradUz
USE MOD_EOS_Vars,       ONLY: mu0
#if PP_VISC==1
USE MOD_Viscosity,      ONLY: muSuth
#elif PP_VISC==2
USE MOD_EOS_Vars,       ONLY: ExpoSuth
#endif
#endif /* PARABOLIC */
USE MOD_Analyze_Vars,   ONLY: NAnalyze,Vdm_GaussN_NAnalyze,wGPVolAnalyze,Vol
USE MOD_EOS_Vars,       ONLY: sKappaM1,Kappa
USE MOD_Mesh_Vars,      ONLY: sJ
USE MOD_ChangeBasis,    ONLY: ChangeBasis3D
USE MOD_Mesh_Vars,      ONLY: nElems
#if USE_MPI
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< simulation time
LOGICAL,INTENT(IN)              :: doFlush                !< indicate that data has to be written
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
REAL                            :: UPrim(PP_nVarPrim)
REAL                            ::  U_NAnalyze(PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: sJ_NAnalyze(      1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: sJ_N(1,0:PP_N,0:PP_N,0:PP_N) ! local array for sJ due to interface of changeBasis
REAL                            :: IntFactor                     ! Integration weights and Jacobian
REAL                            :: T_mean                        ! Mean temperature
REAL                            :: Entropy                       ! Mean entropy in the domain
REAL                            :: uPrime
REAL                            :: Ekin,Ekin_comp                ! Integrated incompr/compr. kin. energy (0.5*u*u and 0.5*rho*u*u)
#if USE_MPI
REAL                            :: Ekin_Glob,Ekin_comp_Glob,T_mean_Glob,Entropy_Glob
#endif
#if PARABOLIC
INTEGER                         :: ii
REAL                            :: GradVel(1:3,1:3)
REAL                            :: GradVelx(1:3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: GradVely(1:3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: GradVelz(1:3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: S( 1:3,1:3)                   ! Strain rate tensor S (symmetric)
REAL                            :: Sd(1:3,1:3)                   ! Deviatoric part of the strain rate tensor S
REAL                            :: divU                          ! Divergence of velocity vector
REAL                            :: Vorticity(1:3)                ! Vorticity = curl(Velocity)
REAL                            :: Vorticity_max                 ! max. vorticity in domain
REAL                            :: Enstr                         ! Enstrophy
REAL                            :: DR_u,DR_S,DR_Sd,DR_p          ! Contributions to dissipation rates
REAL                            :: ED_S,ED_D                     ! Solenoidal and dilitational dissipation rates
#if USE_MPI
REAL                            :: DR_u_Glob,DR_S_glob,DR_Sd_Glob,DR_p_Glob,ED_S_Glob,ED_D_Glob,Enstr_Glob,Vorticity_max_glob
#endif
#endif /* PARABOLIC */
!==================================================================================================================================
! The general workflow of the routine is as follows:
!
! 1. Nullify quantities
! 2. Interpolate solution and gradients to analyze mesh
! 3. Precompute primitive and gradient-based quantities
! 4. Integrate
! 5. Communicate via MPI
! 6. Normalize
! 7. Write to File


!----------------------------------------------------------------------------------------------------------------------------------
! 1. Nullify Quantities
!----------------------------------------------------------------------------------------------------------------------------------
Ekin=0.;Ekin_comp=0.;T_mean=0.;Entropy=0.
#if PARABOLIC
DR_u=0.;DR_S=0.;DR_Sd=0.;DR_p=0.;ED_S=0.;ED_D=0.;Enstr=0.;Vorticity_max=0.
#endif

DO iElem=1,nElems
  !----------------------------------------------------------------------------------------------------------------------------------
  ! 2. Interpolate solution and gradients to analyze mesh
  sJ_N(1,:,:,:)=sJ(:,:,:,iElem,0)
  CALL ChangeBasis3D(      1, PP_N, NAnalyze, Vdm_GaussN_NAnalyze,          sJ_N, sJ_NAnalyze)
  ! Interpolate the solution to the analyze grid
  CALL ChangeBasis3D(PP_nVar, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, U(:,:,:,:,iElem),  U_NAnalyze)
#if PARABOLIC
  ! Interpolate the gradient of the velocity to the analyze grid
  CALL ChangeBasis3D(3, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, GradUx(LIFT_VELV,:,:,:,iElem), GradVelx)
  CALL ChangeBasis3D(3, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, GradUy(LIFT_VELV,:,:,:,iElem), GradVely)
  CALL ChangeBasis3D(3, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, GradUz(LIFT_VELV,:,:,:,iElem), GradVelz)
#endif

  DO k=0,NAnalyze;DO j=0,NAnalyze;DO i=0,NAnalyze
    !----------------------------------------------------------------------------------------------------------------------------------
    ! 3. Precompute primitive and gradient-based quantities
    CALL ConsToPrim(UPrim,U_NAnalyze(:,i,j,k))
#if PARABOLIC
    ! compute velocity gradient tensor GradVel
    GradVel(:,1)=GradVelx(:,i,j,k)
    GradVel(:,2)=GradVely(:,i,j,k)
    GradVel(:,3)=GradVelz(:,i,j,k)
    ! compute divergence of velocity
    divU=GradVel(1,1)+GradVel(2,2)+GradVel(3,3)
    ! compute rate-of-strain tensor S
    S=0.5*(Gradvel+TRANSPOSE(GradVel))
    ! deviatoric part of strain tensor Sd
    Sd=S
    DO ii=1,3
      Sd(ii,ii)=Sd(ii,ii)-1./3.*divU
    END DO
    ! compute vorticity (curl of velocity) and max(vorticity)
    Vorticity(1) = GradVel(3,2) - GradVel(2,3)
    Vorticity(2) = GradVel(1,3) - GradVel(3,1)
    Vorticity(3) = GradVel(2,1) - GradVel(1,2)
    Vorticity_max=MAX(Vorticity_max,NORM2(Vorticity))
#endif

    !----------------------------------------------------------------------------------------------------------------------------------
    ! 4. Integrate
    IntFactor=wGPVolAnalyze(i,j,k)/sJ_NAnalyze(1,i,j,k) ! Integration weight and Jacobian
      ! Temperature
    T_mean   =T_mean   +IntFactor*UPrim(TEMP)
      ! Entropy
    Entropy  =Entropy  +IntFactor*(-sKappaM1)*UPrim(DENS)*(LOG(UPrim(PRES))-Kappa*LOG(UPrim(DENS)))
      ! Kinetic Energy incompressible
    Ekin     =Ekin     +IntFactor*0.5            *DOT_PRODUCT(UPrim(VELV),UPrim(VELV))
      ! Kinetic Energy compressible
    Ekin_comp=Ekin_comp+IntFactor*0.5*UPrim(DENS)*DOT_PRODUCT(UPrim(VELV),UPrim(VELV))
#if PARABOLIC
      ! Enstrophy compressible
    Enstr=Enstr+IntFactor*0.5*UPrim(DENS)*DOT_PRODUCT(Vorticity,Vorticity)
      ! dissipation rate epsilon incomp from velocity gradient tensor (Diss Fauconnier)
    DR_u =DR_u +IntFactor*SUM(GradVel(:,:)*GradVel(:,:))
      ! dissipation rate epsilon incomp from strain rate tensor (incomp) (Sagaut)
    DR_S =DR_S +IntFactor*SUM( S(:,:)* S(:,:))
      ! dissipation rate epsilon 1 from deviatoric part of strain rate tensor Sd (compressible)
    DR_Sd=DR_Sd+IntFactor*SUM(Sd(:,:)*Sd(:,:))
      ! dissipation rate epsilon 3 from compressibility effects (compressible)
    DR_p =DR_p +IntFactor*UPrim(PRES)*divU
      ! solenoidal part of the kinetic energy dissipation
    ED_S =ED_S +IntFactor*VISCOSITY_TEMPERATURE(UPrim(TEMP))*DOT_PRODUCT(Vorticity,Vorticity)
      ! dilatational component of the kinetic energy dissipation
    ED_D =ED_D +IntFactor*VISCOSITY_TEMPERATURE(UPrim(TEMP))*divU**2
#endif
  END DO;END DO;END DO
END DO

!----------------------------------------------------------------------------------------------------------------------------------
! 5. Communicate
#if USE_MPI
CALL MPI_REDUCE(T_mean   ,   T_mean_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(Entropy  ,  Entropy_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(Ekin     ,     Ekin_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(Ekin_comp,Ekin_comp_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
T_mean   =T_mean_Glob
Entropy  =Entropy_Glob
Ekin     =Ekin_Glob
Ekin_comp=Ekin_comp_Glob
#if PARABOLIC
CALL MPI_REDUCE(DR_u   , DR_u_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(DR_S   , DR_S_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(DR_Sd  ,DR_Sd_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(DR_p   , DR_p_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(ED_S   , ED_S_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(ED_D   , ED_D_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(Enstr  ,Enstr_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(Vorticity_max,Vorticity_max_Glob,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_FLEXI,iError)
DR_u =DR_u_Glob
DR_S =DR_S_Glob
DR_Sd=DR_Sd_Glob
DR_p =DR_p_Glob
ED_S =ED_S_Glob
ED_D =ED_D_Glob
Enstr=Enstr_Glob
Vorticity_max=Vorticity_max_Glob
#endif /* PARABOLIC */
#endif /* USE_MPI */

IF(MPIRoot) THEN
  !----------------------------------------------------------------------------------------------------------------------------------
  ! 6. Normalize integrals
  Ekin     =     Ekin/Vol
  Ekin_comp=Ekin_comp/Vol/rho0
  T_mean   =   T_Mean/Vol
  Entropy  =  Entropy/Vol
#if PARABOLIC
  Enstr= Enstr/(rho0*Vol)
  DR_u = DR_u *mu0/Vol
  DR_S = DR_S *2.*mu0/(rho0*Vol)
  DR_Sd= DR_Sd*2.*mu0/(rho0*Vol)
  DR_p =-DR_p /(rho0*Vol)
  ED_S = ED_S /Vol
  ED_D = ED_D *4./3./Vol
#endif

  uPrime=SQRT(2./3.*Ekin)

  !----------------------------------------------------------------------------------------------------------------------------------
  ! 7. Write to file
  ioCounter       = ioCounter+1
  Time(ioCounter) = t
#if PARABOLIC
  writeBuf(1:nTGVvars,ioCounter) = (/DR_S,DR_Sd+DR_p,Ekin,Ekin_comp,Enstr,DR_u,DR_S,DR_Sd,DR_p,&
                                     Vorticity_max,T_mean,uPrime,Entropy,ED_S,ED_D/)
#else
  writeBuf(1:nTGVvars,ioCounter) = (/Ekin,Ekin_comp,T_mean,uPrime,Entropy/)
#endif
  IF((ioCounter.EQ.nWriteStats) .OR. doFlush) THEN
    CALL WriteStats()
    ioCounter=0
  END IF
END IF

END SUBROUTINE AnalyzeTestcase


!==================================================================================================================================
!> Write TGV Analysis Data to File
!==================================================================================================================================
SUBROUTINE WriteStats()
! MODULES
USE MOD_TestCase_Vars
USE MOD_Output,        ONLY: OutputToFile
IMPLICIT NONE
!==================================================================================================================================
CALL OutputToFile(FileName,Time(1:ioCounter),(/nTGVvars,ioCounter/),&
                           RESHAPE(writeBuf(:,1:ioCounter),(/nTGVVars*ioCounter/)))
END SUBROUTINE WriteStats


!==================================================================================================================================
!> Flushes Buffer to File, Deallocation
!==================================================================================================================================
SUBROUTINE FinalizeTestcase()
! MODULES
USE MOD_Globals      ,ONLY:MPIRoot
USE MOD_TestCase_Vars,ONLY:writeBuf,Time
IMPLICIT NONE
!==================================================================================================================================
IF(MPIRoot)THEN
  SDEALLOCATE(Time)
  SDEALLOCATE(writeBuf)
END IF
END SUBROUTINE


SUBROUTINE DO_NOTHING(optionalREAL,optionalREAL2)
IMPLICIT NONE
REAL, OPTIONAL,INTENT(IN) :: optionalREAL,optionalREAL2
END SUBROUTINE DO_NOTHING


SUBROUTINE GetBoundaryFluxTestcase(SideID,t,Nloc,Flux,UPrim_master,                   &
#if PARABOLIC
                           gradUx_master,gradUy_master,gradUz_master,&
#endif
                           NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: SideID  !< ID of current side
REAL,INTENT(IN)      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)   :: Nloc    !< polynomial degree
REAL,INTENT(IN)      :: UPrim_master( PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))    !< inner surface solution
#if PARABOLIC
REAL,INTENT(IN)      :: gradUx_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !> inner surface solution gradients in x-direction
REAL,INTENT(IN)      :: gradUy_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !> inner surface solution gradients in y-direction
REAL,INTENT(IN)      :: gradUz_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !> inner surface solution gradients in z-direction
#endif /*PARABOLIC*/
REAL,INTENT(IN)      :: NormVec (  3,0:Nloc,0:ZDIM(Nloc))  !< normal vectors on surfaces
REAL,INTENT(IN)      :: TangVec1(  3,0:Nloc,0:ZDIM(Nloc))  !< tangential1 vectors on surfaces
REAL,INTENT(IN)      :: TangVec2(  3,0:Nloc,0:ZDIM(Nloc))  !< tangential2 vectors on surfaces
REAL,INTENT(IN)      :: Face_xGP(  3,0:Nloc,0:ZDIM(Nloc))  !< positions of surface flux points
REAL,INTENT(OUT)     :: Flux(PP_nVar,0:Nloc,0:ZDIM(Nloc))  !< resulting boundary fluxes
!==================================================================================================================================
END SUBROUTINE GetBoundaryFluxTestcase


SUBROUTINE GetBoundaryFVgradientTestcase(SideID,t,gradU,UPrim_master)
USE MOD_PreProc
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                                   !< ID of current side
REAL,INTENT(IN)    :: t                                        !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive solution from the inside
REAL,INTENT(OUT)   :: gradU       (PP_nVarPrim,0:PP_N,0:PP_NZ) !< FV boundary gradient
!==================================================================================================================================
END SUBROUTINE GetBoundaryFVgradientTestcase


SUBROUTINE Lifting_GetBoundaryFluxTestcase(SideID,t,UPrim_master,Flux)
USE MOD_PreProc
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                                   !< ID of current side
REAL,INTENT(IN)    :: t                                        !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive solution from the inside
REAL,INTENT(OUT)   :: Flux(     PP_nVarLifting,0:PP_N,0:PP_NZ) !< lifting boundary flux
!==================================================================================================================================
END SUBROUTINE Lifting_GetBoundaryFluxTestcase

END MODULE MOD_TestCase
