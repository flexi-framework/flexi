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
#include "eos.h"

#if FV_ENABLED
#error "This testcase is not tested with FV"
#endif

!==================================================================================================================================
!> The channel case is a setup according to the Moser channel:
!>  - Moser, Robert D., John Kim, and Nagi N. Mansour. "Direct numerical simulation of turbulent channel flow up to Re_tau= 590."
!>    Physics of fluids 11.4 (1999): 943-945.
!>  - Lee, Myoungkyu, and Robert D. Moser. "Direct numerical simulation of turbulent channel flow up to Re_tau=5200."
!>    Journal of Fluid Mechanics 774 (2015): 395-415.
!> The channel halfwidth is set to 1 and the Reynolds number is thus set with mu0 = 1/Re_tau. Further, rho=1 and the pressure is
!> computed to obtain the specified Bulk Mach number (Mach=0.1 for the Moser case). Hence, u_tau = tau = -dp/dx = 1 .
!==================================================================================================================================
MODULE MOD_TestCase
! MODULES
USE MOD_TestCase_Vars
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE DefineParametersTestcase
  MODULE PROCEDURE DefineParametersTestcase
END INTERFACE

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
  MODULE PROCEDURE CalcForcing
END INTERFACE

INTERFACE TestcaseSource
  MODULE PROCEDURE TestcaseSource
END INTERFACE

INTERFACE AnalyzeTestCase
  MODULE PROCEDURE AnalyzeTestCase
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
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersTestcase()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Testcase")
CALL prms%CreateRealOption('ChannelMach', "Bulk mach number used in the channel testcase."                      , '0.1')
CALL prms%CreateIntOption('nWriteStats', "Write testcase statistics to file at every n-th AnalyzeTestcase step.", '100')
CALL prms%CreateIntOption('nAnalyzeTestCase', "Call testcase specific analysis routines every n-th timestep. "//&
                                              "(Note: always called at global analyze level)"                   , '1000')
END SUBROUTINE DefineParametersTestcase

!==================================================================================================================================
!> Initializes the Channel testcase. The initial pressure is set to match the specified Bulk Mach number. For this, the initial
!> Bulk Velocity has to be estimated. Here, an analytical approximation is used to estimate the bulk velocity depending on the
!> specific Reynolds number.
!> TODO: Find source of formula for bulk velocity.
!==================================================================================================================================
SUBROUTINE InitTestcase()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools,        ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_Equation_Vars,      ONLY: RefStatePrim,IniRefState,RefStateCons
USE MOD_EOS,                ONLY: PrimToCons
USE MOD_EOS_Vars,           ONLY: kappa,mu0,R
USE MOD_Output,             ONLY: InitOutputToFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER           :: c1 = 2.4390244 ! Empirical parameter for estimation of bulkVel
REAL                     :: bulkMach,pressure
REAL                     :: UE(PP_2Var)
CHARACTER(LEN=7)         :: varnames(2)
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TESTCASE CHANNEL...'

#if FV_ENABLED
CALL CollectiveStop(__STAMP__,'The testcase has not been implemented for FV yet!')
#endif

nWriteStats      = GETINT('nWriteStats')
nAnalyzeTestCase = GETINT( 'nAnalyzeTestCase')

! Compute initial guess for bulk velocity for given Re_tau to compute background pressure
Re_tau  = 1./mu0
bulkVel = (Re_tau+c1)*LOG(Re_tau+c1) + 1.3064019*(Re_tau + 29.627395*EXP(-1./11.*Re_tau) + 0.66762137*(Re_tau+3)*EXP(-Re_tau/3.))
bulkVel = 1./Re_tau * (c1*bulkVel - 97.4857927165)

! Set the background pressure according to chosen bulk Mach number
bulkMach                       = GETREAL('ChannelMach')
pressure                       = (bulkVel/bulkMach)**2*RefStatePrim(DENS,IniRefState)/kappa
RefStatePrim(PRES,IniRefState) = pressure
! TODO: ATTENTION only sRho and Pressure of UE filled!!!
UE(EXT_SRHO)                   = 1./RefStatePrim(DENS,IniRefState)
UE(EXT_PRES)                   = RefStatePrim(PRES,IniRefState)
RefStatePrim(TEMP,IniRefState) = TEMPERATURE_HE(UE)
CALL PrimToCons(RefStatePrim(:,IniRefState),RefStateCons(:,IniRefState))

IF(MPIRoot) THEN
  WRITE(UNIT_stdOut,*) 'Bulk velocity based on initial velocity Profile =',bulkVel
  WRITE(UNIT_stdOut,*) 'Associated Pressure for Mach = ',bulkMach,' is', pressure

  ! Initialize output of statistics to file
  ALLOCATE(writeBuf(3,nWriteStats))
  Filename = TRIM(ProjectName)//'_Stats'
  varnames(1) = 'dpdx'
  varnames(2) = 'bulkVel'
  CALL InitOutputToFile(Filename,'Statistics',2,varnames)
END IF

SWRITE(UNIT_stdOut,'(A)')' INIT TESTCASE CHANNEL DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitTestcase



!==================================================================================================================================
!> Initial conditions for the channel testcase. Initializes a velocity profile in the streamwise direction and superimposes
!> velocity disturbances to accelerate the development of turbulence.
!==================================================================================================================================
SUBROUTINE ExactFuncTestcase(tIn,x,Resu,Resu_t,Resu_tt)
! MODULES
USE MOD_PreProc,      ONLY: PP_PI
USE MOD_Equation_Vars,ONLY: RefStatePrim,IniRefState
USE MOD_EOS,          ONLY: PrimToCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: x(3),tIn
REAL,INTENT(OUT)                :: Resu(PP_nVar),Resu_t(PP_nVar),Resu_tt(PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Prim(PP_nVarPrim)
REAL                            :: yPlus,Amplitude
!==================================================================================================================================
Prim(:) = RefStatePrim(:,IniRefState)

! Initialize mean turbulent velocity profile in x
IF(x(2).LE.0) THEN
  yPlus = (x(2)+1.)*Re_tau ! Lower half
ELSE
  yPlus = (1.-x(2))*Re_tau ! Upper half
END IF
Prim(VEL1) = 1./0.41*LOG(1+0.41*yPlus)+7.8*(1-EXP(-yPlus/11.)-yPlus/11.*EXP(-yPlus/3.))

! Superimpose sinusoidal disturbances to accelerate development of turbulence
Amplitude = 0.1*Prim(VEL1)
#if EQNSYSNR == 2
Prim(VEL1)=Prim(VEL1)+SIN(20.0*PP_PI*(x(2)/(2.0)))*SIN(20.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL1)=Prim(VEL1)+SIN(30.0*PP_PI*(x(2)/(2.0)))*SIN(30.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL1)=Prim(VEL1)+SIN(35.0*PP_PI*(x(2)/(2.0)))*SIN(35.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL1)=Prim(VEL1)+SIN(40.0*PP_PI*(x(2)/(2.0)))*SIN(40.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL1)=Prim(VEL1)+SIN(45.0*PP_PI*(x(2)/(2.0)))*SIN(45.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL1)=Prim(VEL1)+SIN(50.0*PP_PI*(x(2)/(2.0)))*SIN(50.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude

Prim(VEL2)=Prim(VEL2)+SIN(30.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(30.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL2)=Prim(VEL2)+SIN(35.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(35.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL2)=Prim(VEL2)+SIN(40.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(40.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL2)=Prim(VEL2)+SIN(45.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(45.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude
Prim(VEL2)=Prim(VEL2)+SIN(50.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(50.0*PP_PI*(x(3)/(2*PP_PI)))*Amplitude

Prim(VEL3)=Prim(VEL3)+SIN(30.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(30.0*PP_PI*(x(2)/(2.0)))*Amplitude
Prim(VEL3)=Prim(VEL3)+SIN(35.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(35.0*PP_PI*(x(2)/(2.0)))*Amplitude
Prim(VEL3)=Prim(VEL3)+SIN(40.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(40.0*PP_PI*(x(2)/(2.0)))*Amplitude
Prim(VEL3)=Prim(VEL3)+SIN(45.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(45.0*PP_PI*(x(2)/(2.0)))*Amplitude
Prim(VEL3)=Prim(VEL3)+SIN(50.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(50.0*PP_PI*(x(2)/(2.0)))*Amplitude
#endif

Prim(TEMP) = 0. ! T does not matter for prim to cons
CALL PrimToCons(prim,Resu)

Resu_t =0.
Resu_tt=0.

END SUBROUTINE ExactFuncTestcase


!==================================================================================================================================
!> Compute bulk velocity for forcing term of the channel.
!==================================================================================================================================
SUBROUTINE CalcForcing(t,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,        ONLY: U
USE MOD_Mesh_Vars,      ONLY: sJ
USE MOD_Analyze_Vars,   ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,      ONLY: nElems
#if USE_MPI
USE MOD_MPI_Vars
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t,dt
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!==================================================================================================================================
BulkVel =0.
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    BulkVel = BulkVel+U(MOM1,i,j,k,iElem)/U(DENS,i,j,k,iElem)*wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
END DO

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,BulkVel,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
#endif
BulkVel = BulkVel/Vol
END SUBROUTINE CalcForcing


!==================================================================================================================================
!> Apply forcing term equivalent to a constant streamwise pressure gradient
!==================================================================================================================================
SUBROUTINE TestcaseSource(Ut)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars, ONLY:sJ,nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!==================================================================================================================================
! Apply forcing with the pressure gradient
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    Ut(MOM1,i,j,k,iElem) = Ut(MOM1,i,j,k,iElem) - dpdx/sJ(i,j,k,iElem,0)
    Ut(ENER,i,j,k,iElem) = Ut(ENER,i,j,k,iElem) - dpdx/sJ(i,j,k,iElem,0)*BulkVel
  END DO; END DO; END DO
END DO
END SUBROUTINE TestcaseSource

!==================================================================================================================================
!> Output testcase statistics
!==================================================================================================================================
SUBROUTINE WriteStats()
! MODULES
USE MOD_Output,       ONLY:OutputToFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL OutputToFile(FileName,writeBuf(1,1:ioCounter),(/2,ioCounter/),RESHAPE(writeBuf(2:3,1:ioCounter),(/2*ioCounter/)))
ioCounter=0
END SUBROUTINE WriteStats

!==================================================================================================================================
!> Specifies periodic hill testcase
!==================================================================================================================================
SUBROUTINE AnalyzeTestcase(Time,doFlush)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time                   !< simulation time
LOGICAL,INTENT(IN)              :: doFlush                !< indicate that data has to be written
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF(MPIRoot)THEN
  ioCounter=ioCounter+1
  writeBuf(:,ioCounter) = (/Time, dpdx, BulkVel/)
  IF(ioCounter.GE.nWriteStats .OR. doFlush) CALL WriteStats()
END IF
END SUBROUTINE AnalyzeTestCase

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE FinalizeTestcase()
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
IF(MPIRoot) THEN
  SDEALLOCATE(writeBuf)
END IF
END SUBROUTINE



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
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                                   !< ID of current side
REAL,INTENT(IN)    :: t                                        !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive solution from the inside
REAL,INTENT(OUT)   :: gradU       (PP_nVarPrim,0:PP_N,0:PP_NZ) !< FV boundary gradient
!==================================================================================================================================
END SUBROUTINE GetBoundaryFVgradientTestcase


SUBROUTINE Lifting_GetBoundaryFluxTestcase(SideID,t,UPrim_master,Flux)
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID                                   !< ID of current side
REAL,INTENT(IN)    :: t                                        !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive solution from the inside
REAL,INTENT(OUT)   :: Flux(     PP_nVarLifting,0:PP_N,0:PP_NZ) !< lifting boundary flux
!==================================================================================================================================
END SUBROUTINE Lifting_GetBoundaryFluxTestcase

END MODULE MOD_TestCase
