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
#error "HIT testcase is not tested with FV"
#endif

#if !(PP_dim == 3)
#error "HIT testcase is only implemented for 3D"
#endif

!==================================================================================================================================
!> Subroutines defining homogenenous isotropic turbulence with linear forcing
!> >> This testcase is a heavily simplified version of the anisotropic linear forcing (ALF) approach. See the corresponding git
!> >> branch if you are interested in the full model
!==================================================================================================================================
MODULE MOD_Testcase
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
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Testcase")
CALL prms%CreateIntOption('nWriteStats'     , "Write testcase statistics to file at every n-th AnalyzeTestcase step.", '100')
CALL prms%CreateIntOption('nAnalyzeTestCase', "Call testcase specific analysis routines every n-th timestep. "//&
                                              "(Note: always called at global analyze level)", '10')

! Parameters for HIT forcing
CALL prms%CreateLogicalOption('HIT_Forcing', 'Flag to perform HIT forcing'                            ,'T')
CALL prms%CreateLogicalOption('HIT_Avg'    , 'Flag to perform spatial averaging of HIT forcing'       ,'T')
CALL prms%CreateLogicalOption('HIT_1st'    , 'Flag to update HIT forcing only in first RK stage'      ,'T')
CALL prms%CreateRealOption(   'HIT_k'      , 'Target turbulent kinetic energy in HIT'                 ,'0.')
CALL prms%CreateRealOption(   'HIT_tFilter', 'Temp filter width of exponential, explicit time filter' ,'0.')
CALL prms%CreateRealOption(   'HIT_tauRMS' , 'Strength of RMS forcing.'                               ,'0.')

END SUBROUTINE DefineParametersTestcase


!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE InitTestcase()
! MODULES
USE MOD_Globals
#if PP_N==N
USE MOD_PreProc,            ONLY: N
#endif
USE MOD_Equation_Vars,      ONLY: RefStatePrim,IniRefState
USE MOD_HDF5_Input,         ONLY: File_ID,OpenDataFile,CloseDataFile,ReadArray,DatasetExists,GetDataSize
USE MOD_IO_HDF5,            ONLY: AddToFieldData,FieldOut,HSIZE,nDims
USE MOD_Mesh_Vars,          ONLY: nElems,offsetElem
USE MOD_Output,             ONLY: InitOutputToFile
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_ReadInTools,        ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_Restart_Vars,       ONLY: doRestart,restartFile,interpolateSolution
USE MOD_TestCase_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                  :: HITDataExists
INTEGER                  :: HSize_proc(5)
REAL,ALLOCATABLE         :: HIT_local(:,:,:,:,:)
CHARACTER(LEN=31)        :: varnames(nHITVars)
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TESTCASE HOMOGENEOUS ISOTROPIC TURBULENCE...'

#if FV_ENABLED
CALL CollectiveStop(__STAMP__,'The testcase has not been implemented for FV yet!')
#endif

HIT_Forcing = GETLOGICAL('HIT_Forcing','.TRUE.')

IF(HIT_Forcing) THEN
  ! Initialize HIT filter
  HIT_tFilter = GETREAL('HIT_tFilter','0.')

  ! Read only rho from IniRefState
  HIT_rho = RefStatePrim(1,IniRefState)

  ! Read target turbulent kinetic energy and relaxation time
  HIT_k      = GETREAL('HIT_k'     ,'0.')
  HIT_tauRMS = GETREAL('HIT_tauRMS','0.')

  IF ((HIT_tFilter.EQ.0)) &
    CALL CollectiveStop(__STAMP__, 'HIT target parameters cannot be zero!')

  ! Allocate array for temporally filtered RMS
  ALLOCATE(HIT_RMS(1:3,0:PP_N,0:PP_N,0:PP_NZ,nElems))

  ! Try to restart from state file.
  IF(DoRestart)THEN
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    CALL DatasetExists(File_ID,'HIT',HITDataExists)
    IF (HITDataExists) THEN
      CALL GetDataSize(File_ID,'HIT',nDims,HSize)
      HSize_proc = INT(HSize)
      HSize_proc(5) = nElems
      IF (HSize_proc(1).EQ.3) THEN
          SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='YES') ' | HDF5 state file containing HIT data with ',HSize_proc(1), &
                                                       ' variables and matching current setup. Continuing run...'
      ELSE
        CALL CollectiveStop(__STAMP__,'HDF5 state file containing HIT data with different number of variables!')
      END IF

      ALLOCATE(HIT_local(1:3,0:HSize(2)-1,0:HSize(3)-1,0:HSize(4)-1,nElems))
      CALL ReadArray('HIT',5,HSize_proc,OffsetElem,5,RealArray=HIT_local)
      ! No interpolation needed, read solution directly from file
      IF(.NOT. InterpolateSolution)THEN
        HIT_RMS = HIT_local
        DEALLOCATE(HIT_local)
      ELSE
        CALL CollectiveStop(__STAMP__,'Interpolation currently not supported for HIT test case!')
      END IF
      ! Indicate that we already initialized HIT_RMS from restartfile
      HIT_RMS_InitDone = .TRUE.
    ELSE
      SWRITE(UNIT_stdOut,'(A)')' | No restart variables for HIT test case found. Starting without u_RMS history!'
      HIT_RMS    = 0.
    END IF
    CALL CloseDataFile()
  ELSE
    HIT_RMS    = 0.
  END IF

  ! HIT_RMS required for restart, so add it as additional field data
  CALL AddToFieldData(FieldOut,(/3,PP_N+1,PP_N+1,PP_NZ+1,nElems/),'HIT',(/'HIT_RMSx','HIT_RMSy','HIT_RMSz'/), &
                      RealArray=HIT_RMS,doSeparateOutput=.TRUE.)
END IF ! HIT_Forcing

! Length of Buffer for TGV output
nWriteStats      = GETINT( 'nWriteStats'     ,'100')
nAnalyzeTestCase = GETINT( 'nAnalyzeTestCase','10')

IF(MPIRoot)THEN
  ALLOCATE(Time(nWriteStats))
  ALLOCATE(writeBuf(nHITVars,nWriteStats))
  Filename = TRIM(ProjectName)//'_HITAnalysis'

  varnames(1) ="Dissipation Rate Incompressible"
  varnames(2) ="Dissipation Rate Compressible"
  varnames(3) ="Ekin incomp"
  varnames(4) ="Ekin comp"
  varnames(5) ="U RMS"
  varnames(6) ="Reynolds number"
  varnames(7) ="A ILF"
  CALL InitOutputToFile(FileName,'Homogeneous Isotropic Turbulence Analysis Data',nHITVars,varnames)
END IF

! Spatial averaging according to de Laage de Meux, 2015
HIT_Avg = GETLOGICAL('HIT_Avg','.TRUE.')
HIT_1st = GETLOGICAL('HIT_1st','.FALSE.')

SWRITE(UNIT_stdOut,'(A)')' INIT TESTCASE HOMOGENEOUS ISOTROPIC TURBULENCE DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitTestcase


!==================================================================================================================================
!> Specifies all the initial conditions.
!==================================================================================================================================
SUBROUTINE ExactFuncTestcase(tIn,x,Resu,Resu_t,Resu_tt)
! MODULES
USE MOD_EOS_Vars,       ONLY: kappa
USE MOD_EOS,            ONLY: PrimToCons
USE MOD_TestCase_Vars
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
REAL                            :: A,Ms,prim(PP_nVarPrim)
!==================================================================================================================================
!> Use initialization of TGV for now
!> TODO: Find good initial values here
!A=SQRT(2./3.*HIT_k)                ! magnitude of speed
A=1.                               ! magnitude of speed
Ms=0.1                             ! maximum Mach number
prim(1)=HIT_rho
prim(2)= A*SIN(x(1))*COS(x(2))*COS(x(3))
prim(3)=-A*COS(x(1))*SIN(x(2))*COS(x(3))
prim(4)=0.
prim(5)=(A/Ms*A/Ms/Kappa*prim(1))  ! scaling to get Ms
prim(5)=prim(5)+1./16.*A*A*prim(1)*(COS(2*x(1))*COS(2.*x(3)) + 2.*COS(2.*x(2)) +2.*COS(2.*x(1)) +COS(2*x(2))*COS(2.*x(3)))
prim(6)=0. ! T does not matter for prim to cons

! Convert to conservative variables and set first,second time derivative to zero
CALL PrimToCons(prim,Resu)
Resu_t =0.
Resu_tt=0.
END SUBROUTINE ExactFuncTestcase


!==================================================================================================================================
!> Add testcases source term to solution time derivative. The implemented forcing approach is based on the following literature:
!> de Laage de Meux mode, case 'D. Restriction to the isotropic case and link with the linear forcing of Lundgren'
!> de Laage de Meux, B.; Audebert, B.; Manceau, R. and Perrin, R., "Anisotropic linear forcing for synthetic turbulence generation
!> in large eddy simulation and hybrid RANS/LES modeling." Physics of Fluids, 27 (2015), 035115
!==================================================================================================================================
SUBROUTINE TestcaseSource(Ut)
! MODULES
USE MOD_Globals
#if PP_N==N
USE MOD_PreProc,        ONLY: N
#endif
USE MOD_Analyze_Vars,   ONLY: wGPVol,ElemVol,Vol
USE MOD_DG_Vars,        ONLY: U,UPrim
USE MOD_Mesh_Vars,      ONLY: nElems,sJ
USE MOD_TestCase_Vars
USE MOD_TimeDisc_Vars,  ONLY: dt,CurrentStage
#if USE_MPI
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iElem,i,j,k
REAL            :: fac,TKE
!==================================================================================================================================
! Small trick to survive the first DG call(s) when dt is not yet set
IF (dt.GE.HUGE(1.)) RETURN

! Return if no HIT forcing required
IF (.NOT.HIT_Forcing) RETURN

! Initialize HIT_RMS with initial solution if it was not read from the restart file
IF (.NOT. HIT_RMS_InitDone) THEN
  HIT_RMS = UPrim(VELV,:,:,:,:)**2
  HIT_RMS_InitDone = .TRUE.
END IF

! Time-average the solution to obtain RMS
! TODO: PROUD FILTER
fac = dt/HIT_tFilter
HIT_RMS(1:3,:,:,:,:) = HIT_RMS(1:3,:,:,:,:) + ((UPrim(VELV,:,:,:,:))**2 - HIT_RMS(1:3,:,:,:,:)) * fac

! Calculate scalar kinetic Energy TKE either globally or elementwise by Gaussian integration.
! Then apply computed forcing coefficient to Ut.
IF (HIT_Avg) THEN

  ! Only update forcing coefficient during first RK stage if required
  IF ( (CurrentStage.EQ.1) .OR. (.NOT.HIT_1st) ) THEN

    ! Integrate kinetic energy TKE=1/2*U^2 over domain. Factor 1/2 is applied after integration
    TKE = 0.
    DO iElem=1,nElems
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        TKE = TKE + wGPVol(i,j,k)/sJ(i,j,k,iElem,0)*SUM(HIT_RMS(1:3,i,j,k,iElem))
      END DO
    END DO; END DO; END DO
    TKE = TKE/(2.*Vol) ! Normalize by volume and account for factor 1/2

#if USE_MPI
    ! Sum TKE over all procs and communicate the summed global TKE to all procs
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,TKE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
#endif

    ! Compute forcing coefficient A_ILF
    A_ILF = 1./(2.*HIT_tauRMS)*(HIT_k/TKE - 1.)
  END IF

  ! Apply the forcing to the time derivative
  Ut(MOMV,:,:,:,:) = Ut(MOMV,:,:,:,:) + A_ILF*U(MOMV,:,:,:,:)

ELSE ! HIT_Avg

  A_ILF = 0.
  ! Apply the forcing to every cell
  DO iElem=1,nElems
    TKE = 0.
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      TKE = TKE + wGPVol(i,j,k)/sJ(i,j,k,iElem,0)*SUM(HIT_RMS(1:3,i,j,k,iElem))
    END DO; END DO; END DO
    TKE = TKE/(2.*ElemVol(iElem))

    ! Apply forcing to the time derivative
    Ut(MOMV,:,:,:,iElem) = Ut(MOMV,:,:,:,iElem) + 1./(2.*HIT_tauRMS)*(HIT_k/TKE - 1.)*U(MOMV,:,:,:,iElem)

    ! Integrate the forcing coefficient
    A_ILF = A_ILF + 1./(2.*HIT_tauRMS)*(HIT_k/TKE - 1.)*ElemVol(iElem)
  END DO

END IF ! HIT_Avg

END SUBROUTINE TestcaseSource


!==================================================================================================================================
!> Perform HIT-specific analysis: compute dissipation rates, kinetic energy and Reynolds number
!==================================================================================================================================
SUBROUTINE AnalyzeTestcase(t,doFlush)
! MODULES
USE MOD_Globals
#if PP_N==N
USE MOD_PreProc,        ONLY: N
#endif
USE MOD_Analyze_Vars,   ONLY: wGPVol,Vol
USE MOD_DG_Vars,        ONLY: U
USE MOD_EOS_Vars,       ONLY: mu0,KappaM1
USE MOD_Lifting_Vars,   ONLY: GradUx,GradUy,GradUz
USE MOD_Mesh_Vars,      ONLY: nElems,sJ
USE MOD_TestCase_Vars
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
INTEGER                         :: i,j,k,p,q,iElem
REAL                            :: Vel(1:3),Pressure,rho,uRMS,nu0
REAL                            :: GradVel(1:3,1:3)
REAL                            :: E,Ekin                        ! Integrands: 0.5*v*v
REAL                            :: E_comp,Ekin_comp              ! Integrands: 0.5*rho*v*v
REAL                            :: Intfactor                     ! Integration weights
REAL                            :: S(1:3,1:3)                    ! Strain rate tensor S (symmetric)
REAL                            :: Sd(1:3,1:3)                   ! Deviatoric part of the strain rate tensor S
REAL                            :: divU                          ! Divergence of velocity vector
REAL                            :: eps3                          ! Integrand: p*(div u)
REAL                            :: u_tens, s_tens, sd_tens       ! Matrix product, integrands of Gradvel, S, Sd
REAL                            :: DR_S,DR_Sd,DR_p !,DR_u        ! Contributions to dissipation rate
REAL                            :: Reynolds,lTaylor              ! Reynolds number and Taylor microscale
REAL                            :: A_ILF_Glob                    ! Global average of Forcing Coefficient
#if USE_MPI
REAL                            :: rho_glob
REAL                            :: Ekin_glob,Ekin_comp_glob
REAL                            :: DR_S_glob,DR_Sd_glob,DR_p_glob !,DR_u_Glob
#endif
!==================================================================================================================================
Ekin      = 0.
Ekin_comp = 0.
rho       = 0.
!DR_u      = 0.
DR_S      = 0.
DR_Sd     = 0.
DR_p      = 0.

DO iElem=1,nElems
  DO k = 0,PP_N
    DO j = 0,PP_N
      DO i = 0,PP_N
      ! Compute primitive gradients (of u,v,w) at each GP
      Vel    (1:3) = U(MOMV,i,j,k,iElem)/U(DENS,i,j,k,iElem)
      GradVel(:,1) = GradUx(LIFT_VELV,i,j,k,iElem)
      GradVel(:,2) = GradUy(LIFT_VELV,i,j,k,iElem)
      GradVel(:,3) = GradUz(LIFT_VELV,i,j,k,iElem)

      ! Pressure and density
      Pressure = KappaM1*(U(ENER,i,j,k,iElem)-0.5*SUM(U(MOMV,i,j,k,iElem)*Vel(1:3)))
      rho      = rho + U(DENS,i,j,k,iElem)

      ! Compute divergence of velocity
      divU = GradVel(1,1) + GradVel(2,2) + GradVel(3,3)
      ! Compute tensor of velocity gradients
      S  = 0.5*(Gradvel+TRANSPOSE(GradVel))
      ! Deviatoric part of strain tensor
      Sd = S
      DO p = 1,3
        Sd(p,p) = Sd(p,p) - 1./3.*divU
      END DO

      ! Compute kinetic energy integrand (incompressible/compressible)
      E      = 0.5 * SUM(Vel(1:3)*Vel(1:3))
      E_comp = U(DENS,i,j,k,iElem) * E

      ! Compute integrand for epsilon3, pressure contribution to dissipation (compressiblity effect)
      eps3 = Pressure * divU
      ! Matrix product for velocity gradient tensor, S:S and Sd:Sd
      u_tens=0.; s_tens=0.; sd_tens=0.
      DO p=1,3
        DO q=1,3
          u_tens  = u_tens  + GradVel(p,q)*GradVel(p,q)
          s_tens  = s_tens  + S(p,q)      *S(p,q)
          sd_tens = sd_tens + Sd(p,q)     *Sd(p,q)
        END DO
      END DO

      ! Compute integration factor
      Intfactor = wGPVol(i,j,k) / sJ(i,j,k,iElem,0)

      ! Compute integrals
      !>> Kinetic energy (incompressible/compressible)
      Ekin      = Ekin      + E     *Intfactor
      Ekin_comp = Ekin_comp + E_comp*IntFactor
      !>> Dissipation rate
  !    DR_u  = DR_u  + u_tens *IntFactor  ! epsilon incomp from velocity gradient tensor (Diss Fauconnier)
      DR_S  = DR_S  + 2./U(1,i,j,k,iElem)*S_tens *IntFactor  ! epsilon incomp from strain rate tensor (incomp) (Sagaut)
      DR_SD = DR_SD + 2./U(1,i,j,k,iElem)*sd_tens*Intfactor  ! epsilon 1 from deviatoric part of strain rate tensor Sd (compressible)
      DR_p  = DR_p  - 1./U(1,i,j,k,iElem)*eps3   *Intfactor  ! epsilon 3 from pressure times div u (compressible)
      END DO
    END DO
  END DO
END DO

#if USE_MPI
! MPI case: Globalize analyze variables
CALL MPI_REDUCE(rho      ,rho_Glob        ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(Ekin     ,Ekin_Glob       ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(Ekin_comp,Ekin_comp_Glob  ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
!CALL MPI_REDUCE(DR_u     ,DR_u_Glob       ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(DR_S     ,DR_S_Glob       ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(DR_Sd    ,DR_Sd_Glob      ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
CALL MPI_REDUCE(DR_p     ,DR_p_Glob       ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
IF (HIT_Avg) THEN ! If HIT_Avg = .TRUE., A_ILF is already averaged across all MPI ranks.
  A_ILF_Glob = A_ILF
ELSE ! If A_ILF is computed element-wise, we have to average across all MPI ranks
  CALL MPI_REDUCE(A_ILF,A_ILF_Glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF

! Overwrite the variables only on root
IF (MPIRoot) THEN
  rho       = rho_Glob
  Ekin      = Ekin_Glob
  Ekin_comp = Ekin_comp_Glob
!  DR_u      = DR_u_Glob
  DR_S      = DR_S_Glob
  DR_Sd     = DR_SD_Glob
  DR_p      = DR_p_Glob
ELSE
  ! Only root continues from here
  RETURN
END IF
#else
A_ILF_Glob = A_ILF
#endif

! Normalize the integrals with the volume and calculate turbulent quantities. The density is already accounted for in the integrands
rho       = rho      /Vol
Ekin      = Ekin     /Vol
Ekin_comp = Ekin_comp/Vol
uRMS      = SQRT(Ekin*2./3.)
!DR_u      = DR_u *mu0/Vol
DR_S      = DR_S *mu0/Vol
DR_SD     = DR_SD*mu0/Vol
DR_p      = DR_p     /Vol
IF (.NOT. HIT_Avg) A_ILF_Glob = A_ILF_Glob/Vol ! Do NOT overwrite A_ILF since this is a global variable!

! Taylor microscale Reynolds number (Petersen, 2010)
nu0      = mu0/rho
!>> Taylor microscale should be calculated with compressible eps, but currently unstable. Use only deviatoric part eps1 for now
!lTaylor  = SQRT(15.*nu0/(DR_SD+DR_p))*uRMS
!>> Taylor microscale calculated using incompressible epsilon (Hillewaert, 2012)
lTaylor  = SQRT(15.*nu0/DR_S)*uRMS
Reynolds = uRMS * lTaylor/nu0

! Increment output counter and fill output array at every time step
ioCounter        = ioCounter + 1
Time(ioCounter)  = t
writeBuf(1:nHITVars,ioCounter) = (/DR_S,DR_Sd+DR_p,Ekin,Ekin_comp,uRMS,Reynolds,A_ILF_Glob/)

! Perform output
IF(ioCounter.EQ.nWriteStats .OR. doFlush)THEN
  CALL WriteStats()
  ioCounter = 0
END IF

END SUBROUTINE AnalyzeTestcase


!==================================================================================================================================
!> Write HIT Analysis Data to File
!==================================================================================================================================
SUBROUTINE WriteStats()
! MODULES
USE MOD_TestCase_Vars
USE MOD_Output,        ONLY: OutputToFile
IMPLICIT NONE
!==================================================================================================================================
CALL OutputToFile(FileName,Time(1:ioCounter),(/nHITVars,ioCounter/),&
                           RESHAPE(writeBuf(:,1:ioCounter),(/nHITVars*ioCounter/)))
END SUBROUTINE WriteStats


!==================================================================================================================================
!> Finalizes all the test case variables
!==================================================================================================================================
SUBROUTINE FinalizeTestcase()
! MODULES
USE MOD_Globals
USE MOD_TestCase_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF(MPIroot) THEN
  SDEALLOCATE(Time)
  SDEALLOCATE(writeBuf)
END IF
SDEALLOCATE(HIT_RMS)
END SUBROUTINE

!==================================================================================================================================
!> Empty placeholder routine
!==================================================================================================================================
SUBROUTINE DO_NOTHING(optionalREAL,optionalREAL2)
IMPLICIT NONE
REAL,OPTIONAL,INTENT(IN)  :: optionalREAL,optionalREAL2
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

END MODULE MOD_Testcase
