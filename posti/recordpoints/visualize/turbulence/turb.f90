#include "flexi.h"
!===================================================================================================================================
!> Module to handle the Recordpoints
!===================================================================================================================================
MODULE MOD_Turbulence
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Turbulence
  MODULE PROCEDURE Turbulence
END INTERFACE

PUBLIC :: Turbulence
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize all necessary information to perform interpolation
!===================================================================================================================================
SUBROUTINE Turbulence()
! MODULES
USE MOD_Globals              
USE MOD_VarnamemappingsRP_vars 
USE MOD_RPData_Vars            ,ONLY: RPTime,RPData
USE MOD_RPSetVisuVisu_Vars     ,ONLY: nRP_global
USE MOD_OutputRPVisu_Vars      ,ONLY: nSamples_out
USE MOD_ParametersVisu         ,ONLY: Mu0,cutoffFreq
USE FFTW3
#ifdef WITHTECPLOT
USE MOD_OutputRPVisu_Vars      ,ONLY: CoordNames
USE MOD_Tecplot                ,ONLY: WriteDataToTecplotBinary,WriteTimeAvgDataToTecplotBinary
USE MOD_ParametersVisu         ,ONLY: ProjectName
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample , iVar,iRP
INTEGER                         :: nSamples_spec
INTEGER(KIND=8)                 :: plan  
COMPLEX                         :: in(nSamples_out),out(nSamples_out) 
REAL                            :: dt_equi , PI , df 
REAL , ALLOCATABLE              :: E_kineticSpec(:,:)
REAL , ALLOCATABLE              :: velAbs(:,:) , velAbs_avg(:) , density_avg(:) , kk(:,:) , disRate(:,:) , epsilonMean(:,:)
REAL , ALLOCATABLE              :: nu0(:) , eta(:) , etaK(:,:)
REAL , ALLOCATABLE              :: vel_spec(:,:,:) , velPrim(:,:,:) , RPData_freq(:)
#ifdef WITHTECPLOT
INTEGER                         :: nVar_turb
CHARACTER(LEN=255),ALLOCATABLE  :: VarNameTurb(:)
CHARACTER(LEN=255)              :: Filename,stroutputfile
REAL , ALLOCATABLE              :: RPData_turb(:,:,:) , RPData_turb_avg(:,:)
#endif
!===================================================================================================================================


PI=ATAN(1.)*4.
ALLOCATE(velPrim(1:3,nRP_global,nSamples_out))

DO iSample=1,nSamples_out
  DO iRP=1,nRP_global
    velPrim(1:3,iRP,iSample) = RPData(2:4,iRP,iSample)/RPData(1,iRP,iSample)
  END DO
END DO

WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')' PERFORMING FFT for turbulence spectra...'
df = 1./(RPTime(nSamples_out)-RPTime(1))
WRITE(UNIT_stdOut,'(A,I8)') '         Samples available for FFT:',nSamples_out
WRITE(UNIT_stdOut,'(A,F16.7)') '           Frequency resolution:',df
nSamples_spec=INT((nSamples_out-1)/2)+1
IF(cutoffFreq.NE.-999)THEN
  WRITE(UNIT_stdOut,'(A,F16.7)') '  User defined Cutoff Frequency:',cutoffFreq
  nSamples_spec=NINT(cutoffFreq/dF)
END IF  
WRITE(UNIT_stdOut,'(A,I8)') '    No. spectrum output samples FFT:',nSamples_spec
ALLOCATE(vel_spec(1:3,nRP_global,nSamples_spec))
ALLOCATE(RPData_freq(nSamples_spec))
DO iSample=1,nSamples_spec
  RPData_freq(iSample)=(iSample-1)*df
END DO

WRITE(UNIT_stdOut,'(A,F16.7)') '            Nyquist frequency:',df*REAL(nSamples_out-1)
WRITE(UNIT_stdOut,'(A,F16.7)') '      Max. resolved frequency:',RPData_freq(nSamples_spec)
CALL DFFTW_PLAN_DFT_1D(plan,nSamples_out,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
DO iRP=1,nRP_global
  WRITE(UNIT_stdOut,*)'   Processing RP ',iRP,' of ',nRP_global
  DO iVar=1,3
    in(:)= velPrim(iVar,iRP,:)
    CALL DFFTW_EXECUTE_DFT(plan, in, out)
    vel_spec(iVar,iRP,:)=2./REAL(nSamples_out)*ABS(out(1:nSamples_spec))
    vel_spec(iVar,iRP,1)=0.5*vel_spec(iVar,iRP,1)
  END DO ! iVar
END DO   ! iRP


ALLOCATE(velAbs(nRP_global,nSamples_out))
ALLOCATE(kk(nRP_global,nSamples_spec))
ALLOCATE(E_kineticSpec(nRP_global,nSamples_spec))
ALLOCATE(epsilonMean(nRP_global,nSamples_spec))
ALLOCATE(disRate(nRP_global,nSamples_spec))
ALLOCATE(etaK(nRP_global,nSamples_spec))
ALLOCATE(velAbs_avg(nRP_global))
ALLOCATE(density_avg(nRP_global))
ALLOCATE(nu0(nRP_global))
ALLOCATE(eta(nRP_global))
 
E_kineticSpec = 0.5*(vel_spec(1,:,:)**2 + vel_spec(2,:,:)**2 + vel_spec(3,:,:)**2 )
!write(*,*)'E_kin',E_kineticSpec
!RPData_spec(Prim%IndGlobal(13),:,:)=E_kineticSpec(:,:)
!uAvg=(RPDataTimeAvg_out(Prim%IndGlobal(1))**2+RPDataTimeAvg_out(Prim%IndGlobal(2))+RPDataTimeAvg_out(Prim%IndGlobal(3)))

!Calculate mean transport velocity, density
dt_equi = (RPTime(nSamples_out)-RPTime(1))/nSamples_out 
velAbs=SQRT(velPrim(1,:,:)**2+velPrim(2,:,:)**2+velPrim(3,:,:)**2)  
velAbs_Avg=0.
density_avg = 0.
DO iSample=2,nSamples_out
  velAbs_avg  = velAbs_avg  + 0.5*dt_equi*(velAbs(:,iSample)+velAbs(:,iSample-1))
  density_avg = density_avg + 0.5*dt_equi*(RPData(1,:,iSample)+RPData(1,:,iSample-1))                                
END DO 
velAbs_Avg  = velAbs_Avg /(RPTime(nSamples_out)-RPTime(1))
density_Avg = density_Avg/(RPTime(nSamples_out)-RPTime(1))

!calculate wavenumber for each RP, using mean tranposrt velocity as wavespeed
DO iSample=1,nSamples_Spec
  kk(:,iSample) = 2*PI*RPData_freq(iSample)/velAbs_avg
END DO
!epsilon as far as resolved for run for each RP, good resolutin for EtaK GE 1,
!overprediction of Eta for underresolved turbulence...
!Dissipation Rate
nu0=Mu0/density_avg

DO iSample=1,nSamples_spec
  disRate(:,iSample) = E_kineticSpec(:,iSample)*kk(:,iSample)**2*2*nu0(:)
END DO 

epsilonMean=0.
DO iSample=2,nSamples_spec
  epsilonMean(:,iSample)=epsilonMean(:,iSample-1) + 0.5*(kk(:,iSample)-kk(:,iSample-1))*(disRate(:,iSample)+disRate(:,iSample-1))
END DO

eta = nu0**(3./4.)/epsilonMean(:,nSamples_spec)**(1./4.)

DO iSample=1,nSamples_spec
  etaK(:,iSample) = eta(:)*kk(:,iSample)
END DO


#ifdef WITHTECPLOT
nVar_turb=5
ALLOCATE(VarNameTurb(nVar_turb))
VarNameTurb(1) = 'KineticEnergy'
VarNameTurb(2) = 'DissipationRate'
VarNameTurb(3) = 'MeanDissipation_cumulated'
VarNameTurb(4) = 'Eta_K'
VarNameTurb(5) = 'Wavenumber K'

ALLOCATE(RPData_turb(1:nVar_turb,nRP_global,nSamples_spec))
RPData_turb(1,:,:) = E_kineticSpec(:,:) 
RPData_turb(2,:,:) = disRate(:,:)
RPData_turb(3,:,:) = epsilonMean(:,:)
RPData_turb(4,:,:) = etaK(:,:)
RPData_turb(5,:,:) = kk(:,:)

! Output Turbulence stuff
CoordNames(1)='Frequency'
WRITE(UNIT_StdOut,'(132("-"))')
Filename=TRIM(ProjectName)
FileName=TRIM(FileName)//'_RP_turb'
strOutputFile=TRIM(FileName)//'.plt'
WRITE(UNIT_stdOut,'(A,A)')' WRITING TURBULENCE DATA TO ',strOutputFile
CALL WriteDataToTecplotBinary(nSamples_spec,nRP_global,nVar_turb,VarNameTurb,RPData_freq,RPData_turb,strOutputFile)
WRITE(UNIT_stdOut,'(A)')' DONE.'
WRITE(UNIT_StdOut,'(132("-"))')


DEALLOCATE(VarNameTurb)
nVar_turb=2
ALLOCATE(VarNameTurb(nVar_turb))
VarNameTurb(1) = 'Eta'
VarNameTurb(2) = 'MeanDissipation'

ALLOCATE(RPData_turb_avg(1:nVar_turb,nRP_global))
RPData_turb_avg(1,:) = eta
RPData_turb_avg(2,:) = epsilonMean(:,nSamples_spec) 

! Output Time Average
WRITE(UNIT_StdOut,'(132("-"))')
Filename=TRIM(ProjectName)
FileName=TRIM(FileName)//'_RP_TurbAvg'
strOutputFile=TRIM(FileName)//'.plt'
WRITE(UNIT_stdOut,'(A,A)')' WRITING TURBULENCE AVERAGE DATA TO ',strOutputFile
CALL WriteTimeAvgDataToTecplotBinary(nRP_global,nVar_turb,VarNameTurb,RPData_turb_avg,strOutputFile)
WRITE(UNIT_stdOut,'(A)')' DONE.'
WRITE(UNIT_StdOut,'(132("-"))')
#endif


END SUBROUTINE Turbulence

END MODULE MOD_Turbulence
