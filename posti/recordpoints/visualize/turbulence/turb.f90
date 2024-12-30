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
!===================================================================================================================================
!> Module containing routines to calculate turbulent quantities using a temporal FFT
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
!> This routine performs a temporal FFT for all RPs. Using the result, turbulent quantities like the turbulent kinetic energy
!> will be calculated and a spectrum will be written.
!===================================================================================================================================
SUBROUTINE Turbulence()
! MODULES
USE MOD_Globals
USE MOD_ParametersVisu         ,ONLY: Mu0,cutoffFreq
USE MOD_OutputRPVisu_Vars      ,ONLY: nSamples_out
USE MOD_RPData_Vars            ,ONLY: RPTime,RPData
USE MOD_RPSetVisuVisu_Vars     ,ONLY: nRP_global
USE MOD_Turbulence_Vars
USE FFTW3
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample , iVar,iRP
INTEGER                         :: nSamples_spec
INTEGER(DP)                     :: plan
COMPLEX                         :: in(nSamples_out),out(nSamples_out)
REAL                            :: dt_equi , PI , df
REAL,ALLOCATABLE                :: vel_spec(:,:,:),velPrim(:,:,:)
REAL,ALLOCATABLE                :: velAbs(:,:),velAbs_avg(:),density_avg(:)
REAL,ALLOCATABLE                :: nu0(:)
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
ALLOCATE(RPData_freqTurb(nSamples_spec))
DO iSample=1,nSamples_spec
  RPData_freqTurb(iSample)=(iSample-1)*df
END DO

WRITE(UNIT_stdOut,'(A,F16.7)') '            Nyquist frequency:',df*REAL(nSamples_out-1)
WRITE(UNIT_stdOut,'(A,F16.7)') '      Max. resolved frequency:',RPData_freqTurb(nSamples_spec)
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

!calculate wavenumber for each RP, using mean transport velocity as wavespeed
DO iSample=1,nSamples_Spec
  kk(:,iSample) = 2*PI*RPData_freqTurb(iSample)/velAbs_avg
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

WRITE(UNIT_stdOut,'(A)')' DONE.'
WRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE Turbulence

END MODULE MOD_Turbulence
