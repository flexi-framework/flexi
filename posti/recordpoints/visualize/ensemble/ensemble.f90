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
!> Module to handle ensemble averaging of the record point time signal
!===================================================================================================================================
MODULE MOD_EnsembleRP
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: EnsembleRP
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
! Initialize and perform ensemble averaging of the time signal stored at the record points.
!===================================================================================================================================
SUBROUTINE EnsembleRP()
! MODULES
USE MOD_Globals
USE MOD_EnsembleRP_Vars      ,ONLY: enSamples,nVar_ensTurb,RPData_ens,RPData_turb,RPData_freqEns
USE MOD_RPInterpolation_Vars ,ONLY: dt_out
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out,RPData_out
USE MOD_ParametersVisu       ,ONLY: EnsemblePeriod!,Kappa
USE MOD_RPData_Vars          ,ONLY: RPTime
USE MOD_ParametersVisu       ,ONLY: nVarVisu
USE MOD_RPSetVisuVisu_Vars   ,ONLY: nRP_global
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: i,iSample,iEns
INTEGER,ALLOCATABLE             :: counter(:)
! REAL                            :: KappaM1
REAL                            :: v_ens(1:3,nRP_global),v_ensTurb(nRP_global)
REAL                            :: v_inf(nRP_global)!,UPrim(1:5)
REAL,ALLOCATABLE                :: delta(:,:)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')' PERFORMING ENSEMBLE AVERAGING ...'

! Abort if the frequency is too low for the period
!IF(EnsemblePeriod.LT.(1./EnsembleFreq)) &
!  CALL Abort(__STAMP__,' EnsemblePeriod less than 1./EnsembleFreq')

! Abort if we do not have a full period
IF (EnsemblePeriod.GT.(RPTime(nSamples_out)-RPTime(1))) &
  CALL Abort(__STAMP__,' EnsemblePeriod longer than available RP data!')

! Number of turbulent variables. Hardcoded for now
nVar_ensTurb = 6
! KappaM1      = Kappa-1.

enSamples = CEILING(EnsemblePeriod/dt_out)
ALLOCATE(RPData_freqEns(enSamples))
DO iSample = 1,enSamples
  RPData_freqEns(iSample) = (iSample-1)*dt_out
END DO

ALLOCATE(RPData_ens (1:nVarVisu    ,nRP_global,enSamples))
ALLOCATE(RPData_turb(1:nVar_ensTurb,nRP_global,enSamples))
ALLOCATE(counter    (1:                        enSamples))
RPData_ens  = 0.
RPData_turb = 0.

WRITE(UNIT_stdOut,'(A,I8,A,E9.2)') ' Number of timesteps for ensemble avg.:   ', enSamples, ', ens. timestep:', dt_out

! Ensemble averaged quantities
ALLOCATE(delta(1:nVarVisu,nRP_global))
delta   = 0.
counter = 0
DO iSample = 1,nSamples_out !nSamples_global
  ! Get interval for ensemble averaging
  iEns = FLOOR(MOD(RPTime(iSample),EnsemblePeriod)/dt_out) + 1

  ! Welford's algorithm for averaging
  counter(iEns)        = counter(iEns) + 1
  delta(     :,:     ) = RPData_out(:,:,iSample) - RPData_ens(:,:,iEns)
  RPData_ens(:,:,iEns) = RPData_ens(:,:,iEns   ) + delta(:,:)/counter(iEns)
END DO
DEALLOCATE(delta)
WRITE(UNIT_stdOut,'(A,I8,A,I10)')   ' > Done ensemble  average. Min. samples:  ', MINVAL(counter), ', max. samples:', MAXVAL(counter)

! Turbulent quantities on the ensemble average
ALLOCATE(delta(1:nVar_ensTurb,nRP_global))
delta   = 0.
counter = 0
DO iSample = 1,nSamples_out !nSamples_global
  ! Get interval for ensemble averaging
  iEns = FLOOR(MOD(RPTime(iSample),EnsemblePeriod)/dt_out) + 1

  ! Pressure and kinetic energy compared to ensemble average
  ! Transform to primitive variables. Can't use the function because it assumes the wrong momentum
!  DO iRP=1,nRP_Global
!    sRho       = 1./U(1,i,j,k,RP_ElemID(iRP))
!    UPrim(1)   = U(1,i,j,k,RP_ElemID(iRP))
!    UPrim(2:4) = U(2:4,i,j,k,RP_ElemID(iRP))*sRho
!    UPrim(5)   = KappaM1*(U(5,  i,j,k,RP_ElemID(iRP)) - 0.5*SUM(                                     &
!                       (U(2:4,i,j,k,RP_ElemID(iRP)) - UPrim(1)*Elem_vGP(1:3,i,j,k,RP_ElemID(iRP))) &
!                        *UPrim(2:4)))
!    p_ens(iRP) = UPrim(5)
!  END DO

  DO i=1,3
    v_ens(i,:)          = (RPData_out(i+1,:,iSample)/RPData_out(1,:,iSample)) - (RPData_ens(i+1,:,iEns)/RPData_ens(1,:,iEns))
  END DO
  v_inf    (:)          = SQRT((RPData_ens(2,:,iEns)**2. + RPData_ens(3,:,iEns)**2. + RPData_ens(4,:,iEns)**2.))/RPData_ens(1,:,iEns)
  v_ensTurb(:)          =  v_ens(1,:)*v_ens(1,:) + v_ens(2,:)*v_ens(2,:) + v_ens(3,:)*v_ens(3,:)

  ! Welford's algorithm for averaging
  counter(iEns)         = counter(iEns) + 1
  delta(      1,:     ) = RPData_ens(1,:,iEns)*v_ensTurb(:)/2. - RPData_turb(1,:,iEns) ! turbulent kinetic energy
  delta(      2,:     ) = SQRT(1./3.*v_ensTurb(:))/v_inf(:)    - RPData_turb(2,:,iEns) ! turbulent intensity
  delta(      3,:     ) = v_ensTurb(:)                         - RPData_turb(3,:,iEns) ! velocity [no R] MS
  delta(      4,:     ) = v_ens(1,:)**2.                       - RPData_turb(4,:,iEns) ! velocityX [no R] MS
  delta(      5,:     ) = v_ens(2,:)**2.                       - RPData_turb(5,:,iEns) ! velocityY [no R] MS
  delta(      6,:     ) = v_ens(3,:)**2.                       - RPData_turb(6,:,iEns) ! velocityZ [no R] MS
  RPData_turb(:,:,iEns) = RPData_turb(:,:,iEns) + delta(:,:)/counter(iEns)
END DO

! Calculate the square root of the mean squares
DO i = 3,6
  RPData_turb(i,:,:) = SQRT(RPData_turb(i,:,:))
END DO

WRITE(UNIT_stdOut,'(A,I8,A,I10)')   ' > Done turbulent average. Min. samples:  ', MINVAL(counter), ', max. samples:', MAXVAL(counter)

DEALLOCATE(delta)

END SUBROUTINE EnsembleRP

END MODULE MOD_EnsembleRP
