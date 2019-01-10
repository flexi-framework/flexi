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
!===================================================================================================================================
!> Module to handle temporal filtering of the record point time signal
!===================================================================================================================================
MODULE MOD_FilterRP
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE FilterRP
  MODULE PROCEDURE FilterRP
END INTERFACE

PUBLIC :: FilterRP
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize and perform the temporal filtering of the time signal stored at the record points. Either a high-pass or a
!> low-pass filter can be used.
!===================================================================================================================================
SUBROUTINE FilterRP()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars          ,ONLY: RPTime
USE MOD_RPSetVisuVisu_Vars   ,ONLY: nRP_global
USE MOD_RPInterpolation_Vars
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out,RPData_out
USE MOD_ParametersVisu       ,ONLY: FilterWidth,nVarVisu,FilterMode
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample,iSample2,iRP
INTEGER                         :: dnSamples,iSample_left,iSample_right
REAL                            :: RPData_tmp(nVarVisu,nSamples_out),RPData_f(nVarVisu)
REAL                            :: snSamples_block
REAL                            :: pi
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')

WRITE(UNIT_stdOut,'(A)')' APPLYING FILTER ...'
IF(FilterMode.EQ.0) THEN
  WRITE(UNIT_stdOut,'(A)')' Filter Mode: Low Pass'
ELSE
  WRITE(UNIT_stdOut,'(A)')' Filter Mode: High Pass'
END IF
dnSamples=FLOOR(nSamples_out*0.5*FilterWidth/(RPTime(nSamples_out)-RPTime(1)))

WRITE(UNIT_stdOut,'(A,F16.7)')' Filter Width: ',FilterWidth
WRITE(UNIT_stdOut,'(A,I8)')   ' nSamples_out: ',nSamples_out
WRITE(UNIT_stdOut,'(A,I8)')   ' dnSamples:    ',dnSamples
pi=acos(-1.)
!DO iSample=1,nSamples_out
!    RPData_out(:,:,iSample) = cos(10.*pi*2.*REAL(iSample-1)/REAL(nSamples_out-1)) &
!                             +0.05*cos(100.*pi*2.*REAL(iSample-1)/REAL(nSamples_out-1)+.1)
!END DO
DO iRP=1,nRP_global 
  RPData_tmp=RPData_out(:,iRP,:)
  ! Average around each Sample
  DO iSample=1,nSamples_out
    ! Average around each Sample in a specified range
    iSample_left = max(iSample-dnSamples,1)
    iSample_right= min(iSample+dnSamples,nSamples_out)
    RPData_f=0
    DO iSample2=iSample_left,iSample_right
      RPData_f=RPData_f+RPData_out(:,iRP,iSample2)
    END DO ! iSample2
    snSamples_block=1/(REAL(MAX(iSample_right-iSample_left+1,1)))
    ! low pass filtered signal
    RPData_tmp(:,iSample)=RPData_f*snSamples_block 
  END DO ! iSample
  IF(FilterMode.EQ.1) THEN
    RPData_out(:,iRP,:)=RPData_out(:,iRP,:)-RPData_tmp(:,:)
  ELSE
    RPData_out(:,iRP,:)=RPData_tmp(:,:)
  END IF
END DO ! iRP=1,nRP 
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE FilterRP

END MODULE MOD_FilterRP
