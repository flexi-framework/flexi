#include "flexi.h"
!===================================================================================================================================
!> Module to handle the Recordpoints
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
!> Initialize all necessary information to perform interpolation
!===================================================================================================================================
SUBROUTINE FilterRP()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars          ,ONLY: RPTime
USE MOD_RPSetVisuVisu_Vars           ,ONLY: nRP_global
USE MOD_RPInterpolation_Vars
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out,RPData_out
USE MOD_ParametersVisu           ,ONLY: FilterWidth,nVarVisu,FilterMode
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
