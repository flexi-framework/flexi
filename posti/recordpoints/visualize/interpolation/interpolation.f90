#include "flexi.h"
!===================================================================================================================================
!> Module to handle the Recordpoints
!===================================================================================================================================
MODULE MOD_RPinterpolation
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitInterpolation
  MODULE PROCEDURE InitInterpolation
END INTERFACE

INTERFACE InterpolateEquiTime
  MODULE PROCEDURE InterpolateEquiTime
END INTERFACE

INTERFACE CalcTimeAvg
  MODULE PROCEDURE CalcTimeAvg
END INTERFACE

INTERFACE CalcFluctuations
  MODULE PROCEDURE CalcFluctuations
END INTERFACE

INTERFACE FinalizeInterpolation
  MODULE PROCEDURE FinalizeInterpolation
END INTERFACE

PUBLIC :: InitInterpolation
PUBLIC :: InterpolateEquiTime
PUBLIC :: CalcTimeAvg             
PUBLIC :: CalcFluctuations        
PUBLIC :: FinalizeInterpolation
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize all necessary information to perform interpolation
!===================================================================================================================================
SUBROUTINE InitInterpolation()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars          ,ONLY: RPTime,nSamples_global
USE MOD_RPInterpolation_Vars
USE MOD_OutputRPVisu_Vars          ,ONLY: nSamples_out
USE MOD_Parameters           ,ONLY: equiTimeSpacing,OutputTimeAverage,doFluctuations,nBlocks,doFFT,fourthDeriv,doPSD
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample
INTEGER                         :: nSamples_min,nSamples_2n
REAL                            :: dtMin,factor
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')' INIT INTERPOLATION...'


ALLOCATE(dt(1:nSamples_global))
DO iSample=2,nSamples_global
 dt(iSample)=RPTime(iSample)-RPtime(iSample-1)
END DO
dt(1)=1.

IF(equiTimeSpacing) THEN
  ! get smallest time step
  dtMin=HUGE(1.)
  DO iSample=2,nSamples_global
   ! ignore small time steps due to state file output
   ! i.e. only allow slow timestep 'drift'
   IF(dt(iSample)/dt(iSample-1).GT.0.9) &
     dtMin=MIN(dt(iSample),dtMin)
  END DO
  ! get number of samples with dtMin*factor
  factor=1.
  nSamples_min=INT((RPTime(nSamples_global)-RPTime(1))/(dtMin*factor))+1
  nSamples_out=nSamples_min
  dt_out = (RPTime(nSamples_global)-RPTime(1))/REAL(nSamples_out-1)
  TEnd=RPTime(nSamples_global)
 
  IF(dt_out.LE.dtMin) THEN
    WRITE(UNIT_StdOut,'(A)')'WARNING: Equidistant sampling time step is lower than or equal to the input time step!!!'
  END IF

  WRITE(UNIT_stdOut,'(A)')'  Database:'
  WRITE(UNIT_StdOut,'(A18,ES16.7)')'         TStart = ',RPTIME(1)
  WRITE(UNIT_StdOut,'(A18,ES16.7)')'           TEnd = ',TEnd
  WRITE(UNIT_StdOut,'(A18,ES16.7)')'          dtMin = ',dtMin
  WRITE(UNIT_StdOut,'(A18,I8    )')'       nSamples = ',nSamples_global
  WRITE(UNIT_stdOut,'(A)')'  Equidistant time grid:'
  WRITE(UNIT_StdOut,'(A18,ES16.7)')'             dt = ',dt_out
  WRITE(UNIT_StdOut,'(A18,I8    )')'       nSamples = ',nSamples_out
END IF
WRITE(UNIT_stdOut,'(A)')' DONE.'
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitInterpolation


!===================================================================================================================================
!> Interpolate the RP data on equidistant time grid            
!===================================================================================================================================
SUBROUTINE InterpolateEquiTime()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars          ,ONLY: RPTime,nVar_HDF5,nSamples_global,RPData
USE MOD_RPSet_Vars           ,ONLY: nRP_global
USE MOD_RPInterpolation_Vars
USE MOD_OutputRPVisu_Vars          ,ONLY: nSamples_out
USE MOD_Basis                ,ONLY: ALMOSTEQUAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iSample,iSample_left,iSample_right
REAL                    :: Time,Time_left,Time_right
REAL                    :: RPTime_tmp(nSamples_global)
REAL                    :: RPData_tmp(nVar_HDF5,nRP_global,nSamples_global)
REAL                    :: RPData_right(nVar_HDF5,nRP_global),RPData_left(nVar_HDF5,nRP_global)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Linear interpolation on equidistant time grid...'

RPTime_tmp(:)     = RPTime(:)
RPData_tmp(:,:,:) = RPData(:,:,:)
DEALLOCATE(RPTime,RPData,dt)
ALLOCATE(RPTime(1:nSamples_out))
ALLOCATE(RPData(1:nVar_HDF5,1:nRP_global,1:nSamples_out))
ALLOCATE(dt(1:nSamples_out))
dt=dt_out

RPData(:,:,1)    = RPData_tmp(:,:,1)
RPTime(1)        = RPTime_tmp(1)
Time             = RPTime_tmp(1)
iSample_left =1
iSample_right=2
DO iSample=2,nSamples_out
  Time=RPTime_tmp(1)+dt_out*REAL(iSample-1)
  ! get two closest surrounding samples
  DO
    Time_right=RPTime_tmp(iSample_right)
    IF((Time_right.GE.Time).OR.ALMOSTEQUAL(Time_right,Time)) EXIT
    iSample_right=MIN(iSample_right+1,nSamples_global)
  END DO
  iSample_left=iSample_right-1
  DO
    Time_left=RPTime_tmp(iSample_left)
    IF(Time_left .LE. Time) EXIT
    iSample_left=iSample_left -1
  END DO
  ! linear interpolation 
  RPData_right(:,:)=RPData_tmp(:,:,iSample_right)
  RPData_left(:,:)=RPData_tmp(:,:,iSample_left)
  RPData(:,:,iSample) = RPData_left+(RPData_right-RPData_left)*(Time-Time_left)/(Time_right-Time_left)
  RPTime(iSample)=Time
END DO
WRITE(UNIT_stdOut,'(A)')'done.'
WRITE(UNIT_stdOut,'(132("-"))')
!WRITE(*,*)RPTime(7050)
END SUBROUTINE InterpolateEquiTime


!===================================================================================================================================
!> Calculate temporal TimeAvg values using 'raw' data (not interpolated on equi grid)
!===================================================================================================================================
SUBROUTINE CalcTimeAvg()
! MODULES
USE MOD_Globals
USE MOD_RPSet_Vars           ,ONLY: nRP_global
USE MOD_RPData_Vars          ,ONLY: RPTime
USE MOD_RPInterpolation_Vars
USE MOD_Parameters       ,ONLY: nVarVisu,EquiTimeSpacing
USE MOD_OutputRPVisu_Vars          ,ONLY: nSamples_out,RPData_out,RPDataTimeAvg_out
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iSample
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Calculating Time Average...'

ALLOCATE(RPDataTimeAvg_out(1:nVarVisu,nRP_global))
IF(nSamples_out.EQ.1) THEN
  RPDataTimeAvg_out(:,:)=RPData_out(:,:,1)
  WRITE(UNIT_stdOut,'(A)')'done.'
  WRITE(UNIT_stdOut,'(132("-"))')
  RETURN
END IF
RPDataTimeAvg_out(:,:)=0.
IF(EquiTimeSpacing) THEN ! use arithmetic mean
!  RPDataTimeAvg_out=SUM(RPData_out(:,:,1:nSamples_out))/nSamples_out
  DO iSample=1,nSamples_out
    RPDataTimeAvg_out = RPDataTimeAvg_out + RPData_out(:,:,iSample)
  END DO
  RPDataTimeAvg_out=RPDataTimeAvg_out/nSamples_out
ELSE
  DO iSample=2,nSamples_out
    RPDataTimeAvg_out = RPDataTimeAvg_out + 0.5*dt(iSample)*(RPData_out(:,:,iSample)+RPData_out(:,:,iSample-1))
  END DO
  RPDataTimeAvg_out=RPDataTimeAvg_out/(RPTime(nSamples_out)-RPTime(1))
END IF
WRITE(UNIT_stdOut,'(A)')'done.'
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE CalcTimeAvg


!===================================================================================================================================
!> Calculate temporal TimeAvg values using 'raw' data (not interpolated on equi grid)
!===================================================================================================================================
SUBROUTINE CalcFluctuations()
! MODULES
USE MOD_Globals
USE MOD_RPInterpolation_Vars
USE MOD_OutputRPVisu_Vars          ,ONLY: nSamples_out,RPData_out,RPDataTimeAvg_out,nSamples_out
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iSample
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Calculating Fluctuations...'

DO iSample=1,nSamples_out
  RPData_out(:,:,iSample)=RPData_out(:,:,iSample)-RPDataTimeAvg_out(:,:)
END DO
WRITE(UNIT_stdOut,'(A)')'done.'
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE CalcFluctuations


!===================================================================================================================================
!> Deallocate global variable for Recordpoints
!===================================================================================================================================
SUBROUTINE FinalizeInterpolation()
! MODULES
USE MOD_RPInterpolation_Vars
IMPLICIT NONE
!===================================================================================================================================
SDEALLOCATE(dt)
END SUBROUTINE FinalizeInterpolation

END MODULE MOD_RPInterpolation
