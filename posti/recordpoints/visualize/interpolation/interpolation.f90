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
!===================================================================================================================================
!> Module to handle temporal interpolation for the recordpoints.
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
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out
USE MOD_ParametersVisu       ,ONLY: equiTimeSpacing
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample
INTEGER                         :: nSamples_min
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
    WRITE(UNIT_stdOut,'(A)')'WARNING: Equidistant sampling time step is lower than or equal to the input time step!!!'
  END IF

  WRITE(UNIT_stdOut,'(A)')'  Database:'
  WRITE(UNIT_stdOut,'(A18,ES16.7)')'         TStart = ',RPTIME(1)
  WRITE(UNIT_stdOut,'(A18,ES16.7)')'           TEnd = ',TEnd
  WRITE(UNIT_stdOut,'(A18,ES16.7)')'          dtMin = ',dtMin
  WRITE(UNIT_stdOut,'(A18,I8    )')'       nSamples = ',nSamples_global
  WRITE(UNIT_stdOut,'(A)')'  Equidistant time grid:'
  WRITE(UNIT_stdOut,'(A18,ES16.7)')'             dt = ',dt_out
  WRITE(UNIT_stdOut,'(A18,I8    )')'       nSamples = ',nSamples_out
END IF

WRITE(UNIT_stdOut,'(A)')' INIT INTERPOLATION DONE'
WRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitInterpolation


!===================================================================================================================================
!> Interpolate the RP data on equidistant time grid
!===================================================================================================================================
SUBROUTINE InterpolateEquiTime()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_RPData_Vars          ,ONLY: RPTime,nVar_HDF5,nSamples_global,RPData
USE MOD_RPSetVisuVisu_Vars   ,ONLY: nRP_global
USE MOD_RPInterpolation_Vars
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out
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
    IF((Time_right.GE.Time).OR.ALMOSTEQUALABSORREL(Time_right,Time,PP_RealTolerance)) EXIT
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
!> Calculate temporal TimeAvg values
!> An additonal file may also be specified that contains the temporal average, e.g. from a 2D averaged solution
!===================================================================================================================================
SUBROUTINE CalcTimeAvg()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools          ,ONLY: GETSTR
USE MOD_RPSetVisuVisu_Vars   ,ONLY: nRP_global
USE MOD_RPData_Vars          ,ONLY: RPTime,nVar_HDF5
USE MOD_RPSetVisuVisu_Vars ,  ONLY: nRP_HDF5,RPOutMap
USE MOD_RPInterpolation_Vars
USE MOD_ParametersVisu       ,ONLY: nVarVisu,EquiTimeSpacing
USE MOD_ParametersVisu       ,ONLY: Line_LocalVel,Plane_LocalVel,Box_LocalVel
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out,RPData_out,RPDataTimeAvg_out
USE MOD_EquationRP           ,ONLY: Line_TransformVel,Plane_TransformVel,Box_TransformVel
USE MOD_HDF5_Input
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iSample
CHARACTER(LEN=255)      :: TimeAvgFile
REAL,ALLOCATABLE        :: temparray(:,:,:),temparray2(:,:,:)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' Calculating Time Average...'

ALLOCATE(RPDataTimeAvg_out(1:nVarVisu,nRP_global))

! Check if an external file should be used for the time averages
TimeAvgFile = GETSTR('TimeAvgFile','NOT_SET')
IF (TRIM(TimeAvgFile).EQ.'NOT_SET') THEN
  ! Average the available samples
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
ELSE
  ! Read in the average from the specified file
  ALLOCATE(temparray(0:nVar_HDF5,1:nRP_HDF5,1)) ! storing complete sample set
  CALL OpenDataFile(TimeAvgFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  ! Safety check for the size of the array
  CALL GetDataSize(File_ID,'RP_Data',nDims,HSize)
  IF(nRP_HDF5 .NE. HSize(2)) THEN
    CALL Abort(__STAMP__,'ERROR - Number of RPs in TimeAvg file does not match RP definition file!')
  END IF
  IF(nVar_HDF5 .NE. INT(HSize(1) -1)) THEN
    CALL Abort(__STAMP__,'ERROR - Wrong number of variables in TimeAvg file!')
  END IF
  DEALLOCATE(HSize)
  CALL ReadArray('RP_Data',3,(/nVar_HDF5+1,nRP_HDF5,1/),0,3,RealArray=temparray)
  CALL CloseDataFile()
  ! Coordinate Transform
  ALLOCATE(temparray2(nVarVisu,nRP_global,1))
  temparray2(:,:,1)=temparray(1:nVarVisu,RPOutMap(:),1) !RPOutMap filters out RPs which are to be visualized
  DEALLOCATE(temparray)
  IF(Line_LocalVel) &
    CALL Line_TransformVel(temparray2,1)
  IF(Plane_LocalVel) &
    CALL Plane_TransformVel(temparray2,1)
  IF(Box_LocalVel) &
    CALL Box_TransformVel(temparray2,1)
  ! Save in the TimeAvg array
  RPDataTimeAvg_out=temparray2(:,:,1)
  DEALLOCATE(temparray2)
END IF

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' Calculating Time Average done'
WRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE CalcTimeAvg


!===================================================================================================================================
!> Calculate temporal fluctuations by substracting the mean from the time-resolved value
!===================================================================================================================================
SUBROUTINE CalcFluctuations()
! MODULES
USE MOD_Globals
USE MOD_OutputRPVisu_Vars          ,ONLY: nSamples_out,RPData_out,RPDataTimeAvg_out,nSamples_out
USE MOD_OutputRPVisu_Vars          ,ONLY: RPDataRMS_out
USE MOD_ParametersVisu             ,ONLY: doFluctuations
USE MOD_RPInterpolation_Vars
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
  IF(doFluctuations) RPDataRMS_out(:,:)=RPData_out(:,:,iSample)**2.+RPDataRMS_out(:,:)
END DO

IF (doFluctuations) RPDataRMS_out(:,:)=SQRT(RPDataRMS_out(:,:)/nSamples_out)

WRITE(UNIT_stdOut,'(A)')'done.'
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE CalcFluctuations


!===================================================================================================================================
!> Deallocate global variable for temporal interpolatio
!===================================================================================================================================
SUBROUTINE FinalizeInterpolation()
! MODULES
USE MOD_RPInterpolation_Vars
IMPLICIT NONE
!===================================================================================================================================
SDEALLOCATE(dt)
END SUBROUTINE FinalizeInterpolation

END MODULE MOD_RPInterpolation
