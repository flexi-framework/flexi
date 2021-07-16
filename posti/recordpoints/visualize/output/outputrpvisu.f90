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
!> Module contains the routines that manage the output of te record point data.
!===================================================================================================================================
MODULE MOD_OutputRPVisu
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE

INTERFACE OutputRP
  MODULE PROCEDURE OutputRP
END INTERFACE

INTERFACE FinalizeOutput
  MODULE PROCEDURE FinalizeOutput
END INTERFACE

PUBLIC:: InitOutput,OutputRP,FinalizeOutput
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize output of the record point evaluation
!===================================================================================================================================
SUBROUTINE InitOutput()
! MODULES
USE MOD_Globals
USE MOD_ParametersVisu     ,ONLY: nVarVisu
USE MOD_ParametersVisu     ,ONLY: Line_LocalCoords,Plane_LocalCoords,equiTimeSpacing
USE MOD_ParametersVisu     ,ONLY: doFluctuations
USE MOD_RPSetVisuVisu_Vars ,ONLY: nRP_global
USE MOD_RPData_Vars        ,ONLY: nSamples_global
USE MOD_OutputRPVisu_Vars  ,ONLY: nSamples_out
USE MOD_OutputRPVisu_Vars  ,ONLY: nCoords,CoordNames
USE MOD_OutputRPVisu_Vars  ,ONLY: RPData_out,RPDataRMS_out
USE MOD_EquationRP_Vars    ,ONLY: EquationRPInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: offset
!===================================================================================================================================
IF (.NOT.EquationRPInitIsDone) THEN
  CALL abort(__STAMP__,'InitEquationRP must be called before InitOutput!')
  ! since initequation defines nVarVisu, VarNameVisu
END IF
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'

IF(.NOT.equiTimeSpacing) nSamples_out  = nSamples_global
ALLOCATE(RPData_out(1:nVarVisu,nRP_global,nSamples_out))
RPData_out(:,:,:) = 0.

IF(doFluctuations) THEN
  ALLOCATE(RPDataRMS_out(1:nVarVisu,nRP_global))
  RPDataRMS_out(:,:) = 0.
END IF

! prepare Coordinate Varnames
nCoords=4
IF(Line_LocalCoords) nCoords=nCoords+1
IF(Plane_LocalCoords) nCoords=nCoords+2
ALLOCATE(CoordNames(nCoords))
CoordNames(1)='Time'
offset=1
IF(Line_LocalCoords) THEN
  CoordNames(2)='CoordinateLine'
  offset=offset+1
END IF
IF(Plane_LocalCoords) THEN
  CoordNames(offset+1)='PlaneX'
  CoordNames(offset+2)='PlaneY'
  offset=offset+2
END IF
CoordNames(offset+1)='CoordinateX'
CoordNames(offset+2)='CoordinateY'
CoordNames(offset+3)='CoordinateZ'

WRITE(UNIT_stdOut,'(A)')' INIT OUTPUT DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOutput


!===================================================================================================================================
!> Call of software specific output routines for the user-specified output type. We can perform output of the time-accurate signal,
!> of the temporal averages or of spectral data. Will prepare file names etc. and then call the specific routines to perform
!> the output in the choosen file format.
!===================================================================================================================================
SUBROUTINE OutputRP()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars        ,ONLY: RPTime
USE MOD_ParametersVisu     ,ONLY: OutputFormat,thirdOct,ProjectName
USE MOD_OutputRPVisu_Vars  ,ONLY: RPDataTimeAvg_out,RPDataRMS_out
USE MOD_OutputRPVisu_VTK   ,ONLY: WriteDataToVTK,WriteTimeAvgDataToVTK,WriteBLPropsToVTK
USE MOD_OutputRPVisu_Vars  ,ONLY: nSamples_out,RPData_out,CoordNames
USE MOD_OutputRPVisu_HDF5  ,ONLY: WriteDataToHDF5,WriteBLPropsToHDF5
USE MOD_Spec_Vars          ,ONLY: nSamples_spec,RPData_freq,RPData_spec
USE MOD_Spec_Vars          ,ONLY: nSamples_Oct,RPData_freqOct,RPData_Oct
USE MOD_ParametersVisu     ,ONLY: nVarVisu,VarNameVisu
USE MOD_ParametersVisu     ,ONLY: OutputTimeAverage,OutputTimeData,doSpec,doFluctuations
USE MOD_ParametersVisu     ,ONLY: Plane_doBLProps
USE MOD_ParametersVisu     ,ONLY: doEnsemble,doTurb
USE MOD_EnsembleRP_Vars    ,ONLY: enSamples,nVar_ensTurb,RPData_ens,RPData_freqEns,RPData_turb
USE MOD_RPData_Vars        ,ONLY: VarNames_HDF5,nVar_HDF5
USE MOD_RPSetVisuVisu_Vars ,ONLY: nRP_global
USE MOD_Turbulence_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: FileName
CHARACTER(LEN=255)            :: strOutputFile
CHARACTER(LEN=255),ALLOCATABLE:: VarNameTurb(:)
REAL,ALLOCATABLE              :: RPData_OutTurb(:,:,:),RPData_turbAvg(:,:)
INTEGER                       :: nVar_turb
REAL,ALLOCATABLE              :: TimeAvg_tmp(:,:,:)
!===================================================================================================================================
! Output Time Signal
IF(OutputTimeData) THEN
WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP'
  IF(doFluctuations) FileName=TRIM(FileName)//'_Fluc'
  SELECT CASE(OutputFormat)
    CASE(0) ! Paraview VTK output
      WRITE(UNIT_stdOut,'(A)')' WRITING TIME SIGNAL TO VTK FILE '
      CALL WriteDataToVTK(nSamples_out,nRP_global,nVarVisu,VarNameVisu,RPTime,RPData_out,FileName)
    CASE(2) ! structured HDF5 output
      strOutputFile=TRIM(FileName)//'_PP.h5'
      CALL WriteDataToHDF5(nSamples_out,nRP_global,nVarVisu,VarNameVisu,RPTime,RPData_out,strOutputFile)
  END SELECT

  IF(doFluctuations) THEN
    Filename=TRIM(ProjectName)
    FileName=TRIM(FileName)//'_RP_RMS'
    SELECT CASE(OutputFormat)
      CASE(0) ! ParaView VTK output
        WRITE(UNIT_stdOut,'(A,A)')' WRITING TIME AVERAGE TO VTK FILE'
        CALL WriteTimeAvgDataToVTK(nRP_global,nVarVisu,VarNameVisu,RPDataRMS_out,FileName)
     CASE(2) ! structured HDF5 output
       strOutputFile=TRIM(FileName)//'_PP.h5'
       ALLOCATE(TimeAvg_tmp(1,nVarVisu,nRP_global))
       ! To use the WriteHDF5 routine, we need to create a temporary array since it expect 3 dimensions
       TimeAvg_tmp(1,:,:) = RPDataRMS_out
       CALL WriteDataToHDF5(1,nRP_global,nVarVisu,VarNameVisu,RPTime,TimeAvg_tmp,strOutputFile)
       DEALLOCATE(TimeAvg_tmp)
    END SELECT
  END IF
WRITE(UNIT_StdOut,'(132("-"))')
END IF !output time data

! Output spectra
IF(doSpec) THEN
CoordNames(1)='Frequency'
WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_spec'
  SELECT CASE(OutputFormat)
    CASE(0) ! Paraview VTK output
      WRITE(UNIT_stdOut,'(A,A)')' WRITING SPECTRA TO VTK FILE '
      CALL WriteDataToVTK(nSamples_spec,nRP_global,nVarVisu,VarNameVisu,RPData_freq,RPData_spec,FileName)
    CASE(2) ! structured HDF5 output
     strOutputFile=TRIM(FileName)//'_PP.h5'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING SPECTRA TO ',strOutputFile
      CALL WriteDataToHDF5(nSamples_spec,nRP_global,nVarVisu,VarNameVisu,RPData_freq,RPData_spec,strOutputFile)
  END SELECT
WRITE(UNIT_StdOut,'(132("-"))')
END IF !output time data

IF(ThirdOct) THEN
CoordNames(1)='Frequency'
WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_Octspec'
  SELECT CASE(OutputFormat)
    CASE(0) ! Paraview VTK output
      strOutputFile=TRIM(FileName)//'.plt'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING THIRD OCTAVE SPECTRA TO VTK '
      CALL WriteDataToVTK(nSamples_Oct,nRP_global,nVarVisu,VarNameVisu,RPData_freqOct,RPData_Oct,FileName)
    CASE(2) ! structured HDF5 output
     strOutputFile=TRIM(FileName)//'_PP.h5'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING THIRD OCTAVE SPECTRA TO ',strOutputFile
      CALL WriteDataToHDF5(nSamples_Oct,nRP_global,nVarVisu,VarNameVisu,RPData_freqOct,RPData_Oct,strOutputFile)
  END SELECT
WRITE(UNIT_StdOut,'(132("-"))')
END IF !output 1/3 Oct Spectra

! Output Time Average
IF(OutputTimeAverage) THEN
  WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_TimeAvg'
  SELECT CASE(OutputFormat)
    CASE(0) ! ParaView VTK output
      WRITE(UNIT_stdOut,'(A,A)')' WRITING TIME AVERAGE TO VTK FILE'
      CALL WriteTimeAvgDataToVTK(nRP_global,nVarVisu,VarNameVisu,RPDataTimeAvg_out,FileName)
   CASE(2) ! structured HDF5 output
     strOutputFile=TRIM(FileName)//'_PP.h5'
     ALLOCATE(TimeAvg_tmp(1,nVarVisu,nRP_global))
     ! To use the WriteHDF5 routine, we need to create a temporary array since it expect 3 dimensions
     TimeAvg_tmp(1,:,:) = RPDataTimeAvg_out
     CALL WriteDataToHDF5(1,nRP_global,nVarVisu,VarNameVisu,RPTime,TimeAvg_tmp,strOutputFile)
     DEALLOCATE(TimeAvg_tmp)
  END SELECT
  WRITE(UNIT_StdOut,'(132("-"))')
END IF

IF(Plane_doBLProps)THEN !output the BL stuff along lines
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_BLProps'
  SELECT CASE(OutputFormat)
    CASE(0) ! ParaView VTK output
    WRITE(UNIT_stdOut,'(A,A)')' WRITING BL PROPS TO VTK'
    CALL WriteBLPropsToVTK(FileName)
   CASE(2) ! structured HDF5 output
     strOutputFile=TRIM(FileName)//'_PP.h5'
     CALL WriteBLPropsToHDF5(strOutputFile)
   END SELECT
END IF

! Output Turbulence stuff
IF(doTurb) THEN
  nVar_turb=5
  ALLOCATE(VarNameTurb(nVar_turb))
  VarNameTurb(1) = 'KineticEnergy'
  VarNameTurb(2) = 'DissipationRate'
  VarNameTurb(3) = 'MeanDissipation_cumulated'
  VarNameTurb(4) = 'Eta_K'
  VarNameTurb(5) = 'Wavenumber K'

  ALLOCATE(RPData_OutTurb(1:nVar_turb,nRP_global,nSamples_spec))
  RPData_OutTurb(1,:,:) = E_kineticSpec(:,:)
  RPData_OutTurb(2,:,:) = disRate(:,:)
  RPData_OutTurb(3,:,:) = epsilonMean(:,:)
  RPData_OutTurb(4,:,:) = etaK(:,:)
  RPData_OutTurb(5,:,:) = kk(:,:)

  CoordNames(1)='Frequency'
  WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_turb'
  SELECT CASE(OutputFormat)
    CASE(0) ! ParaView VTK output
      WRITE(UNIT_stdOut,'(A,A)')' WRITING TURBULENCE DATA TO ',TRIM(FileName)
      CALL WriteDataToVTK (nSamples_spec,nRP_global,nVar_turb,VarNameTurb,RPData_freqTurb,RPData_OutTurb,FileName)
    CASE(2) ! structured HDF5 output
      strOutputFile=TRIM(FileName)//'_PP.h5'
      CALL WriteDataToHDF5(nSamples_spec,nRP_global,nVar_turb,VarNameTurb,RPData_freqTurb,RPData_OutTurb,strOutputFile)
  END SELECT
  WRITE(UNIT_StdOut,'(132("-"))')

  DEALLOCATE(VarNameTurb)
  nVar_turb=2
  ALLOCATE(VarNameTurb(nVar_turb))
  VarNameTurb(1) = 'Eta'
  VarNameTurb(2) = 'MeanDissipation'

  ALLOCATE(RPData_turbAvg(1:nVar_turb,nRP_global))
  RPData_turbAvg(1,:) = eta
  RPData_turbAvg(2,:) = epsilonMean(:,nSamples_spec)

  ! Output Time Average
  WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_turbAvg'
  SELECT CASE(OutputFormat)
    CASE(0) ! ParaView VTK output
      WRITE(UNIT_stdOut,'(A,A)')' WRITING TURBULENCE AVERAGE DATA TO ',TRIM(FileName)
      CALL WriteTimeAvgDataToVTK(nRP_global,nVar_turb,VarNameTurb,       RPData_turbAvg,FileName)
    CASE(2) ! structured HDF5 output
      strOutputFile=TRIM(FileName)//'_PP.h5'
     ALLOCATE(TimeAvg_tmp(1,nVar_turb,nRP_global))
     ! To use the WriteHDF5 routine, we need to create a temporary array since it expect 3 dimensions
     TimeAvg_tmp(1,:,:) = RPData_turbAvg
      CALL WriteDataToHDF5(    1,nRP_global,nVar_turb,VarNameTurb,RPTime,TimeAvg_tmp,strOutputFile)
     DEALLOCATE(TimeAvg_tmp)
  END SELECT
  WRITE(UNIT_StdOut,'(132("-"))')
END IF

IF(doEnsemble)THEN
  ! Output of ensemble averages
  WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_ensAvg'
  SELECT CASE(OutputFormat)
    CASE(0) ! ParaView VTK output
      WRITE(UNIT_stdOut,'(A,A)')' WRITING ESEMBLE AVERAGES TO ',TRIM(FileName)
      CALL WriteDataToVTK(enSamples,nRP_global,nVarVisu,VarNameVisu,RPData_freqEns,RPData_ens,FileName)
    CASE(2) ! structured HDF5 output
      strOutputFile=TRIM(FileName)//'_PP.h5'
      CALL WriteDataToHDF5(enSamples,nRP_global,nVarVisu,VarNameVisu,RPData_freqEns,RPData_ens,strOutputFile)
  END SELECT
  WRITE(UNIT_StdOut,'(132("-"))')

  ! Output of turbulent quantities on the ensemble average
  ALLOCATE(VarNameTurb(nVar_ensTurb))
  VarNameTurb(1) = 'KineticEnergy'
  VarNameTurb(2) = 'TurbulentIntensity'
  VarNameTurb(3) = 'VelocityMagnitudeRMS'
  VarNameTurb(4) = 'VelocityXRMS'
  VarNameTurb(5) = 'VelocityYRMS'
  VarNameTurb(6) = 'VelocityZRMS'

  WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_ensTurb'
  SELECT CASE(OutputFormat)
    CASE(0) ! ParaView VTK output
      WRITE(UNIT_stdOut,'(A,A)')' WRITING TURBULENT QUANTITIES ON ESEMBLE AVERAGES TO ',TRIM(FileName)
      CALL WriteDataToVTK( enSamples,nRP_global,nVar_ensTurb,VarNameTurb,RPData_freqEns,RPData_turb,FileName)
    CASE(2) ! structured HDF5 output
      strOutputFile=TRIM(FileName)//'_PP.h5'
      CALL WriteDataToHDF5(enSamples,nRP_global,nVar_ensTurb,VarNameTurb,RPData_freqEns,RPData_turb,strOutputFile)
  END SELECT

  DEALLOCATE(VarNameTurb)
  WRITE(UNIT_StdOut,'(132("-"))')
END IF

END SUBROUTINE OutputRP


!===================================================================================================================================
!> Deallocate global variables
!===================================================================================================================================
SUBROUTINE FinalizeOutput()
! MODULES
USE MOD_Globals
USE MOD_OutputRPVisu_Vars
IMPLICIT NONE
!===================================================================================================================================
OutputInitIsDone = .FALSE.
END SUBROUTINE FinalizeOutput

END MODULE MOD_OutputRPVisu
