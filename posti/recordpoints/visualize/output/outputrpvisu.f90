#include "flexi.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_OutputRPVisu
! MODULES
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
USE MOD_Parameters      ,ONLY: nVar_visu
USE MOD_Parameters      ,ONLY: Line_LocalCoords,Plane_LocalCoords,equiTimeSpacing
USE MOD_RPSet_Vars      ,ONLY: nRP_global
USE MOD_RPData_Vars     ,ONLY: nSamples_global
USE MOD_OutputRPVisu_Vars     ,ONLY: nSamples_out
USE MOD_OutputRPVisu_Vars     ,ONLY: nCoords,CoordNames
USE MOD_OutputRPVisu_Vars     ,ONLY: RPData_out
USE MOD_Equation_Vars   ,ONLY: EquationInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: offset
!===================================================================================================================================
IF (.NOT.EquationInitIsDone) THEN
  CALL abort(__STAMP__,'InitEquation must be called before InitOutput!')
  ! since initequation defines nVarVisu, VarNameVisu
END IF
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'

IF(.NOT.equiTimeSpacing) nSamples_out  = nSamples_global
ALLOCATE(RPData_out(1:nVar_visu,nRP_global,nSamples_out))
RPData_out(:,:,:) = 0.

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
!> Call of software specific output routines
!===================================================================================================================================
SUBROUTINE OutputRP()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars     ,ONLY: RPTime
USE MOD_Parameters      ,ONLY: OutputFormat,thirdOct,ProjectName
#ifdef WITHTECPLOT
USE MOD_Tecplot   
#endif
USE MOD_OutputRPVisu_Vars     ,ONLY: nSamples_out,RPData_out,RPDataTimeAvg_out,CoordNames 
USE MOD_OutputRPVisu_HDF5
USE MOD_spec_Vars       ,ONLY: nSamples_spec,RPData_freq,RPData_spec
USE MOD_spec_Vars       ,ONLY: nSamples_Oct,RPData_freqOct,RPData_Oct
USE MOD_Parameters      ,ONLY: nVar_visu,VarNameVisu
USE MOD_Parameters      ,ONLY: OutputTimeAverage,OutputTimeData,doSpec,doFluctuations
USE MOD_Parameters      ,ONLY: Plane_doBLProps
USE MOD_RPSet_Vars      ,ONLY: nRP_global
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: FileName
CHARACTER(LEN=255)            :: strOutputFile
!===================================================================================================================================
! Output Time Signal
IF(OutputTimeData) THEN
WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP'
  IF(doFluctuations) FileName=TRIM(FileName)//'_Fluc'
  SELECT CASE(OutputFormat)
#ifdef WITHTECPLOT
    CASE(0) ! Tecplot Binary Output, in one file, mesh+solution
      strOutputFile=TRIM(FileName)//'.plt'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING TIME SIGNAL TO ',strOutputFile
      IF(nSamples_out.GT.1) THEN
        CALL WriteDataToTecplotBinary(nSamples_out,nRP_global,nVar_visu,VarNameVisu,RPTime,RPData_out,strOutputFile)
      ELSE 
        ! use time avg routine if only one sample is present.
        CALL WriteTimeAvgDataToTecplotBinary(nRP_global,nVar_visu,VarNameVisu,RPData_out(:,:,1),strOutputFile)
      END IF
#endif
!    CASE(1) ! Tecplot ASCII Output
!      strOutputFile=TRIM(strOutputFile)//'.dat'
!      CALL WriteDataToTecplot(nSamples_out,nVar_visu,VarNameVisu,RPData_out,strOutputFile)
     CASE(2) ! structured HDF5 output
       strOutputFile=TRIM(FileName)//'_PP.h5'
       CALL WriteDataToHDF5(nSamples_out,nRP_global,nVar_visu,VarNameVisu,RPTime,RPData_out,strOutputFile)
  END SELECT
WRITE(UNIT_StdOut,'(132("-"))')
END IF !output time data

! Output spectra
IF(doSpec) THEN
CoordNames(1)='Frequency'
WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_spec'
  SELECT CASE(OutputFormat)
#ifdef WITHTECPLOT
    CASE(0) ! Tecplot Binary Output, in one file, mesh+solution
      strOutputFile=TRIM(FileName)//'.plt'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING SPECTRA TO ',strOutputFile
      CALL WriteDataToTecplotBinary(nSamples_spec,nRP_global,nVar_visu,VarNameVisu,RPData_freq,RPData_spec,strOutputFile)
#endif
    CASE(2) ! structured HDF5 output
     strOutputFile=TRIM(FileName)//'_PP.h5'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING SPECTRA TO ',strOutputFile
      CALL WriteDataToHDF5(nSamples_spec,nRP_global,nVar_visu,VarNameVisu,RPData_freq,RPData_spec,strOutputFile)
  END SELECT
WRITE(UNIT_StdOut,'(132("-"))')
END IF !output time data

IF(ThirdOct) THEN
CoordNames(1)='Frequency'
WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_Octspec'
  SELECT CASE(OutputFormat)
#ifdef WITHTECPLOT
    CASE(0) ! Tecplot Binary Output, in one file, mesh+solution
      strOutputFile=TRIM(FileName)//'.plt'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING THIRD OCTAVE SPECTRA TO ',strOutputFile
      CALL WriteDataToTecplotBinary(nSamples_Oct,nRP_global,nVar_visu,VarNameVisu,RPData_freqOct,RPData_Oct,strOutputFile)
#endif
    CASE(2) ! structured HDF5 output
     strOutputFile=TRIM(FileName)//'_PP.h5'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING THIRD OCTAVE SPECTRA TO ',strOutputFile
      CALL WriteDataToHDF5(nSamples_Oct,nRP_global,nVar_visu,VarNameVisu,RPData_freqOct,RPData_Oct,strOutputFile)
  END SELECT
WRITE(UNIT_StdOut,'(132("-"))')
END IF !output 1/3 Oct Spectra

! Output Time Average
IF(OutputTimeAverage) THEN
  WRITE(UNIT_StdOut,'(132("-"))')
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_TimeAvg'
  SELECT CASE(OutputFormat)
#ifdef WITHTECPLOT
    CASE(0) ! Tecplot Binary Output, in one file, mesh+solution
      strOutputFile=TRIM(FileName)//'.plt'
      WRITE(UNIT_stdOut,'(A,A)')' WRITING TIME AVERAGE TO ',strOutputFile
      CALL WriteTimeAvgDataToTecplotBinary(nRP_global,nVar_visu,VarNameVisu,RPDataTimeAvg_out,strOutputFile)
#endif
   CASE(2) ! structured HDF5 output
     strOutputFile=TRIM(FileName)//'_PP.h5'
     CALL WriteDataToHDF5(1,nRP_global,nVar_visu,VarNameVisu,RPTime,RPData_out,strOutputFile)
  WRITE(UNIT_StdOut,'(132("-"))')
  END SELECT
END IF

#ifdef WITHBLPROPS
IF(Plane_doBLProps)THEN !output the BL stuff along lines
  Filename=TRIM(ProjectName)
  FileName=TRIM(FileName)//'_RP_BLProps'
  SELECT CASE(OutputFormat)
#ifdef WITHTECPLOT
    CASE(0) ! Tecplot Binary Output, in one file, mesh+solution
    strOutputFile=TRIM(FileName)//'.plt'
    WRITE(UNIT_stdOut,'(A,A)')' WRITING BL PROPS TO ',strOutputFile
    CALL WriteBLPropsToTecplotBinary(strOutputFile)
#endif
   CASE(2) ! structured HDF5 output
     strOutputFile=TRIM(FileName)//'_PP.h5'
     CALL WriteBLPropsToHDF5(strOutputFile)
   END SELECT
END IF
#endif

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
