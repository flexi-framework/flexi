#include "flexi.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_EquationDMD
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitEquationDMD
  MODULE PROCEDURE InitEquationDMD
END INTERFACE

INTERFACE CalcEquationDMD
  MODULE PROCEDURE CalcEquationDMD
END INTERFACE

INTERFACE FinalizeEquationDMD
  MODULE PROCEDURE FinalizeEquationDMD
END INTERFACE

PUBLIC::InitEquationDMD,CalcEquationDMD,FinalizeEquationDMD
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize the visualization and map the variable names to classify these in conservative and derived quantities.
!===================================================================================================================================
SUBROUTINE InitEquationDMD()
! MODULES
USE MOD_Globals
USE MOD_DMDData_Vars       ,ONLY:VarNames_HDF5,nVar_HDF5
USE MOD_ParametersVisu        
USE MOD_EquationDMD_Vars
USE MOD_EOS               ,ONLY: InitEOS
USE MOD_EOS_Posti_Vars    ,ONLY: nVarDepEOS,DepTableEOS,DepNames
USE MOD_Readintools       ,ONLY: CountOption,GETSTR
USE MOD_StringTools       ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: iVar,iVar2,strlen,countCons
INTEGER                               :: mapCand(3)
CHARACTER(LEN=255)                    :: tmp255,tmp255_2
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT EquationDMD ...'

justVisualizeState=.FALSE.

! In case no output variables specified: take instead the variable names out of the hDF5 file (used for timeavg-files)
nVarVisu=CountOption("VarName")
IF(nVarVisu .LT. 1) THEN
  justVisualizeState=.TRUE.
  WRITE(*,*) 'No output variables specified, using existing variables from state file.'
  nVarVisu   = nVar_HDF5
  nVarDep    = nVar_HDF5
  ALLOCATE(VarNameVisu(nVar_HDF5))
  ALLOCATE(VarNamesAll(nVar_HDF5))
  VarNameVisu = VarNames_HDF5
  VarNamesAll = VarNames_HDF5
ELSE
  ALLOCATE(VarNameVisu(nVarVisu))
  DO iVar=1,nVarVisu
    VarNameVisu(iVar) = GETSTR("VarName")
  END DO
  !check if Varnames in HDF5 file are conservatives
  ! TODO: also generate mapping in case conservatives are there but not in the correct order
  countCons=0
  IF(nVar_HDF5.GE.PP_nVar)THEN
    DO iVar=1,PP_nVar
      IF (STRICMP(VarNames_HDF5(iVar), DepNames(iVar))) THEN
        countCons=countCons+1
      END IF
    END DO
  END IF
  IF(countCons.NE.PP_nVar) THEN
    CALL CollectiveStop(__STAMP__,'Not all necessary variables are present in HDF5 files')
  END IF
  nVarDep=nVarDepEOS
  ALLOCATE(VarNamesAll(nVarDep))
  VarNamesAll=DepNames
  ALLOCATE(DepTable(nVarDep,0:nVarDep))
  DepTable=DepTableEOS

  ! generate mappings
  CALL Build_mapCalc_mapVisu()

  ! initialize EOS
  CALL InitEOS()
END IF


WRITE(UNIT_stdOut,'(A)')' INIT EquationDMD DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
EquationDMDInitIsDone=.TRUE.
END SUBROUTINE InitEquationDMD



!===================================================================================================================================
!> This routine computes the state on the visualization grid 
!===================================================================================================================================
SUBROUTINE CalcEquationDMD()
! MODULES
USE MOD_Globals
USE MOD_DMD_Vars            ,ONLY: DMDData       
USE MOD_DMDData_Vars        ,ONLY: nVar_HDF5,VarNames_HDF5
USE MOD_DMDSetVisuVisu_Vars ,ONLY: nDMD_global        
USE MOD_OutputDMDVisu_Vars  ,ONLY: nSamples_out
USE MOD_ParametersVisu     ,ONLY: nVarDep,nVarCalc,mapCalc,mapVisu,VarNamesAll,justVisualizeState
USE MOD_ParametersVisu     ,ONLY: Line_LocalVel,Plane_LocalVel
USE MOD_OutputDMDVisu_Vars  ,ONLY: DMDData_out 
USE MOD_EOS_Posti          ,ONLY: CalcQuantities
USE MOD_StringTools        ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nVal(2)
INTEGER            :: iVarOut,iVarIn,iVar,iVarCalc,iVarVisu
REAL,ALLOCATABLE   :: UCalc(:,:,:)
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')" CONVERT DERIVED QUANTITIES..."
! CALCULATE DERIVED QUATITIES -----------------------------------------------------------------------------------------------------!
IF(justVisualizeState)THEN
  DMDData_out=DMDData
ELSE
  maskCalc=1
  nVal=(/nDMD_global,nSamples_out/)
  ! Copy existing variables from solution array
  ! Attention: nVarCalc must be last dimension (needed for CalcQuantities from flexilib!)
  ALLOCATE(UCalc(nDMD_global,nSamples_out,nVarCalc))
  
  DO iVarOut=1,nVarDep ! iterate over all out variables
    IF (mapCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
    DO iVarIn=1,nVar_HDF5 ! iterate over all in variables
      IF( STRICMP(VarNamesAll(iVarOut),VarNames_HDF5(iVarIn))) THEN
        UCalc(:,:,mapCalc(iVarOut))=DMDData(iVarIn,:,:)
        maskCalc(iVarOut)=0 ! remove variable from maskCalc, since they now got copied and must not be calculated.
      END IF
    END DO
  END DO
  
  ! calculate all quantities
  CALL CalcQuantities(nVarCalc,nVal,(/1/),mapCalc,UCalc,maskCalc)
  
  ! fill output array
  DO iVar=1,nVarDep
    IF (mapVisu(iVar).GT.0) THEN
      iVarCalc = mapCalc(iVar) 
      iVarVisu = mapVisu(iVar) 
      DMDData_out(iVarVisu,:,:)=UCalc(:,:,iVarCalc)
    END IF
  END DO 
  DEALLOCATE(UCalc)
END IF

WRITE(UNIT_stdOut,'(A)')" CONVERT DERIVED QUANTITIES DONE!"
END SUBROUTINE CalcEquationDMD


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE FinalizeEquationDMD()
! MODULES
USE MOD_Globals
USE MOD_EquationDMD_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DEALLOCATE(TransMap,is2D) 
WRITE(UNIT_stdOut,'(A)') '  EquationDMD FINALIZED'
END SUBROUTINE FinalizeEquationDMD

END MODULE MOD_EquationDMD

