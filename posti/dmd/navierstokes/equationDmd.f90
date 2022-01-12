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
USE MOD_DMD_Vars          ,ONLY:nVar_State,VarNames_State
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
INTEGER                               :: iVar,countCons
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT EquationDMD ...'

!check if Varnames in HDF5 file are conservatives
! TODO: also generate mapping in case conservatives are there but not in the correct order
countCons=0
IF(nVar_State.GE.PP_nVar)THEN
  DO iVar=1,PP_nVar
    IF (STRICMP(VarNames_State(iVar), DepNames(iVar))) THEN
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


WRITE(UNIT_stdOut,'(A)')' INIT EquationDMD DONE!'
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitEquationDMD



!===================================================================================================================================
!> This routine computes the state on the visualization grid
!===================================================================================================================================
SUBROUTINE CalcEquationDMD(DMDData,DMDData_out)
! MODULES
USE MOD_Globals
USE MOD_DMD_Vars            ,ONLY: nVar_State,VarNames_State,N_State,N_StateZ,nDoFs,nElems_State,nVarDMD,use2D
USE MOD_EquationDMD_Vars
USE MOD_EOS_Posti           ,ONLY: CalcQuantities
USE MOD_StringTools         ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)    :: DMDData(nVar_State,N_State+1,N_State+1,N_StateZ+1,nElems_State)
REAL,INTENT(OUT)   :: DMDData_out(nDoFs*nVarDMD)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nVal(4)
INTEGER            :: iVarOut,iVarIn,iVar
REAL,ALLOCATABLE   :: UCalc(:,:,:,:,:)
REAL,ALLOCATABLE   :: DMDData_tmp(:,:,:,:,:)
REAL,ALLOCATABLE   :: DMDData_tmp2(:,:)
REAL,ALLOCATABLE   :: DMDData_tmp3(:,:)
!INTEGER            :: OutputVarsIndex(MAXVAL(mapVisu))
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')" CONVERT DERIVED QUANTITIES..."
! CALCULATE DERIVED QUATITIES -----------------------------------------------------------------------------------------------------!
maskCalc=1
nVal=(/N_State+1,N_State+1,N_StateZ+1,nElems_State/)
! Copy existing variables from solution array
! Attention: nVarCalc must be last dimension (needed for CalcQuantities from flexilib!)
IF(.NOT. use2D) THEN
  ALLOCATE(UCalc(N_State+1,N_State+1,N_State+1,nElems_State,nVarCalc))
  ALLOCATE(DMDData_tmp(nVarCalc,N_State+1,N_State+1,N_State+1,nElems_State))
  ALLOCATE(DMDData_tmp2(nElems_State*(N_State+1)**3,nVarCalc))
  ALLOCATE(DMDData_tmp3(nElems_State*(N_State+1)**3,nVarVisuTotal))
ELSE
  ALLOCATE(UCalc(N_State+1,N_State+1,1        ,nElems_State,nVarCalc))
  ALLOCATE(DMDData_tmp(nVarCalc,N_State+1,N_State+1,1,nElems_State))
  ALLOCATE(DMDData_tmp2(nElems_State*(N_State+1)**2,nVarCalc))
  ALLOCATE(DMDData_tmp3(nElems_State*(N_State+1)**2,nVarVisuTotal))
END IF

DO iVarOut=1,nVarDep ! iterate over all out variables
  IF (mapCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
  DO iVarIn=1,nVar_State ! iterate over all in variables
    IF( STRICMP(VarNamesAll(iVarOut),VarNames_State(iVarIn))) THEN
      UCalc(:,:,:,:,mapCalc(iVarOut))=DMDData(iVarIn,:,:,:,:)
      maskCalc(iVarOut)=0 ! remove variable from maskCalc, since they now got copied and must not be calculated.
    END IF
  END DO
eND DO

! calculate all quantities
CALL CalcQuantities(nVarCalc,nVal,(/1/),mapCalc,UCalc,maskCalc)

! fill output array
DMDData_tmp2=RESHAPE(UCalc,(/nDofs,nVarCalc/))

DO iVar=1,nVarDep
  IF (mapVisu(iVar) .GT. 0 ) THEN
    DMDData_tmp3(:,mapVisu(ivar))=DMDData_tmp2(:,mapCalc(iVar))
  END IF
END Do

DMDData_out(:)=RESHAPE(TRANSPOSE(DMDData_tmp3(:,:)), (/nDoFs*nVarDMD/))

DEALLOCATE(UCalc)
DEALLOCATE(DMDData_tmp)

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

!===================================================================================================================================
!> This routine builds the mappings from the total number of variables available for visualization to number of calculation
!> and visualization variables.
!>  1. Read 'VarName' options from the parameter file. This are the quantities that will be visualized.
!>  2. Initialize the dependecy table
!>  3. check wether gradients are needed for any quantity. If this is the case, remove the conservative quantities from the
!>     dependecies of the primitive quantities (the primitive quantities are available directly, since the DGTimeDerivative_weakForm
!>     will be executed.
!>  4. build the 'mapCalc' that holds for each quantity that will be calculated the index in 'UCalc' array (0 if not calculated)
!>  5. build the 'mapVisu' that holds for each quantity that will be visualized the index in 'UVisu' array (0 if not visualized)
!===================================================================================================================================
SUBROUTINE Build_mapCalc_mapVisu()
USE MOD_DMD_Vars        ,ONLY: VarNameDMD,nVarDMD
USE MOD_EquationDMD_Vars
USE MOD_ReadInTools     ,ONLY: GETSTR,CountOption
USE MOD_StringTools     ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLESIABLES
INTEGER             :: iVar,iVar2
CHARACTER(LEN=20)   :: format
!===================================================================================================================================
! Read Varnames from parameter file and fill
!   mapVisu = map, which stores at position x the position/index of the x.th quantity in the UVisu array
!             if a quantity is not visualized it is zero
ALLOCATE(mapVisu(1:nVarDep))
mapVisu = 0
nVarVisuTotal = 0
! Compare varnames that should be visualized with availabe varnames
DO iVar=1,nVarDMD
  DO iVar2=1,nVarDep
    IF (STRICMP(VarNameDMD(iVar), VarNamesAll(iVar2))) THEN
      mapVisu(iVar2) = nVarVisuTotal+1
      nVarVisuTotal = nVarVisuTotal + 1
    END IF
  END DO
END DO

! Calculate all dependencies:
! For each quantity copy from all quantities that this quantity depends on the dependencies.
DO iVar=1,nVarDep
  DepTable(iVar,iVar) = 1
  DO iVar2=1,iVar-1
    IF (DepTable(iVar,iVar2).EQ.1) &
      DepTable(iVar,:) = MAX(DepTable(iVar,:), DepTable(iVar2,:))
  END DO
END DO

! Build :
!   mapCalc = map, which stores at position x the position/index of the x.th quantity in the UCalc array
!             if a quantity is not calculated it is zero
ALLOCATE(mapCalc(1:nVarDep))
mapCalc = 0
DO iVar=1,nVarDep
  IF (mapVisu(iVar).GT.0) THEN
    mapCalc = MAX(mapCalc,DepTable(iVar,1:nVarDep))
  END IF
END DO
! enumerate mapCalc
nVarCalc = 0
DO iVar=1,nVarDep
  IF (mapCalc(iVar).GT.0) THEN
    nVarCalc = nVarCalc + 1
    mapCalc(iVar) = nVarCalc
  END IF
END DO

! print the dependecy table
WRITE(format,'(I2)') SIZE(DepTable,2)
DO iVar=1,nVarDep
  WRITE (*,'('//format//'I2,A)') DepTable(iVar,:), " "//TRIM(VarNamesAll(iVar))
END DO

! print the mappings
WRITE(format,'(I2)') nVarDep
WRITE (*,'(A,'//format//'I3)') "mapCalc ",mapCalc
WRITE (*,'(A,'//format//'I3)') "mapVisu ",mapVisu

END SUBROUTINE Build_mapCalc_mapVisu
