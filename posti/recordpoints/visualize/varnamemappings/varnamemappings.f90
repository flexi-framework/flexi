#include "flexi.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_VarNameMappingsRP
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE CreateVarMappings
  MODULE PROCEDURE CreateVarMappings
END INTERFACE

INTERFACE CreateStateMappings
  MODULE PROCEDURE CreateStateMappings
END INTERFACE

PUBLIC::CreateVarMappings,CreateStateMappings
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Create Mappings between VarName_visu List and the derived quantity to be visualized 
!===================================================================================================================================
SUBROUTINE CreateVarMappings(nVar,VarName,DQ)
! MODULES
USE MOD_Globals
USE MOD_VarNameMappingsRP_Vars    ,ONLY: tDerivedQ,max_nVarVisu
USE MOD_Parameters                ,ONLY: nVarVisu,VarNameVisu
USE MOD_RPData_Vars               ,ONLY: nVar_HDF5,VarNames_HDF5
USE MOD_StringTools               ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: nVar
CHARACTER(LEN=255),INTENT(IN)    :: VarName(1:nVar)
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(tDerivedQ),INTENT(OUT)      :: DQ
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: ReadFromFile(1:nVar)
INTEGER                           :: iVar,iDQ
INTEGER                           :: Map(1:nVar)
!===================================================================================================================================
IF(.NOT. ALLOCATED(DQ%VarName)  )ALLOCATE(DQ%VarName(1:nVar))
IF(.NOT. ALLOCATED(DQ%IndGlobal))ALLOCATE(DQ%IndGlobal(1:nVar))
DQ%nVar      =nVar
DQ%Varname(1:nVar)=VarName(1:nVar)
ReadFromFile(:)    =.FALSE.
Map  =0
DQ%nVarVisu=0
DO iDQ=1,nVar
  DO iVar=1,nVar_HDF5
    IF(STRICMP(VarNames_HDF5(iVar),VarName(iDQ)))THEN
      ReadFromFile(iDQ)  =.TRUE.
    END IF
  END DO
END DO
DO iDQ=1,nVar
  DO iVar=1,nVarVisu
    IF(STRICMP(VarNameVisu(iVar),VarName(iDQ)) .AND. .NOT.ReadFromFile(iDQ))THEN
      WRITE(UNIT_StdOut,*) '  ',TRIM(VarNameVisu(iVar))
      DQ%nVarVisu           =DQ%nVarVisu+1
      Map(iDQ)               =DQ%nVarVisu
      DQ%IndGlobal(iDQ)      =iVar
    END IF
  END DO
END DO
IF(DQ%nVarVisu .GT. 0)THEN
  max_nVarVisu=MAX(max_nVarVisu,DQ%nVarVisu)
  ! Only allocate once per run!
  IF(.NOT. ALLOCATED(DQ%Ind))ALLOCATE(DQ%Ind(DQ%nVarVisu))
  DO iDQ=1,DQ%nVar
    IF(Map(iDQ) .GT. 0) DQ%Ind(Map(iDQ))=iDQ
  END DO
END IF

END SUBROUTINE CreateVarMappings


!===================================================================================================================================
!> Create Mappings for non-derived quanitities ( variables already available in the state file)
!===================================================================================================================================
SUBROUTINE CreateStateMappings(nVar,VarName,DQ)
! MODULES
USE MOD_Globals
USE MOD_VarNameMappingsRP_Vars,ONLY: tDerivedQ,max_nVarVisu
USE MOD_Parameters            ,ONLY: nVarVisu,VarNameVisu
USE MOD_StringTools           ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: nVar
CHARACTER(LEN=255),INTENT(IN)    :: VarName(1:nVar)
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(tDerivedQ),INTENT(OUT)      :: DQ
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iVar,iDQ
INTEGER                           :: Map(1:nVar)
!===================================================================================================================================
IF(.NOT. ALLOCATED(DQ%VarName)  )ALLOCATE(DQ%VarName(1:nVar))
IF(.NOT. ALLOCATED(DQ%IndGlobal))ALLOCATE(DQ%IndGlobal(1:nVar))
DQ%nVar      =nVar
DQ%Varname(1:nVar)=VarName(1:nVar)
Map  =0
DQ%nVarVisu=0
DO iDQ=1,nVar
  DO iVar=1,nVarVisu
    IF(STRICMP(VarNameVisu(iVar),VarName(iDQ)))THEN
      WRITE(UNIT_StdOut,*) '  ',TRIM(VarNameVisu(iVar))
      DQ%nVarVisu           =DQ%nVarVisu+1
      Map(iDQ)               =DQ%nVarVisu
      DQ%IndGlobal(iDQ)      =iVar
    END IF
  END DO
END DO
IF(DQ%nVarVisu .GT. 0)THEN
  max_nVarVisu=MAX(max_nVarVisu,DQ%nVarVisu)
  ! Only allocate once per run!
  IF(.NOT. ALLOCATED(DQ%Ind))ALLOCATE(DQ%Ind(DQ%nVarVisu))
  DO iDQ=1,DQ%nVar
    IF(Map(iDQ) .GT. 0) DQ%Ind(Map(iDQ))=iDQ
  END DO
END IF

END SUBROUTINE CreateStateMappings



END MODULE MOD_VarNameMappingsRP
