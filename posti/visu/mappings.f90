!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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
!> Routines to create several mappings:
!> * Mapping that contains the distribution of DG and FV elements
!> * Mappings from all available variables to the ones that should be calculated and visualized
!===================================================================================================================================
MODULE MOD_Posti_Mappings
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE Build_FV_DG_distribution
  MODULE PROCEDURE Build_FV_DG_distribution
END INTERFACE

INTERFACE Build_mapCalc_mapVisu
  MODULE PROCEDURE Build_mapCalc_mapVisu
END INTERFACE

PUBLIC:: Build_FV_DG_distribution
PUBLIC:: Build_mapCalc_mapVisu

CONTAINS

!===================================================================================================================================
!> This routine determines the distribution of DG and FV elements in the state file.
!>  1. reads the 'ElemData' array in the state-file and fills the 'FV_Elems_loc' array.
!>  2. count the number of DG/FV elements in FV_Elems_loc.
!>  3. build the mappings (mapElems_FV/DG) that hold the global indices of the FV/DG element indices.
!>  4. check wether the distribution of FV elements has changed
!===================================================================================================================================
SUBROUTINE Build_FV_DG_distribution(statefile)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_HDF5_Input  ,ONLY: GetArrayAndName,OpenDataFile,CloseDataFile
USE MOD_ReadInTools ,ONLY: GETSTR,CountOption
USE MOD_StringTools ,ONLY: STRICMP
USE MOD_Mesh_ReadIn ,ONLY: BuildPartition
USE MOD_Mesh_Vars   ,ONLY: nElems
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: statefile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem
INTEGER                        :: nElems_FV_glob
#if FV_ENABLED
INTEGER                        :: iElem2,iVar
INTEGER                        :: nVal(15)
REAL,ALLOCATABLE               :: ElemData_loc(:,:),tmp(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesElemData_loc(:)
#endif
!===================================================================================================================================
! Build partition to get nElems
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL BuildPartition() 
CALL CloseDataFile()

SDEALLOCATE(FV_Elems_loc)
ALLOCATE(FV_Elems_loc(1:nElems))
#if FV_ENABLED
IF (.NOT.DGonly) THEN
  NVisu_FV = (PP_N+1)*2-1

  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL GetArrayAndName('ElemData','VarNamesAdd',nVal,tmp,VarNamesElemData_loc)
  CALL CloseDataFile()
  IF (ALLOCATED(VarNamesElemData_loc)) THEN
    ALLOCATE(ElemData_loc(nVal(1),nVal(2)))
    ElemData_loc = RESHAPE(tmp,(/nVal(1),nVal(2)/))
    ! search for FV_Elems 
    FV_Elems_loc = 0
    DO iVar=1,nVal(1)
      IF (STRICMP(VarNamesElemData_loc(iVar),"FV_Elems")) THEN
        FV_Elems_loc = INT(ElemData_loc(iVar,:))
      END IF
    END DO
    DEALLOCATE(ElemData_loc,VarNamesElemData_loc,tmp)
  END IF

  nElems_FV = SUM(FV_Elems_loc)
  nElems_DG = nElems - nElems_FV

  ! build the mapping, that holds the global indices of all FV elements
  SDEALLOCATE(mapElems_FV)
  ALLOCATE(mapElems_FV(1:nElems_FV))
  iElem2 =1
  DO iElem=1,nElems
    IF (FV_Elems_loc(iElem).EQ.1) THEN
      mapElems_FV(iElem2) = iElem
      iElem2 = iElem2 + 1
    END IF
  END DO ! iElem

  ! build the mapping, that holds the global indices of all DG elements
  SDEALLOCATE(mapElems_DG)
  ALLOCATE(mapElems_DG(1:nElems_DG))
  iElem2 =1
  DO iElem=1,nElems
    IF (FV_Elems_loc(iElem).EQ.0) THEN
      mapElems_DG(iElem2) = iElem
      iElem2 = iElem2 + 1
    END IF
  END DO ! iElem
ELSE
#endif
  FV_Elems_loc = 0
  nElems_DG = nElems
  nElems_FV = 0
  NVisu_FV = 0

  ! build the mapping, that holds the global indices of all DG elements
  SDEALLOCATE(mapElems_DG)
  ALLOCATE(mapElems_DG(1:nElems_DG))
  DO iElem=1,nElems
    mapElems_DG(iElem) = iElem
  END DO
#if FV_ENABLED
END IF
#endif


! check if the distribution of DG/FV elements has changed
changedFV_Elems=.TRUE.
IF (ALLOCATED(FV_Elems_old).AND.(SIZE(FV_Elems_loc).EQ.SIZE(FV_Elems_old))) THEN
  changedFV_Elems = .NOT.ALL(FV_Elems_loc.EQ.FV_Elems_old) 
END IF
SDEALLOCATE(FV_Elems_old)
ALLOCATE(FV_Elems_old(1:nElems))
FV_Elems_old = FV_Elems_loc

nElems_FV_glob = SUM(FV_Elems_loc)
#if USE_MPI
! check if any processor has FV elements
CALL MPI_ALLREDUCE(MPI_IN_PLACE,nElems_FV_glob,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,iError)
#endif
hasFV_Elems = (nElems_FV_glob.GT.0)


END SUBROUTINE Build_FV_DG_distribution

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
USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_ReadInTools     ,ONLY: GETSTR,CountOption
USE MOD_StringTools     ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar,iVar2
CHARACTER(LEN=255)  :: VarName
CHARACTER(LEN=255)  :: BoundaryName
CHARACTER(LEN=20)   :: format
!===================================================================================================================================
! Read Varnames from parameter file and fill
!   mapVisu = map, which stores at position x the position/index of the x.th quantity in the UVisu array
!             if a quantity is not visualized it is zero
SDEALLOCATE(mapVisu)
ALLOCATE(mapVisu(1:nVarTotal))
mapVisu = 0
nVarVisuTotal = 0
! Compare varnames that should be visualized with availabe varnames
DO iVar=1,CountOption("VarName")
  VarName = GETSTR("VarName")
  DO iVar2=1,nVarTotal
    IF (STRICMP(VarName, VarNamesTotal(iVar2))) THEN
      mapVisu(iVar2) = nVarVisuTotal+1
      nVarVisuTotal = nVarVisuTotal + 1
    END IF
  END DO
END DO

! check whether gradients are needed for any quantity
DO iVar=1,nVarDep
  IF (mapVisu(iVar).GT.0) THEN
    withDGOperator = withDGOperator .OR. (DepTable(iVar,0).GT.0)
  END IF
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

! print the dependecy table
SWRITE(*,*) "Dependencies: ", withDGOperator
WRITE(format,'(I2)') SIZE(DepTable,2)
DO iVar=1,nVarDep
  SWRITE (*,'('//format//'I2,A)') DepTable(iVar,:), " "//TRIM(VarNamesTotal(iVar))
END DO

! Build :
!   mapCalc = map, which stores at position x the position/index of the x.th quantity in the UCalc array
!             if a quantity is not calculated it is zero
SDEALLOCATE(mapCalc)
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

! check if any varnames changed
changedVarNames = .TRUE.
IF (ALLOCATED(mapVisu_old).AND.(SIZE(mapVisu).EQ.SIZE(mapVisu_old))) THEN
  changedVarNames = .NOT.ALL(mapVisu.EQ.mapVisu_old) 
END IF
changedVarNames = changedVarNames
SDEALLOCATE(mapVisu_old)
ALLOCATE(mapVisu_old(1:nVarTotal))
mapVisu_old = mapVisu

! print the mappings
WRITE(format,'(I2)') nVarTotal
SWRITE (*,'(A,'//format//'I3)') "mapCalc ",mapCalc
SWRITE (*,'(A,'//format//'I3)') "mapVisu ",mapVisu

SDEALLOCATE(mapBCNames)
ALLOCATE(mapBCNames(1:nBCNamesTotal))
mapBCNames = 0
nBCNamesVisu = 0
! Compare varnames that should be visualized with availabe varnames
DO iVar=1,CountOption("BoundaryName")
  BoundaryName = GETSTR("BoundaryName")
  DO iVar2=1,nBCNamesTotal
    IF (STRICMP(BoundaryName, BoundaryNamesTotal(iVar2))) THEN
      mapBCNames(iVar2) = nBCNamesVisu+1
      nBCNamesVisu = nBCNamesVisu + 1
    END IF
  END DO
END DO

WRITE(format,'(I2)') nBCNamesTotal
SWRITE (*,'(A,'//format//'I3)') "mapBCNames ",mapBCNames

END SUBROUTINE Build_mapCalc_mapVisu

END MODULE MOD_Posti_Mappings
