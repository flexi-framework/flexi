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

MODULE MOD_Posti_Mappings
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
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
USE MOD_ReadInTools ,ONLY: GETSTR,CountOption
USE MOD_StringTools ,ONLY: STRICMP
USE MOD_Mesh_Vars   ,ONLY: nElems
#if FV_ENABLED
USE MOD_Restart_Vars,ONLY: nVarElemData,ElemData,VarNamesElemData
#endif
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: statefile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem,iElem2,iVar
INTEGER           :: nElems_FV_glob
!===================================================================================================================================
SDEALLOCATE(FV_Elems_loc)
ALLOCATE(FV_Elems_loc(1:nElems))
#if FV_ENABLED
NVisu_FV = (PP_N+1)*2-1

FV_Elems_loc = 0
DO iVar=1,nVarElemData
  IF (STRICMP(VarNamesElemData(iVar),"FV_Elems")) THEN
    FV_Elems_loc = INT(ElemData(iVar,:))
  END IF
END DO

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

! build the mapping, that holds the global indices of all FV elements
SDEALLOCATE(mapElems_DG)
ALLOCATE(mapElems_DG(1:nElems_DG))
iElem2 =1
DO iElem=1,nElems
  IF (FV_Elems_loc(iElem).EQ.0) THEN
    mapElems_DG(iElem2) = iElem
    iElem2 = iElem2 + 1
  END IF
END DO ! iElem

#else 
FV_Elems_loc = 0
nElems_DG = nElems  
nElems_FV = 0
NVisu_FV = 0

! build the mapping, that holds the global indices of all FV elements
SDEALLOCATE(mapElems_DG)
ALLOCATE(mapElems_DG(1:nElems_DG))
DO iElem=1,nElems
  mapElems_DG(iElem) = iElem
END DO
#endif


! check if the distribution of FV elements has changed
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
USE MOD_EOS_Posti_Vars  ,ONLY: nVarTotal,DepNames,DepTable
USE MOD_EOS_Posti       ,ONLY: InitDepTable,FillDepTable
USE MOD_Restart_Vars    ,ONLY: nVarElemData,VarNamesElemData
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar,iVar2
CHARACTER(LEN=255)  :: VarName
CHARACTER(LEN=20)   :: format
LOGICAL             :: changedVarNames_ElemData
!===================================================================================================================================
! initialize the dependency table
CALL InitDepTable()

! Read Varnames from parameter file and fill
!   mapVisu = map, which stores at position x the position/index of the x.th quantity in the UVisu array
!             if a quantity is not visualized it is zero
SDEALLOCATE(mapVisu)
ALLOCATE(mapVisu(1:nVarTotal))
SDEALLOCATE(VarNamesVisu_ElemData)
ALLOCATE(VarNamesVisu_ElemData(1:CountOption("VarName")))
mapVisu = 0
nVarVisu = 0
nVarVisu_ElemData = 0
DO iVar=1,CountOption("VarName")
  VarName = GETSTR("VarName")
  DO iVar2=1,nVarTotal
    IF (STRICMP(VarName, DepNames(iVar2))) THEN
      mapVisu(iVar2) = nVarVisu+1
      nVarVisu = nVarVisu + 1
    END IF
  END DO
  DO iVar2=1,nVarElemData
    IF (STRICMP(VarName, VarNamesElemData(iVar2))) THEN
      nVarVisu_ElemData = nVarVisu_ElemData + 1
      VarNamesVisu_ElemData(nVarVisu_ElemData) = VarName
    END IF
  END DO
END DO

changedVarNames_ElemData = (nVarVisu_ElemData.NE.nVarVisu_ElemData_old)
IF (.NOT.changedVarNames_ElemData) THEN
  DO iVar=1,nVarVisu_ElemData

    changedVarNames_ElemData = changedVarNames_ElemData &
        .OR.  (.NOT.STRICMP(VarNamesVisu_ElemData(iVar),VarNamesVisu_ElemData_old(iVar)))
  END DO
END IF
SDEALLOCATE(VarNamesVisu_ElemData_old)
ALLOCATE(VarNamesVisu_ElemData_old(SIZE(VarNamesVisu_ElemData)))
VarNamesVisu_ElemData_old = VarNamesVisu_ElemData

! check whether gradients are needed for any quantity
DO iVar=1,nVarTotal
  IF (mapVisu(iVar).GT.0) THEN
    withGradients = withGradients .OR. (DepTable(iVar,0).GT.0)
  END IF
END DO

! calculate all dependencies 
! if withGradients==True then the DGTimeDerivative_weakForm routine is called and the primitive
! quantities are computed in this operator, which can be used for calculation of other quantities
! directly. Therefore the dependecies of the primitive variables on the conservative variables 
! can be removed.
CALL FillDepTable(withGradients)

! print the dependecy table
SWRITE(*,*) "Dependencies: ", withGradients
WRITE(format,'(I2)') SIZE(DepTable,2)
DO iVar=1,nVarTotal
  SWRITE (*,'('//format//'I3,A)') DepTable(iVar,:), TRIM(DepNames(iVar))
END DO

! Build :
!   mapCalc = map, which stores at position x the position/index of the x.th quantity in the UCalc array
!             if a quantity is not calculated it is zero
SDEALLOCATE(mapCalc)
ALLOCATE(mapCalc(1:nVarTotal))
mapCalc = 0
DO iVar=1,nVarTotal
  IF (mapVisu(iVar).GT.0) THEN
    mapCalc = MAX(mapCalc,DepTable(iVar,1:nVarTotal))
  END IF
END DO
! enumerate mapCalc
nVarCalc = 0
DO iVar=1,nVarTotal
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
changedVarNames = changedVarNames .OR. changedVarNames_ElemData
SDEALLOCATE(mapVisu_old)
ALLOCATE(mapVisu_old(1:nVarTotal))
mapVisu_old = mapVisu

! print the mappings
WRITE(format,'(I2)') nVarTotal
SWRITE (*,'(A,'//format//'I3)') "mapCalc ",mapCalc
SWRITE (*,'(A,'//format//'I3)') "mapVisu ",mapVisu

END SUBROUTINE Build_mapCalc_mapVisu

END MODULE MOD_Posti_Mappings
