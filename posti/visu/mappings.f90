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
!> * Mappings from all available variables to the ones that should be calculated and visualized for volume and surface
!> * Mappings from the available boundaries to the ones that should be visualized
!===================================================================================================================================
MODULE MOD_Posti_Mappings
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE Build_FV_DG_distribution
  MODULE PROCEDURE Build_FV_DG_distribution
END INTERFACE

INTERFACE Build_mapDepToCalc_mapAllVarsToVisuVars
  MODULE PROCEDURE Build_mapDepToCalc_mapAllVarsToVisuVars
END INTERFACE

INTERFACE Build_mapBCSides
  MODULE PROCEDURE Build_mapBCSides
END INTERFACE

PUBLIC:: Build_FV_DG_distribution
PUBLIC:: Build_mapDepToCalc_mapAllVarsToVisuVars
PUBLIC:: Build_mapBCSides

CONTAINS

!===================================================================================================================================
!> This routine determines the distribution of DG and FV elements in the state file.
!>  1. reads the 'ElemData' array in the state-file and fills the 'FV_Elems_loc' array.
!>  2. count the number of DG/FV elements in FV_Elems_loc.
!>  3. build the mappings (mapFVElemsToAllElems/DG) that hold the global indices of the FV/DG element indices.
!>  4. check wether the distribution of FV elements has changed
!>  5. store the current distribution in FV_Elems_loc, is needed for some mappings later
!===================================================================================================================================
SUBROUTINE Build_FV_DG_distribution(&
#if FV_ENABLED
    statefile&
#endif
    )
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_HDF5_Input  ,ONLY: GetArrayAndName,OpenDataFile,CloseDataFile
USE MOD_ReadInTools ,ONLY: GETSTR,CountOption
USE MOD_StringTools ,ONLY: STRICMP
USE MOD_Mesh_ReadIn ,ONLY: BuildPartition
USE MOD_Mesh_Vars   ,ONLY: nElems
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
#if FV_ENABLED
CHARACTER(LEN=255),INTENT(IN)  :: statefile
#endif
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
#if FV_RECONSTRUCT
  NCalc_FV = NVisu_FV
#else
  NCalc_FV = PP_N
#endif

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

  nElems_FV = COUNT(FV_Elems_loc.EQ.1)
  nElems_DG = nElems - nElems_FV

  ! build the mapping, that holds the global indices of all FV elements
  SDEALLOCATE(mapFVElemsToAllElems)
  ALLOCATE(mapFVElemsToAllElems(1:nElems_FV))
  iElem2 =1
  DO iElem=1,nElems
    IF (FV_Elems_loc(iElem).EQ.1) THEN
      mapFVElemsToAllElems(iElem2) = iElem
      iElem2 = iElem2 + 1
    END IF
  END DO ! iElem

  ! build the mapping, that holds the global indices of all DG elements
  SDEALLOCATE(mapDGElemsToAllElems)
  ALLOCATE(mapDGElemsToAllElems(1:nElems_DG))
  iElem2 =1
  DO iElem=1,nElems
    IF (FV_Elems_loc(iElem).EQ.0) THEN
      mapDGElemsToAllElems(iElem2) = iElem
      iElem2 = iElem2 + 1
    END IF
  END DO ! iElem
ELSE
#endif
  FV_Elems_loc = 0
  nElems_DG = nElems
  nElems_FV = 0
  NVisu_FV = 1
#if FV_RECONSTRUCT
  NCalc_FV = NVisu_FV
#else
  NCalc_FV = 0
#endif

  ! build the mapping, that holds the global indices of all DG elements
  SDEALLOCATE(mapDGElemsToAllElems)
  ALLOCATE(mapDGElemsToAllElems(1:nElems_DG))
  DO iElem=1,nElems
    mapDGElemsToAllElems(iElem) = iElem
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
!>  2. Initialize the dependency table
!>  3. check wether gradients are needed for any quantity. If this is the case, remove the conservative quantities from the
!>     dependecies of the primitive quantities (the primitive quantities are available directly, since the DGTimeDerivative_weakForm
!>     will be executed.
!>  4. build the 'mapDepToCalc' that holds for each quantity that will be calculated the index in 'UCalc' array (0 if not calculated)
!>  5. build the 'mapAllVarsToVisuVars' that holds for each quantity that will be visualized the index in 'UVisu' array (0 if not visualized)
!>  6. build the 'mapAllBCNamesToVisuBCNames' that holds for each available boundary the visualization index
!===================================================================================================================================
SUBROUTINE Build_mapDepToCalc_mapAllVarsToVisuVars()
USE MOD_Globals
USE MOD_Visu_Vars
USE MOD_Restart_Vars    ,ONLY: RestartMode
USE MOD_ReadInTools     ,ONLY: GETSTR,GETLOGICAL,CountOption
USE MOD_StringTools     ,ONLY: STRICMP,set_formatting,clear_formatting
#if FV_RECONSTRUCT
USE MOD_EOS_Posti       ,ONLY: AppendNeededPrims
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar,iVar2
CHARACTER(LEN=255)  :: VarName
CHARACTER(LEN=255)  :: BoundaryName
CHARACTER(LEN=20)   :: format
LOGICAL             :: UseVarNamesHDF5=.FALSE.,UseVarNamesAll=.FALSE.
!===================================================================================================================================
! Read Varnames from parameter file and fill
!   mapAllVarsToVisuVars = map, which stores at position x the position/index of the x.th quantity in the UVisu array
!             if a quantity is not visualized it is zero
SDEALLOCATE(mapAllVarsToVisuVars)
SDEALLOCATE(mapAllVarsToSurfVisuVars)
ALLOCATE(mapAllVarsToVisuVars    (1:nVarAll))
ALLOCATE(mapAllVarsToSurfVisuVars(1:nVarAll))
mapAllVarsToVisuVars     = 0
mapAllVarsToSurfVisuVars = 0
nVarVisu        = 0
nVarSurfVisuAll = 0

! If no variable names are given in prm file, take the variables given in the HDF5 "VarNames" attribute (if present) or all found
! variables (else). This default can be suppressed via the "noVisuVars" flag (used e.g. in paraview plugin prm files)
IF ((nVarIni.EQ.0).AND..NOT.GETLOGICAL("noVisuVars")) THEN
  IF (ALLOCATED(VarNamesHDF5)) THEN
    IF (RestartMode.LE.1) THEN
      nVarIni         = SIZE(VarNamesHDF5)
      UseVarNamesHDF5 = .TRUE.
    ELSE
      nVarIni         = nVarAll
      UseVarNamesAll  = .TRUE.
    END IF
  ELSE
    nVarIni         = nVarAll
    UseVarNamesAll  = .TRUE.
  END IF
END IF

! Compare varnames that should be visualized with available varnames
DO iVar=1,nVarIni
  ! get var name from prm file or resort to the default vars
  IF(UseVarNamesHDF5) THEN
    VarName=VarNamesHDF5(iVar)
  ELSEIF(UseVarNamesAll)THEN
    VarName=VarNamesAll(iVar)
    IF(TRIM(VarName).EQ."DG_Solution:DUMMY_DO_NOT_VISUALIZE") CYCLE
  ELSE
    VarName = GETSTR("VarName")
  END IF
  DO iVar2=1,nVarAll
    IF (STRICMP(VarName, VarnamesAll(iVar2))) THEN
      nVarSurfVisuAll = nVarSurfVisuAll + 1
      mapAllVarsToSurfVisuVars(iVar2) = nVarSurfVisuAll
      IF (iVar2.LE.nVarDep) THEN
        IF(DepSurfaceOnly(iVar2).EQ.1) CYCLE
      END IF
      nVarVisu = nVarVisu + 1
      mapAllVarsToVisuVars(iVar2) = nVarVisu
    END IF
  END DO
END DO

! Warn the user when restarting from a time-averaged file and no mapping was found
IF (nVarVisu.EQ.0 .AND. RestartMode.EQ.0) THEN
  CALL PrintWarning('Trying to visualize a time-averaged file but no corresponding variable was found.\n'//&
                    'Derived quantities are not available if file is missing conservative/primite variables.\n'//&
                    'To access mean quantites, please prepend VarName with "Mean:VarName", "MeanSquare:VarName" or "Fluc:VarName"!')
END IF

! check whether gradients are needed for any quantity
DO iVar=1,nVarDep
  IF (mapAllVarsToSurfVisuVars(iVar).GT.0) THEN
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

! print the dependency table
IF (.NOT.DependenciesOutputDone) THEN
  SWRITE(UNIT_stdOut,'(A,L1)') "Dependencies: ", withDGOperator
  WRITE(format,'(I2)') SIZE(DepTable,2)
  DO iVar=1,nVarDep
    SWRITE (UNIT_stdOut,'('//format//'I2,A)') DepTable(iVar,:), " "//TRIM(VarnamesAll(iVar))
  END DO
END IF

! Build :
!   mapDepToCalc = map, which stores at position x the position/index of the x.th quantity in the UCalc array
!             if a quantity is not calculated it is zero
SDEALLOCATE(mapDepToCalc)
ALLOCATE(mapDepToCalc(1:nVarDep))
mapDepToCalc = 0
DO iVar=1,nVarDep
  IF (mapAllVarsToSurfVisuVars(iVar).GT.0) THEN
    mapDepToCalc = MAX(mapDepToCalc,DepTable(iVar,1:nVarDep))
  END IF
END DO
! enumerate mapDepToCalc
nVarCalc = 0
DO iVar=1,nVarDep
  IF (mapDepToCalc(iVar).GT.0) THEN
    nVarCalc = nVarCalc + 1
    mapDepToCalc(iVar) = nVarCalc
  END IF
END DO

! check if any varnames changed
changedVarNames = .TRUE.
IF (ALLOCATED(mapAllVarsToSurfVisuVars_old).AND.(SIZE(mapAllVarsToSurfVisuVars).EQ.SIZE(mapAllVarsToSurfVisuVars_old))) THEN
  changedVarNames = .NOT.ALL(mapAllVarsToSurfVisuVars.EQ.mapAllVarsToSurfVisuVars_old)
END IF
SDEALLOCATE(mapAllVarsToSurfVisuVars_old)
ALLOCATE(mapAllVarsToSurfVisuVars_old(1:nVarAll))
mapAllVarsToSurfVisuVars_old = mapAllVarsToSurfVisuVars


SDEALLOCATE(mapDepToCalc_FV)
ALLOCATE(mapDepToCalc_FV(1:nVarDep))
mapDepToCalc_FV = mapDepToCalc
nVarCalc_FV = nVarCalc
#if FV_RECONSTRUCT
! generate a new mapDepToCalc_FV, which is a copy of the original mapDepToCalc but is extended in the following way.
! Since the reconstruction is performed in primitive quantities, the calculation of conservative quantities from them
! introduce for the conservatives dependcies from the primitive ones. Therefore all primitive quantities that
! are needed to build the requested conservatives must be added to the mapDepToCalc_FV.
IF (StateFileMode) CALL AppendNeededPrims(mapDepToCalc,mapDepToCalc_FV,nVarCalc_FV)
#endif

!---------------------- Surface visualization ----------------------------!

! Build the mapping for the surface visualization
! mapAllBCNamesToVisuBCNames(iBC) stores the ascending visualization index of the all boundaries. 0 means no visualization.
! nBCNamesVisu is the number of boundaries to be visualized.
SDEALLOCATE(mapAllBCNamesToVisuBCNames)
ALLOCATE(mapAllBCNamesToVisuBCNames(1:nBCNamesAll))
mapAllBCNamesToVisuBCNames = 0
nBCNamesVisu = 0
! Compare boundary names that should be visualized with available varnames
DO iVar=1,CountOption("BoundaryName")
  BoundaryName = GETSTR("BoundaryName")
  DO iVar2=1,nBCNamesAll
    IF (STRICMP(BoundaryName, BCNamesAll(iVar2))) THEN
      mapAllBCNamesToVisuBCNames(iVar2) = nBCNamesVisu+1
      nBCNamesVisu = nBCNamesVisu + 1
    END IF
  END DO
END DO

! Set flag indicating if surface visualization is needed
doSurfVisu = nBCNamesVisu.GT.0

! check if any boundary changed
changedBCnames = .TRUE.
IF (ALLOCATED(mapAllBCNamesToVisuBCNames_old).AND.(SIZE(mapAllBCNamesToVisuBCNames).EQ.SIZE(mapAllBCNamesToVisuBCNames_old))) THEN
  changedBCnames = .NOT.ALL(mapAllBCNamesToVisuBCNames.EQ.mapAllBCNamesToVisuBCNames_old)
  IF (ALL(mapAllBCNamesToVisuBCNames_old.EQ.0) .AND. changedBCNames) changedStateFile = .TRUE.
END IF
SDEALLOCATE(mapAllBCNamesToVisuBCNames_old)
ALLOCATE(mapAllBCNamesToVisuBCNames_old(1:nBCNamesAll))
mapAllBCNamesToVisuBCNames_old = mapAllBCNamesToVisuBCNames

! print the mappings
IF (.NOT.DependenciesOutputDone) THEN
  WRITE(format,'(I0)') nVarAll
  SWRITE (*,'(A,'//format//'I0)') "mapDepToCalc             ",mapDepToCalc
  SWRITE (*,'(A,'//format//'I0)') "mapDepToCalc_FV          ",mapDepToCalc_FV
  SWRITE (*,'(A,'//format//'I0)') "mapAllVarsToVisuVars     ",mapAllVarsToVisuVars
  SWRITE (*,'(A,'//format//'I0)') "mapAllVarsToSurfVisuVars ",mapAllVarsToSurfVisuVars
  SWRITE (*,'(A,'//format//'I0)') "mapAllBCNamesToVisuBCNames ",mapAllBCNamesToVisuBCNames
  DependenciesOutputDone = .TRUE.
END IF

END SUBROUTINE Build_mapDepToCalc_mapAllVarsToVisuVars

!===================================================================================================================================
!> This routine builds mappings that give for each BC side the index of the visualization side, seperate by FV and DG.
!> We also store the number of visualization sides per BC name.
!===================================================================================================================================
SUBROUTINE Build_mapBCSides()
USE MOD_Visu_Vars
USE MOD_Mesh_Vars   ,ONLY: nBCSides,BC,SideToElem,BoundaryName
USE MOD_StringTools ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iBC,iElem,iSide
!===================================================================================================================================

! Build surface visualization mappings.
! mapAllBCSidesToDGVisuBCSides/FV(iBCSide) contains the ascending index of the visualization boundary sides. They are sorted after the boundary name.
! 0 means no visualization of this boundary side.
! nSidesPerBCNameVisu_DG/FV(iBCNamesVisu) contains how many boundary sides belong to each boundary that should be visualized.
SDEALLOCATE(mapAllBCSidesToDGVisuBCSides)
SDEALLOCATE(mapAllBCSidesToFVVisuBCSides)
SDEALLOCATE(nSidesPerBCNameVisu_DG)
SDEALLOCATE(nSidesPerBCNameVisu_FV)
ALLOCATE(mapAllBCSidesToDGVisuBCSides(1:nBCSides))
ALLOCATE(mapAllBCSidesToFVVisuBCSides(1:nBCSides))
ALLOCATE(nSidesPerBCNameVisu_DG(1:nBCNamesVisu))
ALLOCATE(nSidesPerBCNameVisu_FV(1:nBCNamesVisu))
mapAllBCSidesToDGVisuBCSides = 0
mapAllBCSidesToFVVisuBCSides = 0
nSidesPerBCNameVisu_DG = 0
nSidesPerBCNameVisu_FV = 0
nBCSidesVisu_DG = 0
nBCSidesVisu_FV = 0
IF (doSurfVisu) THEN
  DO iBC=1,nBCNamesAll    ! iterate over all bc names
    IF (mapAllBCNamesToVisuBCNames(iBC).GT.0) THEN
      DO iSide=1,nBCSides   ! iterate over all bc sides
        IF (STRICMP(BoundaryName(BC(iSide)),BCNamesAll(iBC))) THEN ! check if side is of specific boundary name
          iElem = SideToElem(S2E_ELEM_ID,iSide)
          IF (FV_Elems_loc(iElem).EQ.0)  THEN ! DG element
            nBCSidesVisu_DG = nBCSidesVisu_DG + 1
            mapAllBCSidesToDGVisuBCSides(iSide) = nBCSidesVisu_DG
            nSidesPerBCNameVisu_DG(mapAllBCNamesToVisuBCNames(iBC)) = nSidesPerBCNameVisu_DG(mapAllBCNamesToVisuBCNames(iBC)) + 1
          ELSE ! FV Element
            nBCSidesVisu_FV = nBCSidesVisu_FV + 1
            mapAllBCSidesToFVVisuBCSides(iSide) = nBCSidesVisu_FV
            nSidesPerBCNameVisu_FV(mapAllBCNamesToVisuBCNames(iBC)) = nSidesPerBCNameVisu_FV(mapAllBCNamesToVisuBCNames(iBC)) + 1
          END IF
        END IF
      END DO
    END IF
  END DO
END IF ! doSurfVisu

END SUBROUTINE Build_mapBCSides

END MODULE MOD_Posti_Mappings
