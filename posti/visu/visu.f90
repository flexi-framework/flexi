!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
!> Module containing the main procedures for the visu tool: visu_requestInformation is called by ParaView to create a
!> list of available variables and visu is the main routine which is either called by ParaView to get the data it visualizes
!> or by the standalone tool.
!===================================================================================================================================
MODULE MOD_Visu
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: visu
PUBLIC:: FinalizeVisu
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Main routine of the visualization tool visu. Called either by the ParaView plugin or by the standalone program version.
!===================================================================================================================================
SUBROUTINE visu(mpi_comm_IN, prmfile, postifile, statefile)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input          ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,OpenDataFile,CloseDataFile
USE MOD_Interpolation_Vars  ,ONLY: NodeType,NodeTypeVISUFVEqui
USE MOD_IO_HDF5             ,ONLY: InitMPIInfo
USE MOD_MPI                 ,ONLY: InitMPI
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_DG,CalcSurfQuantities_DG
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG,ConvertToSurfVisu_DG,ConvertToVisu_GenericData
USE MOD_Posti_ReadMesh      ,ONLY: VisualizeMesh
USE MOD_Posti_ReadState     ,ONLY: ReadState
USE MOD_Posti_Mappings      ,ONLY: Build_mapBCSides
USE MOD_Posti_VisuMesh      ,ONLY: BuildVisuCoords,BuildSurfVisuCoords
USE MOD_ReadInTools         ,ONLY: prms,addStrListEntry
USE MOD_ReadInTools         ,ONLY: FinalizeParameters,ExtractParameterFile,PrintDefaultParameterFile
USE MOD_Restart_Vars        ,ONLY: RestartMode
USE MOD_StringTools         ,ONLY: STRICMP,set_formatting,clear_formatting
USE MOD_Visu_Avg2D          ,ONLY: Average2D,WriteAverageToHDF5
USE MOD_Visu_Init           ,ONLY: visu_getVarNamesAndFileType,visu_InitFile
USE MOD_Visu_Vars
USE MOD_EOS_Posti_Vars      ,ONLY: DepNames,nVarDepEOS
#if FV_ENABLED
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_FV,CalcSurfQuantities_FV
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV,ConvertToSurfVisu_FV
#endif /*FV_ENABLED*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
TYPE(MPI_Comm),INTENT(IN)        :: mpi_comm_IN
CHARACTER(LEN=255),INTENT(INOUT) :: prmfile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
CHARACTER(LEN=255),INTENT(IN)    :: statefile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: changedPrmFile
INTEGER                          :: iVar
CHARACTER(LEN=2047)              :: str
!===================================================================================================================================

!**********************************************************************************************
! General workflow / principles of the visu ParaView-plugin
!
! * all arrays are SDEALLOCATEd just before they are allocated. This is done to keep their
!   content during successive calls of the visu during a ParaView session. They are only
!   deallocated and reallocated, if there content should change. For example the coords of the
!   mesh file only change if the mesh, NVisu or the distribution of DG/FV elements changes.
!
! VISUALIZE MESH: Call 'VisualizeMesh' routine, which reads the mesh, interpolates it to
!   the visu grid and writes it to VTK 3D arrays.
!
! VISUALIZE STATE:
! * There are two different modes:
!   - without gradients: All quantity that should be visualized can be computed without any
!                        gradient computation (mainly the conservative and primitive quantities)
!                        If FV_RECONSTRUCT is enabled and there are FV elements present in
!                        the state, then this mode is not available.
!                        U is read from the state file directly and the EOS is initialized to
!                        perform ConsToPrim conversions based only on the conservative state U.
!
!   - with gradiens: There are quantities that require the computation of gradients or there
!                    are FV elements with FV_RECONSTRUCT enabled. In this case the DG operator
!                    'DGTimeDerivative_weakForm' is called once to fill the gradients and the
!                    reconstruction of the FV subcell method.
!                    This requires the initialization of several modules of the FLEXI.
!                    U is read via a call of 'Restart'. In the DGTimeDerivative_weakForm the
!                    primitive quantities U_Prim and gradUx/y/z as well as gradUxi/eta/zeta are
!                    filled. These are used to calculate the visu-quantities.
!
! * The calculation of derived quantities is performed on a arbitrary polynomial degree
!   NCalc and afterwards interpolated to NVisu. Default is PP_N.
!   This is not the case for FV elements with FV_RECONSTRUCT enabled.
!   These require to reconstruct the solution first to the visu grid and afterwards can
!   calculate the derived quantities on the NVisu_FV grid.
!
! * The dependencies of the visu-quantities on the state-quantities is stored in a dependency
!   integer table 'DepTable' and it corresponding row/col names 'DepNames' (see eos_vars.f90).
!   The string vector 'DepNames' contains all available quantities available for visualization
!   and the 'DepTable' contains in a row the dependencies of a quantity on other quantities.
!
! * The calculation of the visu-quantities is done in two steps:
!   1. calculate all quantities needed for the visu-quantities (stored in UCalc)
!   2. pick the visu-quantities from UCalc and interpolate them to NVisu (stored in UVisu)
!
! * Therefore two mappings from all available quantities to the calc-quantities and the
!   visu-quantities exist:
!   - 'mapAllVarsToVisuVars' is a integer array of size (1:nVarAll), where nVarAll is the total amount
!     of available quantities. This map contains a zero for all not-to-visu-quantities and for
!     all quantities the index where it is stored in 'UVisu'.
!     This mapping is filled from the 'VarName' entries in the parameter file.
!   - 'mapDepToCalc' is the same as mapAllVarsToVisuVars, but for all (intermediate) quantities stored in 'UCalc'.
!     This mapping is filled from the DepTable.
!
! CHANGED system:
! * There are different logical changedXXX variables, which indicated if XXX changed during
!   successive calls of the visu. These variables control the general workflow of the visu.
!   - changedStateFile:     new state file
!   - changedMeshFile:      new mesh file (only possible if changedStateFile==TRUE)
!   - changedVarNames:      new set of variables to visualize
!   - changedNVisu:         new NVisu, new Nodetype
!   - changedFV_Elems:      new distribution of FV/DG elements (only if changedStateFile==TRUE)
!   - changedWithDGOperator: different mode, with/without gradients
!   - changedDGonly:        the visualization of FV elements as DG elements was set or unset
!   - changedNCalc:         the polynomial degree used for calculations changed
!
! WORKFLOW:
! * The main steps are:
!   1. call InitFile for FV/DG distribution and mappings
!   2. read solution          (if changedStateFile or changedWithDGOperator or changedDGonly)
!   3. build mapping for BC sides that should be visualized (done after read solution since some
!      mesh infos are needed)
!   4. read Mesh              (if changedMeshFile)
!   5. compute UCalc          (if changedStateFile or changedVarNames or changedDGonly or changedNCalc)
!   6. convert to UVisu       (if changedStateFile or changedVarNames or changedNVisu or changedDGonly or changedNCalc)
!   7. build visu mesh        (if changedMeshFile  or changedNVisu or changedFV_Elems or changedDGonly)
!   5. - 7. are done seperately for surface variables if surface visualization is turned on
!   8. write VTK arrays       (always!)
!
!**********************************************************************************************

CALL SetStackSizeUnlimited()
postiMode = .TRUE. ! Flag used in FLEXI routines to do things only for POSTI usage
CALL InitMPI(mpi_comm_IN)
CALL InitMPIInfo()

CALL FinalizeParameters()
! Read Varnames to visualize and build calc and visu dependencies
CALL prms%SetSection("posti")
CALL prms%CreateStringOption(       "MeshFile"        , "Custom mesh file ")
CALL prms%CreateStringOption(       "OutputDirectory" , "Custom output directory")
CALL prms%CreateIntFromStringOption("OutputFormat"    , "File format for visualization: None,ParaView, HDF5."                           &
                                                      , 'paraview')
CALL addStrListEntry('OutputFormat','none'            , OUTPUTFORMAT_NONE)
! CALL addStrListEntry('OutputFormat','tecplot'         , OUTPUTFORMAT_TECPLOT)
! CALL addStrListEntry('OutputFormat','tecplotascii'    , OUTPUTFORMAT_TECPLOTASCII)
CALL addStrListEntry('OutputFormat','paraview'        , OUTPUTFORMAT_PARAVIEW)
CALL addStrListEntry('OutputFormat','hdf5'            , OUTPUTFORMAT_HDF5)
! Build list of all available variables for visualization
str = "Variables that should be visualized. Available are:" // NEW_LINE('A')
DO iVar=1, nVarDepEOS
  str = TRIM(str)//" - "//TRIM(DepNames(iVar))//NEW_LINE('A')
END DO
CALL prms%CreateStringOption(       "VarName"         , str , multiple=.TRUE.)
CALL prms%CreateLogicalOption(      "noVisuVars"      , "If no VarNames are given, this flags supresses visu of standard variables"     &
                                                      , '.FALSE.')
CALL prms%CreateIntOption(          "NVisu"           , "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateIntOption(          "NCalc"           , "Polynomial degree at which calculations are done.")
CALL prms%CreateLogicalOption(      "Avg2D"           , "Average solution in z-direction"                                               &
                                                      , '.FALSE.')
CALL prms%CreateStringOption(       "NodeTypeVisu"    , "NodeType for visualization. Visu, Gauss,Gauss-Lobatto,Visu_inner"              &
                                                      , 'VISU')
CALL prms%CreateLogicalOption(      "DGonly"          , "Visualize FV elements as DG elements."                                         &
                                                      , '.FALSE.')
CALL prms%CreateStringOption(       "BoundaryName"    , "Names of boundaries for surfaces, which should be visualized.", multiple=.TRUE.)
CALL prms%CreateLogicalOption(      "HighOrder"       , "Write high-order element representation"                                       &
                                                      , '.FALSE.')

IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2,statefile) !statefile string conatains --help etc!
  STOP
END IF

SWRITE(UNIT_stdOut,'(A,A)') " READING FROM: ", TRIM(statefile)

changedStateFile      = .FALSE.
changedMeshFile       = .FALSE.
changedNVisu          = .FALSE.
changedNCalc          = .FALSE.
changedVarNames       = .FALSE.
changedFV_Elems       = .FALSE.
changedWithDGOperator = .FALSE.
changedDGonly         = .FALSE.

IF (ISVALIDMESHFILE(statefile)) THEN ! visualize mesh
  SWRITE(UNIT_stdOut,'(A3,A30,A3,A33,A13)')' | ','                   Mode ',' | ','Mesh',' | HDF5    |'
  MeshFileMode = .TRUE.
  MeshFile      = statefile
  nVar_State    = 0
  withDGOperator = .FALSE.
  doSurfVisu     = .FALSE.
  CALL visu_getVarNamesAndFileType(MeshFile,'',VarNamesAll,BCNamesAll)
  CALL VisualizeMesh(postifile,MeshFile)
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! visualize state file
  SWRITE(UNIT_stdOut,'(A3,A30,A3,A33,A13)')' | ','                   Mode ',' | ','State',' | HDF5    |'
  MeshFileMode = .FALSE.
  ! initialize state file
  CALL visu_InitFile(statefile,postifile)

  ! read solution from state file (either direct or including a evaluation of the DG operator)
  IF (LEN_TRIM(prmfile).EQ.0) THEN
    changedPrmFile = .NOT.STRICMP(prmfile_old, ".flexi.ini")
  ELSE
    changedPrmFile = (prmfile .NE. prmfile_old)
  END IF

  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') " DETECTING REQUIRED PARAMETERS..."
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | doSurfVisu              "
  CALL set_formatting(MERGE("blue ","green",doSurfVisu))             ; SWRITE(UNIT_stdOut,'(L1)') doSurfVisu             ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedStateFile        "
  CALL set_formatting(MERGE("blue ","green",changedStateFile))       ; SWRITE(UNIT_stdOut,'(L1)') changedStateFile       ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedMeshFile         "
  CALL set_formatting(MERGE("blue ","green",changedMeshFile))        ; SWRITE(UNIT_stdOut,'(L1)') changedMeshFile        ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedNodeType         "
  CALL set_formatting(MERGE("blue ","green",changedNodeType))        ; SWRITE(UNIT_stdOut,'(L1)') changedNodeType        ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedNVisu            "
  CALL set_formatting(MERGE("blue ","green",changedNVisu))           ; SWRITE(UNIT_stdOut,'(L1)') changedNVisu           ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedNCalc            "
  CALL set_formatting(MERGE("blue ","green",changedNCalc))           ; SWRITE(UNIT_stdOut,'(L1)') changedNCalc           ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedVarNames         "
  CALL set_formatting(MERGE("blue ","green",changedVarNames))        ; SWRITE(UNIT_stdOut,'(L1)') changedVarNames        ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedFV_Elems         "
  CALL set_formatting(MERGE("blue ","green",changedFV_Elems))        ; SWRITE(UNIT_stdOut,'(L1)') changedFV_Elems        ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedWithDGOperator   "
  CALL set_formatting(MERGE("blue ","green",changedWithDGOperator))  ; SWRITE(UNIT_stdOut,'(L1)') changedWithDGOperator  ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedDGonly           "
  CALL set_formatting(MERGE("blue ","green",changedDGonly))          ; SWRITE(UNIT_stdOut,'(L1)') changedDGonly          ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedAvg2D            "
  CALL set_formatting(MERGE("blue ","green",changedAvg2D))           ; SWRITE(UNIT_stdOut,'(L1)') changedAvg2D           ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedPrmFile          "
  CALL set_formatting(MERGE("blue ","green",changedPrmFile))         ; SWRITE(UNIT_stdOut,'(L1)') changedPrmFile         ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedBCNames          "
  CALL set_formatting(MERGE("blue ","green",changedBCNames))         ; SWRITE(UNIT_stdOut,'(L1)') changedBCNames         ; CALL clear_formatting()
  SWRITE(UNIT_StdOut,'(132("-"))')

  IF (changedStateFile.OR.changedWithDGOperator.OR.changedPrmFile.OR.changedDGonly) THEN
    CALL ReadState(prmfile,statefile)
  END IF

  ! build mappings of BC sides for surface visualization
  CALL Build_mapBCSides()

  ! ===== calc solution =====
  IF (changedStateFile.OR.changedWithDGOperator.OR.changedVarNames.OR.changedDGonly.OR.changedNCalc) THEN
    CALL CalcQuantities_DG()
#if FV_ENABLED
    CALL CalcQuantities_FV()
#endif
  END IF
  IF (doSurfVisu) THEN
    ! calc surface solution
    IF (changedStateFile.OR.changedWithDGOperator.OR.changedVarNames.OR.changedDGonly.OR.changedNCalc.OR.changedBCnames) THEN
      CALL CalcSurfQuantities_DG()
#if FV_ENABLED
      CALL CalcSurfQuantities_FV()
#endif
    END IF
  END IF

  ! ===== convert solution to visu grid =====
  IF (changedStateFile.OR.changedWithDGOperator.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedNCalc.OR.changedAvg2D) THEN
    ! ===== Avg2d =====
    IF (Avg2d) THEN
      SDEALLOCATE(UVisu_DG)
      SDEALLOCATE(UVisu_FV)
      ALLOCATE(UVisu_DG(0:NVisu   ,0:NVisu   ,0:0,nElemsAvg2D_DG,nVarVisu))
      ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV,nVarVisu))
      CALL Average2D(nVarCalc    ,nVarCalc_FV  &
                    ,NCalc       ,NCalc_FV     &
                    ,nElems_DG   ,nElems_FV    &
                    ,NodeType                  &
                    ,UCalc_DG    ,UCalc_FV     &
                    ,Vdm_DGToFV  ,Vdm_FVToDG   &
                    ,Vdm_DGToVisu,Vdm_FVToVisu &
                    ,1,nVarDep   ,mapDepToCalc &
                    ,UVisu_DG    ,UVisu_FV)

      SELECT CASE(OutputFormat)
        CASE(OUTPUTFORMAT_HDF5)
          CALL WriteAverageToHDF5(nVarVisu,NVisu,NodeType,OutputTime,MeshFile_state,UVisu_DG &
#if FV_ENABLED
                                 ,NVisu_FV,UVisu_FV                                          &
#endif /* FV_ENABLED */
    )
      END SELECT

    ! .NOT. Avg2D
    ELSE
      CALL ConvertToVisu_DG()
#if FV_ENABLED
      CALL ConvertToVisu_FV()
#endif
    END IF ! Avg2D
  END IF ! changedStateFile.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedNCalc.OR.changedAvg2D

  IF (doSurfVisu) THEN
    ! convert Surface DG solution to visu grid
    IF (changedStateFile.OR.changedWithDGOperator.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedNCalc.OR.changedBCnames) THEN
      CALL ConvertToSurfVisu_DG()
#if FV_ENABLED
      CALL ConvertToSurfVisu_FV()
#endif
    END IF
  END IF ! doSurfVisu

  ! convert generic data to visu grid
  IF (changedStateFile.OR.changedWithDGOperator.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedBCnames.OR.changedAvg2D) THEN
    CALL ConvertToVisu_GenericData(statefile)
  END IF

#if USE_MPI
   IF ((.NOT.MPIRoot).AND.(Avg2d)) THEN
     ! For parallel averaging, all data is gathered on the root. Disable output for other procs.
     nElemsAvg2D_DG = 0
     nElemsAvg2D_FV = 0
     SDEALLOCATE(UVisu_DG)
     SDEALLOCATE(UVisu_FV)
     ALLOCATE(UVisu_DG(0:NVisu   ,0:NVisu   ,0:0,nElemsAvg2D_DG,nVarVisu))
     ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV,nVarVisu))
   END IF
#endif

  ! Convert coordinates to visu grid
  IF (changedMeshFile.OR.changedNodeType.OR.changedNVisu.OR.changedFV_Elems.OR.changedDGonly.OR.changedAvg2D)   &
    CALL BuildVisuCoords()

  IF (doSurfVisu .AND. &
    ! Convert surface coordinates to visu grid
    (changedMeshFile.OR.changedNodeType.OR.changedNVisu.OR.changedFV_Elems.OR.changedDGonly.OR.changedBCnames)) &
      CALL BuildSurfVisuCoords()
END IF

MeshFile_old          = MeshFile
prmfile_old           = prmfile
statefile_old         = statefile
NVisu_old             = NVisu
NCalc_old             = NCalc
nVar_State_old        = nVar_State
withDGOperator_old    = withDGOperator
DGonly_old            = DGonly
Avg2D_old             = Avg2D
NodeTypeVisuPosti_old = NodeTypeVisuPosti
NodeType_State_old    = NodeType_State
NState_old            = PP_N
RestartMode           = -1

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,A)') " Visu finished for state file: ", TRIM(statefile)
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE visu


!===================================================================================================================================
!> Deallocate arrays used by visu.
!===================================================================================================================================
SUBROUTINE FinalizeVisu()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars         ,ONLY: StartTime
USE MOD_Commandline_Arguments,ONLY: FinalizeCommandlineArguments
USE MOD_DG                   ,ONLY: FinalizeDG
USE MOD_DG_Vars
USE MOD_Equation             ,ONLY: FinalizeEquation
USE MOD_Filter               ,ONLY: FinalizeFilter
USE MOD_Interpolation        ,ONLY: FinalizeInterpolation
USE MOD_IO_HDF5              ,ONLY: FinalizeIOHDF5
USE MOD_Mesh                 ,ONLY: FinalizeMesh
USE MOD_Mesh_Vars            ,ONLY: Elem_xGP
USE MOD_Mortar               ,ONLY: FinalizeMortar
USE MOD_Overintegration      ,ONLY: FinalizeOverintegration
USE MOD_ReadInTools          ,ONLY: FinalizeParameters
USE MOD_Restart              ,ONLY: FinalizeRestart
USE MOD_Visu_Vars
#if PARABOLIC
USE MOD_Lifting              ,ONLY: FinalizeLifting
#endif
#if FV_ENABLED
USE MOD_Indicator            ,ONLY: FinalizeIndicator
USE MOD_FV_Basis             ,ONLY: FinalizeFV_Basis
#endif /* FV_ENABLED */
#if USE_MPI
USE MOD_MPI                  ,ONLY: FinalizeMPI
#endif /* USE_MPI */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: Time,SimulationTime,mins,secs,hours,days
!===================================================================================================================================

IF(MPIRoot)THEN
  IF(FILEEXISTS('.posti.ini'))THEN
    OPEN(UNIT=31, FILE='.posti.ini', STATUS='old')
    CLOSE(31, STATUS='delete')
  END IF
  IF(FILEEXISTS('.flexi.ini'))THEN
    OPEN(UNIT=31, FILE='.flexi.ini', STATUS='old')
    CLOSE(31, STATUS='delete')
  END IF
END IF

! Reset all strings and variables
prmfile_old       = ''
statefile_old     = ''
MeshFile          = ''
MeshFile_old      = ''
NodeTypeVisuPosti = 'VISU'
NodeTypeVisuPosti_old = ''
NVisu     = -1
NVisu_old = -1
nVar_State_old = -1
NState_old = -1
withDGOperator_old = .FALSE.
hasFV_Elems = .FALSE.

SDEALLOCATE(mapDepToCalc)
#if FV_ENABLED && FV_RECONSTRUCT
SDEALLOCATE(mapDepToCalc_FV)
#endif
SDEALLOCATE(mapAllVarsToVisuVars)
SDEALLOCATE(mapAllVarsToSurfVisuVars)
SDEALLOCATE(mapAllVarsToSurfVisuVars_old)
SDEALLOCATE(UCalc_DG)
SDEALLOCATE(UCalc_FV)
SDEALLOCATE(FV_Elems_loc)

SDEALLOCATE(mapDGElemsToAllElems)
SDEALLOCATE(mapFVElemsToAllElems)

SDEALLOCATE(CoordsVisu_DG)
SDEALLOCATE(UVisu_DG)
SDEALLOCATE(CoordsVisu_FV)
SDEALLOCATE(UVisu_FV)
SDEALLOCATE(U)
SDEALLOCATE(Elem_xGP)

CALL FinalizeRestart()
CALL FinalizeEquation()
CALL FinalizeDG()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
CALL FinalizeCommandlineArguments()
CALL FinalizeParameters()
CALL FinalizeInterpolation()
CALL FinalizeMesh()
CALL FinalizeIOHDF5()
#if FV_ENABLED
CALL FinalizeIndicator()
CALL FinalizeFV_Basis()
#endif /* FV_ENABLED */
CALL FinalizeMortar()
#if PARABOLIC
CALL FinalizeLifting()
#endif /*PARABOLIC*/

SDEALLOCATE(VarNamesHDF5)
SDEALLOCATE(VarnamesAll)
SDEALLOCATE(BCNamesAll)
SDEALLOCATE(DepTable)
SDEALLOCATE(DepSurfaceOnly)
SDEALLOCATE(DepVolumeOnly)

SDEALLOCATE(FV_Elems_old)
SDEALLOCATE(mapDepToCalc_FV)
SDEALLOCATE(mapAllBCSidesToDGVisuBCSides)
SDEALLOCATE(mapAllBCSidesToFVVisuBCSides)
SDEALLOCATE(mapAllBCNamesToVisuBCNames_old)
SDEALLOCATE(mapAllBCNamesToVisuBCNames)
SDEALLOCATE(nSidesPerBCNameVisu_DG)
SDEALLOCATE(nSidesPerBCNameVisu_FV)

! Calculate simulation time
Time = FLEXITIME()
SimulationTime = Time-StartTime

! Get secs, mins, hours and days
secs = MOD(SimulationTime,60.)
SimulationTime = SimulationTime / 60.
mins = MOD(SimulationTime,60.)
SimulationTime = SimulationTime / 60.
hours = MOD(SimulationTime,24.)
SimulationTime = SimulationTime / 24.
!days = MOD(SimulationTime,365.) ! Use this if years are also to be displayed
days = SimulationTime

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F16.2,A)',ADVANCE='NO')  ' VISU  FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(A2,I6,A1,I0.2,A1,I0.2,A1,I0.2,A1)') ' [',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),']'
SWRITE(UNIT_stdOut,'(132("="))')

#if USE_MPI
CALL FinalizeMPI()
#endif /* USE_MPI */

END SUBROUTINE FinalizeVisu

END MODULE MOD_Visu
