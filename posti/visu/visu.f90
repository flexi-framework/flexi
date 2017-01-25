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
!> Module containing the main procedures for the visu tool: visu_requestInformation is called by ParaView to create a
!> list of available variables and visu is the main routine which is either called by ParaView to get the data it visualizes
!> or by the standalone tool.
!===================================================================================================================================
MODULE MOD_Visu
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE visu_getVarNamesAndFileType
  MODULE PROCEDURE visu_getVarNamesAndFileType
END INTERFACE

INTERFACE Visu_InitFile
  MODULE PROCEDURE Visu_InitFile
END INTERFACE

INTERFACE visu
  MODULE PROCEDURE visu
END INTERFACE

INTERFACE FinalizeVisu
  MODULE PROCEDURE FinalizeVisu
END INTERFACE

PUBLIC:: visu_getVarNamesAndFileType
PUBLIC:: visu_InitFile
PUBLIC:: visu
PUBLIC:: FinalizeVisu

CONTAINS

!===================================================================================================================================
!> Create a list of available variables for ParaView. This list contains the conservative, primitve and derived quantities
!> that are available in the current equation system as well as the additional variables read from the state file.
!> The additional variables are stored in the datasets 'ElemData' (elementwise data) and 'FieldData' (pointwise data).
!> Also a list of all available boundary names is created for surface visualization.
!===================================================================================================================================
SUBROUTINE visu_getVarNamesAndFileType(statefile,varnames_loc, bcnames_loc) 
USE MOD_Globals
USE MOD_Visu_Vars      ,ONLY: FileType,VarNamesHDF5,nBCNamesAll
USE MOD_HDF5_Input     ,ONLY: OpenDataFile,CloseDataFile,GetDataSize,GetVarNames,ISVALIDMESHFILE,ISVALIDHDF5FILE,ReadAttribute
USE MOD_HDF5_Input     ,ONLY: DatasetExists,HSize,nDims,ReadArray
USE MOD_IO_HDF5        ,ONLY: GetDatasetNamesInGroup,File_ID
USE MOD_StringTools    ,ONLY: STRICMP
USE MOD_EOS_Posti_Vars ,ONLY: DepNames,nVarDepEOS
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=255),INTENT(IN)                       :: statefile
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: varnames_loc(:)
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: bcnames_loc(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: i,j,nVar,dims
LOGICAL                                             :: varnames_found,readDGsolutionVars,sameVars,VarNamesExist
CHARACTER(LEN=255),ALLOCATABLE                      :: datasetNames(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: varnames_tmp(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: tmp(:)
CHARACTER(LEN=255)                                  :: MeshFile_loc
INTEGER                                             :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================

IF (ISVALIDMESHFILE(statefile)) THEN      ! MESH
  SDEALLOCATE(varnames_loc)
  FileType='Mesh'
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! other file
  SDEALLOCATE(varnames_loc)
  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadAttribute(File_ID,'File_Type',   1,StrScalar =FileType)

  ! check if variables in state file are the same as in the EQNSYS
  ! and set FileType to 'Generic' if not
  IF (STRICMP(FileType,'State')) THEN
    SDEALLOCATE(VarNamesHDF5)
    CALL GetVarNames("VarNames",VarNamesHDF5,VarNamesExist)
    IF (VarNamesExist) THEN
      sameVars=.TRUE.
      DO i=1,SIZE(VarNamesHDF5)
        sameVars = sameVars.AND.(STRICMP(VarNamesHDF5(i), DepNames(i)))
      END DO
    ELSE 
      sameVars=.FALSE.
    END IF
    IF (.NOT.sameVars) FileType='Generic'
  END IF

  IF (STRICMP(FileType,'State')) THEN
    nVar = nVarDepEOS
    ALLOCATE(varnames_loc(nVar))
    varnames_loc(1:nVar)=DepNames
    readDGsolutionVars = .FALSE.
  ELSE
    nVar=0
    readDGsolutionVars = .TRUE.
  END IF

  CALL GetDatasetNamesInGroup("/",datasetNames)

  DO i=1,SIZE(datasetNames)
    SDEALLOCATE(varnames_tmp)
    VarNamesExist=.FALSE.
    CALL DatasetExists(File_ID,"VarNames_"//TRIM(datasetNames(i)),varnames_found,attrib=.TRUE.)
    IF (varnames_found) THEN
          CALL GetVarNames("VarNames_"//TRIM(datasetNames(i)),varnames_tmp,VarNamesExist)
    ELSE
      IF (STRICMP(datasetNames(i), "DG_Solution")) THEN
        IF (readDGsolutionVars) THEN
          CALL GetVarNames("VarNames",varnames_tmp,VarNamesExist)
        END IF
      ELSE IF(STRICMP(datasetNames(i), "ElemData")) THEN
        CALL GetVarNames("VarNamesAdd",varnames_tmp,VarNamesExist)
      ELSE IF(STRICMP(datasetNames(i), "FieldData")) THEN
        CALL GetVarNames("VarNamesAddField",varnames_tmp,VarNamesExist)
      ELSE
        CALL GetDataSize(File_ID,TRIM(datasetNames(i)),dims,HSize)
        ALLOCATE(varnames_tmp(INT(HSize(1))))
        DO j=1,INT(HSize(1))
          WRITE(varnames_tmp(j),'(I0)') j
        END DO
        VarNamesExist=.TRUE.
      END IF
    END IF
    IF (.NOT.VarNamesExist) CYCLE
    
    ! increase array 'varnames_loc'  
    IF (nVar.GT.0) THEN
      ALLOCATE(tmp(nVar))
      tmp = varnames_loc
    END IF
    SDEALLOCATE(varnames_loc)
    ALLOCATE(varnames_loc(nVar+SIZE(varnames_tmp)))
    IF (nVar.GT.0) varnames_loc(1:nVar) = tmp(1:nVar)
    SDEALLOCATE(tmp)

    ! copy new varnames from varnames_tmp to varnames_loc 
    DO j=1,SIZE(varnames_tmp) 
      varnames_loc(nVar+j) = TRIM(datasetNames(i))//":"//TRIM(varnames_tmp(j))
    END DO
    nVar = nVar + SIZE(varnames_tmp)
  END DO

  ! Save mesh file to get boundary names later
  CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar =MeshFile_loc)

  CALL CloseDataFile()

  ! Open the mesh file and read all boundary names for surface visualization
  CALL OpenDataFile(MeshFile_loc,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
  CHECKSAFEINT(HSize(1),4)
  nBCNamesAll=INT(HSize(1),4)
  DEALLOCATE(HSize)
  SDEALLOCATE(bcnames_loc)
  ALLOCATE(bcnames_loc(nBCNamesAll))
  CALL ReadArray('BCNames',1,(/nBCNamesAll/),Offset,1,StrArray=bcnames_loc)
  CALL CloseDataFile()

  SDEALLOCATE(datasetNames)
END IF
END SUBROUTINE visu_getVarNamesAndFileType

!===================================================================================================================================
!> This routine is used to prepare everything we need to visualize data from a statefile. 
!> This includes:
!> * Get the mesh file
!> * Read the desired visualization polynomial degree, the visualization dimennsion, the node type we want to visualize on and the 
!>   Dg only option
!> * Decide whether the state file, the mesh file, the visualization polynomial degree or the dg only option changed. This is
!>   needed to decide what parts of the visualization routines should be called.
!> * Call routines that build the distribution between FV and DG elements and the mappings needed to calculate and visualize the
!>   desired variables.
!===================================================================================================================================
SUBROUTINE visu_InitFile(statefile,postifile)
! MODULES
USE HDF5
USE MOD_Preproc
USE MOD_Globals
USE MOD_Visu_Vars
USE MOD_EOS_Posti_Vars
USE MOD_MPI                ,ONLY: InitMPI
USE MOD_HDF5_Input         ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,GetArrayAndName
USE MOD_HDF5_Input         ,ONLY: ReadAttribute,File_ID,OpenDataFile,GetDataProps,CloseDataFile,ReadArray,DatasetExists
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_Output_Vars        ,ONLY: ProjectName
USE MOD_StringTools        ,ONLY: STRICMP
USE MOD_ReadInTools        ,ONLY: prms,GETINT,GETLOGICAL,addStrListEntry,GETSTR,CountOption,FinalizeParameters
USE MOD_Posti_Mappings     ,ONLY: Build_FV_DG_distribution,Build_mapDepToCalc_mapAllVarsToVisuVars
USE MOD_Mesh_Vars          ,ONLY: nElems,offsetElem

IMPLICIT NONE
CHARACTER(LEN=255),INTENT(IN)    :: statefile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
INTEGER                          :: nElems_State,iElem
CHARACTER(LEN=255)               :: NodeType_State 
LOGICAL                          :: exists
INTEGER                          :: ii,jj
!===================================================================================================================================
CALL visu_getVarNamesAndFileType(statefile,VarnamesAll,BCNamesAll)
IF (STRICMP(statefile,'Mesh')) THEN
    CALL CollectiveStop(__STAMP__, &
        "FileType==Mesh, but we try to initialize a state file!")
END IF

! open state file to be able to read attributes
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! read the meshfile attribute from statefile
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar =MeshFile_state)

! read options from posti parameter file
CALL prms%read_options(postifile)
NVisu             = GETINT("NVisu")
! again read MeshFile from posti prm file (this overwrites the MeshFile read from the state file)

!!!!!!
! WARNING: GETCWD is a GNU extension to the Fortran standard and will probably not work on other compilers
CALL GETCWD(MeshFile)
!!!!!!
Meshfile          =  TRIM(Meshfile) // "/" // GETSTR("MeshFile",MeshFile_state) 
Avg2D             = GETLOGICAL("Avg2D")
NodeTypeVisuPosti = GETSTR('NodeTypeVisu')
DGonly            = GETLOGICAL('DGonly')

! check if state, mesh, NVisu or DGonly changed 
changedStateFile = .NOT.STRICMP(statefile,statefile_old)
changedMeshFile  = .NOT.(STRICMP(MeshFile,MeshFile_old))
changedNVisu     = ((NVisu.NE.NVisu_old) .OR. (NodeTypeVisuPosti.NE.NodeTypeVisuPosti_old))
changedDGonly    = (DGonly.NEQV.DGonly_old)
changedAvg2D     = (Avg2D.NEQV.Avg2D_old)

SWRITE(*,*) "state file old -> new: ", TRIM(statefile_old), " -> ",TRIM(statefile)
SWRITE(*,*) " mesh file old -> new: ", TRIM(MeshFile_old) , " -> ",TRIM(MeshFile)

! if Mesh or State changed readin some more attributes/parameters
IF (changedStateFile.OR.changedMeshFile) THEN
  CALL GetDataProps(nVar_State,PP_N,nElems_State,NodeType_State)
  IF (.NOT.STRICMP(NodeType_State, NodeType)) THEN
    CALL CollectiveStop(__STAMP__, &
        "NodeType of state does not match with NodeType the visu-posti is compiled with!")
  END IF
  CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =ProjectName)
  CALL ReadAttribute(File_ID,'Time',        1,RealScalar=OutputTime)
END IF

CALL CloseDataFile()

! set number of dependent and raw variables 
SDEALLOCATE(DepTable)
SDEALLOCATE(DepSurfaceOnly)
nVarAll=SIZE(VarnamesAll)
IF (STRICMP(FileType,'State')) THEN
  nVarDep = nVarDepEOS
  ALLOCATE(DepTable(nVarDep,0:nVarDep))
  ALLOCATE(DepSurfaceOnly(nVarDep))
  DepTable = DepTableEOS
  DepSurfaceOnly = DepSurfaceOnlyEOS
ELSE 
  nVarDep = 0
  ALLOCATE(DepTable(nVarDep,0:nVarDep))
  ALLOCATE(DepSurfaceOnly(nVarDep))
  DepTable = 0
  DepSurfaceOnly = 0
END IF

! build distribution of FV and DG elements, which is stored in FV_Elems_loc
IF (changedStateFile.OR.changedMeshFile.OR.changedDGonly) THEN
  CALL Build_FV_DG_distribution(statefile) 
END IF 

! reset withDGOperator flag and check if it is needed due to existing FV elements
withDGOperator = .FALSE.
#if FV_RECONSTRUCT
! force calculation of metrics if there are any FV elements
IF (hasFV_Elems) withDGOperator = .TRUE.
#endif

! build mappings of variables which must be calculated/visualized
! also set withDGOperator flag if a dependent variable requires the evaluation of the DG operator
CALL Build_mapDepToCalc_mapAllVarsToVisuVars()

IF (Avg2D) THEN
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL DatasetExists(File_ID,'nElems_IJK',exists)
  IF (exists) THEN
    SDEALLOCATE(Elem_IJK)
    ALLOCATE(Elem_IJK(3,nElems))
    CALL ReadArray('nElems_IJK',1,(/3/),0,1,IntegerArray=nElems_IJK)
    CALL ReadArray('Elem_IJK',2,(/3,nElems/),offsetElem,2,IntegerArray=Elem_IJK)
  ELSE 
    CALL CollectiveStop(__STAMP__,&
        "No Elem_IJK sorting found in mesh file. Required for Avg2D!")
  END IF
  CALL CloseDataFile()

  SDEALLOCATE(FVAmountAvg2D)
  ALLOCATE(FVAmountAvg2D(nElems_IJK(1),nElems_IJK(2)))
  FVAmountAvg2D = 0.
  DO iElem=1,nElems
    ii = Elem_IJK(1,iElem)
    jj = Elem_IJK(2,iElem)
    FVAmountAvg2D(ii,jj) = FVAmountAvg2D(ii,jj) + FV_Elems_loc(iElem)
  END DO
  FVAmountAvg2D = FVAmountAvg2D / REAL(nElems_IJK(3))

  nElemsAvg2D_DG = 0
  nElemsAvg2D_FV = 0
  SDEALLOCATE(mapElemIJToDGElemAvg2D)
  SDEALLOCATE(mapElemIJToFVElemAvg2D)
  ALLOCATE(mapElemIJToDGElemAvg2D(nElems_IJK(1),nElems_IJK(2)))
  ALLOCATE(mapElemIJToFVElemAvg2D(nElems_IJK(1),nElems_IJK(2)))
  mapElemIJToDGElemAvg2D = 0
  mapElemIJToDGElemAvg2D = 0
  DO iElem=1,nElems
    IF (Elem_IJK(3,iElem).EQ.1) THEN
      ii = Elem_IJK(1,iElem)
      jj = Elem_IJK(2,iElem)
      IF (FVAmountAvg2D(ii,jj).LE.0.5) THEN
        nElemsAvg2D_DG = nElemsAvg2D_DG + 1
        mapElemIJToDGElemAvg2D(ii,jj) = nElemsAvg2D_DG
      ELSE
        nElemsAvg2D_FV = nElemsAvg2D_FV + 1
        mapElemIJToFVElemAvg2D(ii,jj) = nElemsAvg2D_FV
      END IF
    END IF
  END DO
END IF

changedWithDGOperator = (withDGOperator.NEQV.withDGOperator_old)
END SUBROUTINE visu_InitFile

!===================================================================================================================================
!> Main routine of the visualization tool visu. Called either by the ParaView plugin or by the standalone program version.
!===================================================================================================================================
SUBROUTINE visu(mpi_comm_IN, prmfile, postifile, statefile)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_MPI                 ,ONLY: InitMPI
USE MOD_HDF5_Input          ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,OpenDataFile,CloseDataFile
USE MOD_Posti_ReadState     ,ONLY: ReadState
USE MOD_Posti_VisuMesh      ,ONLY: VisualizeMesh
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_DG,CalcSurfQuantities_DG
#if FV_ENABLED
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_FV,CalcSurfQuantities_FV
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV,ConvertToSurfVisu_FV
USE MOD_FV_Vars             ,ONLY: FV_Vdm,FV_sVdm
#endif
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG,ConvertToSurfVisu_DG,ConvertToVisu_GenericData
USE MOD_ReadInTools         ,ONLY: prms,FinalizeParameters,ExtractParameterFile
USE MOD_StringTools         ,ONLY: STRICMP
USE MOD_Posti_VisuMesh      ,ONLY: BuildVisuCoords,BuildSurfVisuCoords
USE MOD_Posti_Mappings      ,ONLY: Build_mapBCSides
USE MOD_Visu_Avg2D          ,ONLY: Average2D
USE MOD_Interpolation_Vars  ,ONLY: NodeType,NodeTypeVISUFVEqui
USE MOD_Interpolation       ,ONLY: GetVandermonde
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: mpi_comm_IN    
CHARACTER(LEN=255),INTENT(INOUT) :: prmfile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
CHARACTER(LEN=255),INTENT(IN)    :: statefile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: changedPrmFile
REAL,ALLOCATABLE :: Vdm_DGToFV  (:,:)
REAL,ALLOCATABLE :: Vdm_FVToDG  (:,:)
REAL,ALLOCATABLE :: Vdm_DGToVisu(:,:)
REAL,ALLOCATABLE :: Vdm_FVToVisu(:,:)
REAL,ALLOCATABLE :: FVdouble(:,:)
INTEGER           :: i
!===================================================================================================================================
CALL InitMPI(mpi_comm_IN) 
SWRITE (*,*) "READING FROM: ", TRIM(statefile)

!**********************************************************************************************
! General workflow / principles of the visu ParaView-plugin
!
! * all arrays are SDEALLOCATEd just before they are allocated. This is done to keep there
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
!                    primitive quantities U_Prim and gradUx/y/z as well as gradUxi/eta/zet are
!                    filled. These are used to calculate the visu-quantities.
!                    
! * The whole calculation of derived quantities is performed on PP_N and afterwards 
!   interpolated to NVisu. This is not the case for FV elements with FV_RECONSTRUCT enabled.
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
!   
! WORKFLOW:
! * The main steps are:
!   1. call InitFile for FV/DG distribution and mappings
!   2. read solution          (if changedStateFile or changedWithDGOperator or changedDGonly)
!   3. build mapping for BC sides that should be visualized (done after read solution since some
!      mesh infos are needed)
!   4. read Mesh              (if changedMeshFile)
!   5. compute UCalc          (if changedStateFile or changedVarNames or changedDGonly) 
!   6. convert to UVisu       (if changedStateFile or changedVarNames or changedNVisu or changedDGonly)
!   7. build visu mesh        (if changedMeshFile  or changedNVisu or changedFV_Elems or changedDGonly)
!   5. - 7. are done seperately for surface variables if surface visualization is turned on
!   8. write VTK arrays       (always!) 
!
!**********************************************************************************************

CALL FinalizeParameters()
! Read Varnames to visualize and build calc and visu dependencies
CALL prms%SetSection("posti")
CALL prms%CreateStringOption( "MeshFile"     , "Custom mesh file ")
CALL prms%CreateStringOption( "VarName"      , "Names of variables, which should be visualized.", multiple=.TRUE.)
CALL prms%CreateIntOption(    "NVisu"        , "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateLogicalOption("Avg2D"        , "Average solution in z-direction",".FALSE.")
CALL prms%CreateStringOption( "NodeTypeVisu" , "NodeType for visualization. Visu, Gauss,Gauss-Lobatto,Visu_inner"    ,"VISU")
CALL prms%CreateLogicalOption("DGonly"       , "Visualize FV elements as DG elements."    ,".FALSE.")
CALL prms%CreateStringOption( "BoundaryName" , "Names of boundaries for surfaces, which should be visualized.", multiple=.TRUE.)

changedStateFile      = .FALSE.
changedMeshFile       = .FALSE.
changedNVisu          = .FALSE.
changedVarNames       = .FALSE.
changedFV_Elems       = .FALSE.
changedWithDGOperator = .FALSE.
changedDGonly         = .FALSE.

IF (ISVALIDMESHFILE(statefile)) THEN ! visualize mesh
  SWRITE(*,*) "MeshFile Mode"
  MeshFileMode = .TRUE.
  MeshFile      = statefile
  nVar_State    = 0
  withDGOperator = .FALSE.
  CALL VisualizeMesh(postifile,MeshFile)
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! visualize state file
  SWRITE(*,*) "State Mode"
  MeshFileMode = .FALSE.
  ! initialize state file
  CALL visu_InitFile(statefile,postifile)

  ! read solution from state file (either direct or including a evaluation of the DG operator)
  IF (LEN_TRIM(prmfile).EQ.0) THEN
    changedPrmFile = .NOT.STRICMP(prmfile_old, ".flexi.ini")
  ELSE
    changedPrmFile = (prmfile .NE. prmfile_old)  
  END IF
  SWRITE (*,*) "changedStateFile     ", changedStateFile     
  SWRITE (*,*) "changedMeshFile      ", changedMeshFile      
  SWRITE (*,*) "changedNVisu         ", changedNVisu         
  SWRITE (*,*) "changedVarNames      ", changedVarNames      
  SWRITE (*,*) "changedFV_Elems      ", changedFV_Elems      
  SWRITE (*,*) "changedWithDGOperator", changedWithDGOperator
  SWRITE (*,*) "changedDGonly        ", changedDGonly
  SWRITE (*,*) "changedAvg2D         ", changedAvg2D
  SWRITE (*,*) "changedPrmFile       ", changedPrmFile, TRIM(prmfile_old), " -> ", TRIM(prmfile)
  SWRITE (*,*) "changedBCnames       ", changedBCnames
  IF (changedStateFile.OR.changedWithDGOperator.OR.changedPrmFile.OR.changedDGonly) THEN
      CALL ReadState(prmfile,statefile)
  END IF

  ! build mappings of BC sides for surface visualization
  CALL Build_mapBCSides()

  ! ===== calc solution =====
  IF (changedStateFile.OR.changedVarNames.OR.changedDGonly) THEN
    CALL CalcQuantities_DG()
#if FV_ENABLED
    CALL CalcQuantities_FV()
#endif
  END IF
  IF (doSurfVisu) THEN
    ! calc surface solution 
    IF (changedStateFile.OR.changedVarNames.OR.changedDGonly.OR.changedBCnames) THEN
      CALL CalcSurfQuantities_DG()
#if FV_ENABLED
      CALL CalcSurfQuantities_FV()
#endif
    END IF
  END IF

  ! ===== convert solution to visu grid =====
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedAvg2D) THEN
    ! ===== Avg2d =====
    IF (Avg2d) THEN
      SDEALLOCATE(UVisu_DG)
      SDEALLOCATE(UVisu_FV)
      ALLOCATE(UVisu_DG(0:NVisu   ,0:NVisu   ,0:0,nElemsAvg2D_DG,nVarVisu))
      ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV,nVarVisu))
#if FV_RECONSTRUCT      
      ALLOCATE(Vdm_DGToFV  (0:NVisu_FV,0:PP_N    ))
      ALLOCATE(Vdm_FVToDG  (0:PP_N    ,0:NVisu_FV))
      ALLOCATE(Vdm_DGToVisu(0:NVisu   ,0:PP_N    ))
      ALLOCATE(Vdm_FVToVisu(0:NVisu_FV,0:NVisu_FV))
      ALLOCATE(FVdouble(0:NVisu_FV,0:PP_N))
      FVdouble = 0.
      DO i = 0, PP_N
        FVdouble(i*2  ,i) = 1. 
        FVdouble(i*2+1,i) = 1.
      END DO ! i = 0, PP_N
      Vdm_DGToFV = MATMUL(FVdouble,FV_Vdm)
      Vdm_FVToDG = MATMUL(FV_sVdm,TRANSPOSE(FVdouble))
      Vdm_FVToVisu = 0.
      DO i = 0, NVisu_FV
        Vdm_FVToVisu(i,i) = 1.
      END DO
      CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVisuPosti,Vdm_DGToVisu,modal=.FALSE.)

      CALL Average2D(nVarCalc,nVarCalc_FV,PP_N,NVisu_FV,nElems_DG,nElems_FV,NodeType,UCalc_DG,UCalc_FV,&
          Vdm_DGToFV,Vdm_FVToDG,Vdm_DGToVisu,Vdm_FVToVisu, &
          1,nVarDep,mapDepToCalc,&
          UVisu_DG,UVisu_FV)
#else
      todo
#endif
      
    ELSE 
      CALL ConvertToVisu_DG()
#if FV_ENABLED
      CALL ConvertToVisu_FV()
#endif
    END IF
  END IF
  IF (doSurfVisu) THEN
    ! convert Surface DG solution to visu grid
    IF (changedStateFile.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedBCnames) THEN
      CALL ConvertToSurfVisu_DG()
#if FV_ENABLED
      CALL ConvertToSurfVisu_FV()
#endif
    END IF
  END IF

  ! convert generic data to visu grid
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedBCnames.OR.changedAvg2D) THEN
    CALL ConvertToVisu_GenericData(statefile)
  END IF


  ! Convert coordinates to visu grid
  IF (changedMeshFile.OR.changedNVisu.OR.changedFV_Elems.OR.changedDGonly.OR.changedAvg2D) THEN
    CALL BuildVisuCoords()
  END IF
  IF (doSurfVisu) THEN
    ! Convert surface coordinates to visu grid
    IF (changedMeshFile.OR.changedNVisu.OR.changedFV_Elems.OR.changedDGonly.OR.changedBCnames) THEN
      CALL BuildSurfVisuCoords()
    END IF
  END IF


END IF

MeshFile_old          = MeshFile
prmfile_old           = prmfile
statefile_old         = statefile
NVisu_old             = NVisu
nVar_State_old        = nVar_State
withDGOperator_old    = withDGOperator
DGonly_old            = DGonly
Avg2D_old             = Avg2D
NodeTypeVisuPosti_old = NodeTypeVisuPosti

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(*,*) "Visu finished for state file: ", TRIM(statefile)
SWRITE(UNIT_StdOut,'(132("="))')
END SUBROUTINE visu

!===================================================================================================================================
!> Deallocate arrays used by visu.
!===================================================================================================================================
SUBROUTINE FinalizeVisu()
USE MOD_Globals,ONLY: MPIRoot
USE MOD_Visu_Vars
IMPLICIT NONE
!===================================================================================================================================
SWRITE (*,*) "VISU FINALIZE"
prmfile_old = ""
statefile_old = ""
MeshFile = ""
MeshFile_old = ""
NodeTypeVisuPosti = "VISU"
NodeTypeVisuPosti_old = ""
NVisu     = -1
NVisu_old = -1
nVar_State_old = -1
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

END SUBROUTINE FinalizeVisu

END MODULE MOD_Visu
