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
!> Module containing the main procedures for the POSTI tool: visu3d_requestInformation is called by ParaView to create a
!> list of available variables and visu3D is the main routine of POSTI.
!===================================================================================================================================
MODULE MOD_Visu3D
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE visu3d_requestInformation
  MODULE PROCEDURE visu3d_requestInformation
END INTERFACE

INTERFACE Visu3D_InitFile
  MODULE PROCEDURE Visu3D_InitFile
END INTERFACE

INTERFACE visu3D
  MODULE PROCEDURE visu3D
END INTERFACE

INTERFACE visu3D_CWrapper
  MODULE PROCEDURE visu3D_CWrapper
END INTERFACE

INTERFACE visu3d_dealloc_nodeids
  MODULE PROCEDURE visu3d_dealloc_nodeids
END INTERFACE

INTERFACE FinalizeVisu3D
  MODULE PROCEDURE FinalizeVisu3D
END INTERFACE

PUBLIC:: visu3d_requestInformation
PUBLIC:: visu3D_InitFile
PUBLIC:: visu3D
PUBLIC:: visu3D_CWrapper
PUBLIC:: visu3d_dealloc_nodeids
PUBLIC:: FinalizeVisu3D

CONTAINS


FUNCTION cstrToChar255(cstr, strlen) 
USE ISO_C_BINDING
! INPUT / OUTPUT VARIABLES 
TYPE(C_PTR),TARGET,INTENT(IN)  :: cstr
INTEGER,INTENT(IN)             :: strlen
CHARACTER(LEN=255)             :: cstrToChar255
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(KIND=C_CHAR),POINTER :: tmp(:)
!===================================================================================================================================
CALL C_F_POINTER(C_LOC(cstr), tmp, [strlen])
cstrToChar255 = TRANSFER(tmp(1:strlen), cstrToChar255)
cstrToChar255(strlen+1:255) = ' ' 
END FUNCTION cstrToChar255

!===================================================================================================================================
!> Wrapper to visu3D_InitFile for Paraview plugin
!===================================================================================================================================
SUBROUTINE visu3d_requestInformation(mpi_comm_IN, strlen_state, statefile_IN, varnames)
! MODULES
USE MOD_Globals
USE MOD_MPI     ,ONLY: InitMPI
USE MOD_Posti_Vars,ONLY: VarNamesTotal
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                    :: mpi_comm_IN    
INTEGER,INTENT(IN)                    :: strlen_state
TYPE(C_PTR),TARGET,INTENT(IN)         :: statefile_IN
TYPE (CARRAY), INTENT(INOUT)          :: varnames
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
CHARACTER(LEN=255)                    :: statefile
CHARACTER(LEN=255),POINTER            :: varnames_pointer(:)
!===================================================================================================================================
statefile = cstrToChar255(statefile_IN, strlen_state)

CALL InitMPI(mpi_comm_IN) 
CALL visu3d_getVarNamesAndFileType(statefile,VarNamesTotal)
IF (ALLOCATED(VarNamesTotal)) THEN
  varnames_pointer => VarNamesTotal
  varnames%len  = SIZE(varnames_pointer)*255 
  varnames%data = C_LOC(varnames_pointer(1))
ELSE
  varnames%len  = 0
  varnames%data = C_NULL_PTR
END IF
END SUBROUTINE visu3d_requestInformation

!===================================================================================================================================
! TODO:
!> Create a list of available variables for ParaView. This list contains the conservative, primitve and derived quantities
!> that are available in the current equation system as well as the additional variables read from the state file.
!> The additional variables are stored in the datasets 'ElemData' (elementwise data) and 'FieldData' (pointwise data).
!===================================================================================================================================
SUBROUTINE visu3d_getVarNamesAndFileType(statefile,varnames_loc) 
USE MOD_Globals
USE MOD_Posti_Vars     ,ONLY: FileType,VarNamesHDF5
USE MOD_HDF5_Input     ,ONLY: OpenDataFile,CloseDataFile,GetDataSize,GetVarNames,ISVALIDMESHFILE,ISVALIDHDF5FILE,ReadAttribute
USE MOD_HDF5_Input     ,ONLY: DatasetExists
USE MOD_IO_HDF5        ,ONLY: GetDatasetNamesInGroup,File_ID,HSize
USE MOD_StringTools    ,ONLY: STRICMP
USE MOD_EOS_Posti_Vars ,ONLY: DepNames,nVarTotalEOS
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=255),INTENT(IN)                       :: statefile
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: varnames_loc(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: i,j,nVar,dims
LOGICAL                                             :: varnames_found,readDGsolutionVars,sameVars
CHARACTER(LEN=255),ALLOCATABLE                      :: datasetNames(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: varnames_tmp(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: tmp(:)
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
    CALL GetVarNames("VarNames",VarNamesHDF5)
    sameVars=.TRUE.
    DO i=1,SIZE(VarNamesHDF5)
      sameVars = sameVars.AND.(STRICMP(VarNamesHDF5(i), DepNames(i)))
    END DO
    IF (.NOT.sameVars) FileType='Generic'
  END IF

  IF (STRICMP(FileType,'State')) THEN
    nVar = nVarTotalEOS
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
    CALL DatasetExists(File_ID,"VarNames_"//TRIM(datasetNames(i)),varnames_found,attrib=.TRUE.)
    IF (varnames_found) THEN
          CALL GetVarNames("VarNames_"//TRIM(datasetNames(i)),varnames_tmp)
    ELSE
      IF (STRICMP(datasetNames(i), "DG_Solution")) THEN
        IF (readDGsolutionVars) THEN
          CALL GetVarNames("VarNames",varnames_tmp)
        END IF
      ELSE IF(STRICMP(datasetNames(i), "ElemData")) THEN
        CALL GetVarNames("VarNamesAdd",varnames_tmp)
      ELSE IF(STRICMP(datasetNames(i), "FieldData")) THEN
        CALL GetVarNames("VarNamesAddField",varnames_tmp)
      ELSE
        CALL GetDataSize(File_ID,TRIM(datasetNames(i)),dims,HSize)
        ALLOCATE(varnames_tmp(INT(HSize(1))))
        DO j=1,INT(HSize(1))
          WRITE(varnames_tmp(j),'(I0)') j
        END DO
      END IF
    END IF
    IF (.NOT.ALLOCATED(varnames_tmp)) CYCLE
    
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
  CALL CloseDataFile()
  SDEALLOCATE(datasetNames)
END IF
END SUBROUTINE visu3d_getVarNamesAndFileType

!===================================================================================================================================
!===================================================================================================================================
SUBROUTINE visu3d_InitFile(statefile,postifile)
! MODULES
USE HDF5
USE MOD_Preproc
USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_EOS_Posti_Vars
USE MOD_MPI                ,ONLY: InitMPI
USE MOD_HDF5_Input         ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,GetArrayAndName
USE MOD_HDF5_Input         ,ONLY: ReadAttribute,File_ID,OpenDataFile,GetDataProps,CloseDataFile
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_Output_Vars        ,ONLY: ProjectName
USE MOD_StringTools        ,ONLY: STRICMP
USE MOD_ReadInTools        ,ONLY: prms,GETINT,GETLOGICAL,addStrListEntry,GETSTR,CountOption,FinalizeParameters
USE MOD_Posti_Mappings     ,ONLY: Build_FV_DG_distribution,Build_mapCalc_mapVisu

IMPLICIT NONE
CHARACTER(LEN=255),INTENT(IN)    :: statefile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
INTEGER                          :: nElems_State,iVar 
CHARACTER(LEN=255)               :: NodeType_State 
INTEGER                          :: postiUnit 
INTEGER                          :: stat
!===================================================================================================================================
CALL visu3d_getVarNamesAndFileType(statefile,VarNamesTotal)
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
Meshfile          = GETSTR("MeshFile",MeshFile_state) 
VisuDimension     = GETINT("VisuDimension")
NodeTypeVisuPosti = GETSTR('NodeTypeVisu')

! check if state, mesh or NVisu changed 
changedStateFile = .NOT.STRICMP(statefile,statefile_old)
changedMeshFile  = .NOT.(STRICMP(MeshFile,MeshFile_old))
changedNVisu     = ((NVisu.NE.NVisu_old) .OR. (NodeTypeVisuPosti.NE.NodeTypeVisuPosti_old))

SWRITE(*,*) "state file old -> new: ", TRIM(statefile_old), " -> ",TRIM(statefile)
SWRITE(*,*) " mesh file old -> new: ", TRIM(MeshFile_old) , " -> ",TRIM(MeshFile)

! if Mesh or State changed readin some more attributes/parameters
IF (changedStateFile.OR.changedMeshFile) THEN
  CALL GetDataProps(nVar_State,PP_N,nElems_State,NodeType_State)
  IF (.NOT.STRICMP(NodeType_State, NodeType)) THEN
    CALL CollectiveStop(__STAMP__, &
        "NodeType of state does not match with NodeType the visu3D-posti is compiled with!")
  END IF
  CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =ProjectName)
  CALL ReadAttribute(File_ID,'Time',        1,RealScalar=OutputTime)

  ! if no explicit posti parameter file is given, then generate a default '.posti.ini' file
  IF (LEN_TRIM(postifile).EQ.0) THEN
    postifile = ".posti.ini"
    IF (MPIRoot) THEN
      OPEN(NEWUNIT=postiUnit,FILE=TRIM(postifile),STATUS='UNKNOWN',ACTION='WRITE',ACCESS='SEQUENTIAL',IOSTAT=stat)
      WRITE(postiUnit,'(A,I3)') "NVisu   = ", PP_N+1
      DO iVar=1,nVar_State
        WRITE(postiUnit,'(A,A)') "VarName = ", TRIM(VarNamesHDF5(iVar))
      END DO
      CLOSE(postiUnit)
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
  END IF
END IF

CALL CloseDataFile()

! set number of dependent and raw variables 
SDEALLOCATE(DepTable)
nVarTotal=SIZE(VarNamesTotal)
IF (STRICMP(FileType,'State')) THEN
  nVarDep = nVarTotalEOS
  nVarRaw = nVarTotal - nVarDep
  calcDeps= .TRUE.
  ALLOCATE(DepTable(nVarDep,nVarDep))
  DepTable = DepTableEOS
ELSE 
  nVarDep = 0
  nVarRaw = nVarTotal
  calcDeps= .FALSE.
  ALLOCATE(DepTable(nVarDep,nVarDep))
  DepTable = 0
END IF

! build distribution of FV and DG elements, which is stored in FV_Elems_loc
IF (changedStateFile.OR.changedMeshFile) THEN
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
CALL Build_mapCalc_mapVisu()

changedWithDGOperator = (withDGOperator.NEQV.withDGOperator_old)
END SUBROUTINE visu3d_InitFile


!===================================================================================================================================
!> C wrapper routine for the visu3D call from ParaView.
!===================================================================================================================================
SUBROUTINE visu3D_CWrapper(mpi_comm_IN, strlen_prm, prmfile_IN, strlen_posti, postifile_IN, strlen_state, statefile_IN,&
        coordsDG_out,valuesDG_out,nodeidsDG_out, &
        coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)
USE ISO_C_BINDING
USE MOD_Globals
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: mpi_comm_IN    
INTEGER,INTENT(IN)            :: strlen_prm    
INTEGER,INTENT(IN)            :: strlen_posti    
INTEGER,INTENT(IN)            :: strlen_state    
TYPE(C_PTR),TARGET,INTENT(IN) :: prmfile_IN
TYPE(C_PTR),TARGET,INTENT(IN) :: postifile_IN
TYPE(C_PTR),TARGET,INTENT(IN) :: statefile_IN
TYPE (CARRAY), INTENT(INOUT)  :: coordsDG_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesDG_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsDG_out
TYPE (CARRAY), INTENT(INOUT)  :: coordsFV_out
TYPE (CARRAY), INTENT(INOUT)  :: valuesFV_out
TYPE (CARRAY), INTENT(INOUT)  :: nodeidsFV_out
TYPE (CARRAY), INTENT(INOUT)  :: varnames_out
TYPE (CARRAY), INTENT(INOUT)  :: components_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: prmfile
CHARACTER(LEN=255)            :: postifile
CHARACTER(LEN=255)            :: statefile
!===================================================================================================================================
prmfile   = cstrToChar255(prmfile_IN,   strlen_prm)
postifile = cstrToChar255(postifile_IN, strlen_posti)
statefile = cstrToChar255(statefile_IN, strlen_state)
CALL visu3D(mpi_comm_IN, prmfile, postifile, statefile, &
        coordsDG_out,valuesDG_out,nodeidsDG_out, &
        coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)
END SUBROUTINE visu3D_CWrapper


!===================================================================================================================================
!> Main routine of the visualization tool POSTI. Called either by the ParaView plugin or by the standalone program version.
!===================================================================================================================================
SUBROUTINE visu3D(mpi_comm_IN, prmfile, postifile, statefile, &
                  coordsDG_out,valuesDG_out,nodeidsDG_out,    &
                  coordsFV_out,valuesFV_out,nodeidsFV_out,    &
                  varnames_out,components_out)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_MPI                 ,ONLY: InitMPI
USE MOD_HDF5_Input          ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,OpenDataFile,CloseDataFile
USE MOD_Posti_ReadState     ,ONLY: ReadState
USE MOD_Posti_VisuMesh      ,ONLY: VisualizeMesh
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_DG
#if FV_ENABLED
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_ConvertToVisu_FV
#endif
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG
USE MOD_ReadInTools         ,ONLY: prms,FinalizeParameters,ExtractParameterFile
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: mpi_comm_IN    
CHARACTER(LEN=255),INTENT(INOUT) :: prmfile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
CHARACTER(LEN=255),INTENT(IN)    :: statefile
TYPE (CARRAY), INTENT(INOUT)     :: coordsDG_out
TYPE (CARRAY), INTENT(INOUT)     :: valuesDG_out
TYPE (CARRAY), INTENT(INOUT)     :: nodeidsDG_out
TYPE (CARRAY), INTENT(INOUT)     :: coordsFV_out
TYPE (CARRAY), INTENT(INOUT)     :: valuesFV_out
TYPE (CARRAY), INTENT(INOUT)     :: nodeidsFV_out
TYPE (CARRAY), INTENT(INOUT)     :: varnames_out
TYPE (CARRAY), INTENT(INOUT)     :: components_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: changedPrmFile
!===================================================================================================================================
CALL InitMPI(mpi_comm_IN) 
SWRITE (*,*) "READING FROM: ", TRIM(statefile)

!**********************************************************************************************
! General workflow / principles of the visu3D ParaView-plugin
!
! * all arrays are SDEALLOCATEd just before they are allocated. This is done to keep there
!   content during successive calls of the visu3D during a ParaView session. They are only
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
!   - 'mapVisu' is a integer array of size (1:nVarTotal), where nVarTotal is the total amount 
!     of available quantities. This map contains a zero for all not-to-visu-quantities and for
!     all quantities the index where it is stored in 'UVisu'.
!     This mapping is filled from the 'VarName' entries in the parameter file.
!   - 'mapCalc' is the same as mapVisu, but for all (intermediate) quantities stored in 'UCalc'.
!     This mapping is filled from the DepTable.
!
! CHANGED system:
! * There are different logical changedXXX variables, which indicated if XXX changed during 
!   successive calls of the visu3D. These variables control the general workflow of the visu3D.
!   - changedStateFile:     new state file 
!   - changedMeshFile:      new mesh file (only possible if changedStateFile==TRUE)
!   - changedVarNames:      new set of variables to visualize
!   - changedNVisu:         new NVisu, new Nodetype
!   - changedFV_Elems:      new distribution of FV/DG elements (only if changedStateFile==TRUE)
!   - changedWithDGOperator: different mode, with/without gradients
!   
! WORKFLOW:
! * The main steps are:
!   1. get nElems             (if changedStateFile)
!   2. get FV/DG distribution (if changedStateFile)
!   3. read solution          (if changedStateFile or changedWithDGOperator)
!   4. read Mesh              (if changedMeshFile)
!   5. compute UCalc          (if changedStateFile or changedVarNames) 
!   6. convert to UVisu       (if changedStateFile or changedVarNames or changedNVisu)
!   6. build visu mesh        (if changedMeshFile  or changedNVisu or changedFV_Elems)
!   7. write VTK arrays       (always!) 
!
!**********************************************************************************************

CALL FinalizeParameters()
! Read Varnames to visualize and build calc and visu dependencies
CALL prms%SetSection("posti")
CALL prms%CreateStringOption("MeshFile"     , "Custom mesh file ")
CALL prms%CreateStringOption("VarName"      , "Names of variables, which should be visualized.", multiple=.TRUE.)
CALL prms%CreateIntOption(   "NVisu"        ,  "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateIntOption(   "VisuDimension", "2 = Slice at first Gauss point in zeta-direction to get 2D solution.","3")
CALL prms%CreateStringOption("NodeTypeVisu" , "NodeType for visualization. Visu, Gauss,Gauss-Lobatto,Visu_inner"    ,"VISU")

changedStateFile     = .FALSE.
changedMeshFile      = .FALSE.
changedNVisu         = .FALSE.
changedVarNames      = .FALSE.
changedFV_Elems      = .FALSE.
changedWithDGOperator = .FALSE.

IF (ISVALIDMESHFILE(statefile)) THEN ! visualize mesh
  SWRITE(*,*) "MeshFile Mode"
  MeshFile      = statefile
  nVar_State    = 0
  withDGOperator = .FALSE.
  CALL VisualizeMesh(postifile,MeshFile,coordsDG_out,valuesDG_out,nodeidsDG_out, &
      coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! visualize state file
  SWRITE(*,*) "State Mode"
  ! initialize state file
  CALL visu3d_InitFile(statefile,postifile)

  ! read solution from state file (either direct or including a evaluation of the DG operator)
  changedPrmFile = (prmfile .NE. prmfile_old)  
  IF (changedStateFile.OR.changedWithDGOperator.OR.changedPrmFile) THEN
      CALL ReadState(prmfile,statefile)
  END IF

  ! calc DG solution 
  IF (changedStateFile.OR.changedVarNames) THEN
    CALL CalcQuantities_DG()
  END IF
  ! convert DG solution to visu grid
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu) THEN
    CALL ConvertToVisu_DG()
  END IF

#if FV_ENABLED
  ! calc FV solution and convert to visu grid
  IF ((changedStateFile.OR.changedVarNames).AND.hasFV_Elems) THEN
    CALL CalcQuantities_ConvertToVisu_FV()
  END IF
#endif /* FV_ENABLED */

  ! convert ElemData and FieldData to visu grid
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu) THEN
  END IF

  CALL Visu3D_Build_VTK(coordsDG_out,valuesDG_out,nodeidsDG_out, &
                        coordsFV_out,valuesFV_out,nodeidsFV_out, &
                        varnames_out,components_out)
END IF

MeshFile_old          = MeshFile
prmfile_old           = prmfile
statefile_old         = statefile
NVisu_old             = NVisu
nVar_State_old        = nVar_State
withDGOperator_old    = withDGOperator
NodeTypeVisuPosti_old = NodeTypeVisuPosti

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(*,*) "Visu3D finished for state file: ", TRIM(statefile)
SWRITE(UNIT_StdOut,'(132("="))')
END SUBROUTINE visu3D

!===================================================================================================================================
!> Build VTK arrays from solution
!===================================================================================================================================
SUBROUTINE Visu3D_Build_VTK(coordsDG_out,valuesDG_out,nodeidsDG_out, &
                            coordsFV_out,valuesFV_out,nodeidsFV_out, &
                            varnames_out,components_out)

USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_VTK           ,ONLY: WriteCoordsToVTK_array,WriteDataToVTK_array,WriteVarnamesToVTK_array
USE MOD_Posti_VisuMesh,ONLY: BuildVisuCoords
USE MOD_Output_Vars   ,ONLY: ProjectName
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
TYPE (CARRAY), INTENT(INOUT)     :: coordsDG_out
TYPE (CARRAY), INTENT(INOUT)     :: valuesDG_out
TYPE (CARRAY), INTENT(INOUT)     :: nodeidsDG_out
TYPE (CARRAY), INTENT(INOUT)     :: coordsFV_out
TYPE (CARRAY), INTENT(INOUT)     :: valuesFV_out
TYPE (CARRAY), INTENT(INOUT)     :: nodeidsFV_out
TYPE (CARRAY), INTENT(INOUT)     :: varnames_out
TYPE (CARRAY), INTENT(INOUT)     :: components_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iounit,i,iElem
CHARACTER(LEN=255) :: strOutputFile
!===================================================================================================================================
! write UVisu to VTK 2D / 3D arrays (must be done always!)
IF (VisuDimension.EQ.3) THEN

  CALL WriteDataToVTK_array(nVarVisuTotal,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG,3)
  CALL WriteDataToVTK_array(nVarVisuTotal,NVisu_FV,nElems_FV,valuesFV_out,UVisu_FV,3)

ELSE IF (VisuDimension.EQ.2) THEN

  ! allocate Visu 2D array and copy from first zeta-slice of 3D array
  SDEALLOCATE(UVisu_DG_2D)
  ALLOCATE(UVisu_DG_2D(0:NVisu,0:NVisu,0:0,1:nElems_DG,1:(nVarVisuTotal)))
  UVisu_DG_2D = UVisu_DG(:,:,0:0,:,:)
  SDEALLOCATE(UVisu_FV_2D)
  ALLOCATE(UVisu_FV_2D(0:NVisu_FV,0:NVisu_FV,0:0,1:nElems_FV,1:(nVarVisuTotal)))
#if FV_ENABLED
  UVisu_FV_2D = UVisu_FV(:,:,0:0,:,:)
#else
  CoordsVisu_FV_2D = 0
#endif

  CALL WriteDataToVTK_array(nVarVisuTotal,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG_2D,2)
  CALL WriteDataToVTK_array(nVarVisuTotal,NVisu_FV,nElems_FV,valuesFV_out,UVisu_FV_2D,2)

END IF

! Convert coordinates to visu grid
IF (changedMeshFile.OR.changedNVisu.OR.changedFV_Elems) THEN
  CALL BuildVisuCoords()
END IF

! write coords, UVisu to VTK  2D / 3D arrays (must be done always!)
IF (VisuDimension.EQ.3) THEN

  CALL WriteCoordsToVTK_array(NVisu   ,nElems_DG,coordsDG_out,nodeidsDG_out,&
      CoordsVisu_DG,nodeids_DG,dim=3,DGFV=0)
  CALL WriteCoordsToVTK_array(NVisu_FV,nElems_FV,coordsFV_out,nodeidsFV_out,&
      CoordsVisu_FV,nodeids_FV,dim=3,DGFV=1)

  CALL WriteVarnamesToVTK_array(nVarTotal,mapVisu,varnames_out,components_out)

ELSE IF (VisuDimension.EQ.2) THEN

  ! allocate Coords 2D array and copy from first zeta-slice of 3D array
  SDEALLOCATE(CoordsVisu_DG_2D)
  ALLOCATE(CoordsVisu_DG_2D(1:3,0:NVisu,0:NVisu,0:0,1:nElems_DG))
  CoordsVisu_DG_2D = CoordsVisu_DG(:,:,:,0:0,:)
  SDEALLOCATE(CoordsVisu_FV_2D)
  ALLOCATE(CoordsVisu_FV_2D(1:3,0:NVisu_FV,0:NVisu_FV,0:0,1:nElems_FV))
#if FV_ENABLED    
  CoordsVisu_FV_2D = CoordsVisu_FV(:,:,:,0:0,:)
#else
  CoordsVisu_FV_2D = 0
#endif

  CALL WriteCoordsToVTK_array(NVisu   ,nElems_DG,coordsDG_out,nodeidsDG_out,&
      CoordsVisu_DG_2D,nodeids_DG_2D,dim=2,DGFV=0)
  CALL WriteCoordsToVTK_array(NVisu_FV,nElems_FV,coordsFV_out,nodeidsFV_out,&
      CoordsVisu_FV_2D,nodeids_FV_2D,dim=2,DGFV=1)

  CALL WriteVarnamesToVTK_array(nVarTotal,mapVisu,varnames_out,components_out)

ELSE IF (VisuDimension.EQ.1) THEN ! CSV along 1d line

  IF (nProcessors.GT.1) &
    CALL CollectiveStop(__STAMP__,"1D csv output along lines only supported for single execution")

  strOutputFile=TRIM(TIMESTAMP(TRIM(ProjectName)//'_extract1D',OutputTime))

  OPEN(NEWUNIT = iounit, STATUS='REPLACE',FILE=TRIM(strOutputFile)//'_DG.csv')
  DO iElem=1,nElems_DG
    DO i=0,NVisu
      WRITE(iounit,*) CoordsVisu_DG(1,i,0,0,iElem), UVisu_DG(i,0,0,iElem,:)
    END DO 
  END DO
  CLOSE(iounit) ! close the file

#if FV_ENABLED
  OPEN(NEWUNIT = iounit, STATUS='REPLACE',FILE=TRIM(strOutputFile)//'_FV.csv')
  DO iElem=1,nElems_FV
    DO i=0,NVisu_FV
      WRITE(iounit,*) CoordsVisu_FV(1,i,0,0,iElem), UVisu_FV(i,0,0,iElem,:)
    END DO 
  END DO
  CLOSE(iounit) ! close the file
#endif

END IF

END SUBROUTINE Visu3d_Build_VTK

!===================================================================================================================================
!> Deallocate the different NodeID arrays.
!===================================================================================================================================
SUBROUTINE visu3d_dealloc_nodeids() 
USE MOD_Posti_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(nodeids_DG)
SDEALLOCATE(nodeids_FV)
SDEALLOCATE(nodeids_DG_2D)
SDEALLOCATE(nodeids_FV_2D)
END SUBROUTINE visu3d_dealloc_nodeids


!===================================================================================================================================
!> Deallocate arrays used by visu3D.
!===================================================================================================================================
SUBROUTINE FinalizeVisu3D()
USE MOD_Globals,ONLY: MPIRoot
USE MOD_Posti_Vars
IMPLICIT NONE
!===================================================================================================================================
SWRITE (*,*) "VISU3D FINALIZE"
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

SDEALLOCATE(mapCalc)
#if FV_ENABLED && FV_RECONSTRUCT
SDEALLOCATE(mapCalc_FV)
#endif
SDEALLOCATE(mapVisu)
SDEALLOCATE(mapVisu_old)
SDEALLOCATE(UCalc_DG)
SDEALLOCATE(UCalc_FV)

SDEALLOCATE(mapElems_DG)
SDEALLOCATE(mapElems_FV)
SDEALLOCATE(FV_Elems_loc)

END SUBROUTINE FinalizeVisu3D

END MODULE MOD_Visu3D
