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
USE MOD_Posti_Vars,ONLY: VarNamesTotal,VarNames_ElemData,VarNames_FieldData,nVarTotal,nVar_ElemData,nVar_FieldData,FileType
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)             :: mpi_comm_IN    
INTEGER,INTENT(IN)             :: strlen_state
TYPE(C_PTR),TARGET,INTENT(IN)  :: statefile_IN
TYPE (CARRAY), INTENT(INOUT)   :: varnames
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
CHARACTER(LEN=255),POINTER     :: varnames_loc(:)
CHARACTER(LEN=255)             :: statefile
INTEGER                        :: s
!===================================================================================================================================
print*, 'requestinfostart'; CALL FLUSH()
statefile = cstrToChar255(statefile_IN, strlen_state)
CALL visu3d_InitFile(mpi_comm_IN, statefile)

print*, 'copystuff'; CALL FLUSH()
SELECT CASE(TRIM(FileType))
CASE('Mesh')
  varnames%len  = 0
  varnames%data = C_NULL_PTR
CASE('State','Generic')
print*, 'dealloc'; CALL FLUSH()
  ALLOCATE(varnames_loc(nVarTotal+nVar_ElemData+nVar_FieldData))
  
print*, 'copy'; CALL FLUSH()
  varnames_loc(1:nVarTotal) = VarNamesTotal(1:nVarTotal)
  
  s=nVarTotal
  IF(nVar_ElemData.GT.0) &
    varnames_loc(s+1:s+nVar_ElemData)  = VarNames_ElemData(1:nVar_ElemData)
  s=s+nVar_ElemData
  IF(nVar_FieldData.GT.0)&
    varnames_loc(s+1:s+nVar_FieldData) = VarNames_FieldData(1:nVar_FieldData)
print*, 'copy2'; CALL FLUSH()
  
  varnames%len  = SIZE(varnames_loc)*255 
print*, 'len'; CALL FLUSH()
  varnames%data = C_LOC(varnames_loc(1))
print*, 'data'; CALL FLUSH()
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,"The type of '"//TRIM(statefile)//"' is unknown.")
END SELECT

print*, 'requestinfodone'; CALL FLUSH()
END SUBROUTINE visu3d_requestInformation



!===================================================================================================================================
!> Create a list of available variables for ParaView. This list contains the conservative, primitve and derived quantities
!> that are available in the current equation system as well as the additional variables read from the state file.
!> The additional variables are stored in the datasets 'ElemData' (elementwise data) and 'FieldData' (pointwise data).
!===================================================================================================================================
SUBROUTINE visu3d_InitFile(mpi_comm_IN, statefile)
! MODULES
USE ISO_C_BINDING      ,ONLY: C_PTR, C_LOC, C_F_POINTER
USE MOD_Preproc
USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_EOS_Posti_Vars
USE MOD_MPI            ,ONLY: InitMPI
USE MOD_HDF5_Input     ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,GetArrayAndName
USE MOD_HDF5_Input     ,ONLY: ReadAttribute,File_ID,OpenDataFile,GetDataProps,CloseDataFile
USE MOD_Interpolation_Vars,ONLY: NodeType
USE MOD_Output_Vars    ,ONLY: ProjectName
USE MOD_StringTools    ,ONLY: STRICMP
USE MOD_Mesh_ReadIn    ,ONLY: BuildPartition
USE MOD_Posti_Mappings ,ONLY: Build_FV_DG_distribution
IMPLICIT NONE
INTEGER,INTENT(IN)             :: mpi_comm_IN    
CHARACTER(LEN=255),INTENT(IN)  :: statefile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
INTEGER                        :: nElems_State,iVar
INTEGER                        :: nVal(15) ! max rank according to F08 standard
CHARACTER(LEN=255)             :: NodeType_State 
REAL,ALLOCATABLE               :: tmp(:)
!===================================================================================================================================
print*, 'initfilestart'; CALL FLUSH()
CALL FLUSH()
CALL InitMPI(mpi_comm_IN) 

changedMeshFile      = .FALSE.
changedStateFile = .NOT.STRICMP(statefile,statefile_old)
IF(changedStateFile)THEN
  SWRITE(*,*) "ChangedFile: ", TRIM(statefile_old), " -> ",TRIM(statefile)
ELSE
  RETURN
END IF

nVarTotal=0
nVar_ElemData=0
nVar_FieldData=0
SDEALLOCATE(DepTable)
SDEALLOCATE(FieldData)
SDEALLOCATE(ElemData)
SDEALLOCATE(VarNamesTotal)
SDEALLOCATE(VarNames_ElemData)
SDEALLOCATE(VarNames_FieldData)

print*, 'checkfilest'; CALL FLUSH()
IF (ISVALIDMESHFILE(statefile)) THEN      ! MESH
  FileType      = 'Mesh'
  MeshFile      = statefile
  nVar_State    = 0
  withGradients = .FALSE.
  calcDeps      = .FALSE.
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! other file
  ! Read attributes
  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

  CALL GetDataProps(nVar_State,PP_N,nElems_State,NodeType_State)
  IF (.NOT.STRICMP(NodeType_State, NodeType)) THEN
    CALL CollectiveStop(__STAMP__, &
        "NodeType of state does not match with NodeType the visu3D-posti is compiled with!")
  END IF
  SDEALLOCATE(VarNamesHDF5); ALLOCATE(VarNamesHDF5(nVar_State))

  CALL ReadAttribute(File_ID,'File_Type',   1,StrScalar =FileType)
  CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =ProjectName)
  CALL ReadAttribute(File_ID,'Time',        1,RealScalar=OutputTime)
  CALL ReadAttribute(File_ID,'VarNames',nVar_State,StrArray=VarNamesHDF5)
  CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar =MeshFile_state)

  CALL CloseDataFile()
print*, 'getattribs'; CALL FLUSH()

  ! Build partition to get nElems
  CALL OpenDataFile(MeshFile_state,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL BuildPartition()
  CALL CloseDataFile()
print*, 'new partition'; CALL FLUSH()

  changedMeshFile = .NOT.(STRICMP(MeshFile_state,MeshFile_old))

  ! TODO: for other specific file types extend comparison by list
  IF(TRIM(FileType).NE.'State') FileType='Generic'

ELSE
  CALL CollectiveStop(__STAMP__,"'"//TRIM(statefile)//"' is not a valid mesh or state file.")
END IF


print*,'readdata'; CALL FLUSH()
SELECT CASE(TRIM(FileType))
CASE('Mesh')
CASE('State')

  ! Read in ElemData and additional Field Data
  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL GetArrayAndName('FieldData','VarNamesAddField',nVal,tmp,VarNames_FieldData)
  IF(ALL(nVal(1:5).GT.0))THEN
    ALLOCATE(FieldData(nVal(1),0:nVal(2)-1,0:nVal(3)-1,0:nVal(4)-1,nVal(5)))
    FieldData=RESHAPE(tmp,nVal(1:5))
    nVar_FieldData = nVal(1)
  END IF
  CALL GetArrayAndName('ElemData' ,'VarNamesAdd'     ,nVal,tmp,VarNames_ElemData)
  IF(ALL(nVal(1:2).GT.0))THEN
    ALLOCATE(ElemData(nVal(1),nVal(2)))
    ElemData=RESHAPE(tmp, nVal(1:2))
    nVar_ElemData  = nVal(1)
  END IF
  SDEALLOCATE(tmp)
  CALL CloseDataFile()

  CALL Build_FV_DG_distribution()


  nVarTotal=nVarTotalEOS
  ALLOCATE(VarNamesTotal(nVarTotal))
  VarNamesTotal=DepNames

  ALLOCATE(DepTable(1:nVarTotal,0:nVarTotal))
  DepTable=DepTableEOS
  calcDeps      = .TRUE.

CASE('Generic')

  CALL Build_FV_DG_distribution()

  nVarTotal=SIZE(VarNamesHDF5,1)
  ALLOCATE(VarNamesTotal(nVarTotal))

  VarNamesTotal=VarNamesHDF5
  ALLOCATE(DepTable(1:nVarTotal,0:nVarTotal))
  DepTable=0
  DO iVar=1,nVarTotal
    DepTable(iVar,iVar)=1
  END DO
  calcDeps      = .FALSE.

CASE DEFAULT
  CALL CollectiveStop(__STAMP__,"The type of '"//TRIM(statefile)//"' is unknown.")
END SELECT

print*,'initfiledone'; CALL FLUSH()
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
print*,'startwrap'; CALL FLUSH()
prmfile   = cstrToChar255(prmfile_IN,   strlen_prm)
postifile = cstrToChar255(postifile_IN, strlen_posti)
statefile = cstrToChar255(statefile_IN, strlen_state)
print*,'startwrapvisu'; CALL FLUSH()
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
USE MOD_Posti_Mappings      ,ONLY: Build_mapCalc_mapVisu
USE MOD_Posti_ReadState     ,ONLY: ReadStateAndGradients,ReadState
USE MOD_Posti_VisuMesh      ,ONLY: VisualizeMesh
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_DG
#if FV_ENABLED
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_ConvertToVisu_FV
#endif
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG,ConvertToVisu_ElemData,ConvertToVisu_FieldData
USE MOD_ReadInTools         ,ONLY: prms,FinalizeParameters
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
!===================================================================================================================================
print*,'start'; CALL FLUSH()
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
!   - changedWithGradients: different mode, with/without gradients
!   
! WORKFLOW:
! * The main steps are:
!   1. get nElems             (if changedStateFile)
!   2. get FV/DG distribution (if changedStateFile)
!   3. read solution          (if changedStateFile or changedWithGradients)
!   4. read Mesh              (if changedMeshFile)
!   5. compute UCalc          (if changedStateFile or changedVarNames) 
!   6. convert to UVisu       (if changedStateFile or changedVarNames or changedNVisu)
!   6. build visu mesh        (if changedMeshFile  or changedNVisu or changedFV_Elems)
!   7. write VTK arrays       (always!) 
!
!**********************************************************************************************

changedVarNames      = .FALSE.
changedNVisu         = .FALSE.
changedFV_Elems      = .FALSE.
changedWithGradients = .FALSE.

changedPrmFile   = (prmfile .NE. prmfile_old)

CALL FinalizeParameters()

! Read Varnames to visualize and build calc and visu dependencies
CALL prms%SetSection("posti")
CALL prms%CreateStringOption("MeshFile"     , "Custom mesh file ")
CALL prms%CreateStringOption("VarName"      , "Names of variables, which should be visualized.", multiple=.TRUE.)
CALL prms%CreateIntOption(   "NVisu"        ,  "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateIntOption(   "VisuDimension", "2 = Slice at first Gauss point in zeta-direction to get 2D solution.","3")
CALL prms%CreateStringOption("NodeTypeVisu" , "NodeType for visualization. Visu, Gauss,Gauss-Lobatto,Visu_inner"    ,"VISU")

SELECT CASE (TRIM(FileType))
CASE('Mesh')  ! visualize mesh

  SWRITE(*,*) "MeshFile Mode"
  CALL VisualizeMesh(postifile,MeshFile,coordsDG_out,valuesDG_out,nodeidsDG_out, &
      coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)

CASE('State','Generic')   ! visualize state (EOS related data)
  print*,'state'; CALL FLUSH()

  SWRITE(*,*) "State Mode"
  CALL Visu3D_Get_Ini(postifile)

  withGradients = .FALSE.
#if FV_RECONSTRUCT
  ! force calculation of metrics if there are any FV elements
  IF (hasFV_Elems) withGradients = .TRUE.
#endif

  print*,withGradients; CALL FLUSH()
  CALL Build_mapCalc_mapVisu()
  print*,'mapdone'; CALL FLUSH()

  changedWithGradients = (withGradients.NEQV.withGradients_old)

  ! get solution
  IF (changedStateFile.OR.changedWithGradients.OR.changedPrmFile) THEN
    CALL Visu3D_GetFlexiParameterFile(statefile,prmfile)
    SWRITE(*,*) "[ALL] get solution. withGradients = ", withGradients
    IF (withGradients) THEN
      CALL ReadStateAndGradients(prmfile,statefile)
    ELSE
      CALL ReadState(prmfile,statefile)
    END IF
  END IF
  print*,'readstatedone'; CALL FLUSH()

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
    CALL ConvertToVisu_ElemData()
    CALL ConvertToVisu_FieldData()
  END IF

  CALL Visu3D_Build_VTK(coordsDG_out,valuesDG_out,nodeidsDG_out, &
                        coordsFV_out,valuesFV_out,nodeidsFV_out, &
                        varnames_out,components_out)

CASE DEFAULT
  CALL CollectiveStop(__STAMP__,"The type of '"//TRIM(statefile)//"' is unknown.")
END SELECT

MeshFile_old          = MeshFile
prmfile_old           = prmfile
statefile_old         = statefile
NVisu_old             = NVisu
nVar_State_old        = nVar_State
withGradients_old     = withGradients
NodeTypeVisuPosti_old = NodeTypeVisuPosti

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(*,*) "Visu3D finished for state file: ", TRIM(statefile)
SWRITE(UNIT_StdOut,'(132("="))')
END SUBROUTINE visu3D

!===================================================================================================================================
!> Extact parameter file from state if non given
!===================================================================================================================================
SUBROUTINE Visu3D_GetFlexiParameterFile(statefile,prmfile)
USE MOD_Globals
USE MOD_ReadInTools         ,ONLY: ExtractParameterFile
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)    :: statefile
CHARACTER(LEN=255),INTENT(INOUT) :: prmfile
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
LOGICAL :: userblockFound
!===================================================================================================================================
IF (LEN_TRIM(prmfile).EQ.0) THEN
   prmfile = ".flexi.ini"
   CALL ExtractParameterFile(statefile,prmfile,userblockFound)
   IF (.NOT.userblockFound) THEN
     CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(statefile)//"'")
   END IF
END IF
END SUBROUTINE Visu3D_GetFlexiParameterFile


!===================================================================================================================================
!> Read Visu3D settings from parameterfile and build generic file if necessary
!===================================================================================================================================
SUBROUTINE Visu3D_Get_Ini(postifile)
USE MOD_Preproc
USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_StringTools         ,ONLY: STRICMP
USE MOD_ReadInTools         ,ONLY: prms,CountOption,GETINT,GETSTR,FinalizeParameters
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iVar,postiUnit,stat
!===================================================================================================================================

IF (LEN_TRIM(postifile).EQ.0) THEN
  postifile = ".posti.ini"
  ! no posti parameter file given => write a default .posti.ini file
  IF (MPIRoot) THEN
    OPEN(NEWUNIT=postiUnit,FILE=TRIM(postifile),STATUS='UNKNOWN',ACTION='WRITE',ACCESS='SEQUENTIAL',IOSTAT=stat)
    WRITE(postiUnit,'(A,I3)') "NVisu   = ", PP_N+1
    DO iVar=1,SIZE(VarNamesHDF5,1)
      WRITE(postiUnit,'(A,A)') "VarName = ", TRIM(VarNamesHDF5(iVar))
    END DO
    CLOSE(postiUnit)
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END IF

! Read Varnames to visualize and build calc and visu dependencies
CALL prms%read_options(postifile)
NVisu         = GETINT("NVisu")
VisuDimension = GETINT("VisuDimension")
Meshfile      = GETSTR("MeshFile",MeshFile_state)
NodeTypeVisuPosti = GETSTR('NodeTypeVisu')

SDEALLOCATE(VarNamesIni)
nVarIni = CountOption("VarName")
ALLOCATE(VarNamesIni(nVarIni))
DO iVar=1,nVarIni
  VarNamesIni(iVar) = GETSTR("VarName")
END DO

changedMeshFile = .NOT.(STRICMP(MeshFile,MeshFile_old))
changedNVisu = ((NVisu.NE.NVisu_old) .OR. (NodeTypeVisuPosti.NE.NodeTypeVisuPosti_old))

CALL FinalizeParameters()

END SUBROUTINE Visu3D_Get_Ini

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

  CALL WriteDataToVTK_array(nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG,3)
  CALL WriteDataToVTK_array(nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData,NVisu_FV,nElems_FV,valuesFV_out,UVisu_FV,3)

ELSE IF (VisuDimension.EQ.2) THEN

  ! allocate Visu 2D array and copy from first zeta-slice of 3D array
  SDEALLOCATE(UVisu_DG_2D)
  ALLOCATE(UVisu_DG_2D(0:NVisu,0:NVisu,0:0,1:nElems_DG,1:(nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData)))
  UVisu_DG_2D = UVisu_DG(:,:,0:0,:,:)
  SDEALLOCATE(UVisu_FV_2D)
  ALLOCATE(UVisu_FV_2D(0:NVisu_FV,0:NVisu_FV,0:0,1:nElems_FV,1:(nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData)))
#if FV_ENABLED
  UVisu_FV_2D = UVisu_FV(:,:,0:0,:,:)
#else
  CoordsVisu_FV_2D = 0
#endif

  CALL WriteDataToVTK_array(nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG_2D,2)
  CALL WriteDataToVTK_array(nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData,NVisu_FV,nElems_FV,valuesFV_out,UVisu_FV_2D,2)

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
withGradients_old = .FALSE.
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

SDEALLOCATE(ElemData)
SDEALLOCATE(VarNames_ElemData)
SDEALLOCATE(VarNamesVisu_ElemData)
SDEALLOCATE(VarNamesVisu_ElemData_old)

SDEALLOCATE(FieldData)
SDEALLOCATE(VarNames_FieldData)
SDEALLOCATE(VarNamesVisu_FieldData)
SDEALLOCATE(VarNamesVisu_FieldData_old)
END SUBROUTINE FinalizeVisu3D

END MODULE MOD_Visu3D
