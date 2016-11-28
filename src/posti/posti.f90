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

MODULE MOD_Visu3D
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
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

INTERFACE visu3D
  MODULE PROCEDURE visu3D
END INTERFACE

INTERFACE visu3D_CWrapper
  MODULE PROCEDURE visu3D_CWrapper
END INTERFACE

INTERFACE visu3d_dealloc_nodeids
  MODULE PROCEDURE visu3d_dealloc_nodeids
END INTERFACE

INTERFACE visu3D_finalize
  MODULE PROCEDURE visu3D_finalize
END INTERFACE

PUBLIC:: visu3d_requestInformation
PUBLIC:: visu3D
PUBLIC:: visu3D_CWrapper
PUBLIC:: visu3d_dealloc_nodeids
PUBLIC:: visu3D_finalize

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


SUBROUTINE visu3d_requestInformation(mpi_comm_IN, strlen_state, statefile_IN, varnames)
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE ISO_C_BINDING      ,ONLY: C_PTR, C_LOC, C_F_POINTER
USE MOD_Globals
USE MOD_MPI            ,ONLY: InitMPI
USE MOD_HDF5_Input     ,ONLY: ISVALIDMESHFILE
USE MOD_EOS_Posti_Vars ,ONLY: nVarTotal
USE MOD_EOS_Posti      ,ONLY: GetVarnames
USE MOD_HDF5_Input     ,ONLY: GetDataSize,DatasetExists,ReadAttribute,File_ID,nDims,HSize,OpenDataFile,CloseDataFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
INTEGER,INTENT(IN)             :: mpi_comm_IN    
INTEGER,INTENT(IN)             :: strlen_state
TYPE(C_PTR),TARGET,INTENT(IN)  :: statefile_IN
TYPE (CARRAY), INTENT(INOUT)   :: varnames
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
CHARACTER(LEN=255),POINTER     :: varnames_loc(:)
CHARACTER(LEN=255)             :: statefile
LOGICAL                        :: elemDataFound
INTEGER                        :: nVarAdd,iVar
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAdd(:)
LOGICAL                        :: FieldDataFound
INTEGER                        :: nVarFieldData
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAddField(:)
!===================================================================================================================================
statefile = cstrToChar255(statefile_IN, strlen_state)

CALL InitMPI(mpi_comm_IN) 

IF (ISVALIDMESHFILE(statefile)) THEN
  varnames%len     = 0
ELSE ! state file


  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  ! check if additional elem data exists
  CALL DatasetExists(File_ID,'ElemData',elemDataFound)
  IF(elemDataFound)THEN
    CALL GetDataSize(File_ID,'ElemData',nDims,HSize)
    nVarAdd=INT(HSize(1),4)
    ALLOCATE(VarNamesAdd(nVarAdd))
    CALL ReadAttribute(File_ID,'VarNamesAdd',nVarAdd,StrArray=VarNamesAdd)
  ELSE
    nVarAdd=0
  END IF
  ! check if additional elem data exists
  CALL DatasetExists(File_ID,'FieldData',FieldDataFound)
  IF(FieldDataFound)THEN
    CALL GetDataSize(File_ID,'FieldData',nDims,HSize)
    nVarFieldData=INT(HSize(1),4)
    ALLOCATE(VarNamesAddField(nVarFieldData))
    CALL ReadAttribute(File_ID,'VarNamesAddField',nVarFieldData,StrArray=VarNamesAddField)
  ELSE
    nVarFieldData=0
  END IF
  CALL CloseDataFile()

  ! Fill varnames
  CALL GetVarnames(varnames_loc,nVarAdd+nVarFieldData)
  ! Add additional element data
  DO iVar=1,nVarAdd
    varnames_loc(nVarTotal+iVar) = VarNamesAdd(iVar)
  END DO
  ! Add additional pointwise data
  DO iVar=1,nVarFieldData
    varnames_loc(nVarTotal+nVarAdd+iVar) = VarNamesAddField(iVar)
  END DO

  varnames%len  = (nVarTotal+nVarAdd+nVarFieldData)*255 
  varnames%data = C_LOC(varnames_loc(1))

  SDEALLOCATE(VarNamesAdd)
  SDEALLOCATE(VarNamesAddField)
END IF 

END SUBROUTINE visu3d_requestInformation


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


SUBROUTINE visu3D(mpi_comm_IN, prmfile, postifile, statefile, &
        coordsDG_out,valuesDG_out,nodeidsDG_out, &
        coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Posti_Vars
USE MOD_EOS_Posti_Vars      ,ONLY: nVarTotal
USE MOD_Posti_Mappings      ,ONLY: Build_FV_DG_distribution,Build_mapCalc_mapVisu
USE MOD_Posti_ReadState     ,ONLY: ReadStateAndGradients,ReadState,ReadFieldData
USE MOD_Posti_VisuMesh      ,ONLY: BuildVisuCoords,VisualizeMesh
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_DG
#if FV_ENABLED
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_ConvertToVisu_FV
#endif
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG,ConvertToVisu_ElemData,ConvertToVisu_FieldData
USE MOD_MPI                 ,ONLY: InitMPI
USE MOD_HDF5_Input          ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE
USE MOD_HDF5_Input          ,ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_Mesh_ReadIn         ,ONLY: BuildPartition
USE MOD_Interpolation_Vars  ,ONLY: NodeTypeVisu,NodeTypeFVEqui
USE MOD_VTK                 ,ONLY: WriteCoordsToVTK_array,WriteDataToVTK_array,WriteVarnamesToVTK_array
USE MOD_StringTools         ,ONLY: STRICMP,LowCase
USE MOD_Interpolation_Vars  ,ONLY: NodeType
USE MOD_ReadInTools         ,ONLY: prms,GETINT,GETLOGICAL
USE MOD_ReadInTools         ,ONLY: FinalizeParameters,ExtractParameterFile
USE MOD_Output_Vars         ,ONLY: ProjectName
USE MOD_Restart             ,ONLY: ReadElemData
! IMPLICIT VARIABLE HANDLING
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
! LOCAL VARIABLES
INTEGER                          :: postiUnit = 104
INTEGER                          :: stat
INTEGER                          :: nElems_State,iVar
CHARACTER(LEN=255)               :: NodeType_State 
CHARACTER(LEN=255),ALLOCATABLE   :: StrVarNames(:)
LOGICAL                          :: userblockFound
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
!   - changedNVisu:         new NVisu
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
changedStateFile     = .FALSE.
changedMeshFile      = .FALSE.
changedNVisu         = .FALSE.
changedVarNames      = .FALSE.
changedFV_Elems      = .FALSE.
changedWithGradients = .FALSE.

IF (ISVALIDMESHFILE(statefile)) THEN ! visualize mesh
  SWRITE(*,*) "MeshFile Mode"
  MeshFile      = statefile
  nVar_State    = 0
  withGradients = .FALSE.
  CALL VisualizeMesh(postifile,MeshFile,coordsDG_out,valuesDG_out,nodeidsDG_out, &
      coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! visualize state file
  changedStateFile = .NOT.STRICMP(statefile,statefile_old)
  SWRITE(*,*) "ChangedStateFile: ", TRIM(statefile_old), " -> ",TRIM(statefile)
  IF (changedStateFile) THEN
    ! Read in parameters from the State file
    CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    CALL GetDataProps(nVar_State,PP_N,nElems_State,NodeType_State)
    IF (.NOT.STRICMP(NodeType_State, NodeType)) THEN
      CALL CollectiveStop(__STAMP__, &
          "NodeType of state does not match with NodeType the visu3D-posti is compiled with!")
    END IF
    CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar =MeshFile)
    CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =ProjectName)
    CALL ReadAttribute(File_ID,'Time',        1,RealScalar=OutputTime)
    changedMeshFile = .NOT.(STRICMP(MeshFile,MeshFile_old))
    SDEALLOCATE(StrVarNames)
    ALLOCATE(StrVarNames(nVar_State))
    CALL ReadAttribute(File_ID,'VarNames',nVar_State,StrArray=StrVarNames)
    CALL CloseDataFile()

    ! Build partition to get nElems
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    CALL BuildPartition() 
    CALL CloseDataFile()

    ! Read in ElemData and additional Field Data - needs nElems and can be called only after BuildPartition
    CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    CALL ReadElemData()
    CALL ReadFieldData()
    CALL CloseDataFile()

    CALL Build_FV_DG_distribution(statefile) 
  END IF ! changedStateFile

  IF (LEN_TRIM(postifile).EQ.0) THEN
    postifile = ".posti.ini"
    ! no posti parameter file given => write a default .posti.ini file
    IF (MPIRoot) THEN
      OPEN(UNIT=postiUnit,FILE=TRIM(postifile),STATUS='UNKNOWN',ACTION='WRITE',ACCESS='SEQUENTIAL',IOSTAT=stat)
      WRITE(postiUnit,'(A,I3)') "NVisu   = ", PP_N+1
      DO iVar=1,nVar_State
        WRITE(postiUnit,'(A,A)') "VarName = ", TRIM(StrVarNames(iVar))
      END DO
      CLOSE(postiUnit)
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
  END IF

  ! Read Varnames to visualize and build calc and visu dependencies
  CALL prms%SetSection("posti")
  CALL prms%CreateStringOption('VarName',"Names of variables, which should be visualized.", multiple=.TRUE.)
  CALL prms%CreateIntOption(   'NVisu',  "Polynomial degree at which solution is sampled for visualization.")
  CALL prms%CreateIntOption(   'VisuDimension',"2 = Slice at first Gauss point in zeta-direction to get 2D solution.","3")
  CALL prms%read_options(postifile)
  NVisu  = GETINT("NVisu")
  VisuDimension = GETINT("VisuDimension")

  changedNVisu = (NVisu.NE.NVisu_old)

  withGradients = .FALSE.
#if FV_ENABLED && FV_RECONSTRUCT
  ! force calculation of metrics if there are any FV elements
  IF (hasFV_Elems) withGradients = .TRUE.
#endif
  CALL Build_mapCalc_mapVisu()
  CALL FinalizeParameters()

  changedWithGradients = (withGradients.NEQV.withGradients_old)

  ! get solution 
  changedPrmFile = (prmfile .NE. prmfile_old)
  IF (changedStateFile.OR.changedWithGradients.OR.changedPrmFile) THEN
    ! extact parameter file if non given
    IF (LEN_TRIM(prmfile).EQ.0) THEN
       prmfile = ".flexi.ini"
       CALL ExtractParameterFile(statefile,prmfile,userblockFound)
       IF (.NOT.userblockFound) THEN
         CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(statefile)//"'")
       END IF
    END IF     
    SWRITE(*,*) "[ALL] get solution. withGradients = ", withGradients
    IF (withGradients) THEN
      CALL ReadStateAndGradients(prmfile,statefile)
    ELSE
      CALL ReadState(prmfile,statefile)
    END IF
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
  IF (changedStateFile.OR.changedVarNames) THEN
    CALL CalcQuantities_ConvertToVisu_FV()
  END IF
#endif /* FV_ENABLED */

  ! convert ElemData and FieldData to visu grid
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu) THEN
    CALL ConvertToVisu_ElemData()
    CALL ConvertToVisu_FieldData()
  END IF

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
    UVisu_FV_2D = UVisu_FV(:,:,0:0,:,:)

    CALL WriteDataToVTK_array(nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData,NVisu   ,nElems_DG,valuesDG_out,UVisu_DG_2D,2)
    CALL WriteDataToVTK_array(nVarVisu+nVarVisu_ElemData+nVarVisu_FieldData,NVisu_FV,nElems_FV,valuesFV_out,UVisu_FV_2D,2)
  END IF

  ! Convert coordinates to visu grid
  changedMeshFile = .NOT.(STRICMP(MeshFile,MeshFile_old))
  IF (changedMeshFile.OR.changedNVisu.OR.changedFV_Elems) THEN
    CALL BuildVisuCoords(withGradients)
  END IF

  ! write coords, UVisu to VTK  2D / 3D arrays (must be done always!)
  IF (VisuDimension.EQ.3) THEN
    CALL WriteCoordsToVTK_array(NVisu   ,nElems_DG,coordsDG_out,nodeidsDG_out,&
        CoordsVisu_DG,nodeids_DG,dim=3,DGFV=0)
    CALL WriteCoordsToVTK_array(NVisu_FV,nElems_FV,coordsFV_out,nodeidsFV_out,&
        CoordsVisu_FV,nodeids_FV,dim=3,DGFV=1)
  ELSE IF (VisuDimension.EQ.2) THEN
    ! allocate Coords 2D array and copy from first zeta-slice of 3D array
    SDEALLOCATE(CoordsVisu_DG_2D)
    ALLOCATE(CoordsVisu_DG_2D(1:3,0:NVisu,0:NVisu,0:0,1:nElems_DG))
    CoordsVisu_DG_2D = CoordsVisu_DG(:,:,:,0:0,:)
    SDEALLOCATE(CoordsVisu_FV_2D)
    ALLOCATE(CoordsVisu_FV_2D(1:3,0:NVisu_FV,0:NVisu_FV,0:0,1:nElems_FV))
    CoordsVisu_FV_2D = CoordsVisu_FV(:,:,:,0:0,:)

    CALL WriteCoordsToVTK_array(NVisu   ,nElems_DG,coordsDG_out,nodeidsDG_out,&
        CoordsVisu_DG_2D,nodeids_DG_2D,dim=2,DGFV=0)
    CALL WriteCoordsToVTK_array(NVisu_FV,nElems_FV,coordsFV_out,nodeidsFV_out,&
        CoordsVisu_FV_2D,nodeids_FV_2D,dim=2,DGFV=1)
  END IF

  ! write varnames to VTK (must be done always!)
  CALL WriteVarnamesToVTK_array(nVarTotal,mapVisu,varnames_out,components_out)
ELSE
  CALL CollectiveStop(__STAMP__,"'"//TRIM(statefile)//"' is not a valid mesh or state file.")
END IF ! visualize mesh/state file

MeshFile_old      = MeshFile
prmfile_old       = prmfile
statefile_old     = statefile
NVisu_old         = NVisu
nVar_State_old    = nVar_State
withGradients_old = withGradients

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(*,*) "Vius3D finished for state file: ", TRIM(statefile)
SWRITE(UNIT_StdOut,'(132("="))')
END SUBROUTINE visu3D


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


SUBROUTINE visu3D_finalize() 
USE MOD_Posti_Vars
IMPLICIT NONE
!===================================================================================================================================
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

SDEALLOCATE(VarNamesVisu_ElemData)
SDEALLOCATE(VarNamesVisu_ElemData_old)
END SUBROUTINE visu3D_finalize

END MODULE MOD_Visu3D




PROGRAM Posti_Visu3D
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_Posti_Vars
USE MOD_Commandline_Arguments
USE MOD_Visu3D
USE MOD_ISO_VARYING_STRING
USE MOD_MPI                   ,ONLY: InitMPI
USE MOD_VTK                   ,ONLY: WriteDataToVTK,WriteVTKMultiBlockDataSet
USE MOD_Output_Vars           ,ONLY: ProjectName
USE MOD_StringTools           ,ONLY: STRICMP,GetFileExtension
impliCIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iArg
CHARACTER(LEN=255),TARGET      :: prmfile
CHARACTER(LEN=255),TARGET      :: postifile
CHARACTER(LEN=255),TARGET      :: statefile
TYPE (CARRAY)                  :: coordsDG_out
TYPE (CARRAY)                  :: valuesDG_out
TYPE (CARRAY)                  :: nodeidsDG_out
TYPE (CARRAY)                  :: coordsFV_out
TYPE (CARRAY)                  :: valuesFV_out
TYPE (CARRAY)                  :: nodeidsFV_out
TYPE (CARRAY)                  :: varnames_out
TYPE (CARRAY)                  :: components_out

INTEGER                        :: skipArgs
INTEGER                        :: nVal
CHARACTER(LEN=255),POINTER     :: VarNames_p(:)
REAL,POINTER                   :: Coords_p(:,:,:,:,:)
REAL,POINTER                   :: Values_p(:,:,:,:,:)
CHARACTER(LEN=255)             :: FileString_DG
CHARACTER(LEN=255)             :: FileString_FV
CHARACTER(LEN=255)             :: FileString_multiblock
#if !USE_MPI
INTEGER                        :: MPI_COMM_WORLD = 0
#endif
INTEGER                        :: NVisu_k,NVisu_k_FV
!==================================================================================================================================
CALL InitMPI()
CALL ParseCommandlineArguments()
IF (nArgs.LT.1) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti [posti-prm-file [flexi-prm-file]] statefile [statefiles]')
END IF

prmfile = ""
! check if parameter file is given
IF(STRICMP(GetFileExtension(Args(1)),'ini')) THEN
  skipArgs = 1 ! first argument is the parameter file
  postifile = Args(1)
  ! check if a second parameter file is given (this is used instead of the parameter file stored in the userblock of a state file)
  IF (nArgs.GT.2) THEN
    IF (STRICMP(GetFileExtension(Args(2)),'ini')) THEN
      prmfile = Args(2)
      skipArgs = 2
    END IF
  END IF
ELSE IF(STRICMP(GetFileExtension(Args(1)),'h5')) THEN
  skipArgs = 0 ! do not skip a argument. first argument is a h5 file
  postifile = ""
ELSE
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti [posti-prm-file [flexi-prm-file]] statefile [statefiles]')
END IF

DO iArg=1+skipArgs,nArgs
  statefile = TRIM(Args(iArg))
  SWRITE(*,*) "Processing state-file: ",TRIM(statefile)
  
  CALL visu3D(MPI_COMM_WORLD, prmfile, postifile, statefile, &
      coordsDG_out,valuesDG_out,nodeidsDG_out, &
      coordsFV_out,valuesFV_out,nodeidsFV_out,varnames_out,components_out)

#if FV_ENABLED                            
  FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DG',OutputTime))//'.vtu'
#else
  FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtu'
#endif
  nVal = varnames_out%len/255
  CALL C_F_POINTER(varnames_out%data, VarNames_p, [nVal])
  NVisu_k = NVisu * (VisuDimension-2)
  CALL C_F_POINTER(coordsDG_out%data, Coords_p, [3,NVisu+1,NVisu+1,NVisu_k+1,nElems_DG])
  CALL C_F_POINTER(valuesDG_out%data, Values_p, [  NVisu+1,NVisu+1,NVisu_k+1,nElems_DG,nVal])
  CALL WriteDataToVTK(nVal,NVisu   ,nElems_DG,VarNames_p,Coords_p,Values_p,FileString_DG,&
      dim=VisuDimension,DGFV=0,nValAtLastDimension=.TRUE.)
#if FV_ENABLED                            
  FileString_FV=TRIM(TIMESTAMP(TRIM(ProjectName)//'_FV',OutputTime))//'.vtu'
  NVisu_k_FV = NVisu_FV * (VisuDimension-2)
  CALL C_F_POINTER(coordsFV_out%data, Coords_p, [3,NVisu_FV+1,NVisu_FV+1,NVisu_k_FV+1,nElems_FV])
  CALL C_F_POINTER(valuesFV_out%data, Values_p, [  NVisu_FV+1,NVisu_FV+1,NVisu_k_FV+1,nElems_FV,nVal])
  CALL WriteDataToVTK(nVal,NVisu_FV,nElems_FV,VarNames_p,Coords_p,Values_p,FileString_FV,&
      dim=VisuDimension,DGFV=1,nValAtLastDimension=.TRUE.)

  IF (MPIRoot) THEN                   
    ! write multiblock file
    FileString_multiblock=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtm'
    CALL WriteVTKMultiBlockDataSet(FileString_multiblock,FileString_DG,FileString_FV)
  ENDIF
#endif
END DO
END PROGRAM 

