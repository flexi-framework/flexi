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
MODULE MOD_Visu_Init
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE visu_getVarNamesAndFileType
  MODULE PROCEDURE visu_getVarNamesAndFileType
END INTERFACE

INTERFACE Visu_InitFile
  MODULE PROCEDURE Visu_InitFile
END INTERFACE

PUBLIC:: visu_getVarNamesAndFileType
PUBLIC:: visu_InitFile
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Create a list of available variables for ParaView. This list contains the conservative, primitive and derived quantities
!> that are available in the current equation system as well as the additional variables read from the state file.
!> The additional variables are stored in the datasets 'ElemData' (elementwise data) and 'FieldData' (pointwise data).
!> Also a list of all available boundary names is created for surface visualization.
!===================================================================================================================================
SUBROUTINE visu_getVarNamesAndFileType(statefile,meshfile,varnames_loc,bcnames_loc)
! MODULES
USE MOD_Globals
USE MOD_EOS_Posti_Vars ,ONLY: DepNames,nVarDepEOS
USE MOD_IO_HDF5        ,ONLY: GetDatasetNamesInGroup,File_ID
USE MOD_HDF5_Input     ,ONLY: OpenDataFile,CloseDataFile,GetDataSize,GetVarNames,ISVALIDMESHFILE,ISVALIDHDF5FILE,ReadAttribute
USE MOD_HDF5_Input     ,ONLY: DatasetExists,HSize,nDims,ReadArray
USE MOD_StringTools    ,ONLY: STRICMP
USE MOD_Restart        ,ONLY: InitRestartFile
USE MOD_Restart_Vars   ,ONLY: RestartMode
USE MOD_Visu_Vars      ,ONLY: FileType,VarNamesHDF5,nBCNamesAll,nVarIni,nVar_State,IJK_exists
USE MOD_Visu_Vars      ,ONLY: statefile_old,changedStateFile
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)                       :: statefile
CHARACTER(LEN=*)  ,INTENT(IN)                       :: meshfile
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: varnames_loc(:)
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: bcnames_loc(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: i,j,nVar,dims
LOGICAL                                             :: varnames_found,readDGsolutionVars,sameVars,VarNamesExist,file_exists
CHARACTER(LEN=255),ALLOCATABLE                      :: datasetNames(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: varnames_tmp(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: tmp(:)
CHARACTER(LEN=255)                                  :: MeshFile_loc
INTEGER                                             :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================

IF (ISVALIDMESHFILE(statefile)) THEN ! MESH
  SDEALLOCATE(varnames_loc)

  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  ! IJK-sorted mesh
  CALL DatasetExists(File_ID,'Elem_IJK'  ,IJK_exists)

  nVar = 3
  IF (IJK_exists)       nVar = nVar + 3            ! IJK sorting
  ALLOCATE(varnames_loc(nVar))

  varnames_loc(1) = 'ScaledJacobian'
  varnames_loc(2) = 'ScaledJacobianElem'
  varnames_loc(3) = 'ElemID'
  nVar            = 3
  ! Add IJK sorting
  IF (IJK_exists) THEN
    varnames_loc(nVar+1) = 'Elem_I'
    varnames_loc(nVar+2) = 'Elem_J'
    varnames_loc(nVar+3) = 'Elem_K'
    nVar                 = nVar + 3
  END IF
  FileType='Mesh'

ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! other file
  SDEALLOCATE(varnames_loc)
  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadAttribute(File_ID,'File_Type',   1,StrScalar =FileType)

  SELECT CASE(TRIM(FileType))
    ! check if variables in state file are the same as in the EQNSYS and set FileType to 'Generic' if not
    CASE('State')
      SDEALLOCATE(VarNamesHDF5)
      CALL GetVarNames("VarNames",VarNamesHDF5,VarNamesExist)

      sameVars = .FALSE.
      IF (VarNamesExist .AND. PP_nVar.EQ.SIZE(VarNamesHDF5)) THEN
        sameVars = .TRUE.
        DO i = 1,SIZE(VarNamesHDF5)
          sameVars = sameVars.AND.(STRICMP(VarNamesHDF5(i),DepNames(i)))
        END DO
      END IF
      IF (.NOT.sameVars) FileType = 'Generic'

    CASE('BaseFlow','TimeAvg')
      IF (nVarIni.EQ.0) THEN
        FileType = 'Generic'
      ELSE
        CALL CloseDataFile()
        ! This routine requires the file to be closed
        CALL InitRestartFile(statefile)
        CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
        SELECT CASE(RestartMode)
          ! TimeAvg with vars missing
          CASE(0)
            SDEALLOCATE(VarNamesHDF5)
            CALL GetVarNames("VarNames_Mean",VarNamesHDF5,VarNamesExist)
            FileType = 'Generic'

          ! BaseFlow
          CASE(1)
            SDEALLOCATE(VarNamesHDF5)
            CALL GetVarNames("VarNames",VarNamesHDF5,VarNamesExist)

            sameVars = .FALSE.
            FileType = 'State'
            IF (VarNamesExist .AND. PP_nVar.EQ.SIZE(VarNamesHDF5)) THEN
              sameVars = .TRUE.
              DO i = 1,SIZE(VarNamesHDF5)
                sameVars = sameVars.AND.(STRICMP(VarNamesHDF5(i),DepNames(i)))
              END DO
            END IF
            IF (.NOT.sameVars) FileType = 'Generic'

          ! TimeAvg with complete vars
          CASE(2,3)
            SDEALLOCATE(VarNamesHDF5)
            CALL GetVarNames("VarNames_Mean",VarNamesHDF5,VarNamesExist)
            FileType = 'State'
            ! When restarting from a time-averaged file, we convert to U array to PP_nVar
            nVar_State = PP_nVar
          CASE DEFAULT
            FileType = 'Generic'
        END SELECT
      END IF
  END SELECT

  IF (STRICMP(FileType,'State')) THEN
    nVar = nVarDepEOS
    ALLOCATE(varnames_loc(nVar))
    varnames_loc(1:nVar) = DepNames
    readDGsolutionVars   = .FALSE.
  ELSE
    nVar=0
    readDGsolutionVars   = .TRUE.
  END IF

  CALL GetDatasetNamesInGroup("/",datasetNames)

  DO i = 1,SIZE(datasetNames)
    SDEALLOCATE(varnames_tmp)
    VarNamesExist=.FALSE.
    CALL DatasetExists(File_ID,"VarNames_"//TRIM(datasetNames(i)),varnames_found,attrib=.TRUE.)
    IF (varnames_found) THEN
      CALL GetVarNames("VarNames_"//TRIM(datasetNames(i)),varnames_tmp,VarNamesExist)
    ELSE
      SELECT CASE(TRIM(datasetNames(i)))
        CASE('DG_Solution')
          IF (readDGsolutionVars) THEN
            CALL GetVarNames("VarNames",varnames_tmp,VarNamesExist)
          END IF
        ! CASE('Mean')
        !   IF (RestartMode.GT.1 .AND. readDGsolutionVars) THEN
        !     CALL GetVarNames("VarNames_Mean",varnames_tmp,VarNamesExist)
        !   END IF
        CASE('ElemTime')
          ! Ignore, already part of ElemData
          CYCLE
        CASE('ElemData')
          CALL GetVarNames("VarNamesAdd",varnames_tmp,VarNamesExist)
        CASE('FieldData')
          CALL GetVarNames("VarNamesAddField",varnames_tmp,VarNamesExist)
        CASE DEFAULT
          CALL GetDataSize(File_ID,TRIM(datasetNames(i)),dims,HSize)
          IF ((dims.NE.5).AND.(dims.NE.2)) CYCLE ! Do not add datasets to the list that can not contain elementwise or field data
          ALLOCATE(varnames_tmp(INT(HSize(1))))
          DO j=1,INT(HSize(1))
            WRITE(varnames_tmp(j),'(I0)') j
          END DO
          VarNamesExist=.TRUE.
          DEALLOCATE(HSize)
      END SELECT
    END IF ! varnames_found
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
  END DO ! i = 1,SIZE(datasetNames)

  IF (LEN_TRIM(meshfile).EQ.0) THEN
    ! Save mesh file to get boundary names later
    CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar =MeshFile_loc)
  ELSE
    MeshFile_loc = meshfile
  END IF

  CALL CloseDataFile()

  INQUIRE(FILE=TRIM(MeshFile_loc), EXIST=file_exists)
  IF (file_exists) THEN
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
  ELSE
    CALL PrintWarning('Mesh file "'//TRIM(MeshFile_loc)//'" specified in the loaded file does not exist! Surface visualization will not be available!')
  END IF

  SDEALLOCATE(datasetNames)
END IF

! check if state changed
changedStateFile = .NOT.STRICMP(statefile,statefile_old)
IF (changedStateFile .AND. MPIRoot) WRITE(*,*) "state file old -> new: ", TRIM(statefile_old), " -> ",TRIM(statefile)

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
USE MOD_Globals
USE MOD_Preproc
USE MOD_EOS_Posti_Vars
USE MOD_HDF5_Input         ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,GetArrayAndName
USE MOD_HDF5_Input         ,ONLY: ReadAttribute,File_ID,OpenDataFile,GetDataProps,CloseDataFile,ReadArray,DatasetExists
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_MPI                ,ONLY: InitMPI
USE MOD_Output_Vars        ,ONLY: ProjectName
USE MOD_Posti_Mappings     ,ONLY: Build_FV_DG_distribution,Build_mapDepToCalc_mapAllVarsToVisuVars
USE MOD_ReadInTools        ,ONLY: prms,addStrListEntry,FinalizeParameters,CountOption
USE MOD_ReadInTools        ,ONLY: GETINT,GETINTFROMSTR,GETLOGICAL,GETSTR
USE MOD_StringTools        ,ONLY: STRICMP,INTTOSTR
USE MOD_Visu_Avg2D         ,ONLY: InitAverage2D,BuildVandermonds_Avg2D
USE MOD_Visu_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)    :: statefile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
LOGICAL                          :: RestartMean
CHARACTER(LEN=255)               :: cwd
INTEGER                          :: nElems_State
#if PP_N!=N
INTEGER                          :: N_State
#endif
!===================================================================================================================================
IF (STRICMP(fileType,'Mesh')) THEN
    CALL CollectiveStop(__STAMP__, &
        "FileType==Mesh, but we try to initialize a state file!")
END IF

! open state file to be able to read attributes
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! read the meshfile attribute from statefile
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar =MeshFile_state)

! read options from posti parameter file
CALL prms%read_options(postifile)

! read output directory
SELECT CASE(CountOption("OutputDirectory"))
  CASE(0) ! Do nothing
  CASE(1)
    OutputDirectory = GETSTR("OutputDirectory")
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'OutputDirectory is not a multiple option!')
END SELECT

! Get number of variables to be visualized
nVarIni = CountOption("VarName")

! get properties
#if EQNSYSNR != 1
! Check if the file is a time-averaged file
CALL DatasetExists(File_ID,'Mean',RestartMean)
! Read in attributes
IF (.NOT.RestartMean .OR. nVarIni.EQ.0) THEN
#endif /* EQNSYSNR != 1 */
#if PP_N==N
  CALL GetDataProps(nVar_State,PP_N,nElems_State,NodeType_State)
#else
  CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State)
#endif /* PP_N==N */
#if EQNSYSNR != 1
! Check if the file is a time-averaged file
ELSE
#if PP_N==N
  CALL GetDataProps(nVar_State,PP_N,nElems_State,NodeType_State,'Mean')
#else
  CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State,'Mean')
#endif /* PP_N==N */
END IF
#endif /* EQNSYSNR != 1 */

! read options from posti parameter file
NVisu             = GETINT("NVisu",INTTOSTR(PP_N))
HighOrder         = GETLOGICAL('HighOrder')

! again read MeshFile from posti prm file (this overwrites the MeshFile read from the state file)
Meshfile          =  GETSTR("MeshFile",MeshFile_state)
IF (.NOT.FILEEXISTS(MeshFile) .OR. ((Meshfile(1:1) .NE. "/") .OR. (Meshfile(1:1) .NE. "~") .OR. (Meshfile(1:1) .NE. "."))) THEN
  !!!!!!
  ! WARNING: GETCWD is a GNU extension to the Fortran standard and will probably not work on other compilers
  CALL GETCWD(cwd)
  !!!!!!
  Meshfile          =  TRIM(cwd) // "/" // TRIM(Meshfile)
END IF
Avg2D             = GETLOGICAL("Avg2D")
#if PP_dim == 2
IF (Avg2D) THEN
  CALL PrintWarning("Avg2D not available for 2D-Posti! Switching it OFF.")
  Avg2D = .FALSE.
END IF
#endif
NodeTypeVisuPosti = GETSTR('NodeTypeVisu')
DGonly            = GETLOGICAL('DGonly')
CALL CloseDataFile()

CALL visu_getVarNamesAndFileType(statefile,'',VarnamesAll,BCNamesAll)

CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Get output format for visualization
OutputFormat = GETINTFROMSTR('OutputFormat')
SELECT CASE(OutputFormat)
  CASE(OUTPUTFORMAT_NONE)
    CALL CollectiveStop(__STAMP__,'Output disabled, exiting.')
  ! CASE(OUTPUTFORMAT_TECPLOT)
  !   CALL CollectiveStop(__STAMP__,'Tecplot output removed due to license issues (possible GPL incompatibility).')
  ! CASE(OUTPUTFORMAT_TECPLOTASCII)
  !   CALL CollectiveStop(__STAMP__,'Tecplot output removed due to license issues (possible GPL incompatibility).')
END SELECT

! check if state, mesh, NVisu, DGonly or Avg2D changed
changedStateFile = .NOT.STRICMP(statefile,statefile_old)
changedMeshFile  = .NOT.STRICMP(MeshFile,MeshFile_old)
changedDGonly    = (DGonly.NEQV.DGonly_old)
changedAvg2D     = (Avg2D .NEQV.Avg2D_old)

IF (changedStateFile .AND. MPIRoot) WRITE(*,*) "state file old -> new: ", TRIM(statefile_old), " -> ",TRIM(statefile)
IF (changedMeshFile  .AND. MPIRoot) WRITE(*,*) " mesh file old -> new: ", TRIM(MeshFile_old) , " -> ",TRIM(MeshFile)

! if Mesh or State changed readin some more attributes/parameters
IF (changedStateFile.OR.changedMeshFile) THEN
  CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =ProjectName)
  CALL ReadAttribute(File_ID,'Time',        1,RealScalar=OutputTime)
  ! If the polynomial degree is changing, we could need new mesh mappings.
  changedMeshFile = (NState_old.NE.PP_N)
  changedNodeType = .NOT.STRICMP(NodeType_State,NodeType_State_old)
END IF

CALL CloseDataFile()

! Polynomial degree for calculations
NCalc                 = GETINT("NCalc",INTTOSTR(PP_N))
IF (NCalc.LE.0) NCalc = PP_N
changedNCalc          = NCalc.NE.NCalc_old

! Output of averaged data is only available for NVisu = PP_N and NodeTypeVisuPosti=NodeType_State
! These settings are enforced here!
IF (Avg2D .AND. OutputFormat.EQ.OUTPUTFORMAT_HDF5) THEN
  NVisu             = PP_N
  NodeTypeVisuPosti = NodeType_State
END IF
! Check for changed visualization basis here to take change done for average output into account
changedNVisu     = ((NVisu.NE.NVisu_old) .OR. (NodeTypeVisuPosti.NE.NodeTypeVisuPosti_old))

! set number of dependent and raw variables
SDEALLOCATE(DepTable)
SDEALLOCATE(DepSurfaceOnly)
SDEALLOCATE(DepVolumeOnly)
nVarAll=SIZE(VarnamesAll)
IF (STRICMP(FileType,'State')) THEN
  StateFileMode = .TRUE.
  nVarDep = nVarDepEOS
  ALLOCATE(DepTable(nVarDep,0:nVarDep))
  ALLOCATE(DepSurfaceOnly(nVarDep))
  ALLOCATE(DepVolumeOnly(nVarDep))
  DepTable = DepTableEOS
  DepSurfaceOnly = DepSurfaceOnlyEOS
  DepVolumeOnly  = DepVolumeOnlyEOS
ELSE
  StateFileMode = .FALSE.
  nVarDep = 0
  ALLOCATE(DepTable(nVarDep,0:nVarDep))
  ALLOCATE(DepSurfaceOnly(nVarDep))
  ALLOCATE(DepVolumeOnly(nVarDep))
  DepTable = 0
  DepSurfaceOnly = 0
  DepVolumeOnly  = 0
END IF

! build distribution of FV and DG elements, which is stored in FV_Elems_loc
IF (changedStateFile.OR.changedMeshFile.OR.changedDGonly) THEN
  CALL Build_FV_DG_distribution(&
#if FV_ENABLED
    statefile&
#endif
    )
END IF

! reset withDGOperator flag and check if it is needed due to existing FV elements
withDGOperator = .FALSE.
#if FV_RECONSTRUCT
! If what we want to visualize is a state and has FV elements, the DG operator needs to be called for reconstruction
IF (StateFileMode) THEN
  IF (hasFV_Elems) withDGOperator = .TRUE.
END IF
#endif

! build mappings of variables which must be calculated/visualized
! also set withDGOperator flag if a dependent variable requires the evaluation of the DG operator
CALL Build_mapDepToCalc_mapAllVarsToVisuVars()

IF (Avg2D) THEN
  CALL InitAverage2D()
  CALL BuildVandermonds_Avg2D(NCalc&
#if FV_ENABLED
    ,NCalc_FV&
#endif
    )
END IF

changedWithDGOperator = (withDGOperator.NEQV.withDGOperator_old)

END SUBROUTINE visu_InitFile

END MODULE MOD_Visu_Init
