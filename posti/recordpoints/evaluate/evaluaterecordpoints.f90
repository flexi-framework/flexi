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
!> Tool to extract recordpoint type signals which are normally generated during the simulation, from state files or from other
!> type of files (e.g. TimeAverages)
!> If state files are used, no parameter files has to be given, it will re-use an extracted parameter file from the userblock.
!===================================================================================================================================
PROGRAM evalrec
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_ReadInTools,       ONLY:prms,IgnoredParameters,PrintDefaultParameterFile,FinalizeParameters,ExtractParameterFile
USE MOD_ReadInTools,       ONLY:GETSTR
USE MOD_StringTools       ,ONLY:STRICMP,GetFileExtension
USE MOD_HDF5_Input,        ONLY:ISVALIDHDF5FILE
!!! These variables are actually modified by the tool itself !!!
USE MOD_MPI,               ONLY:DefineParametersMPI,InitMPI
#if USE_MPI
USE MOD_MPI,               ONLY:InitMPIvars,FinalizeMPI
#endif
USE MOD_IO_HDF5
USE MOD_HDF5_Input,        ONLY:GetDataSize,ReadArray,OpenDataFile,CloseDataFile,DatasetExists,GetVarNames,ReadAttribute
USE MOD_Output,            ONLY:DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Mesh,              ONLY:DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars,         ONLY:nElems,nGlobalElems,OffsetElem
USE MOD_Equation_Vars,     ONLY:StrVarNames
USE MOD_Interpolation,     ONLY:DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Interpolation_Vars,ONLY:NodeType
USE MOD_RecordPoints
USE MOD_RecordPoints_Vars, ONLY:RP_onProc
USE MOD_Restart,           ONLY:DefineParametersRestart,InitRestart,Restart,FinalizeRestart
USE MOD_Restart_Vars
!!! These variables are actually modified by the tool itself !!!
USE MOD_DG_Vars,           ONLY:U
USE MOD_Analyze_Vars,      ONLY:WriteData_dt,tWriteData
USE MOD_Timedisc_Vars,     ONLY:dt
#if FV_ENABLED
USE MOD_FV_Vars,           ONLY:FV_Elems
USE MOD_Indicator_Vars,    ONLY:IndValue
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg,start
CHARACTER(LEN=255)                 :: FileType
LOGICAL                            :: isValid,userblockFound
LOGICAL                            :: stateFileMode = .TRUE.
CHARACTER(LEN=255)                 :: DataSetName
INTEGER                            :: HSize_proc(5),nVar,HSize_tmp(5)
CHARACTER(LEN=255),ALLOCATABLE     :: StrVarNames_loc(:)
LOGICAL                            :: VarNamesExist,varnames_found
INTEGER                            :: j
CHARACTER(LEN=255)                 :: NodeType_HDF5
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL InitMPIInfo()
CALL ParseCommandlineArguments()
SWRITE(UNIT_stdOut,'(A)') " ||=============================||"
SWRITE(UNIT_stdOut,'(A)') " || Recordpoints Evaluation Tool||"
SWRITE(UNIT_stdOut,'(A)') " ||=============================||"


CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersRestart()
CALL DefineParametersOutput()
CALL DefineParametersMesh()
CALL DefineParametersRecordPoints()
! Parameters for the tool itself
CALL prms%SetSection("EvaluateRecordPoints")
CALL prms%CreateStringOption('RecordpointsDataSetName', "If no state files are given to evaluate, specify the data set name to be&
                                                         & used here.")

IF (doPrintHelp.LE.0) THEN
  ! Check if the number of arguments is correct
  IF (nArgs.LT.2) THEN
    ! Print out error message containing valid syntax
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: evalrec [parameter.ini] statefile1.h5...statefileN.h5 or evalrec --help'// &
    '[option/section name] to print help for a single parameter, parameter sections or all parameters.')
  END IF

  ! Check if a parameter file was given
  iArg=1
  ParameterFile='NOT SET'
  IF(STRICMP(GetFileExtension(Args(iArg)), "ini"))THEN
    ParameterFile = Args(iArg)
    iArg=iArg+1
  END IF

  ! Two options: either valid state files are used, or another valid HDF5 file - otherwise, abort
  isValid = ISVALIDHDF5FILE(Args(iArg),FileType=FileType)
  RestartFile=Args(iArg)
  IF(isValid.AND.STRICMP(FileType,'State'))THEN
    IF(STRICMP(ParameterFile,'NOT SET'))THEN
      ParameterFile = ".flexi.ini"
      CALL ExtractParameterFile(Args(iArg), ParameterFile, userblockFound)
      IF (.NOT.userblockFound)&
        CALL CollectiveStop(__STAMP__, "No userblock provided either by user or state file "//TRIM(RestartFile))
    END IF
    stateFileMode = .TRUE.
  ELSE IF (isValid) THEN
    ! Not a state, but a valid HDF5 file (e.g. TimeAvg)
    stateFileMode = .FALSE.
  ELSE
    CALL CollectiveStop(__STAMP__,'ERROR - Not a valid HDF5 file: '//TRIM(RestartFile))
  END IF

! check for command line argument --help or --markdown
ELSE
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
CALL prms%read_options(ParameterFile)

IF (stateFileMode) THEN
  ! The variable names are taken from the equation system
  ALLOCATE(StrVarNames_loc(PP_nVar))
  StrVarNames_loc = StrVarNames
ELSE
  ! We need to specify the name of the dataset to evaluate for non-state files.
  DataSetName = GETSTR('RecordpointsDataSetName','DG_Solution')
  ! Determine the size of the data set
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL GetDataSize(File_ID,TRIM(DataSetName),nDims,HSize)
  HSize_proc = INT(HSize)
  ! Try to get the names of the variables from the attributes of the file
  IF (STRICMP(TRIM(DataSetName),'DG_Solution')) THEN
    CALL GetVarNames("VarNames",StrVarNames_loc,VarNamesExist)
  ELSE
    CALL DatasetExists(File_ID,"VarNames_"//TRIM(DataSetName),varnames_found,attrib=.TRUE.)
    IF (varnames_found) THEN
      CALL GetVarNames("VarNames_"//TRIM(DataSetName),StrVarNames_loc,VarNamesExist)
    END IF
  END IF
  CALL CloseDataFile()
  ! If that was unsuccessfull, simply number the variables
  IF (.NOT.VarNamesExist) THEN
    ALLOCATE(StrVarNames_loc(INT(HSize_proc(1))))
    DO j=1,INT(HSize_proc(1))
      WRITE(StrVarNames_loc(j),'(I0)') j
    END DO
  END IF
END IF

CALL InitIOHDF5()
IF (stateFileMode) THEN
  ! The polynomial degree is read from the parameter file extracted from the user block
  CALL InitInterpolation()
ELSE
  ! The polynomial degree is determined from the size of the input data set
  CALL InitInterpolation(INT(HSize_proc(2))-1)
END IF
prms%removeAfterRead=.TRUE.
CALL InitOutput()
CALL InitMesh(meshMode=0)
#if USE_MPI
CALL InitMPIvars()
#endif
CALL InitRecordpoints()
CALL IgnoredParameters()

IF (stateFileMode) THEN
  ! As the file is a state, we can just allocate the normal U and perform a restart later on
  ALLOCATE(U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
  nVar = PP_nVar
ELSE
  ! This might be a TimeAvg file or similar, allocate depending on the size of the data set
  HSize_proc(5) = nElems
  ALLOCATE(U(HSize_proc(1),0:HSize_proc(2)-1,0:HSize_proc(3)-1,0:HSize_proc(4)-1,nElems))
  nVar = HSize_proc(1)
END IF
#if FV_ENABLED
ALLOCATE(FV_Elems(nElems))
ALLOCATE(IndValue(nElems))
#endif
! used to define RP_Buffersize
dt=1.
WriteData_dt=nArgs*dt
tWriteData=HUGE(1.)

prms%removeAfterRead=.FALSE.
IF (stateFileMode) THEN
  CALL InitRestart(RestartFile)
  CALL FinalizeRestart()
END IF

! make a loop over all files
start=iArg
DO iArg=start,nArgs
  RestartFile=Args(iArg)

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iArg-1,' of ',nArgs-1,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(RestartFile),'" )'

  isValid = ISVALIDHDF5FILE(RestartFile,FileType=FileType)
  IF(.NOT.isValid) &
    CALL CollectiveStop(__STAMP__,'ERROR - Not a valid file: '//TRIM(RestartFile))

  ! Read the file
  IF (stateFileMode) THEN
    CALL InitRestart(RestartFile)
    IF((PP_nVar.NE.NVar_Restart) .OR. (PP_N.NE.N_Restart) .OR. (nGlobalElems.NE.nElems_Restart) .OR. (NodeType.NE.NodeType_Restart)) THEN
      CALL Abort(__STAMP__,'All state files must be identical regarding the number of variables, elements, interpolation points, the mesh file and node type!')
    END IF

    CALL Restart(doFlushFiles=.FALSE.)
  ELSE
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    ! Safety check if the data size stays the same
    CALL GetDataSize(File_ID,TRIM(DataSetName),nDims,HSize)
    HSize_tmp = INT(HSize)
    HSize_tmp(5) = nElems
    IF (.NOT.ALL(HSize_proc.EQ.HSize_tmp)) &
      CALL Abort(__STAMP__,'The size of the data array changed!')
    ! Safety check if the node type is the correct one
    CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType_HDF5)
    IF (TRIM(NodeType_HDF5).NE.TRIM(NodeType)) &
      CALL Abort(__STAMP__,'Wrong node type!')
    CALL ReadArray(TRIM(DataSetName),5,HSize_proc,OffsetElem,5,RealArray=U)
    ! Try to get the time stamp from the file
    CALL DatasetExists(File_ID,"Time",varnames_found,attrib=.TRUE.)
    IF (varnames_found) THEN
      CALL ReadAttribute(File_ID,"Time",1,RealScalar=RestartTime)
    ELSE
      RestartTime = 0.
    END IF
    CALL CloseDataFile()
  END IF

  ! Evaluate the recordpoints
  IF(RP_onProc) CALL RecordPoints(nVar,StrVarNames_loc,0_8,RestartTime,.TRUE.)

  IF (stateFileMode) CALL FinalizeRestart()
END DO
IF(RP_onProc) CALL WriteRP(nVar,StrVarNames_loc,RestartTime,.TRUE.)

DEALLOCATE(U)
#if FV_ENABLED
DEALLOCATE(FV_Elems,IndValue)
#endif

CALL FinalizeMesh()
CALL FinalizeInterpolation()
CALL FinalizeOutput()
CALL FinalizeParameters()
CALL FinalizeCommandlineArguments()

#if USE_MPI
CALL FinalizeMPI()
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' EVALREC TOOL FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM evalrec
