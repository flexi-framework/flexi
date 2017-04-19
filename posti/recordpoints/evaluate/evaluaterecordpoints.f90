#include "flexi.h"

!===================================================================================================================================
!> Tool to extract recordpoint type signals which are normally generated during the simulation, from state files 
!===================================================================================================================================
PROGRAM evalrec
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_ReadInTools,       ONLY:prms,IgnoredParameters,PrintDefaultParameterFile,FinalizeParameters,ExtractParameterFile
USE MOD_StringTools       ,ONLY:STRICMP,GetFileExtension
USE MOD_HDF5_Input,        ONLY:ISVALIDHDF5FILE
!!! These variables are actually modified by the tool itself !!!
USE MOD_MPI,               ONLY:DefineParametersMPI,InitMPI
#if USE_MPI
USE MOD_MPI,               ONLY:InitMPIvars,FinalizeMPI
#endif
USE MOD_IO_HDF5,           ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Output,            ONLY:DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Mesh,              ONLY:DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars,         ONLY:nElems,nGlobalElems
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
CHARACTER(LEN=255)                 :: RPFile,FileType
LOGICAL                            :: isValid,userblockFound
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
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

IF(doPrintHelp.LE.0)THEN
  ! Check if the number of arguments is correct
  IF (nArgs.LT.2) THEN
    ! Print out error message containing valid syntax
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: evalrec [parameter.ini] RPDefFile.h5 statefile1.h5...statefileN.h5 or evalrec --help'// &
    '[option/section name] to print help for a single parameter, parameter sections or all parameters.')
  END IF
  
  iArg=1
  ParameterFile='NOT SET'
  IF(STRICMP(GetFileExtension(Args(iArg)), "ini"))THEN
    ParameterFile = Args(iArg)
    iArg=iArg+1
  END IF
  
  RPFile=Args(iArg)
  isValid = ISVALIDHDF5FILE(RPFile,FileType=FileType)
  IF(isValid.AND.STRICMP(FileType,'RecordPoints'))THEN
    iArg=iArg+1
  ELSE
    CALL CollectiveStop(__STAMP__,'ERROR: No record point definition file provided.')
  END IF
  
  isValid = ISVALIDHDF5FILE(Args(iArg),FileType=FileType)
  IF(isValid.AND.STRICMP(FileType,'State'))THEN
    RestartFile=Args(iArg)
    IF(STRICMP(ParameterFile,'NOT SET'))THEN
      ParameterFile = ".flexi.ini" 
      CALL ExtractParameterFile(Args(iArg), ParameterFile, userblockFound)
      IF (.NOT.userblockFound)&
        CALL CollectiveStop(__STAMP__, "No userblock provided either by user or state file "//TRIM(RestartFile))
    END IF
  ELSE
    CALL CollectiveStop(__STAMP__,'ERROR - Not a valid state file: '//TRIM(RestartFile))
  END IF
END IF ! doPrintHelp.LE.0

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
CALL prms%read_options(ParameterFile)

CALL InitIOHDF5()
CALL InitInterpolation()
prms%removeAfterRead=.TRUE.
CALL InitOutput()
CALL InitMesh(meshMode=0)
#if USE_MPI
CALL InitMPIvars()
#endif
CALL InitRecordpoints(RPFile)
CALL IgnoredParameters()

ALLOCATE(U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
#if FV_ENABLED
ALLOCATE(FV_Elems(nElems))
ALLOCATE(IndValue(nElems))
#endif
! used to define RP_Buffersize
dt=1.
WriteData_dt=nArgs*dt
tWriteData=HUGE(1.)

prms%removeAfterRead=.FALSE.
CALL InitRestart(RestartFile)
CALL FinalizeRestart()

! make a loop over all files
start=iArg
DO iArg=start,nArgs
  RestartFile=Args(iArg)

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iArg-1,' of ',nArgs-1,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(RestartFile),'" )'

  isValid = ISVALIDHDF5FILE(RestartFile,FileType=FileType)
  isValid = STRICMP(FileType, "State")
  IF(.NOT.isValid) &
    CALL CollectiveStop(__STAMP__,'ERROR - Not a valid state file: '//TRIM(RestartFile))

  ! Read the file
  CALL InitRestart(RestartFile)
  IF((PP_nVar.NE.NVar_Restart) .OR. (PP_N.NE.N_Restart) .OR. (nGlobalElems.NE.nElems_Restart) .OR. (NodeType.NE.NodeType_Restart)) THEN
    CALL Abort(__STAMP__,'All state files must be identical regarding the number of variables, elements, interpolation points, the mesh file and node type!')
  END IF

  CALL Restart(doFlushFiles=.FALSE.)

  ! Evaluate the recordpoints
  IF(RP_onProc) CALL RecordPoints(0_8,RestartTime,.TRUE.)

  CALL FinalizeRestart()
END DO
IF(RP_onProc) CALL WriteRP(RestartTime,.TRUE.)

DEALLOCATE(U)
#if FV_ENABLED
DEALLOCATE(FV_Elems,IndValue)
#endif

CALL FinalizeMesh()
CALL FinalizeInterpolation()
CALL FinalizeOutput()
CALL FinalizeParameters()
CALL FinalizeCommandlineArguments()

#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
  CALL FinalizeMPI()
#endif

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' EVALREC TOOL FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM evalrec
