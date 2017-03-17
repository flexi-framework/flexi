#include "flexi.h"

!===================================================================================================================================
!> Add comments please!
!===================================================================================================================================
PROGRAM PrepareRecordPoints
! MODULES
! General and support
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_StringTools,        ONLY:STRICMP,GetFileExtension
USE MOD_ReadInTools,        ONLY:prms,IgnoredParameters,PrintDefaultParameterFile,FinalizeParameters
! Flexilib initialization
USE MOD_MPI,                ONLY:DefineParametersMPI,InitMPI
USE MOD_IO_HDF5,            ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Interpolation,      ONLY:DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Output,             ONLY:DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Output_Vars,        ONLY:ProjectName
USE MOD_Mesh,               ONLY:DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars,          ONLY:MeshFile
USE MOD_Mortar,             ONLY:InitMortar,FinalizeMortar
! Recordpoints
USE MOD_RPSet
USE MOD_VisuRP,             ONLY:VisuRP
USE MOD_RPParametricCoords, ONLY:GetRecordPoints
USE MOD_HDF5_OutputRP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                            :: success=.TRUE.
!===================================================================================================================================
CALL InitMPI() ! NO PARALLELIZATION, ONLY FOR COMPILING WITH MPI FLAGS ON SOME MACHINES OR USING MPI-DEPENDANT HDF5
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
  'This tool is designed for single execution only!')
WRITE(UNIT_stdOut,'(A)') " ||=================================||" 
WRITE(UNIT_stdOut,'(A)') " || Prepare record points for Flexi ||"
WRITE(UNIT_stdOut,'(A)') " ||=================================||"

CALL ParseCommandlineArguments()
IF(nArgs .NE. 1) success=.FALSE.

IF(success)THEN
  ParameterFile = Args(1)
  IF (.NOT.STRICMP(GetFileExtension(ParameterFile), "ini")) success=.FALSE.
END IF
IF(.NOT.success.AND.doPrintHelp.LE.0)THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: preparerec parameter.ini or preparerec --help'// &
  '[option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF

CALL DefineParameters()
CALL DefineParametersRPSet()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersInterpolation()
CALL DefineParametersOutput()
CALL DefineParametersMesh()

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
CALL prms%read_options(ParameterFile)

CALL InitParameters()
CALL InitIOHDF5()
CALL InitInterpolation()
CALL InitMortar()
CALL InitOutput()
CALL InitMesh(meshMode=2)

CALL InitRPSet()
CALL GetRecordPoints()
CALL WriteRecordPointstoHDF5(ProjectName,MeshFile)
CALL visuRP()

CALL FinalizeOutput()
CALL FinalizeInterpolation()
CALL FinalizeMesh()
CALL FinalizeMortar()
CALL FinalizeParameters()
CALL FinalizeCommandlineArguments()

#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif

WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' PREPARE RECORDPOINTS TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM PrepareRecordPoints


!===================================================================================================================================
!> Initialize parameter variables of Posti tool
!===================================================================================================================================
SUBROUTINE DefineParameters()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!===================================================================================================================================
CALL prms%SetSection('Prepare Record Points')
CALL prms%CreateIntOption( 'NSuper',"Number of Newton start values per element per direction.",'0')
CALL prms%CreateRealOption('maxTolerance',"Tolerance in parameter space at the element "//&
                           "boundaries, required to mark a recordpoint as found.",'1.E-3')
CALL prms%CreateLogicalOption('doVisuRP',"Visualize recordpoints.",".TRUE.")
END SUBROUTINE DefineParameters


!===================================================================================================================================
!> Read parameters of Posti tool
!===================================================================================================================================
SUBROUTINE InitParameters()
! MODULES
USE MOD_Parameters
USE MOD_Readintools   ,ONLY:GETINT,GETREAL,GETLOGICAL
IMPLICIT NONE
!===================================================================================================================================
NSuper     =GETINT('NSuper','0')
maxTol     =1.+ABS(GETREAL('maxTolerance','1.e-3'))
doVisuRP   =GETLOGICAL('doVisuRP','T')
END SUBROUTINE InitParameters

