!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
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
!> Tool used to create record points. The points are sorted in groups and can take different shapes, e.g. simple lines, circles,
!> planes etc. For all points, the reference coordinates and respective element IDs are found in the mesh and stored as a .h5 file
!> which can later be read by FLEXI or POSTI tools.
!> In this way, we can e.g. produce a time-series of flow variables at the RPs in a very fine time intervall.
!===================================================================================================================================
PROGRAM PrepareRecordPoints
! MODULES
! General and support
USE MOD_Globals
USE MOD_PreProc
USE MOD_Commandline_Arguments
USE MOD_StringTools,        ONLY:STRICMP,GetFileExtension
USE MOD_ReadInTools,        ONLY:prms,IgnoredParameters,PrintDefaultParameterFile,FinalizeParameters,GETSTR
! Flexilib initialization
USE MOD_MPI,                ONLY:DefineParametersMPI,InitMPI
USE MOD_IO_HDF5,            ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_HDF5_Input
USE MOD_Interpolation,      ONLY:InitInterpolation,FinalizeInterpolation
USE MOD_Output,             ONLY:DefineParametersOutput,InitOutput,FinalizeOutput
USE MOD_Output_Vars,        ONLY:ProjectName
USE MOD_Mesh,               ONLY:DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars,          ONLY:MeshFile
USE MOD_Mortar,             ONLY:InitMortar,FinalizeMortar
#if FV_ENABLED
USE MOD_FV_Basis,           ONLY:InitFV_Basis,FinalizeFV_Basis
#endif
! Recordpoints
USE MOD_RPSet
USE MOD_VisuRP,             ONLY:VisuRP
USE MOD_RPParametricCoords, ONLY:GetRecordPoints
USE MOD_HDF5_OutputRP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                            :: success=.TRUE.
INTEGER                            :: Ntmp
!===================================================================================================================================
#if PP_dim ==2
STOP 'Please compile with 3D to use the recordpoints tool!'
#endif
CALL SetStackSizeUnlimited()
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
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti_preparerecordpoints parameter.ini or&
   & posti_preparerecordpoints --help [option/section name] to print help for a single parameter, parameter sections&
   & or all parameters.')
END IF

CALL DefineParameters()
CALL DefineParametersRPSet()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
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
! Open mesh and read NGeo. The polynomial basis will be initialized with 2*NGeo.
MeshFile = GETSTR('MeshFile')
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=Ntmp)
CALL CloseDataFile()
Ntmp = Ntmp*2
CALL InitInterpolation(Ntmp)
#if FV_ENABLED
CALL InitFV_Basis()
#endif
CALL InitMortar()
CALL InitOutput()
CALL InitMesh(meshMode=2,MeshFile_IN=MeshFile)

CALL InitRPSet()
CALL GetRecordPoints()
CALL WriteRecordPointstoHDF5(ProjectName,MeshFile)
CALL visuRP()

CALL FinalizeOutput()
CALL FinalizeInterpolation()
CALL FinalizeMesh()
CALL FinalizeMortar()
#if FV_ENABLED
CALL FinalizeFV_Basis()
#endif
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

