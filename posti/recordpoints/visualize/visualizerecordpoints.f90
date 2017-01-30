#include "flexi.h"

!===================================================================================================================================
!> Add comments please!
!===================================================================================================================================
PROGRAM postrec
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_StringTools                 ,ONLY:STRICMP, GetFileExtension
USE MOD_ReadInTools                 ,ONLY: prms
USE MOD_Parameters                  ,ONLY:equiTimeSpacing,doSpec,doFluctuations,doTurb,doFilter,Plane_doBLProps
USE MOD_RPSet                       ,ONLY:InitRPSet,FinalizeRPSet
USE MOD_RPData                      ,ONLY:ReadRPData,AssembleRPData,FinalizeRPData
USE MOD_OutputRPVisu                      
USE MOD_RPInterpolation
USE MOD_RPInterpolation_Vars        ,ONLY:CalcTimeAverage
USE MOD_EquationRP                
USE MOD_FilterRP                    ,ONLY:FilterRP
USE MOD_Spec                        ,ONLY:InitSpec,spec,FinalizeSpec
USE MOD_Turbulence
USE MOD_MPI                         ,ONLY:InitMPI
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg,iExt,nDataFiles
LOGICAL                            :: success=.TRUE.
CHARACTER(LEN=255)                 :: InputIniFile
CHARACTER(LEN=255)                 :: InputDataFile
CHARACTER(LEN=255),ALLOCATABLE     :: DataFiles(:)
!===================================================================================================================================
CALL InitMPI() ! NO PARALLELIZATION, ONLY FOR COMPILING WITH MPI FLAGS ON SOME MACHINES OR USING MPI-DEPENDANT HDF5
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
  'This tool is designed for single execution only!')
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)') " || Postprocessing for Flexi Recordpoints||"
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)')

CALL ParseCommandlineArguments()
IF(nArgs .LT. 2) success=.FALSE.

IF(success)THEN
  ParameterFile = Args(1)
  IF (.NOT.STRICMP(GetFileExtension(ParameterFile), "ini")) success=.FALSE.
END IF
IF(.NOT.success)THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: postrec parameter.ini RPdatafiles.h5')
END IF

CALL DefineParameters()

CALL prms%read_options(ParameterFile)

CALL InitParameters()

nDataFiles=nArgs-1
ALLOCATE(DataFiles(1:nDataFiles))

! get list of input files
DO iArg=2,nArgs
  CALL GET_COMMAND_ARGUMENT(iArg,DataFiles(iArg-1))
END DO

! readin RP Data from all input files
DO iArg=1,nDataFiles
  InputDataFile=DataFiles(iArg)
  WRITE(UNIT_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iArg,' of ',nDataFiles,' FILES.'
  WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputDataFile),'" )'
  WRITE(UNIT_stdOut,'(132("="))')

  ! Get start index of file extension to check if it is a h5 file
  iExt=INDEX(InputDataFile,'.',BACK = .TRUE.)
  IF(InputDataFile(iExt+1:iExt+2) .NE. 'h5') &
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid file extension!')
  ! Read in main attributes from given HDF5 State File
  WRITE(UNIT_stdOUT,*) "READING DATA FROM RP FILE """,TRIM(InputDataFile), """"
  IF(iArg.EQ.1) THEN
    CALL ReadRPData(InputDataFile,firstFile=.TRUE.)
  ELSE 
    CALL ReadRPData(InputDataFile)
  END IF
END DO

! assemble data to one global array
CALL AssembleRPData()

CALL InitEquationRP()
CALL InitInterpolation()
IF(doSpec)             CALL InitSpec()
CALL InitOutput()
IF(equiTimeSpacing)    CALL InterpolateEquiTime()
CALL CalcEquationRP()
IF(calcTimeAverage)    CALL CalcTimeAvg() 
IF(doFluctuations)     CALL CalcFluctuations()
IF(doFilter)           CALL FilterRP()
IF(doSpec)             CALL spec()
#ifdef WITHBLPROPS
IF(Plane_doBLProps)    CALL Plane_BLProps() 
#endif
CALL OutputRP()
IF(doTurb)             CALL Turbulence()
CALL FinalizeInterpolation()
CALL FinalizeEquationRP()
CALL FinalizeOutput()
CALL FinalizeRPSet()
CALL FinalizeRPData()
CALL FinalizeSpec()
CALL FinalizeCommandlineArguments()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL CollectiveStop(__STAMP__,'MPI finalize error',iError)
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' RECORDPOINTS POSTPROC FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM postrec


!===================================================================================================================================
!> Initialize parameter variables of Posti tool
!===================================================================================================================================
SUBROUTINE DefineParameters()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!===================================================================================================================================
CALL prms%SetSection('Visualize_Record_Points')

CALL prms%CreateStringOption( 'ProjectName'        ,"TODO")
CALL prms%CreateStringOption( 'GroupName'          ,"TODO",multiple=.TRUE.)
CALL prms%CreateStringOption( 'VarName'            ,"TODO",multiple=.TRUE.)
CALL prms%CreateStringOption( 'RP_DefFile'         ,"TODO")

CALL prms%CreateLogicalOption('usePrims'           ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputTimeData'     ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputTimeAverage'  ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doFluctuations'     ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doFilter'           ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doFFT'              ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doPSD'              ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doTurb'             ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('hanning'            ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('fourthDeriv'        ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('ThirdOct'           ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('equiTimeSpacing'    ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Plane_LocalCoords'  ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Plane_LocalVel'     ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Plane_doBLProps'    ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Line_LocalCoords'   ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Line_LocalVel'      ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputPoints'       ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputLines'        ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputPlanes'       ,"TODO",".FALSE.")

CALL prms%CreateRealArrayOption('Line_LocalVel_vec',"TODO",".FALSE.")
CALL prms%CreateRealOption   ('FilterWidth'        ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('SamplingFreq'       ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('u_inf'              ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('chord'              ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('mu0'                ,"TODO",".FALSE.") 

CALL prms%CreateIntOption    ('FilterMode'         ,"TODO",".FALSE.")
CALL prms%CreateIntOption    ('nBlocks'            ,"TODO",".FALSE.")
CALL prms%CreateIntOption    ('BlockSize'          ,"TODO",".FALSE.")
CALL prms%CreateIntOption    ('Plane_BLvelScaling' ,"TODO",".FALSE.")
CALL prms%CreateIntOption    ('SkipSample'         ,"TODO",".FALSE.")
CALL prms%CreateIntOption    ('OutputFormat'       ,"TODO",".FALSE.")
END SUBROUTINE DefineParameters

!===================================================================================================================================
!> Read parameters of Posti tool
!===================================================================================================================================
SUBROUTINE InitParameters()
! MODULES
USE MOD_Globals
USE MOD_Readintools         ,ONLY:GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY,CountOption
USE MOD_Parameters
USE MOD_RPInterpolation_Vars,ONLY:calcTimeAverage
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: iGroup,iVar
!===================================================================================================================================
Projectname=GETSTR('ProjectName')

! =============================================================================== !
! RP INFO
! =============================================================================== !
nGroups_visu=CountOption('GroupName')
ALLOCATE(GroupNames_visu(nGroups_visu))
DO iGroup=1,nGroups_visu
 GroupNames_visu(iGroup)=GETSTR('Groupname','')
END DO
RP_SET_defined=.FALSE.
RP_DefFile=GETSTR('RP_DefFile','none')
IF(TRIM(RP_defFile).NE.'none') THEN
  RP_SET_defined=.TRUE.
END IF

! use primitive variables for derived quantities if they exist in the state file
usePrims=GETLOGICAL('usePrims','.FALSE.')

! =============================================================================== !
! TIME INTERVAL
! =============================================================================== !
OutputTimeData   =GETLOGICAL('OutputTimeData','.FALSE.')
OutputTimeAverage=GETLOGICAL('OutputTimeAverage','.FALSE.')
doFluctuations   =GETLOGICAL('doFluctuations','.FALSE.')
IF(OutputTimeAverage .OR. doFluctuations) CalcTimeAverage=.TRUE.

doFilter=GETLOGICAL('doFilter','.FALSE.')
IF(doFilter) THEN
  FilterWidth=GETREAL('FilterWidth')
  FilterMode=GETINT('FilterMode','0')
END IF
! =============================================================================== !
! FOURIER TRANSFORM
! =============================================================================== !
doFFT=GETLOGICAL('doFFT','.FALSE.')
doPSD=GETLOGICAL('doPSD','.FALSE.')
IF(doFFT.OR.doPSD) doSpec=.TRUE.
IF(doSpec) THEN
  ! two readin "modes" for spectrum averaging:
  ! 1. Prescription of number of blocks
  nBlocks=GETINT('nBlocks','1')
  ! 2. Prescription of Sampling Frequency and Blocksize
  samplingFreq=GETREAL('SamplingFreq','-999')
  IF(samplingFreq.GT.0.) THEN
    BlockSize=GETINT('BlockSize') 
  END IF
  doHanning=GETLOGICAL('hanning','.FALSE.')
  fourthDeriv=GETLOGICAL('fourthDeriv','.FALSE.')
  ThirdOct=GETLOGICAL('ThirdOct','.FALSE.')
  IF (ThirdOct) THEN
    u_infPhys   = GETREAL('u_inf') !velocity for re-dimensionalization of frequency
    chordPhys = GETREAL('chord')   !length for re-dimensionalization of frequency
  END IF
END IF

doTurb=GETLOGICAL('doTurb','.FALSE.')

IF(doSpec .OR. doTurb) cutoffFreq=GETREAL('cutoffFreq','-999.')

! for any FFT operation we need equidistant time spacing
equiTimeSpacing=.FALSE.
IF(OutputTimeData) equiTimeSpacing=GETLOGICAL('equiTimeSpacing','.FALSE.')
IF(doTurb.OR.doSpec) equiTimeSpacing=.TRUE.

! =============================================================================== !
! PLANE OPTIONS 
! =============================================================================== !
OutputPlanes      =GETLOGICAL('OutputPlanes','.TRUE.')
Plane_LocalCoords =GETLOGICAL('Plane_LocalCoords','.FALSE.')
Plane_LocalVel    =GETLOGICAL('Plane_LocalVel','.FALSE.')
Plane_doBLProps   =GETLOGICAL('Plane_doBLProps','.FALSE.')
IF(Plane_doBLProps) THEN ! for BL properties we need local coords and velocities
  WRITE(UNIT_StdOut,'(A)')' BL properties depend on local velocities and coordinates'
  WRITE(UNIT_StdOut,'(A)')' and are calculated based on time-averaged data.'
  WRITE(UNIT_StdOut,'(A)')' Setting Plane_localCoords=.TRUE. and Plane_localVel=.TRUE..'
  CalcTimeAverage  =.TRUE.
  OutputTimeAverage=.TRUE.
  Plane_LocalCoords=.TRUE.
  Plane_LocalVel   =.TRUE.
  Plane_BLvelScaling  =GETINT('Plane_BLvelScaling','0') ! 0 - no scaling. 
  ! 1 - "laminar scaling": scale velocity with u_delta and PlaneY with delta99
  ! 2 - "turbulent scaling:" calculate u+ and y+
END IF
! =============================================================================== !
! LINE OPTIONS 
! =============================================================================== !
OutputLines      =GETLOGICAL('OutputLines',     '.TRUE.')
Line_LocalCoords =GETLOGICAL('Line_LocalCoords','.FALSE.')
Line_LocalVel    =GETLOGICAL('Line_LocalVel',   '.FALSE.')
IF(Line_LocalVel) THEN 
  Line_LocalVel_vec=GETREALARRAY('Line_LocalVel_vec',3)
  ! normalize the vector
  Line_LocalVel_vec = Line_LocalVel_vec/NORM2(Line_LocalVel_vec)
END IF
! =============================================================================== !
! POINT OPTIONS 
! =============================================================================== !
OutputPoints     =GETLOGICAL('OutputPoints','.TRUE.')

skip = GETINT('SkipSample','1')

nVarVisu=CountOption('VarName')
ALLOCATE(VarNameVisu(nVarVisu))
DO iVar=1,nVarVisu
 VarNameVisu(iVar)=GETSTR('VarName')
END DO

IF(doTurb.OR.Plane_doBLProps) Mu0=GETREAL('Mu0')

OutputFormat = GETINT('OutputFormat','0')

! when no time averaging and no fft shall be performed, set the output signal to true
IF(.NOT.ANY( (/doSpec, OutputTimeAverage , doTurb /)) ) OutputTimeData=.TRUE.
END SUBROUTINE InitParameters


!===================================================================================================================================
!> This routine builds the mappings from the total number of variables available for visualization to number of calculation
!> and visualization variables.
!>  1. Read 'VarName' options from the parameter file. This are the quantities that will be visualized.
!>  2. Initialize the dependecy table
!>  3. check wether gradients are needed for any quantity. If this is the case, remove the conservative quantities from the 
!>     dependecies of the primitive quantities (the primitive quantities are available directly, since the DGTimeDerivative_weakForm
!>     will be executed.
!>  4. build the 'mapCalc' that holds for each quantity that will be calculated the index in 'UCalc' array (0 if not calculated)
!>  5. build the 'mapVisu' that holds for each quantity that will be visualized the index in 'UVisu' array (0 if not visualized)
!===================================================================================================================================
SUBROUTINE Build_mapCalc_mapVisu()
USE MOD_Parameters
USE MOD_ReadInTools     ,ONLY: GETSTR,CountOption
USE MOD_StringTools     ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLESIABLES
INTEGER             :: iVar,iVar2
CHARACTER(LEN=20)   :: format
!===================================================================================================================================
! Read Varnames from parameter file and fill
!   mapVisu = map, which stores at position x the position/index of the x.th quantity in the UVisu array
!             if a quantity is not visualized it is zero
ALLOCATE(mapVisu(1:nVarDep))
mapVisu = 0
nVarVisuTotal = 0
! Compare varnames that should be visualized with availabe varnames
DO iVar=1,nVarVisu
  DO iVar2=1,nVarDep
    IF (STRICMP(VarNameVisu(iVar), VarNamesAll(iVar2))) THEN
      mapVisu(iVar2) = nVarVisuTotal+1
      nVarVisuTotal = nVarVisuTotal + 1
    END IF
  END DO
END DO

! Calculate all dependencies:
! For each quantity copy from all quantities that this quantity depends on the dependencies.
DO iVar=1,nVarDep
  DepTable(iVar,iVar) = 1
  DO iVar2=1,iVar-1
    IF (DepTable(iVar,iVar2).EQ.1) &
      DepTable(iVar,:) = MAX(DepTable(iVar,:), DepTable(iVar2,:))
  END DO
END DO

! Build :
!   mapCalc = map, which stores at position x the position/index of the x.th quantity in the UCalc array
!             if a quantity is not calculated it is zero
ALLOCATE(mapCalc(1:nVarDep))
mapCalc = 0
DO iVar=1,nVarDep
  IF (mapVisu(iVar).GT.0) THEN
    mapCalc = MAX(mapCalc,DepTable(iVar,1:nVarDep))
  END IF
END DO
! enumerate mapCalc
nVarCalc = 0
DO iVar=1,nVarDep
  IF (mapCalc(iVar).GT.0) THEN
    nVarCalc = nVarCalc + 1
    mapCalc(iVar) = nVarCalc
  END IF
END DO

! print the dependecy table
WRITE(format,'(I2)') SIZE(DepTable,2)
DO iVar=1,nVarDep
  WRITE (*,'('//format//'I2,A)') DepTable(iVar,:), " "//TRIM(VarNamesAll(iVar))
END DO

! print the mappings
WRITE(format,'(I2)') nVarDep
WRITE (*,'(A,'//format//'I3)') "mapCalc ",mapCalc
WRITE (*,'(A,'//format//'I3)') "mapVisu ",mapVisu

END SUBROUTINE Build_mapCalc_mapVisu

