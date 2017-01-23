#include "flexi.h"

!===================================================================================================================================
!> Add comments please!
!===================================================================================================================================
PROGRAM postrec
! MODULES
USE MOD_Globals
USE MOD_Parameters                  ,ONLY:equiTimeSpacing,doSpec,doFluctuations,doTurb,doFilter,Plane_doBLProps
USE MOD_Parameters                  ,ONLY:ProjectName
USE MOD_RPSet                       ,ONLY:InitRPSet,FinalizeRPSet
USE MOD_RPData                      ,ONLY:ReadRPData,AssembleRPData,FinalizeRPData
USE MOD_OutputRPVisu                      
USE MOD_RPInterpolation
USE MOD_RPInterpolation_Vars        ,ONLY:CalcTimeAverage
USE MOD_Equation                
USE MOD_FilterRP                    ,ONLY:FilterRP
USE MOD_Spec                        ,ONLY:InitSpec,spec,FinalizeSpec
USE MOD_Turbulence
#ifdef MPI
USE MOD_MPI                         ,ONLY:InitMPI
#endif /* MPI */
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg,nArgs,iExt,nDataFiles
CHARACTER(LEN=255)                 :: InputIniFile
CHARACTER(LEN=255)                 :: InputDataFile
CHARACTER(LEN=255),ALLOCATABLE     :: DataFiles(:)
!===================================================================================================================================
#ifdef MPI
CALL InitMPI() ! NO PARALLELIZATION, ONLY FOR COMPILING WITH MPI FLAGS ON SOME MACHINES OR USING MPI-DEPENDANT HDF5
#endif /* MPI */
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)') " || Postprocessing for Flexi Recordpoints||"
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)')


nArgs=COMMAND_ARGUMENT_COUNT()
IF(nArgs .LT. 2) CALL Abort(__STAMP__,'Missing argument')
CALL GET_COMMAND_ARGUMENT(1,InputIniFile)
! Get start index of file extension
iExt=INDEX(InputIniFile,'.',BACK = .TRUE.)
! check if first file is a .ini file
IF(InputIniFile(iExt+1:iExt+3) .NE. 'ini') &
  CALL Abort(__STAMP__,'ERROR - No / invalid parameter file given.')
CALL InitParameters()
!CALL InitRPSet(RP_DefFile)

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
    CALL Abort(__STAMP__,'ERROR - Invalid file extension!')
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

CALL InitEquation()
CALL InitInterpolation()
IF(doSpec)             CALL InitSpec()
CALL InitOutput()
IF(equiTimeSpacing)    CALL InterpolateEquiTime()
CALL CalcEquation()
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
CALL FinalizeEquation()
CALL FinalizeOutput()
CALL FinalizeRPSet()
CALL FinalizeRPData()
CALL FinalizeSpec()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
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
CALL prms%SetSection('Visualize Record Points')

CALL prms%CreateStringOption( 'GroupName'          ,"TODO",multiple=.TRUE.)
CALL prms%CreateStringOption( 'VarName'            ,"TODO",multiple=.TRUE.)
CALL prms%CreateStringOption( 'RP_DefFile'         ,"TODO")

CALL prms%CreateLogicalOption('usePrims'           ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputTimeData'     ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputTimeAverage'  ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doFluctuations'     ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doFilter'           ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doFFT'              ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('doTurb'             ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('hanning'            ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('fourthDeriv'        ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('ThirdOct'           ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('equiTimeSpacing'    ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputPlanes'       ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Plane_LocalCoords'  ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Plane_LocalVel'     ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Plane_doBLProps'    ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputLines'        ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Line_LocalCoords'   ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('Line_LocalVel'      ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption('OutputPlanes'       ,"TODO",".FALSE.")

CALL prms%CreateRealArrayOption('Line_LocalVel_vec',"TODO",".FALSE.")
CALL prms%CreateRealOption   ('FilterWidth'        ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('SamplingFreq'       ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('u_inf'              ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('chord'              ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('mu0'                ,"TODO",".FALSE.") ! probably not needed

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

Projectname=GETSTR('ProjectName','')
! =============================================================================== !
! RP INFO
! =============================================================================== !
nGroups_visu=CountOption('GroupName')
ALLOCATE(GroupNames_visu(nGroups_visu))
DO iGroup=1,nGroups_visu
 GroupNames_visu(iGroup)=GETSTR('Groupname','')
END DO
RP_SET_defined=.FALSE.
RP_DefFile=GETSTR('RP_DefFile','')
IF(TRIM(RP_defFile).NE.'') THEN
  RP_SET_defined=.TRUE.
END IF

! use primitive variables for derived quantities if they exist in the state file
usePrims=GETLOGICAL('usePrims','.FALSE.')

! =============================================================================== !
! TIME INTERVAL
! =============================================================================== !
OutputTimeData   =GETLOGICAL('OutputTimeData','F')
OutputTimeAverage=GETLOGICAL('OutputTimeAverage','F')
doFluctuations   =GETLOGICAL('doFluctuations','F')
IF(OutputTimeAverage .OR. doFluctuations) CalcTimeAverage=.TRUE.

doFilter=GETLOGICAL('doFilter','F')
IF(doFilter) THEN
  FilterWidth=GETREAL('FilterWidth')
  FilterMode=GETINT('FilterMode','0')
END IF
! =============================================================================== !
! FOURIER TRANSFORM
! =============================================================================== !
doFFT=GETLOGICAL('doFFT','F')
doPSD=GETLOGICAL('doPSD','F')
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
  doHanning=GETLOGICAL('hanning','F')
  fourthDeriv=GETLOGICAL('fourthDeriv','F')
  ThirdOct=GETLOGICAL('ThirdOct','F')
  IF (ThirdOct) THEN
    u_infPhys   = GETREAL('u_inf') !velocity for re-dimensionalization of frequency
    chordPhys = GETREAL('chord')   !length for re-dimensionalization of frequency
  END IF
END IF

doTurb=GETLOGICAL('doTurb','F')

IF(doSpec .OR. doTurb) cutoffFreq=GETREAL('cutoffFreq','-999.')

! for any FFT operation we need equidistant time spacing
equiTimeSpacing=.FALSE.
IF(OutputTimeData) equiTimeSpacing=GETLOGICAL('equiTimeSpacing','F')
IF(doTurb.OR.doSpec) equiTimeSpacing=.TRUE.

! =============================================================================== !
! PLANE OPTIONS 
! =============================================================================== !
OutputPlanes      =GETLOGICAL('OutputPlanes','T')
Plane_LocalCoords =GETLOGICAL('Plane_LocalCoords','F')
Plane_LocalVel    =GETLOGICAL('Plane_LocalVel','F')
Plane_doBLProps   =GETLOGICAL('Plane_doBLProps','F')
IF(Plane_doBLProps) THEN ! for BL properties we need local coords and velocities
  WRITE(UNIT_StdOut,'(A)')' BL properties depend on local velocities and coordinates'
  WRITE(UNIT_StdOut,'(A)')' and are calculated based on time-averaged data.'
  WRITE(UNIT_StdOut,'(A)')' Setting Plane_localCoords=T and Plane_localVel=T.'
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
OutputLines      =GETLOGICAL('OutputLines',     'T')
Line_LocalCoords =GETLOGICAL('Line_LocalCoords','F')
Line_LocalVel    =GETLOGICAL('Line_LocalVel',   'F')
IF(Line_LocalVel) THEN 
  Line_LocalVel_vec=GETREALARRAY('Line_LocalVel_vec',3)
  ! normalize the vector
  Line_LocalVel_vec = Line_LocalVel_vec/NORM2(Line_LocalVel_vec)
END IF
! =============================================================================== !
! POINT OPTIONS 
! =============================================================================== !
OutputPoints     =GETLOGICAL('OutputPoints','T')

skip = GETINT('SkipSample','1')

nVar_visu=CountOption('VarName')
ALLOCATE(VarNameVisu(nVar_visu))
DO iVar=1,nVar_visu
 VarNameVisu(iVar)=GETSTR('VarName')
END DO

IF(doTurb.OR.Plane_doBLProps) Mu0=GETREAL('Mu0')

OutputFormat = GETINT('OutputFormat','0')

! when no time averaging and no fft shall be performed, set the output signal to true
IF(.NOT.ANY( (/doSpec, OutputTimeAverage , doTurb /)) ) OutputTimeData=.TRUE.
END SUBROUTINE InitParameters
