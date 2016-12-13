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
!==================================================================================================================================
!> Module providing routines for reading Flexi parameter files. 
!> 
!> The whole structure to read options from the parameter file is as follows:
!> 
!> All the options are stored in a linked list, which is defined as a class and has a single global instance'prms'.
!> 
!> The options are appended to this list via the DefineParametersXXX() routines, which exist for all the modules that 
!> have an option, which can be specified in the parameter file. This is done at the beginning of the execution. After calling
!> all DefineParametersXXX() routines the prms list contains all possible options (with name, description, default value (optional)).
!> 
!> After that the prms\%read_options() routine is called, which actually reads the options from the parameter file. Therefore the
!> parameter file is read line by line and each line is parsed for an option.
!> By this the values of the options, that are already in the linked list 'prms' are set.
!> 
!> Now all the options are filled with the data from the parameter file and can be accessed via the functions GETINT(ARRAY), 
!> GETREAL(ARRAY), ...
!> A call of these functions then removes the specific option from the linked list, such that 
!> every option can only be read once. This is necessary for options with the same name, that occure multiple times in the parameter
!> file.
!==================================================================================================================================
MODULE MOD_ReadInTools

USE MOD_Globals
USE MOD_ISO_VARYING_STRING
USE MOD_Options
USE MOD_StringTools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

!================================================
!> Link for linked List 
!================================================
TYPE,PUBLIC :: LINK
  CLASS(OPTION), POINTER :: opt => null()
  CLASS(LINK), POINTER   :: next => null()
END TYPE LINK

!==================================================================================================================================
!> Class to store all options.
!> This is basically a linked list of options.
!==================================================================================================================================
TYPE,PUBLIC :: Parameters
  CLASS(LINK), POINTER :: firstLink => null() !< first option in the list
  CLASS(LINK), POINTER :: lastLink  => null() !< last option in the list
  INTEGER              :: maxNameLen          !< maximal string length of the name of an option in the list
  INTEGER              :: maxValueLen         !< maximal string length of the value of an option in the list
  CHARACTER(LEN=255)   :: actualSection = ""  !< actual section, to set section of an option, when inserted into list
CONTAINS 
  PROCEDURE :: SetSection                 !< routine to set 'actualSection'
  PROCEDURE :: CreateOption               !< general routine to create a option and insert it into the linked list
  PROCEDURE :: CreateIntOption            !< routine to generate an integer option
  PROCEDURE :: CreateIntFromStringOption  !< routine to generate an integer option with a optional string representation
  PROCEDURE :: CreateLogicalOption        !< routine to generate an logical option 
  PROCEDURE :: CreateRealOption           !< routine to generate an real option
  PROCEDURE :: CreateStringOption         !< routine to generate an string option
  PROCEDURE :: CreateIntArrayOption       !< routine to generate an integer array option
  PROCEDURE :: CreateLogicalArrayOption   !< routine to generate an logical array option 
  PROCEDURE :: CreateRealArrayOption      !< routine to generate an real array option
  !PROCEDURE :: CreateStringArrayOption    !< routine to generate an string array option
  PROCEDURE :: CountOption_               !< function to count the number of options of a given name
  PROCEDURE :: read_options               !< routine that loops over the lines of a parameter files 
                                          !< and calls read_option for every option. Outputs all unknow options
  PROCEDURE :: read_option                !< routine that parses a single line from the parameter file.
END TYPE Parameters

INTERFACE IgnoredParameters
  MODULE PROCEDURE IgnoredParameters
END INTERFACE

INTERFACE PrintDefaultParameterFile
  MODULE PROCEDURE PrintDefaultParameterFile
END INTERFACE

INTERFACE CountOption
  MODULE PROCEDURE CountOption
END INTERFACE

INTERFACE GETINT
  MODULE PROCEDURE GETINT
END INTERFACE

INTERFACE GETLOGICAL
  MODULE PROCEDURE GETLOGICAL
END INTERFACE

INTERFACE GETREAL
  MODULE PROCEDURE GETREAL
END INTERFACE

INTERFACE GETSTR
  MODULE PROCEDURE GETSTR
END INTERFACE

INTERFACE GETINTARRAY
  MODULE PROCEDURE GETINTARRAY
END INTERFACE

INTERFACE GETLOGICALARRAY
  MODULE PROCEDURE GETLOGICALARRAY
END INTERFACE

INTERFACE GETREALARRAY
  MODULE PROCEDURE GETREALARRAY
END INTERFACE

INTERFACE GETSTRARRAY
  MODULE PROCEDURE GETSTRARRAY
END INTERFACE

INTERFACE GETINTFROMSTR
  MODULE PROCEDURE GETINTFROMSTR
END INTERFACE

INTERFACE GETDESCRIPTION
  MODULE PROCEDURE GETDESCRIPTION
END INTERFACE

INTERFACE addStrListEntry
  MODULE PROCEDURE addStrListEntry
END INTERFACE

INTERFACE ExtractParameterFile
  MODULE PROCEDURE ExtractParameterFile
END INTERFACE

INTERFACE FinalizeParameters
  MODULE PROCEDURE FinalizeParameters
END INTERFACE

PUBLIC :: IgnoredParameters
PUBLIC :: PrintDefaultParameterFile
PUBLIC :: CountOption
PUBLIC :: GETINT
PUBLIC :: GETLOGICAL
PUBLIC :: GETREAL
PUBLIC :: GETSTR
PUBLIC :: GETINTARRAY
PUBLIC :: GETLOGICALARRAY
PUBLIC :: GETREALARRAY
PUBLIC :: GETSTRARRAY
PUBLIC :: GETDESCRIPTION
PUBLIC :: GETINTFROMSTR
PUBLIC :: addStrListEntry
PUBLIC :: FinalizeParameters
PUBLIC :: ExtractParameterFile

TYPE(Parameters) :: prms
PUBLIC :: prms

  type, public :: STR255
     private
     character(LEN=255) :: chars
  end type STR255
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Set actual section. All options created after calling this subroutine are in this 'section'. The section property is only
!> used to get nicer looking parameter files when using --help or --markdown.
!==================================================================================================================================
SUBROUTINE SetSection(this, section)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT) :: this                   !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)     :: section                !< section to set
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
this%actualSection = section
END SUBROUTINE SetSection

!==================================================================================================================================
!> General routine to create an option. 
!> Fills all fields of the option. Since the prms\%parse function is used to set the value, this routine can be abstract for all 
!> types of options. 
!==================================================================================================================================
SUBROUTINE CreateOption(this, opt, name, description, value, multiple)
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)       :: this             !< CLASS(Parameters)
CLASS(OPTION),INTENT(INOUT)           :: opt              !< option class
CHARACTER(LEN=*),INTENT(IN)           :: name             !< option name
CHARACTER(LEN=*),INTENT(IN)           :: description      !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL  :: value            !< option value
LOGICAL,INTENT(IN),OPTIONAL           :: multiple         !< marker if multiple option
! LOCAL VARIABLES
CLASS(link), POINTER :: newLink
!==================================================================================================================================
opt%hasDefault = PRESENT(value)
IF (opt%hasDefault) THEN
  CALL opt%parse(value)
END IF

opt%multiple   = .FALSE.
IF (PRESENT(multiple)) opt%multiple = multiple
IF (opt%multiple.AND.opt%hasDefault) THEN
  CALL Abort(__STAMP__, &
      "A default value can not be given, when multiple=.TRUE. in creation of option: '"//TRIM(name)//"'")
END IF

opt%isSet = .FALSE.
opt%name = name
opt%description = description
opt%section = this%actualSection
opt%isRemoved = .FALSE.

! insert option into linked list
IF (.not. associated(this%firstLink)) then
  this%firstLink => constructor_Link(opt, this%firstLink)
  this%lastLink => this%firstLink
ELSE
  newLink => constructor_Link(opt, this%lastLink%next)
  this%lastLink%next => newLink
  this%lastLink => newLink
END IF
END SUBROUTINE CreateOption

!==================================================================================================================================
!> Create a new integer option. Only calls the general prms\%createoption routine.
!==================================================================================================================================
SUBROUTINE CreateIntOption(this, name, description, value, multiple) 
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(IntOption),ALLOCATABLE,TARGET :: intopt
!==================================================================================================================================
ALLOCATE(intopt)
CALL this%CreateOption(intopt, name, description, value=value, multiple=multiple)
END SUBROUTINE CreateIntOption

!==================================================================================================================================
!> Create a new integer option with a optional string representation. Only calls the general prms\%createoption routine.
!==================================================================================================================================
SUBROUTINE CreateIntFromStringOption(this, name, description, value, multiple) 
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(IntFromStringOption),ALLOCATABLE,TARGET :: intfromstropt
!==================================================================================================================================
ALLOCATE(intfromstropt)
CALL this%CreateOption(intfromstropt, name, description, value=value, multiple=multiple)
END SUBROUTINE CreateIntFromStringOption

!==================================================================================================================================
!> Create a new logical option. Only calls the general prms\%createoption routine.
!==================================================================================================================================
SUBROUTINE CreateLogicalOption(this, name, description, value, multiple) 
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(LogicalOption),ALLOCATABLE,TARGET :: logicalopt
!==================================================================================================================================
ALLOCATE(logicalopt)
CALL this%CreateOption(logicalopt, name, description, value=value, multiple=multiple)
END SUBROUTINE CreateLogicalOption

!==================================================================================================================================
!> Create a new real option. Only calls the general prms\%createoption routine.
!==================================================================================================================================
SUBROUTINE CreateRealOption(this, name, description, value, multiple) 
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(RealOption),ALLOCATABLE,TARGET :: realopt
!==================================================================================================================================
ALLOCATE(realopt)
CALL this%CreateOption(realopt, name, description, value=value, multiple=multiple)
END SUBROUTINE CreateRealOption

!==================================================================================================================================
!> Create a new string option. Only calls the general prms\%createoption routine.
!==================================================================================================================================
SUBROUTINE CreateStringOption(this, name, description, value, multiple) 
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(StringOption),ALLOCATABLE,TARGET :: stringopt
!==================================================================================================================================
ALLOCATE(stringopt)
CALL this%CreateOption(stringopt, name, description, value=value, multiple=multiple)
END SUBROUTINE CreateStringOption

!==================================================================================================================================
!> Create a new integer array option. Only calls the general prms\%createoption routine.
!==================================================================================================================================
SUBROUTINE CreateIntArrayOption(this, name, description, value, multiple) 
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(IntArrayOption),ALLOCATABLE,TARGET :: intopt
!==================================================================================================================================
ALLOCATE(intopt)
CALL this%CreateOption(intopt, name, description, value=value, multiple=multiple)
END SUBROUTINE CreateIntArrayOption

!==================================================================================================================================
!> Create a new logical array option. Only calls the general prms\%createoption routine.
!==================================================================================================================================
SUBROUTINE CreateLogicalArrayOption(this, name, description, value, multiple) 
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(LogicalArrayOption),ALLOCATABLE,TARGET :: logicalopt
!==================================================================================================================================
ALLOCATE(logicalopt)
CALL this%CreateOption(logicalopt, name, description, value=value, multiple=multiple)
END SUBROUTINE CreateLogicalArrayOption

!==================================================================================================================================
!> Create a new real array option. Only calls the general prms\%createoption routine.
!==================================================================================================================================
SUBROUTINE CreateRealArrayOption(this, name, description, value, multiple) 
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(RealArrayOption),ALLOCATABLE,TARGET :: realopt
!==================================================================================================================================
ALLOCATE(realopt)
CALL this%CreateOption(realopt, name, description, value=value, multiple=multiple)
END SUBROUTINE CreateRealArrayOption

!==================================================================================================================================
!> Create a new string array option. Only calls the general prms\%createoption routine.
!==================================================================================================================================
!SUBROUTINE CreateStringArrayOption(this, name, description, value, multiple) 
!! INPUT/OUTPUT VARIABLES
!CLASS(Parameters),INTENT(INOUT)      :: this           !< CLASS(Parameters)
!CHARACTER(LEN=*),INTENT(IN)          :: name           !< option name
!CHARACTER(LEN=*),INTENT(IN)          :: description    !< option description
!CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: value          !< option value
!LOGICAL,INTENT(IN),OPTIONAL          :: multiple       !< marker if multiple option
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES 
!CLASS(StringArrayOption),ALLOCATABLE,TARGET :: stringopt
!!==================================================================================================================================
!ALLOCATE(stringopt)
!CALL this%CreateOption(stringopt, name, description, value=value, multiple=multiple)
!END SUBROUTINE CreateStringArrayOption

!==================================================================================================================================
!> Count number of occurrence of option with given name.
!==================================================================================================================================
FUNCTION CountOption_(this, name) result(count)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT) :: this  !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)     :: name  !< Search for this keyword in ini file
INTEGER                         :: count !< number of found occurences of keyword
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(link),POINTER :: current
!==================================================================================================================================
count = 0
! iterate over all options and compare names
current => this%firstLink
DO WHILE (associated(current))
  IF (current%opt%NAMEEQUALS(name)) THEN
    IF (current%opt%isSet) count = count + 1
  END IF
  current => current%next
END DO
END FUNCTION  CountOption_

!==================================================================================================================================
!> Insert a option in front of option with same name in the 'prms' linked list. 
!==================================================================================================================================
SUBROUTINE insertOption(first, opt)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(link),POINTER,INTENT(IN) :: first !< first item in linked list
CLASS(OPTION),INTENT(IN)       :: opt   !< option to be inserted
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(link),POINTER :: newLink
CLASS(link),POINTER :: current
!==================================================================================================================================
current =>  first
DO WHILE (associated(current%next))
  IF (.NOT.current%next%opt%NAMEEQUALS(opt%name)) THEN
    EXIT
  END IF
  current => current%next
END DO
newLink => constructor_Link(opt, current%next)
current%next => newLink
END SUBROUTINE insertOption    

!==================================================================================================================================
!> Read options from parameter file. 
!>
!> Therefore the file is read line by line. After removing comments and all white spaces each line is parsed in the 
!> prms\%read_option() routine. Outputs all unknown options.
!==================================================================================================================================
SUBROUTINE read_options(this, filename) 
! MODULES
USE MOD_StringTools ,ONLY: STRICMP,GetFileExtension
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(INOUT) :: this     !< CLASS(Parameters)
CHARACTER(LEN=255),INTENT(IN)   :: filename !< name of file to be read
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(link), POINTER  :: current
INTEGER               :: stat,iniUnit
TYPE(Varying_String)  :: aStr,bStr
CHARACTER(LEN=255)    :: HelpStr
LOGICAL               :: firstWarn=.TRUE.,file_exists
!==================================================================================================================================
CALL this%CreateLogicalOption('ColoredOutput','Colorize stdout, included for compatibility with FLEXI', '.TRUE.')

! Check if first argument is the ini-file 
IF(.NOT.(STRICMP(GetFileExtension(filename),'ini'))) THEN
  SWRITE(*,*) "Usage: flexi parameter.ini [restart.h5] [keyword arguments]"
  SWRITE(*,*) "   or: flexi restart.h5 [keyword arguments]"
  CALL CollectiveStop(__STAMP__,&
    'ERROR - Not an parameter file (file-extension must be .ini) or restart file (*.h5): '//TRIM(filename))
END IF
INQUIRE(FILE=TRIM(filename), EXIST=file_exists)  
IF (.NOT.file_exists) THEN
  CALL CollectiveStop(__STAMP__,&
      "Ini file does not exist.")
END IF

SWRITE(UNIT_StdOut,*)'| Reading from parameter file from "',TRIM(filename),'":'

! Open parameter file for reading
iniUnit= 100 !GETFREEUNIT()
OPEN(UNIT=iniUnit,FILE=TRIM(filename),STATUS='OLD',ACTION='READ',ACCESS='SEQUENTIAL',IOSTAT=stat)
IF(stat.NE.0) THEN
  CALL Abort(__STAMP__,&
      "Could not open ini file.")
END IF

! infinte loop. Exit at EOF
DO
  ! read a line into 'aStr'
  CALL Get(iniUnit,aStr,iostat=stat)
  ! exit loop if EOF
  IF(IS_IOSTAT_END(stat)) EXIT 
  IF(.NOT.IS_IOSTAT_EOR(stat)) THEN
    CALL Abort(__STAMP__,&
       'Error during ini file read')
  END IF

  ! Remove comments with "!"
  CALL Split(aStr,bStr,"!")
  ! Remove comments with "#"
  CALL Split(bStr,aStr,"#")
  ! aStr may hold an option 

  ! Remove blanks
  aStr=Replace(aStr," ","",Every=.true.)
  ! Replace brackets
  aStr=Replace(aStr,"(/"," ",Every=.true.)
  aStr=Replace(aStr,"/)"," ",Every=.true.)
  ! Lower case
  CALL LowCase(CHAR(aStr),HelpStr)
  ! If something remaind, this should be an option
  IF (LEN_TRIM(HelpStr).GT.2) THEN
    ! read the option
    IF (.NOT.this%read_option(HelpStr)) THEN
      IF (firstWarn) THEN
        firstWarn=.FALSE.
        SWRITE(UNIT_StdOut,'(100("!"))')
        SWRITE(UNIT_StdOut, *) "WARNING: The following options are unknown!"
      END IF
      SWRITE(UNIT_StdOut,*) "   ", TRIM(HelpStr)
    END IF
  END IF
END DO
IF (.NOT.firstWarn) THEN
  SWRITE(UNIT_StdOut,'(100("!"))')
END IF
CLOSE(iniUnit)

! calculate the maximal string lenght of all option-names and option-values
this%maxNameLen  = 0
this%maxValueLen = 0
current => prms%firstLink
DO WHILE (associated(current))
  this%maxNameLen = MAX(this%maxNameLen, current%opt%GETNAMELEN())
  this%maxValueLen = MAX(this%maxValueLen, current%opt%GETVALUELEN())
  current => current%next
END DO

! check for colored output 
use_escape_codes = GETLOGICAL("ColoredOutput")
END SUBROUTINE read_options

!==================================================================================================================================
!> Parses one line of parameter file and sets the value of the specific option in the 'prms' linked list.
!> Therefore it iterate over all entries of the linked list and compares the names.
!==================================================================================================================================
FUNCTION read_option(this, line) result(found)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(Parameters),INTENT(IN) :: this  !< CLASS(Parameters)
CHARACTER(LEN=*),INTENT(IN)  :: line  !< line to be parsed
LOGICAL                      :: found !< marker if option found
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)           :: name
CHARACTER(LEN=255)           :: rest
CLASS(link), POINTER         :: current
CLASS(OPTION),ALLOCATABLE    :: newopt
INTEGER                      :: i
!==================================================================================================================================
found = .FALSE.

! split at '='
i = index(line, '=')
IF (i==0) return
name = line(1:i-1)
rest = line(i+1:)

! iterate over all options and compare names
current => this%firstLink
DO WHILE (associated(current))
  ! compare name
  IF (current%opt%NAMEEQUALS(name)) THEN
    found = .TRUE.
    IF (current%opt%isSet) THEN
      IF (.NOT.(current%opt%multiple)) THEN
        ! option already set, but is not a multiple option 
        SWRITE(UNIT_StdOut,*) 'Option "', TRIM(name), '" is already set, but is not a multiple option!'
        STOP
      ELSE
        ! create new instance of multiple option
        ALLOCATE(newopt, source=current%opt)
        CALL newopt%parse(rest)
        newopt%isSet = .TRUE.
        ! insert option
        CALL insertOption(current, newopt)
        RETURN
      END IF
    END IF
    ! parse option
    IF(LEN_TRIM(rest).NE.0)THEN
      CALL current%opt%parse(rest)
      current%opt%isSet = .TRUE.
    ELSE
      CALL set_formatting("bright red")
      SWRITE(UNIT_StdOut,*) 'WARNING: Option "', TRIM(name), '" is specified in file but is empty!'
      CALL clear_formatting()
    END IF
    RETURN
  END IF
  current => current%next
END DO
END FUNCTION read_option

!==================================================================================================================================
!> Output all parameters, which are defined but NOT set in the parameter file. 
!==================================================================================================================================
SUBROUTINE IgnoredParameters() 
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(link), POINTER :: current
!==================================================================================================================================
current => prms%firstLink
CALL set_formatting("bright red")
SWRITE(UNIT_StdOut,'(100("!"))')
SWRITE(UNIT_StdOut,'(A)') "WARNING: The following options are defined, but NOT set in parameter-file or readin:"
DO WHILE (associated(current))
  IF (.NOT.current%opt%isRemoved) THEN
    SWRITE(UNIT_StdOut,*) "   ", TRIM(current%opt%name)
  END IF
  current => current%next
END DO
SWRITE(UNIT_StdOut,'(100("!"))')
CALL clear_formatting()
END SUBROUTINE IgnoredParameters 

!==================================================================================================================================
!> Print a default parameter file. The command line argument --help prints it in the format, that is used for reading the parameter
!> file. With --markdown one can print a default parameter file in markdown format.
!> Also prints the descriptions of a single parameter or parameter sections, if name corresponds to one of them.
!==================================================================================================================================
SUBROUTINE PrintDefaultParameterFile(markdown,name)
! MODULES
USE MOD_StringTools ,ONLY: STRICMP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN)   :: markdown  !< marker whether markdown format is used for output
CHARACTER(LEN=255)   :: name      !< for this parameter help is printed. If empty print all.
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CLASS(link), POINTER   :: current
CLASS(OPTION), POINTER :: currentOpt
INTEGER                :: maxNameLen
INTEGER                :: maxValueLen
INTEGER                :: lineLen
INTEGER                :: spaceNameLen
INTEGER                :: spaceValueLen
INTEGER                :: mode
CHARACTER(LEN=255)     :: section = "-"
CHARACTER(LEN=255)     :: singlesection = ""
CHARACTER(LEN=255)     :: singleoption = ""
CHARACTER(LEN=20)      :: fmtLineLen
CHARACTER(LEN=20)      :: fmtName
CHARACTER(LEN=20)      :: fmtValue
CHARACTER(LEN=20)      :: fmtComment
CHARACTER(LEN=20)      :: fmtNamespace
CHARACTER(LEN=20)      :: fmtValuespace
INTEGER                :: i
CHARACTER(LEN=255)     :: intFromStringOutput
CHARACTER(LEN=255)     :: fmtIntFromStringLength
CHARACTER(LEN=255)     :: fmtStringIntFromString
!==================================================================================================================================


maxNameLen  = 0
maxValueLen = 0
current => prms%firstLink
! check if name is a section or a option
DO WHILE (associated(current))
  IF (STRICMP(current%opt%section,name)) THEN
    singlesection = TRIM(name)
    EXIT
  END IF
  IF (current%opt%NAMEEQUALS(name)) THEN
    singleoption = TRIM(name)
    singlesection = TRIM(current%opt%section)
    EXIT
  END IF
  current => current%next
END DO

! if name is not specified, the complete parameter files needs to be printed
IF ((.NOT.markdown).AND.(LEN_TRIM(name).EQ.0)) THEN
  SWRITE(UNIT_StdOut,'(A80)')  "!==============================================================================="
  SWRITE(UNIT_StdOut,'(A)')    "! Default Parameter File generated using 'flexi --help' "
  SWRITE(UNIT_StdOut,'(4A)')   "!   compiled at : ", __DATE__," ", __TIME__ 
  SWRITE(UNIT_StdOut,'(A80)')  "!==============================================================================="
END IF

mode = 1
IF (markdown) THEN
  mode = 2
  SWRITE(UNIT_StdOut,'(A)') "## Parameterfile"
  SWRITE(UNIT_StdOut,'(A)') ""
END IF

! Find longest parameter name and length of the standard values
current => prms%firstLink
DO WHILE (associated(current))
  maxNameLen = MAX(maxNameLen, current%opt%GETNAMELEN())
  maxValueLen = MAX(maxValueLen, current%opt%GETVALUELEN())
  current => current%next
END DO
lineLen = maxNameLen + maxValueLen + 4 + 50
spaceNameLen = maxNameLen - 9
spaceValueLen = maxValueLen - 10
WRITE(fmtLineLen,*) lineLen
WRITE(fmtName,*)    maxNameLen
WRITE(fmtValue,*)   maxValueLen
WRITE(fmtComment,*) 50
WRITE(fmtNamespace,*) spaceNameLen
WRITE(fmtValuespace,*) spaceValueLen
current => prms%firstLink
DO WHILE (associated(current))
  IF (STRICMP(current%opt%section,'RecordPoints')) THEN
    current => current%next
    CYCLE
  END IF
  IF ((LEN_TRIM(singlesection).EQ.0).OR.(STRICMP(singlesection,current%opt%section))) THEN
    IF (.NOT.STRICMP(section,current%opt%section)) THEN
      section = current%opt%section
      IF (markdown) THEN
        SWRITE(UNIT_StdOut,'('//fmtLineLen//'("-"))')
        SWRITE(UNIT_StdOut,'(A2,A,A2)')                                 "**",TRIM(section),"**"
        SWRITE(UNIT_StdOut,'('//fmtName//'("-")"--"A1)', ADVANCE='NO')  " "
        SWRITE(UNIT_StdOut,'('//fmtValue//'("-")A1)', ADVANCE='NO')     " "
        SWRITE(UNIT_StdOut,'('//fmtComment//'("-"))')
        SWRITE(UNIT_StdOut,'(A)', ADVANCE='NO')                         "**Variable**"
        SWRITE(UNIT_StdOut,'('//fmtNamespace//'(" "))', ADVANCE='NO') 
        SWRITE(UNIT_StdOut,'(A)', ADVANCE='NO')                         "**Default**"
        SWRITE(UNIT_StdOut,'('//fmtValuespace//'(" "))', ADVANCE='NO') 
        SWRITE(UNIT_StdOut,'(A)')                                       "**Description**"
        SWRITE(UNIT_StdOut,'(A80)')                                     ""
      ELSE
        SWRITE(UNIT_StdOut,'(A1,'//fmtLineLen//'("="))') "!"
        SWRITE(UNIT_StdOut,'(A2,A)') "! ", TRIM(section)
        SWRITE(UNIT_StdOut,'(A1,'//fmtLineLen//'("="))') "!"
      END IF
    END IF

    IF ((LEN_TRIM(singleoption).EQ.0).OR.(current%opt%NAMEEQUALS(singleoption))) THEN
      CALL current%opt%print(maxNameLen, maxValueLen,mode)
    END IF

    ! If help is called for a single IntFromStringOption, print the possible values of this parameter
    IF (current%opt%NAMEEQUALS(singleoption)) THEN
      currentOpt => current%opt
      SELECT TYPE(currentOpt)
      CLASS IS (IntFromStringOption)
        SWRITE(UNIT_StdOut,'(A)') 'Possible options for this parameter are:'
        WRITE(fmtIntFromStringLength,*) currentOpt%maxLength   ! The biggest lenght of a named option
        WRITE(fmtStringIntFromString,*) "(A"//TRIM(fmtIntFromStringLength)//",A,I0,A)"
        DO i=1,SIZE(currentOpt%strList)
          ! Output is in the format STRING (INTEGER)
          WRITE(intFromStringOutput,TRIM(fmtStringIntFromString)) TRIM(currentOpt%strList(i)), ' (', currentOpt%intList(i), ')'
          SWRITE(UNIT_StdOut,'(A)') TRIM(intFromStringOutput)
        END DO
      END SELECT
    END IF
 
    ! print ------ line at the end of a section in markdown mode
    IF (associated(current%next).AND.markdown) THEN 
      IF (.NOT.STRICMP(section,current%next%opt%section)) THEN
        SWRITE(UNIT_StdOut,'('//fmtLineLen//'("-"))')
        SWRITE(UNIT_StdOut,*) ''    
      END IF
    END IF
  END IF
  current => current%next
END DO
END SUBROUTINE

!==================================================================================================================================
!> Creates a new link to a option-object, 'next' is the following link in the linked list
!==================================================================================================================================
FUNCTION constructor_Link(opt, next)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(link),POINTER            :: constructor_Link  !< new link
CLASS(OPTION),INTENT(IN)       :: opt               !< option to be linked
CLASS(link),INTENT(IN),POINTER :: next              !< next link
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
ALLOCATE(constructor_Link)
constructor_Link%next => next
ALLOCATE(constructor_Link%opt, SOURCE=opt)
END FUNCTION constructor_Link

!==================================================================================================================================
!> Count number of times a parameter is used within a file in case of multiple parameters. This only calls the internal
!> function countoption_ of the parameters class.
!==================================================================================================================================
FUNCTION CountOption(name) result(no)
! MODULES
USE MOD_Options
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN) :: name  !< parameter name
INTEGER                     :: no    !< number of parameters
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
no = prms%CountOption_(name)
END FUNCTION CountOption

!==================================================================================================================================
!> General routine to get an option. This routine is called from GETINT,GETREAL,GETLOGICAL,GETSTR to get the value a non-array
!> option. 
!==================================================================================================================================
SUBROUTINE GetGeneralOption(value, name, proposal)
USE MOD_Options
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name     !< parameter name
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal !< reference value
CLASS(*)                             :: value    !< parameter value 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CLASS(link),POINTER   :: current
CLASS(Option),POINTER :: opt
CHARACTER(LEN=255)    :: proposal_loc
!==================================================================================================================================

! iterate over all options
current => prms%firstLink
DO WHILE (associated(current))
  ! if name matches option
  IF (current%opt%NAMEEQUALS(name).AND.(.NOT.current%opt%isRemoved)) THEN
    opt => current%opt
    ! if proposal is present and the option is not set due to the parameter file, then return the proposal
    IF ((PRESENT(proposal)).AND.(.NOT.opt%isSet)) THEN
      proposal_loc = TRIM(proposal)
      CALL opt%parse(proposal_loc)
    ELSE
      ! no proposal, no default and also not set in parameter file => abort
      IF ((.NOT.opt%hasDefault).AND.(.NOT.opt%isSet)) THEN
        CALL ABORT(__STAMP__, &
            "Required option '"//TRIM(name)//"' not set in parameter file and has no default value.")
        RETURN
      END IF
    END IF
    ! copy value from option to result variable
    SELECT TYPE (opt)
    CLASS IS (IntOption)
      SELECT TYPE(value)
      TYPE IS (INTEGER)
        value = opt%value
      END SELECT
    CLASS IS (RealOption)
      SELECT TYPE(value)
      TYPE IS (REAL)
        value = opt%value
      END SELECT
    CLASS IS (LogicalOption)
      SELECT TYPE(value)
      TYPE IS (LOGICAL)
        value = opt%value
      END SELECT
    CLASS IS (StringOption)
      SELECT TYPE(value)
      TYPE IS (STR255)
        value%chars = opt%value
      END SELECT
    END SELECT
    ! print option and value to stdout
    CALL opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
    ! remove the option from the linked list of all parameters
    current%opt%isRemoved = .TRUE.
    RETURN
  END IF 
  current => current%next
END DO
CALL ABORT(__STAMP__, &
    'Option "'//TRIM(name)//'" is not defined in any DefineParameters... routine '//&
    'or already read (use GET... routine only for multiple options more than once).')
END SUBROUTINE GetGeneralOption
    
!==================================================================================================================================
!> General routine to get an array option. This routine is called from GETINTARRAY,GETREALARRAY,GETLOGICALARRAY,GETSTRARRAY to get
!> the value an array option.
!==================================================================================================================================
SUBROUTINE GetGeneralArrayOption(value, name, no, proposal)
USE MOD_Options
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name      !< parameter name
INTEGER,INTENT(IN)                   :: no        !< size of array
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal  !< reference value
CLASS(*)                             :: value(no) !< parameter value 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CLASS(link),POINTER   :: current
CLASS(Option),POINTER :: opt
CHARACTER(LEN=255)    :: proposal_loc
!INTEGER               :: i
!==================================================================================================================================

! iterate over all options
current => prms%firstLink
DO WHILE (associated(current))
  ! if name matches option
  IF (current%opt%NAMEEQUALS(name).AND.(.NOT.current%opt%isRemoved)) THEN
    opt => current%opt
    ! if proposal is present and the option is not set due to the parameter file, then return the proposal
    IF ((PRESENT(proposal)).AND.(.NOT.opt%isSet)) THEN
      proposal_loc = TRIM(proposal)
      CALL opt%parse(proposal_loc)
    ELSE
      ! no proposal, no default and also not set in parameter file => abort
      IF ((.NOT.opt%hasDefault).AND.(.NOT.opt%isSet)) THEN
        CALL ABORT(__STAMP__, &
            "Required option '"//TRIM(name)//"' not set in parameter file and has no default value.")
        RETURN
      END IF
    END IF
    ! copy value from option to result variable
    SELECT TYPE (opt)
    CLASS IS (IntArrayOption)
      IF (SIZE(opt%value).NE.no) CALL Abort(__STAMP__,"Array size of option '"//TRIM(name)//"' is not correct!")
      SELECT TYPE(value)
      TYPE IS (INTEGER)
        value = opt%value
      END SELECT
    CLASS IS (RealArrayOption)
      IF (SIZE(opt%value).NE.no) CALL Abort(__STAMP__,"Array size of option '"//TRIM(name)//"' is not correct!")
      SELECT TYPE(value)
      TYPE IS (REAL)
        value = opt%value
      END SELECT
    CLASS IS (LogicalArrayOption)
      IF (SIZE(opt%value).NE.no) CALL Abort(__STAMP__,"Array size of option '"//TRIM(name)//"' is not correct!")
      SELECT TYPE(value)
      TYPE IS (LOGICAL)
        value = opt%value
      END SELECT
    !CLASS IS (StringArrayOption)
      !IF (SIZE(opt%value).NE.no) CALL Abort(__STAMP__,"Array size of option '"//TRIM(name)//"' is not correct!")
      !SELECT TYPE(value)
      !TYPE IS (STR255)
        !DO i=1,no
          !value(i)%chars = opt%value(i)
        !END DO
      !END SELECT
    END SELECT
    ! print option and value to stdout
    CALL opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
    ! remove the option from the linked list of all parameters
    current%opt%isRemoved = .TRUE.
    RETURN
  END IF 
  current => current%next
END DO
CALL ABORT(__STAMP__, &
    'Option "'//TRIM(name)//'" is not defined in any DefineParameters... routine '//&
    'or already read (use GET... routine only for multiple options more than once).')
END SUBROUTINE GetGeneralArrayOption

!==================================================================================================================================
!> Get integer, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETINT(name, proposal) result(value)
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN) :: name              !< parameter name
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal !< reference value 
INTEGER                     :: value             !< parameter value 
!==================================================================================================================================
value = -1
CALL GetGeneralOption(value, name, proposal)
END FUNCTION GETINT

!==================================================================================================================================
!> Get logical, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETLOGICAL(name, proposal) result(value)
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN) :: name              !< parameter name
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal !< reference value 
LOGICAL                     :: value             !< parameter value 
!==================================================================================================================================
value = .FALSE.
CALL GetGeneralOption(value, name, proposal)
END FUNCTION GETLOGICAL

!==================================================================================================================================
!> Get real, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETREAL(name, proposal) result(value)
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name     !< parameter name
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal !< reference value 
REAL                                 :: value    !< parameter value 
!==================================================================================================================================
value = -1.0
CALL GetGeneralOption(value, name, proposal)
END FUNCTION GETREAL

!==================================================================================================================================
!> Get string, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETSTR(name, proposal) result(value)
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name     !< parameter name
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal !< reference value
CHARACTER(LEN=255)                   :: value    !< parameter value
! LOCAL VARIABLES
TYPE(STR255) :: tmp ! compiler bug workaround (gfortran 4.8.4)
!==================================================================================================================================
CALL GetGeneralOption(tmp, name, proposal)
value = tmp%chars
END FUNCTION GETSTR

!==================================================================================================================================
!> Get integer array, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETINTARRAY(name, no, proposal) result(value)
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name      !< parameter name
INTEGER,INTENT(IN)                   :: no        !< size of array
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal  !< reference value 
INTEGER                              :: value(no) !< array of integers
!==================================================================================================================================
value = -1
CALL GetGeneralArrayOption(value, name, no, proposal)
END FUNCTION GETINTARRAY

!==================================================================================================================================
!> Get logical array, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETLOGICALARRAY(name, no, proposal) result(value)
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name      !< parameter name
INTEGER,INTENT(IN)                   :: no        !< size of array
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal  !< reference value 
LOGICAL                              :: value(no) !< array of logicals
!==================================================================================================================================
value = .FALSE.
CALL GetGeneralArrayOption(value, name, no, proposal)
END FUNCTION GETLOGICALARRAY

!==================================================================================================================================
!> Get real array, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETREALARRAY(name, no, proposal) RESULT(value)
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name      !< parameter name
INTEGER,INTENT(IN)                   :: no        !< size of array
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal  !< reference value 
REAL                                 :: value(no) !< array of reals
!==================================================================================================================================
value = -1.
CALL GetGeneralArrayOption(value, name, no, proposal)
END FUNCTION GETREALARRAY

!==================================================================================================================================
!> Get string array, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETSTRARRAY(name, no, proposal) result(value)
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name      !< parameter name
INTEGER,INTENT(IN)                   :: no        !< size of array
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: proposal  !< reference value 
CHARACTER(LEN=255)                   :: value(no) !< array of strings
! LOCAL VARIABLES
TYPE(STR255) :: tmp(no) ! compiler bug workaround (gfortran 4.8.4)
INTEGER      :: i
!==================================================================================================================================
CALL GetGeneralArrayOption(tmp, name, no, proposal)
DO i = 1, no
  value(i)=tmp(i)%chars
END DO ! i = 1, no
END FUNCTION GETSTRARRAY

!==================================================================================================================================
!> Get string array, where proposal is used as default value, if the option was not set in parameter file 
!==================================================================================================================================
FUNCTION GETDESCRIPTION(name) result(description)
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: name        !< parameter name
CHARACTER(LEN=1000)                  :: description !< description
! LOCAL VARIABLES
CLASS(link),POINTER :: current
!==================================================================================================================================
! iterate over all options and compare names
current => prms%firstLink
DO WHILE (associated(current))
  IF (current%opt%NAMEEQUALS(name)) THEN
    description = current%opt%description
  END IF
  current => current%next
END DO
END FUNCTION GETDESCRIPTION

!==================================================================================================================================
!> GETINT for options with string values. Requires a map that provides the link between the 
!> possible integer values and the corresponding named values. This map is set using the addStrListEntry routine during 
!> parameter definition. If there is no named value to an option passed as int a warning is returned.
!==================================================================================================================================
FUNCTION GETINTFROMSTR(name) result(value)
USE MOD_StringTools ,ONLY: ISINT, STRICMP
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)   :: name        !< parameter name
INTEGER                       :: value       !< return value
! LOCAL VARIABLES
CLASS(link),POINTER           :: current
CLASS(Option),POINTER         :: opt
INTEGER                       :: i
LOGICAL                       :: found
INTEGER                       :: listSize         ! current size of list
!==================================================================================================================================
! iterate over all options and compare names
current => prms%firstLink
DO WHILE (associated(current))
  IF (current%opt%NAMEEQUALS(name)) THEN
    opt => current%opt
    SELECT TYPE (opt)
    CLASS IS (IntFromStringOption)
      ! Set flag indicating the given option has an entry in the mapping
      opt%foundInList = .TRUE.
      ! Size of list with string-integer pairs
      listSize = SIZE(opt%strList) 
      ! Check if an integer has been specfied directly
      IF (ISINT(opt%value)) THEN
        READ(opt%value,*) value
        found=.FALSE.
        ! Check if the integer is present in the list of possible integers
        DO i=1,listSize
          IF (opt%intList(i).EQ.value)THEN
            found=.TRUE.
            opt%listIndex = i ! Store index of the mapping
            EXIT
          END IF
        END DO
        ! If it is not found, print a warning and set the flag to later use the correct output format
        IF(.NOT.found)THEN
          CALL PrintWarning("No named option for parameter " //TRIM(name)// " exists for this number, please ensure your input is correct.")
          opt%foundInList = .FALSE.
        END IF
        CALL opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
        RETURN
      END IF
      ! If a string has been supplied, check if this string exists in the list and set it's integer representation according to the
      ! mapping
      DO i=1,listSize
        IF (STRICMP(opt%strList(i), opt%value)) THEN
          value = opt%intList(i)
          opt%listIndex = i ! Store index of the mapping
          CALL opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
          RETURN
        END IF
      END DO
      CALL Abort(__STAMP__,&
          "Unknown value for option: "//TRIM(name))
    END SELECT
  END IF
  current => current%next
END DO
CALL Abort(__STAMP__,&
    "Unknown option: "//TRIM(name))
END FUNCTION GETINTFROMSTR

!===================================================================================================================================
!> Add an entry to the mapping of string and integer values for the StringToInt option.
!===================================================================================================================================
SUBROUTINE addStrListEntry(name,string_in,int_in) 
USE MOD_Globals,     ONLY: abort
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)    :: name      !< parameter name
CHARACTER(LEN=*),INTENT(IN)    :: string_in !< (IN) string used for the option value
INTEGER         ,INTENT(IN)    :: int_in    !< (IN) integer used internally for the option value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CLASS(link), POINTER           :: current
CLASS(OPTION),POINTER          :: opt
INTEGER                        :: listSize         ! current size of list
CHARACTER(LEN=255),ALLOCATABLE :: strListTmp(:)    ! temporary string list
INTEGER           ,ALLOCATABLE :: intListTmp(:)    ! temporary integer list
!===================================================================================================================================
! iterate over all options and compare names
current => prms%firstLink
DO WHILE (associated(current))
  IF (current%opt%NAMEEQUALS(name)) THEN
    opt => current%opt
    SELECT TYPE (opt)
    CLASS IS (IntFromStringOption)
      ! Check if the arrays containing the string and integer values are already allocated
      IF (.NOT.(ALLOCATED(opt%strList))) THEN
        ! This is the first call to addEntry, allocate the arrays with dimension one
        ALLOCATE(opt%strList(1))
        ALLOCATE(opt%intList(1))
        ! Store the values in the lists
        opt%strList(1) = TRIM(string_in)
        opt%intList(1) = int_in
        ! Save biggest length of string entry
        opt%maxLength = LEN_TRIM(string_in)+4+INT(LOG10(REAL(int_in)))
      ELSE
        ! Subsequent call to addEntry, re-allocate the lists with one additional entry
        listSize = SIZE(opt%strList)    ! opt size of the list
        ! store opt values in temporary arrays
        ALLOCATE(strListTmp(listSize))
        ALLOCATE(intListTmp(listSize))
        strListTmp = opt%strList
        intListTmp = opt%intList
        ! Deallocate and re-allocate the list arrays
        SDEALLOCATE(opt%strList)
        SDEALLOCATE(opt%intList)
        ALLOCATE(opt%strList(listSize+1))
        ALLOCATE(opt%intList(listSize+1))
        ! Re-write the old values
        opt%strList(1:listSize) = strListTmp
        opt%intList(1:listSize) = intListTmp
        ! Deallocate temp arrays
        SDEALLOCATE(strListTmp)
        SDEALLOCATE(intListTmp)
        ! Now save the actual new entry in the list
        opt%strList(listSize+1) = TRIM(string_in)
        opt%intList(listSize+1) = int_in
        ! Save biggest length of string entry
        opt%maxLength = MAX(opt%maxLength,LEN_TRIM(string_in)+4+INT(LOG10(REAL(int_in))))
      END IF
      RETURN
    CLASS DEFAULT
      CALL Abort(__STAMP__,&
        "Option is not of type IntFromString: "//TRIM(name))
    END SELECT
  END IF
  current => current%next
END DO
CALL Abort(__STAMP__,&
    "Option not yet set: "//TRIM(name))

END SUBROUTINE addStrListEntry

SUBROUTINE ExtractParameterFile(filename,prmfile,userblockFound) 
! MODULES
USE MOD_StringTools ,ONLY: STRICMP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: filename !< name of file to be read
CHARACTER(LEN=*),INTENT(IN)   :: prmfile  !< name of file to be written
LOGICAL,INTENT(OUT)           :: userblockFound
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER               :: stat,iniUnit,fileUnit
TYPE(Varying_String)  :: aStr
CHARACTER(LEN=3)      :: tmp
LOGICAL               :: file_exists,iniFound
!==================================================================================================================================

IF (MPIRoot) THEN
  INQUIRE(FILE=TRIM(filename), EXIST=file_exists)  
  IF (.NOT.file_exists) THEN
    CALL CollectiveStop(__STAMP__,&
        "File '"//TRIM(filename)//"' does not exist.")
  END IF

  SWRITE(UNIT_StdOut,*)'| Extract parameter file from "',TRIM(filename),'" to "',TRIM(prmfile),'"'

  ! Open parameter file for reading
  fileUnit= 100 !GETFREEUNIT()
  OPEN(UNIT=fileUnit,FILE=TRIM(filename),STATUS='OLD',ACTION='READ',ACCESS='SEQUENTIAL',IOSTAT=stat)
  IF(stat.NE.0) THEN
    CALL Abort(__STAMP__,&
        "Could not open '"//TRIM(filename)//"'")
  END IF

  iniUnit= 101 !GETFREEUNIT()
  OPEN(UNIT=iniUnit,FILE=TRIM(prmfile),STATUS='UNKNOWN',ACTION='WRITE',ACCESS='SEQUENTIAL',IOSTAT=stat)
  IF(stat.NE.0) THEN
    CALL Abort(__STAMP__,&
        "Could not open '"//TRIM(prmfile)//"'")
  END IF

  iniFound = .FALSE.
  userblockFound = .FALSE.
  ! infinte loop. Exit at EOF
  DO
    ! read a line into 'aStr'
    CALL Get(fileUnit,aStr,iostat=stat)
    ! exit loop if EOF
    IF(IS_IOSTAT_END(stat)) EXIT 
    IF(.NOT.IS_IOSTAT_EOR(stat)) THEN
      CALL Abort(__STAMP__,&
          'Error during ini file read')
    END IF
    ! check if file starts "{[(" and therewith has a userblock
    IF (.NOT.userblockFound) THEN
      tmp = CHAR(extract(aStr,1,3))
      userblockFound = STRICMP(tmp,"{[(")
    END IF
    IF (.NOT.userblockFound) THEN
      SWRITE(*,*) "No Userblock found!"
      EXIT
    END IF

    ! search for begin of inifile
    IF (STRICMP(CHAR(aStr),"{[( INIFILE )]}")) THEN
      iniFound = .TRUE.
      CYCLE
    END IF
    IF (.NOT.iniFound) THEN
      ! if not found cycle (other userblock stuff)
      CYCLE
    ELSE
      ! if found and string starts with {[(, than this is the beginning of another userblock entry
      ! => finish reading of inifile
      tmp = CHAR(extract(aStr,1,3))
      IF (STRICMP(tmp, "{[(")) THEN
        EXIT
      END IF
    END IF
    WRITE(iniUnit,'(A)') CHAR(aStr)
  END DO

  CLOSE(fileUnit)
  CLOSE(iniUnit)
END IF
CALL MPI_BCAST(userblockFound,1,MPI_LOGICAL,0,MPI_COMM_WORLD,iError)

END SUBROUTINE ExtractParameterFile


!===================================================================================================================================
!> Clear parameters list 'prms'.
!===================================================================================================================================
SUBROUTINE FinalizeParameters() 
! LOCAL VARIABLES
CLASS(link), POINTER         :: current, tmp
!===================================================================================================================================

if(associated(prms%firstlink))then
  current => prms%firstLink
  DO WHILE (associated(current%next))
    DEALLOCATE(current%opt)
    NULLIFY(current%opt)
    tmp => current%next
    DEALLOCATE(current)
    NULLIFY(current)
    current => tmp
  END DO
end if
prms%firstLink => null()
prms%lastLink  => null()
END SUBROUTINE FinalizeParameters

END MODULE MOD_ReadInTools
