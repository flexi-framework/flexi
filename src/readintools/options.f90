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
!> Option-classes for values that are read from the parameter file (integer,logical, real, string; each single or array).
!==================================================================================================================================
MODULE MOD_Options
  USE MOD_Globals ,ONLY:MPIRoot,UNIT_stdOut
  IMPLICIT NONE
  PRIVATE
!================================================
!> Genereal, abstract OPTION
!================================================
  TYPE,PUBLIC  :: OPTION
    CLASS(OPTION),POINTER :: next         !< pointer to next option, used for a linked list of options
    CHARACTER(LEN=255)    :: name         !< name of the option, case-insensitive (part before '=' in parameter file)
    CHARACTER(LEN=1000)   :: description  !< comment in parameter file, after '!' character
    CHARACTER(LEN=255)    :: section      !< section to which the option belongs. Not mandatory.
    LOGICAL               :: isSet        !< default false. Becomes true, if set in parameter file
    LOGICAL               :: hasDefault   !< default false. True if a default value is given in CreateXXXOption routine
    LOGICAL               :: multiple     !< default false. Indicates if an option can occur multiple times in parameter file
    LOGICAL               :: isRemoved    !< default false. Indicates if the option is already used (GET... call) and therefore is
                                          !< no longer available in the list of parameters

  CONTAINS
    PROCEDURE :: print                    !< function used to print option for a default parameter file
    PROCEDURE :: printValue               !< function used to print the value
    PROCEDURE :: parse                    !< function that parses a string from the parameter file to fill the value of the option
    PROCEDURE :: parseReal                !< function that parses a string from the parameter file to fill the value of the option
    PROCEDURE :: NAMEEQUALS               !< function to compare case-insensitive a string with the name of this option
    PROCEDURE :: GETNAMELEN               !< function that returns the string length of the name
    PROCEDURE :: GETVALUELEN              !< function that returns the string length required to print the value
  END TYPE OPTION

!================================================
!> Integer Option
!================================================
  TYPE,PUBLIC,EXTENDS(OPTION) :: IntOption
    INTEGER :: value
  END TYPE IntOption

!================================================
!> \brief Integer from String Option
!> Used for options that are integer values, but also have a string representation.
!>
!> Many options are set as integers in the code for historical reasons and to use a short notation, but the user
!> can also specify a telling name for them in the parameter file which is more intuitive.
!================================================
  TYPE,PUBLIC,EXTENDS(OPTION) :: IntFromStringOption
    CHARACTER(LEN=255)              :: value
    INTEGER,ALLOCATABLE             :: intList(:)
    CHARACTER(LEN=255),ALLOCATABLE  :: strList(:)
    INTEGER                         :: listIndex
    LOGICAL                         :: foundInList = .FALSE.
    INTEGER                         :: maxLength=0
  END TYPE IntFromStringOption

!================================================
!> Integer Array Option
!================================================
  TYPE,PUBLIC,EXTENDS(OPTION) :: IntArrayOption
    INTEGER,ALLOCATABLE :: value(:)
  END TYPE IntArrayOption

!================================================
!> Logical Option
!================================================
  TYPE,PUBLIC,EXTENDS(OPTION) :: LogicalOption
    LOGICAL :: value
  END TYPE LogicalOption

!================================================
!> Logical Array Option
!================================================
  TYPE,PUBLIC,EXTENDS(OPTION) :: LogicalArrayOption
    LOGICAL,ALLOCATABLE :: value(:)
  END TYPE LogicalArrayOption

!================================================
!> Real Option
!================================================
  TYPE,PUBLIC,EXTENDS(OPTION) :: RealOption
    REAL    :: value
    INTEGER :: digits = 0 !< number of digits, the value has in parameter file
                          !< negative: -number of digits in exponential representation
                          !< 0 means not given
  END TYPE RealOption

!================================================
!> Real Array Option
!================================================
  TYPE,PUBLIC,EXTENDS(OPTION) :: RealArrayOption
    REAL,ALLOCATABLE    :: value(:)
    INTEGER,ALLOCATABLE :: digits(:) !< number of digits, the value has in parameter file
                                     !< negative: -number of digits in exponential representation
                                     !< 0 means not given
  END TYPE RealArrayOption

!================================================
!> String Option
!================================================
  TYPE,PUBLIC,EXTENDS(OPTION) :: StringOption
    CHARACTER(LEN=255) :: value
  END TYPE StringOption

!================================================
!> String Array Option
!================================================
  !TYPE,PUBLIC,EXTENDS(OPTION) :: StringArrayOption
    !CHARACTER(LEN=255),ALLOCATABLE  :: value(:)
  !END TYPE StringArrayOption

  INTERFACE GETSTRLENREAL
    MODULE PROCEDURE GETSTRLENREAL
  END INTERFACE
  PUBLIC :: GETSTRLENREAL

CONTAINS


!==================================================================================================================================
!> Compares name with the name of the option (case-insensitive)
!==================================================================================================================================
FUNCTION NAMEEQUALS(this, name)
! MODULES
USE MOD_StringTools ,ONLY: STRICMP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(OPTION),INTENT(IN)    :: this !< CLASS(OPTION)
CHARACTER(LEN=*),INTENT(IN) :: name !< incoming name, which is compared with the name of this option
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: NAMEEQUALS
!==================================================================================================================================
NAMEEQUALS = STRICMP(this%name, name)
END FUNCTION NAMEEQUALS

!==================================================================================================================================
!> return string-length of name
!==================================================================================================================================
FUNCTION GETNAMELEN(this)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(OPTION),INTENT(IN) :: this         !< CLASS(OPTION)
INTEGER                  :: GETNAMELEN   !< length of option name
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
GETNAMELEN = LEN_TRIM(this%name)
END FUNCTION GETNAMELEN

!==================================================================================================================================
!> return string-length required to print the value
!==================================================================================================================================
FUNCTION GETVALUELEN(this)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(OPTION),INTENT(IN)    :: this         !< CLASS(OPTION)
INTEGER                     :: GETVALUELEN  !< string length
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
CHARACTER(LEN=50) :: tmp
!==================================================================================================================================

GETVALUELEN = 0  ! default return value

! only if option is set of has default value, we have to find the length of this value
IF ((this%isSet).OR.(this%hasDefault)) THEN
  ! each class has a different type, which requires different commands to get the string-length of the value
  SELECT TYPE (this)
  CLASS IS (IntOption)
    WRITE(tmp,"(I0)") this%value
    GETVALUELEN = LEN_TRIM(tmp)
  CLASS IS (LogicalOption)
    GETVALUELEN = 1
  CLASS IS (RealOption)
    GETVALUELEN = GETSTRLENREAL(this%value,this%digits)
  CLASS IS (StringOption)
    GETVALUELEN = LEN_TRIM(this%value)
  CLASS IS (IntFromStringOption)
    GETVALUELEN = this%maxLength
  CLASS IS (IntArrayOption)
    GETVALUELEN = 3 ! '(/ '
    DO i=1,SIZE(this%value)
      WRITE(tmp,"(I0)") this%value(i)
      GETVALUELEN = GETVALUELEN + LEN_TRIM(tmp)
    END DO
    GETVALUELEN = GETVALUELEN + 2*(SIZE(this%value)-1) ! ', ' between array elements
    GETVALUELEN = GETVALUELEN + 3 ! ' /)'
  CLASS IS (LogicalArrayOption)
    GETVALUELEN = 3 ! '(/ '
    GETVALUELEN = GETVALUELEN + SIZE(this%value) ! each value needs only one character
    GETVALUELEN = GETVALUELEN + 2*(SIZE(this%value)-1) ! ', ' between array elements
    GETVALUELEN = GETVALUELEN + 3 ! ' /)'
  CLASS IS (RealArrayOption)
    GETVALUELEN = 3 ! '(/ '
    DO i=1,SIZE(this%value)
      GETVALUELEN = GETVALUELEN + GETSTRLENREAL(this%value(i), this%digits(i))
    END DO
    GETVALUELEN = GETVALUELEN + 2*(SIZE(this%value)-1) ! ', ' between array elements
    GETVALUELEN = GETVALUELEN + 3 ! ' /)'
  !CLASS IS (StringArrayOption)
    !GETVALUELEN = 3 ! '(/ '
    !DO i=1,SIZE(this%value)
      !GETVALUELEN = GETVALUELEN + LEN_TRIM(this%value(i))
    !END DO
    !GETVALUELEN = GETVALUELEN + 2*(SIZE(this%value)-1) ! ', ' between array elements
    !GETVALUELEN = GETVALUELEN + 3 ! ' /)'
  CLASS DEFAULT
    STOP 'Unknown TYPE'
  END SELECT
END IF
END FUNCTION GETVALUELEN

!===================================================================================================================================
!> Returns length of a real represented as string with a given number of digits
!===================================================================================================================================
FUNCTION GETSTRLENREAL(value,digits)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)         :: value         !< real value to print
INTEGER,INTENT(IN)      :: digits        !< number of digits (if < 1, then print as scientific)
INTEGER                 :: GETSTRLENREAL !< length of real as string
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=20)       :: fmtDigits
CHARACTER(LEN=50)       :: tmp
!===================================================================================================================================
IF (digits.GE.1) THEN ! floating point representation
  WRITE(fmtDigits,*) digits
  WRITE(tmp,'(F0.'//fmtDigits//')') value
  IF (index(tmp,'.').EQ.1) tmp = '0'//tmp(1:49)
ELSE IF (digits.LE.-1) THEN ! scientific (exponential) representation
  WRITE(fmtDigits,*) -digits
  WRITE(tmp,'(E24.'//fmtDigits//')') value
ELSE ! digits not given
  WRITE(tmp,'(E24.19)') value
END IF
GETSTRLENREAL = LEN(TRIM(ADJUSTL(tmp)))
END FUNCTION GETSTRLENREAL

!==================================================================================================================================
!> print option
!==================================================================================================================================
SUBROUTINE print(this, maxNameLen, maxValueLen, mode)
! MODULES
USE MOD_StringTools
USE MOD_ISO_VARYING_STRING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(OPTION),INTENT(IN)    :: this         !< option to print
INTEGER,INTENT(IN)          :: maxNameLen   !< max string length of name
INTEGER,INTENT(IN)          :: maxValueLen  !< max string length of value
INTEGER,INTENT(IN)          :: mode         !< 0: during readin, 1: default parameter file, 2: markdown
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=20)    :: fmtName
CHARACTER(LEN=20)    :: fmtValue
TYPE(VARYING_STRING) :: comment,headNewline,headSpace
INTEGER              :: length
INTEGER              :: commentLen
!==================================================================================================================================
IF(mode.EQ.1) commentLen=80 !--help
IF(mode.EQ.2) commentLen=50 !--markdown
WRITE(fmtName,*) maxNameLen
WRITE(fmtValue,*) maxValueLen
! print name
IF (mode.EQ.0) THEN
  WRITE(fmtName,*) maxNameLen
  SWRITE(UNIT_stdOut,'(a3)', ADVANCE='NO')  " | "
  CALL set_formatting("blue")
  SWRITE(UNIT_stdOut,"(a"//fmtName//")", ADVANCE='NO') TRIM(this%name)
  CALL clear_formatting()
ELSE
  SWRITE(UNIT_stdOut,"(A" // ADJUSTL(fmtName) // ")",ADVANCE='NO') this%name(:maxNameLen)
END IF

! print delimiter between name and value
SELECT CASE(mode)
CASE(0)
  SWRITE(UNIT_stdOut,'(a3)', ADVANCE='NO')  " | "
CASE(1)
  SWRITE(UNIT_stdOut,"(A3)",ADVANCE='NO') " = "
CASE(2)
  SWRITE(UNIT_stdOut,"(A3)",ADVANCE='NO') "   "
END SELECT

! print value
IF ((mode.EQ.0).OR.(this%hasDefault)) THEN
  CALL this%printValue(maxValueLen)
ELSE
  SWRITE(UNIT_stdOut, "(A"//fmtValue//")", ADVANCE='NO') ""
END IF


! print DEFAULT/CUSTOM or print comment
IF (mode.EQ.0) THEN
  ! print DEFAULT/CUSTOM
  IF (this%isSet) THEN
    SWRITE(UNIT_stdOut,"(a3)", ADVANCE='NO') ' | '
    CALL set_formatting("green")
    SWRITE(UNIT_stdOut,'(a7)', ADVANCE='NO')  "*CUSTOM"
    CALL clear_formatting()
    SWRITE(UNIT_stdOut,"(a3)") ' | '
  ELSE
    SWRITE(UNIT_stdOut,"(a3)", ADVANCE='NO') ' | '
    CALL set_formatting("red")
    SWRITE(UNIT_stdOut,'(a7)', ADVANCE='NO')  "DEFAULT"
    CALL clear_formatting()
    SWRITE(UNIT_stdOut,"(a3)") ' | '
  END IF
ELSE
  ! print comment: this is complicated, since it includes line breaks for long comments
  ! line breaks are inserted at spaces and at '\n' characters
  comment = TRIM(this%description)
  WRITE(fmtValue,*) maxNameLen + maxValueLen + 4
  length = 0
  SWRITE (UNIT_stdOut,'(A)',ADVANCE='NO') " " ! insert space after value
  ! loop until comment is empty (split by newline)
  DO WHILE (LEN_TRIM(comment) .GT. 0)
    ! split comment at first newline. After the split:
    !   - comment contains remaining part after newline
    !   - headNewline contains part before newline or ==comment, if there is no newline character in the comment anymore
    CALL SPLIT(comment, headNewline, "\n")
    ! loop until headNewline is empty (split by words)
    DO WHILE(LEN_TRIM(headNewline) .GT. 0)
      ! split comment at first space. After the split:
      !   - headNewline contains remaining part after the space
      !   - headSpace contains part before space or ==headNewline, if there is no space character in the headNewline anymore
      CALL SPLIT(headNewline, headSpace, " ")
      ! if word in headSpace does not fit into actual line -> insert newline
      IF (length+LEN_TRIM(headSpace).GT.commentLen) THEN
        SWRITE (UNIT_stdOut,*) ''
        SWRITE(UNIT_stdOut, "(A"//fmtValue//")", ADVANCE='NO') ""
        length = 0 ! reset length of line
      END IF
      ! insert word in headSpace and increase length of actual line
      IF ((length.EQ.0).AND.(mode.EQ.1)) THEN
        SWRITE(UNIT_stdOut,'(A2)',ADVANCE='NO') "! "
      END IF
      SWRITE (UNIT_stdOut,'(A)',ADVANCE='NO') CHAR(headSpace)//" "
      length = length + LEN_TRIM(headSpace)+1
    END DO
    ! insert linebreak due to newline character in comment
    SWRITE(UNIT_stdOut,*) ''
    IF ((LEN_TRIM(comment).GT.0).OR.(mode.EQ.2)) THEN
      SWRITE(UNIT_stdOut, "(A"//fmtValue//")", ADVANCE='NO') ""
    END IF
    length = 0
  END DO
  ! insert empty line after each option
  IF (mode.EQ.2) THEN
    SWRITE(UNIT_stdOut,*) ''
  END IF
END IF
END SUBROUTINE print

!==================================================================================================================================
!> print value of an option
!==================================================================================================================================
SUBROUTINE printValue(this,maxValueLen)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(OPTION),INTENT(IN)    :: this         !< option to print
INTEGER,INTENT(IN)          :: maxValueLen  !< max string length of name
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=20)            :: fmtValue
CHARACTER(LEN=20)            :: fmtDigits
CHARACTER(LEN=maxValueLen)   :: intFromStringOutput
INTEGER                      :: i,length
!===================================================================================================================================
WRITE(fmtValue,*) maxValueLen
SELECT TYPE (this)
CLASS IS (IntOption)
  SWRITE(UNIT_stdOut,"(I"//fmtValue//")",ADVANCE='NO') this%value
CLASS IS (LogicalOption)
  SWRITE(UNIT_stdOut,"(L"//fmtValue//")",ADVANCE='NO') this%value
CLASS IS (RealOption)
  IF (this%digits.GE.1) THEN ! floating point representation
    WRITE(fmtDigits,*) this%digits
    SWRITE(UNIT_stdOut,'(F'//fmtValue//'.'//fmtDigits//')',ADVANCE='NO') this%value
  ELSE IF (this%digits.LE.-1) THEN ! scientific (exponential) representation
    WRITE(fmtDigits,*) -this%digits
    SWRITE(UNIT_stdOut,'(E'//fmtValue//'.'//fmtDigits//')',ADVANCE='NO') this%value
  ELSE ! digits not given
    SWRITE(UNIT_stdOut,'(E'//fmtValue//'.19)',ADVANCE='NO') this%value
  END IF
CLASS IS (StringOption)
  IF (TRIM(this%value).EQ."") THEN
    SWRITE(UNIT_stdOut,"(A"//fmtValue//")",ADVANCE='NO') '-'
  ELSE
    SWRITE(UNIT_stdOut,"(A"//fmtValue//")",ADVANCE='NO') TRIM(this%value)
  END IF
CLASS IS (IntFromStringOption)
  IF (TRIM(this%value).EQ."") THEN
    SWRITE(UNIT_stdOut,"(A"//fmtValue//")",ADVANCE='NO') '-'
  ELSE
    ! IntFromStringOption: Print in the format STRING (INTEGER) if the given value is found in the mapping,
    ! otherwise print just the integer
    IF (this%foundInList) THEN
      WRITE(intFromStringOutput,"(A,A,I0,A)") TRIM(this%strList(this%listIndex)), ' (', this%intList(this%listIndex), ')'
      SWRITE(UNIT_stdOut,"(A"//fmtValue//")",ADVANCE='NO') TRIM(intFromStringOutput)
    ELSE
      SWRITE(UNIT_stdOut,"(A"//fmtValue//")",ADVANCE='NO') TRIM(this%value)
    END IF
  END IF
CLASS IS (IntArrayOption)
  length=this%GETVALUELEN()
  IF (maxValueLen - length.GT.0) THEN
    WRITE(fmtValue,*) (maxValueLen - length)
    SWRITE(UNIT_stdOut,'('//fmtValue//'(" "))',ADVANCE='NO')
  END IF
  SWRITE(UNIT_stdOut,"(A3)",ADVANCE='NO') "(/ "
  DO i=1,SIZE(this%value)
    WRITE(fmtValue,'(I0)') this%value(i)
    WRITE(fmtValue,*) LEN_TRIM(fmtValue)
    SWRITE(UNIT_stdOut,"(I"//fmtValue//")",ADVANCE='NO') this%value(i)
    IF (i.NE.SIZE(this%value)) THEN
      SWRITE(UNIT_stdOut,"(A2)",ADVANCE='NO') ", "
    END IF
  END DO
  SWRITE(UNIT_stdOut,"(A3)",ADVANCE='NO') " /)"
CLASS IS (LogicalArrayOption)
  length=this%GETVALUELEN()
  IF (maxValueLen - length.GT.0) THEN
    WRITE(fmtValue,*) (maxValueLen - length)
    SWRITE(UNIT_stdOut,'('//fmtValue//'(" "))',ADVANCE='NO')
  END IF
  SWRITE(UNIT_stdOut,"(A3)",ADVANCE='NO') "(/ "
  DO i=1,SIZE(this%value)
    SWRITE(UNIT_stdOut,"(L1)",ADVANCE='NO') this%value(i)
    IF (i.NE.SIZE(this%value)) THEN
      SWRITE(UNIT_stdOut,"(A2)",ADVANCE='NO') ", "
    END IF
  END DO
  SWRITE(UNIT_stdOut,"(A3)",ADVANCE='NO') " /)"
CLASS IS (RealArrayOption)
  length=this%GETVALUELEN()
  IF (maxValueLen - length.GT.0) THEN
    WRITE(fmtValue,*) (maxValueLen - length)
    SWRITE(UNIT_stdOut,'('//fmtValue//'(" "))',ADVANCE='NO')
  END IF
  SWRITE(UNIT_stdOut,'(a3)',ADVANCE='NO') '(/ '
  DO i=1,SIZE(this%value)
    WRITE(fmtValue,*) GETSTRLENREAL(this%value(i), this%digits(i))
      IF (this%digits(i).GE.1) THEN ! floating point representation
      WRITE(fmtDigits,*) this%digits(i)
      SWRITE(UNIT_stdOut,'(F'//fmtValue//'.'//fmtDigits//',A3)',ADVANCE='NO') this%value(i)
    ELSE IF (this%digits(i).LE.-1) THEN ! scientific (exponential) representation
      WRITE(fmtDigits,*) -this%digits(i)
      SWRITE(UNIT_stdOut,'(E'//fmtValue//'.'//fmtDigits//')',ADVANCE='NO') this%value(i)
    ELSE ! digits not given
      SWRITE(UNIT_stdOut,'(E'//fmtValue//'.19,A3)',ADVANCE='NO') this%value(i)
    END IF
    IF (i.NE.SIZE(this%value)) THEN
      SWRITE(UNIT_stdOut,'(A2)',ADVANCE='NO') ', '
    END IF
  END DO
  SWRITE(UNIT_stdOut,'(a3)',ADVANCE='NO') ' /)'
!###
! TODO: Causes internal compiler error with GNU 6+ due to compiler bug (older GNU and Intel,Cray work). Uncomment as unused.
!###
!CLASS IS (StringArrayOption)
  !length=this%GETVALUELEN()
  !IF (maxValueLen - length.GT.0) THEN
    !WRITE(fmtValue,*) (maxValueLen - length)
    !SWRITE(UNIT_stdOut,'('//fmtValue//'(" "))',ADVANCE='NO')
  !END IF
  !SWRITE(UNIT_stdOut,"(A3)",ADVANCE='NO') "(/ "
  !DO i=1,SIZE(this%value)
    !WRITE(fmtValue,*) LEN_TRIM(this%value(i))
    !SWRITE(UNIT_stdOut,"(A"//fmtValue//")",ADVANCE='NO') this%value(i)
    !IF (i.NE.SIZE(this%value)) THEN
      !SWRITE(UNIT_stdOut,"(A2)",ADVANCE='NO') ", "
    !END IF
  !END DO
  !SWRITE(UNIT_stdOut,"(A3)",ADVANCE='NO') " /)"
CLASS DEFAULT
  STOP
END SELECT
END SUBROUTINE printValue

!==================================================================================================================================
!> parse value from string 'rest_in'. This subroutine is used to readin values from the parameter file.
!==================================================================================================================================
SUBROUTINE parse(this, rest_in)
! MODULES
USE MOD_Globals, ONLY:abort
USE MOD_ISO_VARYING_STRING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CLASS(OPTION)               :: this     !< CLASS(OPTION)
CHARACTER(LEN=*),INTENT(IN) :: rest_in  !< string to parse
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)  :: tmp,tmp2,rest
INTEGER             :: count,i,stat
INTEGER,ALLOCATABLE :: inttmp(:)
LOGICAL,ALLOCATABLE :: logtmp(:)
REAL,ALLOCATABLE    :: realtmp(:)
!CHARACTER(LEN=255),ALLOCATABLE :: strtmp(:)
!==================================================================================================================================
stat=0
! Replace brackets
rest=Replace(rest_in,"(/"," ",Every=.true.)
rest=Replace(rest,   "/)"," ",Every=.true.)
rest = TRIM(rest)
IF(LEN_TRIM(rest).EQ.0)THEN
  CALL Abort(__STAMP__,'Variable '//TRIM(this%name)//' is empty!')
END IF

SELECT TYPE (this)
CLASS IS (IntOption)
  READ(rest, *,IOSTAT=stat) this%value
CLASS IS (LogicalOption)
  READ(rest, *,IOSTAT=stat) this%value
CLASS IS (RealOption)
  CALL this%parseReal(rest, this%value, this%digits)
CLASS IS (StringOption)
  READ(rest, "(A)",IOSTAT=stat) this%value
CLASS IS (IntFromStringOption)
  READ(rest, "(A)",IOSTAT=stat) this%value
CLASS IS (IntArrayOption)
  ! Array options are complicated, since we do not know a priori, how long the array will be.
  ! Therefore we must read entry by entry and always increase the array size by reallocation.
  ! This is done using the temporary arrays inttmp/realtmp and with the MOVE_ALLOC routine.
  tmp2 = TRIM(ADJUSTL(rest))
  IF (ALLOCATED(this%value)) DEALLOCATE(this%value)
  count = 0
  ALLOCATE(this%value(count))
  DO ! loop until no more entry
    i = index(TRIM(tmp2), ',')
    ! store text of tmp2 until next , in tmp
    IF (i.GT.0) THEN
      tmp = tmp2(1:i-1)
    ELSE
      tmp = tmp2
    END IF
    ! remove entry and trim the result
    tmp2 = tmp2(i+1:)
    tmp2 = TRIM(ADJUSTL(tmp2))
    ! increase the number of entries
    count = count + 1
    ! allocate temporary array and copy content of this%value to the first entries
    ! aftwards use MOVE_ALLOC to move it back to this%value (this deallocates the temporary array)
    ALLOCATE(inttmp(count))
    inttmp(:size(this%value)) = this%value
    CALL MOVE_ALLOC(inttmp, this%value) ! inttmp gets deallocated
    ! finally this%value has the correct size and we can readin the entry
    READ(tmp, *,IOSTAT=stat) this%value(count)
    IF (i.EQ.0) EXIT
  END DO
CLASS IS (LogicalArrayOption)
  ! comments see IntArrayOption!
  tmp2 = TRIM(ADJUSTL(rest))
  IF (ALLOCATED(this%value)) DEALLOCATE(this%value)
  count = 0
  ALLOCATE(this%value(count))
  DO ! loop until no more entry
    i = index(TRIM(tmp2), ',')
    ! store text of tmp2 until next , in tmp
    IF (i.GT.0) THEN
      tmp = tmp2(1:i-1)
    ELSE
      tmp = tmp2
    END IF
    ! remove entry and trim the result
    tmp2 = tmp2(i+1:)
    tmp2 = TRIM(ADJUSTL(tmp2))
    ! increase the number of entries
    count = count + 1
    ! allocate temporary array and copy content of this%value to the first entries
    ! aftwards use MOVE_ALLOC to move it back to this%value (this deallocates the temporary array)
    ALLOCATE(logtmp(count))
    logtmp(:size(this%value)) = this%value
    CALL MOVE_ALLOC(logtmp, this%value) ! logtmp gets deallocated
    ! finally this%value has the correct size and we can readin the entry
    READ(tmp, *,IOSTAT=stat) this%value(count)
    IF (i.EQ.0) EXIT
  END DO
CLASS IS (RealArrayOption)
  ! comments see IntArrayOption!
  tmp2 = TRIM(ADJUSTL(rest))
  IF (ALLOCATED(this%value)) DEALLOCATE(this%value)
  IF (ALLOCATED(this%digits)) DEALLOCATE(this%digits)
  count = 0
  count = 0
  ALLOCATE(this%value(count))
  ALLOCATE(this%digits(count))
  DO
    i = index(TRIM(tmp2), ',')
    IF (i.GT.0) THEN
      tmp = tmp2(1:i-1)
    ELSE
      tmp = tmp2
    END IF
    tmp2 = tmp2(i+1:)
    tmp2 = TRIM(ADJUSTL(tmp2))
    count = count + 1
    ALLOCATE(realtmp(count))
    ALLOCATE(inttmp(count))
    realtmp(:size(this%value)) = this%value
    inttmp(:size(this%digits)) = this%digits
    CALL MOVE_ALLOC(realtmp, this%value)  ! realtmp gets deallocated
    CALL MOVE_ALLOC(inttmp,  this%digits) ! realtmp gets deallocated
    this%digits(count) = 0

    CALL this%parseReal(tmp, this%value(count), this%digits(count))

    IF (i.EQ.0) EXIT
  END DO
!CLASS IS (StringArrayOption)
  !! comments see IntArrayOption!
  !tmp2 = TRIM(ADJUSTL(rest))
  !IF (ALLOCATED(this%value)) DEALLOCATE(this%value)
  !count = 0
  !ALLOCATE(this%value(count))
  !DO ! loop until no more entry
    !i = index(TRIM(tmp2), ',')
    !! store text of tmp2 until next , in tmp
    !IF (i.GT.0) THEN
      !tmp = tmp2(1:i-1)
    !ELSE
      !tmp = tmp2
    !END IF
    !! remove entry and trim the result
    !tmp2 = tmp2(i+1:)
    !tmp2 = TRIM(ADJUSTL(tmp2))
    !! increase the number of entries
    !count = count + 1
    !! allocate temporary array and copy content of this%value to the first entries
    !! aftwards use MOVE_ALLOC to move it back to this%value (this deallocates the temporary array)
    !ALLOCATE(strtmp(count))
    !strtmp(:size(this%value)) = this%value
    !CALL MOVE_ALLOC(strtmp, this%value) ! strtmp gets deallocated
    !! finally this%value has the correct size and we can readin the entry
    !READ(tmp, *,IOSTAT=stat) this%value(count)
    !IF (i.EQ.0) EXIT
  !END DO
CLASS DEFAULT
  STOP
END SELECT
IF(stat.GT.0)THEN
  CALL Abort(__STAMP__,&
    "Failed to parse: "//TRIM(this%name))
END IF

END SUBROUTINE parse

!===================================================================================================================================
!> parse string to real and get the format of the number (floating,scientific)
!===================================================================================================================================
SUBROUTINE parseReal(this,string_in, value, digits)
USE MOD_Globals, ONLY:abort
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CLASS(OPTION)                 :: this      !< CLASS(OPTION)
CHARACTER(LEN=255),INTENT(IN) :: string_in !< (IN) string containing a real number
REAL,INTENT(OUT)              :: value     !< (OUT) the converted real
INTEGER,INTENT(OUT)           :: digits    !< (OUT) the number of digits if floating representation, or -1 if scientific
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: pos,posE,posMinus,stat
CHARACTER(LEN=255) :: string  ! left adjusted string
!===================================================================================================================================
string = ADJUSTL(string_in)
posE=MAX(index(string, 'E'),index(string, 'e'))
IF (posE.GT.0) THEN
  ! exponential representation, will be printed as '0.xxEyy'. The number of digits (here: the length of 'xx') is needed.
  ! This number is the number of significant digits in the string, before and after the dot.
  pos = index(string, '.')
  IF (pos.EQ.0) THEN
    ! no dot in string, number of digits is length of string before the 'e' == posE-1
    digits = -(posE-1)
  ELSE
    ! dot in string, number of digits is length of string before the 'e' -1 for the dot == posE-2
    digits = -(posE-2)
  END IF
  posMinus = index(string,'-')
  IF (posMinus.EQ.1) THEN
    ! string starts with a minus => decrease number of digits by one (here +1, since negative values for scientific format)
    digits = digits+1
    ! string is -0.xxx, meaning no significant stuff before the dot => remove the zero form the number of digits
    IF (index(string,'0.').EQ.2) digits = digits+1
  ELSE
    ! string is 0.xxx, meaning no significant stuff before the dot => remove the zero form the number of digits
    IF (index(string,'0.').EQ.1) digits = digits+1
  END IF

ELSE
  ! floating representation
  pos = index(string, '.')
  IF (pos.EQ.0) THEN
    digits = 1
  ELSE
    digits = LEN_TRIM(string) - pos
    IF (digits.EQ.0) digits = 1
  END IF
END IF
READ(string, *,IOSTAT=stat) value
IF(stat.GT.0)THEN
  CALL Abort(__STAMP__,&
    "Failed to parse: "//TRIM(this%name))
END IF

END SUBROUTINE parseReal

END module
