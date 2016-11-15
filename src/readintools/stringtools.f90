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
!
! ATTENTION:
! The routines 'clear_formatting', 'set_formatting', 'get_escape_sequence' and 'split_string' are copied from the fortran output 
! library (foul).
! Copyright and license see below.
! Full version of the foul-library can be found here:
!   http://foul.sourceforge.net
!
!------------------------------------------------------------
! foul - The Fortran Output Library
!------------------------------------------------------------
! Provides routines enabling Fortran programs to:
!
! - Write formatted console output using ANSI escape codes
!   (see http://en.wikipedia.org/wiki/ANSI_escape_code
!    for more information)
! - Convert numbers to strings in a variety of
!   finely-controlled ways, including fully trimmed
!   of whitespace characters
! - Time program execution and other processes with the
!   highest accuracy provided by the system
!------------------------------------------------------------
! Copyright (C) 2010-2011 by Philipp Emanuel Weidmann
! E-Mail: philipp.weidmann@gmx.de
!------------------------------------------------------------
! This library is free software; you can redistribute it
! and/or modify it under the terms of the GNU General Public
! License as published by the Free Software Foundation;
! version 3 of the License.
!
! This library is distributed in the hope that it will be
! useful, but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
! PURPOSE.
!
! See the GNU General Public License for more details.
!------------------------------------------------------------

#include "flexi.h"

!==================================================================================================================================
!> Routines for performing operations on strings, which are not covered by ISO_VARYING_STRING.
!> The routines 'clear_formatting', 'set_formatting', 'get_escape_sequence' and 'split_string' are copied from the fortran output 
!> library (foul).
!> Full version of the foul-library can be found here:
!>   http://foul.sourceforge.net
!==================================================================================================================================
MODULE MOD_StringTools

USE MOD_Globals
USE MOD_ISO_VARYING_STRING
IMPLICIT NONE
PRIVATE

INTERFACE LowCase
  MODULE PROCEDURE LowCase
  MODULE PROCEDURE LowCase_overwrite
END INTERFACE

INTERFACE STRICMP
  MODULE PROCEDURE STRICMP
END INTERFACE

INTERFACE StripSpaces
  MODULE PROCEDURE StripSpaces
END INTERFACE

INTERFACE INTTOSTR
  MODULE PROCEDURE INTTOSTR
END INTERFACE

INTERFACE ISINT
  MODULE PROCEDURE ISINT
END INTERFACE

INTERFACE set_formatting
  MODULE PROCEDURE set_formatting
END INTERFACE

INTERFACE clear_formatting
  MODULE PROCEDURE clear_formatting
END INTERFACE

PUBLIC :: LowCase
PUBLIC :: STRICMP
PUBLIC :: StripSpaces
PUBLIC :: INTTOSTR
PUBLIC :: ISINT
PUBLIC :: set_formatting
PUBLIC :: clear_formatting

LOGICAL :: use_escape_codes = .TRUE.  !< If set to .FALSE., output will consist only of standard text, allowing the 
                                      !< escape characters to be switched off in environments which don't support them.
PUBLIC :: use_escape_codes                                      
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Transform upper case letters in "Str1" into lower case letters, result is "Str2" (in place version)
!==================================================================================================================================
SUBROUTINE LowCase_overwrite(Str1)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(INOUT) :: Str1  !< Input/output string
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                    :: iLen,nLen,Upper
CHARACTER(LEN=*),PARAMETER :: lc='abcdefghijklmnopqrstuvwxyz'
CHARACTER(LEN=*),PARAMETER :: UC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
LOGICAL                    :: HasEq
!==================================================================================================================================
HasEq=.FALSE.
nLen=LEN_TRIM(Str1)
DO iLen=1,nLen
  ! Transformation stops at "="
  IF(Str1(iLen:iLen).EQ.'=') HasEq=.TRUE.
  Upper=INDEX(UC,Str1(iLen:iLen))
  IF ((Upper > 0).AND. .NOT. HasEq) THEN
    Str1(iLen:iLen) = lc(Upper:Upper)
  END IF  
END DO
END SUBROUTINE LowCase_overwrite

!==================================================================================================================================
!> Transform upper case letters in "Str1" into lower case letters, result is "Str2"
!==================================================================================================================================
SUBROUTINE LowCase(Str1,Str2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: Str1  !< Input string
CHARACTER(LEN=*),INTENT(OUT) :: Str2  !< Output string, lower case letters only
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                    :: iLen,nLen,Upper
CHARACTER(LEN=*),PARAMETER :: lc='abcdefghijklmnopqrstuvwxyz'
CHARACTER(LEN=*),PARAMETER :: UC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
LOGICAL                    :: HasEq
!==================================================================================================================================
HasEq=.FALSE.
Str2=Str1
nLen=LEN_TRIM(Str1)
DO iLen=1,nLen
  ! Transformation stops at "="
  IF(Str1(iLen:iLen).EQ.'=') HasEq=.TRUE.
  Upper=INDEX(UC,Str1(iLen:iLen))
  IF ((Upper > 0).AND. .NOT. HasEq) Str2(iLen:iLen) = lc(Upper:Upper)
END DO
END SUBROUTINE LowCase

!==================================================================================================================================
!> Case insensitive string comparison 
!==================================================================================================================================
FUNCTION STRICMP(a, b)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: a,b !< strings to compare with each other
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
LOGICAL            :: STRICMP
CHARACTER(LEN=255) :: alow
CHARACTER(LEN=255) :: blow
!==================================================================================================================================
CALL LowCase(a, alow)
CALL LowCase(b, blow)
STRICMP = (TRIM(alow).EQ.TRIM(blow))
END FUNCTION STRICMP

!==================================================================================================================================
!> Removes ALL whitespace from a string
!==================================================================================================================================
SUBROUTINE StripSpaces(string)
! MODULES
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(INOUT) :: string  !< input string
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: stringLen 
INTEGER :: last, actual
!==================================================================================================================================
stringLen = LEN(string)
last = 1
actual = 1
DO WHILE (actual < stringLen)
  IF(string(last:last) == ' ') then
    actual = actual + 1
    string(last:last) = string(actual:actual)
    string(actual:actual) = ' '
  ELSE
    last = last + 1
    IF (actual < last) &
        actual = last
  ENDIF
END DO
END SUBROUTINE

!==================================================================================================================================
!> Converts integer to string
!==================================================================================================================================
PURE FUNCTION INTTOSTR(value) 
INTEGER,INTENT(IN)  :: value
CHARACTER(LEN=255)  :: INTTOSTR
WRITE(INTTOSTR,"(I20)") value
END FUNCTION INTTOSTR

!==================================================================================================================================
!> Checks if a string is an integer
!==================================================================================================================================
PURE FUNCTION ISINT(value) 
CHARACTER(LEN=255),INTENT(IN)  :: value
LOGICAL                        :: ISINT
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: i,stat
!==================================================================================================================================
READ(value, *, IOSTAT=stat) i
ISINT=(stat.EQ.0)
END FUNCTION ISINT

!==================================================================================================================================
!> Splits the supplied string along a delimiter.
!> This function is copied from the fortran output library (foul). For the full version of this library see:
!>   http://foul.sourceforge.net
!==================================================================================================================================
SUBROUTINE split_string(string, delimiter, substrings, substring_count)
! MODULES
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER (LEN = *), INTENT(IN)  :: string          !< Variable-length character string that is to be split
CHARACTER,           INTENT(IN)  :: delimiter       !< Character along which to split
CHARACTER (LEN = *), INTENT(OUT) :: substrings(*)   !< Array of substrings generated by split operation
INTEGER,             INTENT(OUT) :: substring_count !< Number of substrings generated
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: start_position, end_position
!==================================================================================================================================
start_position  = 1
substring_count = 0

DO
  end_position = INDEX(string(start_position:), delimiter)
  substring_count = substring_count + 1

  IF (end_position == 0) THEN
    substrings(substring_count) = string(start_position:)
    EXIT
  ELSE
    substrings(substring_count) = string(start_position : start_position + end_position - 2)
    start_position = start_position + end_position
  END IF
END DO
END SUBROUTINE split_string

!==================================================================================================================================
!> Generates an ANSI escape sequence from the supplied style string.
!> This function is copied from the fortran output library (foul). For the full version of this library see:
!>   http://foul.sourceforge.net
!==================================================================================================================================
SUBROUTINE get_escape_sequence(style_string, escape_sequence)
! MODULES
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)   :: style_string    !< String describing which styles to set (separated by space)
                                                 !< see source code for supported styles
CHARACTER(LEN=16),INTENT(OUT) :: escape_sequence !< escape_sequence: ANSI escape sequence generated from the specified styles
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
CHARACTER(LEN=32) :: style_substrings(16)
INTEGER           :: style_substring_count
!==================================================================================================================================
! Start sequence with command to clear any previous attributes
escape_sequence = CHAR(27) // '[0'

CALL split_string(TRIM(style_string), ' ', style_substrings, style_substring_count)

DO i = 1, style_substring_count
  CALL LowCase(style_substrings(i))

  SELECT CASE (TRIM(style_substrings(i)))
  CASE ('bright')
    escape_sequence = TRIM(escape_sequence) // ';1'
  CASE ('faint')
    escape_sequence = TRIM(escape_sequence) // ';2'
  CASE ('italic')
    escape_sequence = TRIM(escape_sequence) // ';3'
  CASE ('underline')
    escape_sequence = TRIM(escape_sequence) // ';4'
  CASE ('blink_slow')
    escape_sequence = TRIM(escape_sequence) // ';5'
  CASE ('blink_fast')
    escape_sequence = TRIM(escape_sequence) // ';6'
  CASE ('black')
    escape_sequence = TRIM(escape_sequence) // ';30'
  CASE ('red')
    escape_sequence = TRIM(escape_sequence) // ';31'
  CASE ('green')
    escape_sequence = TRIM(escape_sequence) // ';32'
  CASE ('yellow')
    escape_sequence = TRIM(escape_sequence) // ';33'
  CASE ('blue')
    escape_sequence = TRIM(escape_sequence) // ';34'
  CASE ('magenta')
    escape_sequence = TRIM(escape_sequence) // ';35'
  CASE ('cyan')
    escape_sequence = TRIM(escape_sequence) // ';36'
  CASE ('white')
    escape_sequence = TRIM(escape_sequence) // ';37'
  CASE ('background_black')
    escape_sequence = TRIM(escape_sequence) // ';40'
  CASE ('background_red')
    escape_sequence = TRIM(escape_sequence) // ';41'
  CASE ('background_green')
    escape_sequence = TRIM(escape_sequence) // ';42'
  CASE ('background_yellow')
    escape_sequence = TRIM(escape_sequence) // ';43'
  CASE ('background_blue')
    escape_sequence = TRIM(escape_sequence) // ';44'
  CASE ('background_magenta')
    escape_sequence = TRIM(escape_sequence) // ';45'
  CASE ('background_cyan')
    escape_sequence = TRIM(escape_sequence) // ';46'
  CASE ('background_white')
    escape_sequence = TRIM(escape_sequence) // ';47'
  END SELECT
END DO

! Append end of sequence marker
escape_sequence = TRIM(escape_sequence) // 'm'
END SUBROUTINE get_escape_sequence

!==================================================================================================================================
!> Sets output formatting to the supplied styles.
!> This function is copied from the fortran output library (foul). For the full version of this library see:
!>   http://foul.sourceforge.net
!==================================================================================================================================
SUBROUTINE set_formatting(style_string)
! MODULES
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
CHARACTER (LEN = *), INTENT(IN) :: style_string !< String describing which styles to set (separated by space).
                                                !< See get_escape_sequence for supported styles.
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=16) :: escape_sequence
CHARACTER         :: escape_sequence_array(16)
EQUIVALENCE         (escape_sequence, escape_sequence_array)
CHARACTER(LEN=16) :: format_string
CHARACTER(LEN=16) :: istring
INTEGER           :: i
!==================================================================================================================================
IF (use_escape_codes) THEN
  CALL get_escape_sequence(TRIM(style_string), escape_sequence)

  WRITE(istring,*) LEN_TRIM(escape_sequence)
  format_string = '(' // TRIM(istring) // 'A1)'

  SWRITE(UNIT_stdOut, TRIM(format_string), ADVANCE="NO") (escape_sequence_array(i), i = 1, LEN_TRIM(escape_sequence))
END IF
END SUBROUTINE set_formatting

!==================================================================================================================================
!> Resets output formatting to normal.
!> This function is copied from the fortran output library (foul). For the full version of this library see:
!>   http://foul.sourceforge.net
!==================================================================================================================================
SUBROUTINE clear_formatting()
IMPLICIT NONE
!==================================================================================================================================
IF (use_escape_codes) THEN
  ! Clear all previously set styles
  SWRITE(UNIT_stdOut, '(3A1)', ADVANCE="NO") (/ CHAR(27), '[', 'm' /)
END IF
END SUBROUTINE clear_formatting

END MODULE MOD_StringTools
