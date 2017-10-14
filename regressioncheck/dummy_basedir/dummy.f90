PROGRAM DUMMY_FLEXI
  INTEGER :: i,stat
  CHARACTER(len=255) :: arg
  CHARACTER(len=255) :: line
WRITE(*,'(80("="))')
WRITE(*,*) "   DUMMY FLEXI!!!"
WRITE(*,'(80("-"))')

DO i = 1, iargc()
  CALL getarg(i, arg)
  WRITE (*,'(A,I2,A,A)') "COMMANDLINE ARGUMENT", i, " = ", arg
END DO

IF (iargc().GT.0) THEN
  CALL getarg(1, arg)
  WRITE (*,*) "READING FROM FILE: ", TRIM(arg)

  open(unit = 100, file = TRIM(arg), status = 'old', action = 'read', ACCESS = 'SEQUENTIAL', iostat=stat)
  IF(stat.NE.0)THEN
    WRITE(*,*) "Could not open ini file."
    stop
  END IF
  DO
    READ(100,"(A)",IOSTAT=stat)line
    WRITE (*,*) line
    IF(stat.NE.0) EXIT
  END DO
  close(100)
END IF
WRITE(*,'(80("-"))')
WRITE(*,*) "   FINISH: DUMMY FLEXI!!!"
WRITE(*,'(80("="))')
END PROGRAM DUMMY_FLEXI
