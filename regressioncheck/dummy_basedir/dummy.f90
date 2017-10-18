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

!#warning FV_ENABLED
#if FAIL
  FLEXI_FAIL=ON must fail
#endif

WRITE(*,*) "Sys date  :    18.10.2017 09:12:07                 "
WRITE(*,*) "CALCULATION TIME PER TSTEP/DOF: [ 4.62091E-06 sec ]"
WRITE(*,*) "Timestep  :    6.7836667E-02                       "
WRITE(*,*) "#Timesteps :    1.7700000E+02                      "
WRITE(*,*) "Sim time  :    1.2000000E+01                       "
#if FV_ENABLED
WRITE(*,*) "L_2       :    2.9689145E+00   2.9650478E+00   3.3355633E-02   1.6525407E-02   1.7087397E-02   1.0025174E-01   1.9094894E-03   3.9832497E+00"
WRITE(*,*) "L_inf     :    3.3300955E+00   3.2441740E+00   1.1862045E-01   5.7398722E-02   6.2248911E-02   1.7292781E-01   5.2882593E-03   4.3510114E+00"
#else
WRITE(*,*) "L_2       :    2.9689145E-99   2.9650478E-99   3.3355633E-99   1.6525407E-99   2.9650478E-99   3.3355633E-99   1.6525407E-99   1.6525407E-99"
WRITE(*,*) "L_inf     :    3.3300955E-99   3.2441740E-99   1.1862045E-99   5.7398722E-99   2.9650478E-99   3.3355633E-99   1.6525407E-99   1.6525407E-99"
#endif

END PROGRAM DUMMY_FLEXI
