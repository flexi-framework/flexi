#include "flexi.h"

!==================================================================================================================================
!> Contains the routines to 
!> - compare the LNorm norm
!> - compare Datasets of H5-Files
!> - reuired io-routines
!==================================================================================================================================
MODULE MOD_RegressionCheck_Compare
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE CompareResults
  MODULE PROCEDURE CompareResults
END INTERFACE

INTERFACE CompareConvergence
  MODULE PROCEDURE CompareConvergence
END INTERFACE

INTERFACE CompareNorm
  MODULE PROCEDURE CompareNorm
END INTERFACE

INTERFACE CompareDataSet
  MODULE PROCEDURE CompareDataSet
END INTERFACE

INTERFACE CompareRuntime
  MODULE PROCEDURE CompareRuntime
END INTERFACE

INTERFACE ReadNorm
  MODULE PROCEDURE ReadNorm
END INTERFACE

PUBLIC::CompareResults,CompareConvergence,CompareNorm,CompareDataSet,CompareRuntime,ReadNorm
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compare the results that were created by the binary execution
!==================================================================================================================================
SUBROUTINE CompareResults(iExample,iSubExample,MPIthreadsStr)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Tools,   ONLY: AddError
USE MOD_RegressionCheck_Vars,    ONLY: Examples
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iSubExample
CHARACTER(LEN=*),INTENT(IN)    :: MPIthreadsStr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: ReferenceNorm(:,:)                !> L2 and Linf norm of the executed example from a reference
                                                                    !> solution
INTEGER                        :: ErrorStatus                       !> Error-code of regressioncheck
!==================================================================================================================================
! -----------------------------------------------------------------------------------------------------------------------
! compare the results and write error messages for the current case
! -----------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A)',ADVANCE='no')  ' Comparing results...'
! check error norms  L2/LInf
ALLOCATE(ReferenceNorm(Examples(iExample)%nVar,2))
IF(Examples(iExample)%ReferenceNormFile.EQ.'')THEN
  ! constant value, should be zero no reference file given
  CALL CompareNorm(ErrorStatus,iExample,iSubExample)
ELSE
  ! read in reference and compare to reference solution
  CALL ReadNorm(iExample,ReferenceNorm)
  CALL CompareNorm(ErrorStatus,iExample,iSubExample,ReferenceNorm)
END IF
DEALLOCATE(ReferenceNorm)
IF(ErrorStatus.EQ.1)THEN
  SWRITE(UNIT_stdOut,'(A)')   ' Error-norm mismatched! Example failed! '
  SWRITE(UNIT_stdOut,'(A)')   ' For more information: '
  SWRITE(UNIT_stdOut,'(A,A)') ' Out-file: ', TRIM(Examples(iExample)%PATH)//'std.out'
  SWRITE(UNIT_stdOut,'(A,A)') ' Errorfile: ', TRIM(Examples(iExample)%PATH)//'err.out'
  CALL AddError(MPIthreadsStr,'Mismatch of error norms',iExample,iSubExample,ErrorStatus=1,ErrorCode=3)
END IF

! ConvergenceTest
IF(Examples(iExample)%ConvergenceTest)THEN
  IF(iSubExample.EQ.MAX(1,Examples(iExample)%SubExampleNumber))THEN ! after subexample 
    ! the subexample must be executed with "N" or "MeshFile": check if the convergence was successful
    CALL CompareConvergence(iExample)
    IF(Examples(iExample)%ErrorStatus.EQ.3)THEN
      CALL AddError(MPIthreadsStr,'Mismatch Order of '//TRIM(Examples(iExample)%ConvergenceTestType)&
                                                      //'-Convergence',iExample,iSubExample,ErrorStatus=3,ErrorCode=3)
    END IF
  END IF
END IF

! diff h5 file
IF(Examples(iExample)%ReferenceStateFile.NE.'')THEN
  CALL CompareDataSet(iExample)
  IF(Examples(iExample)%ErrorStatus.EQ.3)THEN
    CALL AddError(MPIthreadsStr,'Mismatch in HDF5-files. Datasets are unequal',iExample,iSubExample,ErrorStatus=3,ErrorCode=4)
    !SWRITE(UNIT_stdOut,'(A)')  ' Mismatch in HDF5-files'
  END IF
END IF

! Integrate over line
IF(Examples(iExample)%IntegrateLine)THEN
  CALL IntegrateLine(ErrorStatus,iExample)
  IF(Examples(iExample)%ErrorStatus.EQ.5)THEN
    CALL AddError(MPIthreadsStr,'Mismatch in LineIntegral',iExample,iSubExample,ErrorStatus=5,ErrorCode=5)
  END IF
END IF

! read a single row from a file and compare each entry
IF(Examples(iExample)%CompareDatafileRow)THEN
  CALL CompareDatafileRow(ErrorStatus,iExample)
  IF(Examples(iExample)%ErrorStatus.EQ.5)THEN
    CALL AddError(MPIthreadsStr,'Mismatch in CompareDatafileRow',iExample,iSubExample,ErrorStatus=5,ErrorCode=5)
  END IF
END IF

! successful execution and comparison
IF(Examples(iExample)%ErrorStatus.EQ.0)THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Example successful! '
END IF

END SUBROUTINE CompareResults


!==================================================================================================================================
!> Compare the results that were created by the binary execution
!==================================================================================================================================
SUBROUTINE CompareConvergence(iExample)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,    ONLY: Examples
USE MOD_RegressionCheck_tools,   ONLY: str2int,CalcOrder
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSTATUS
INTEGER                        :: I,J
INTEGER                        :: NumberOfCellsInteger
INTEGER                        :: iSubExample,p
REAL,ALLOCATABLE               :: Order(:,:),OrderAveraged(:)
INTEGER,ALLOCATABLE            :: OrderIncrease(:,:)
LOGICAL,ALLOCATABLE            :: OrderReached(:)
REAL                           :: DummyReal
LOGICAL                        :: DoDebugOutput
!==================================================================================================================================
DoDebugOutput=.TRUE. ! change to ".TRUE." if problems with this routine occur for info printed to screen
SWRITE(UNIT_stdOut,'(A)')''
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' ConvergenceTest: '
! 
IF(DoDebugOutput)THEN
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)')' L2 Error for nVar=[',Examples(iExample)%nVar,&
                                 '] and SubExampleNumber=[',Examples(iExample)%SubExampleNumber,']'
  SWRITE(UNIT_stdOut,'(A5)', ADVANCE="NO") ''
  DO J=1,Examples(iExample)%nVar
    SWRITE(UNIT_stdOut, '(A10,I3,A1)',ADVANCE="NO") 'nVar=[',J,']'
  END DO
  SWRITE(UNIT_stdOut,'(A)')''
  DO I=1,Examples(iExample)%SubExampleNumber
    SWRITE(UNIT_stdOut,'(I5)', ADVANCE="NO") I
    DO J=1,Examples(iExample)%nVar
        SWRITE(UNIT_stdOut, '(E14.6)',ADVANCE="NO") Examples(iExample)%ConvergenceTestError(I,J)
    END DO
    SWRITE(UNIT_stdOut,'(A)')''
  END DO
  SWRITE(UNIT_stdOut,'(A)')''
END IF

! Calculate the approximate distance between the DG DOF
! -----------------------------------------------------------------------------------------------------------------------
! p-convergence
IF(TRIM(Examples(iExample)%ConvergenceTestType).EQ.'p')THEN
  ! for p-convergence, the number of cells is constant: convert type from CHARACTER to INTEGER
  CALL str2int(ADJUSTL(TRIM(Examples(iExample)%NumberOfCellsStr(1))),NumberOfCellsInteger,iSTATUS) ! NumberOfCellsStr -> Int
  SWRITE(UNIT_stdOut,'(A,I4,A)')' Selecting p-convergence: Number of cells in one direction=[',NumberOfCellsInteger,'] (const.)'
  SWRITE(UNIT_stdOut,'(A)')''
  ! Calculate the approximate distance between the DG DOF
  DO iSubExample=1,Examples(iExample)%SubExampleNumber
    CALL str2int(ADJUSTL(TRIM(Examples(iExample)%SubExampleOption(iSubExample))),p,iSTATUS) ! SubExampleOption -> Int
    Examples(iExample)%ConvergenceTestGridSize(iSubExample)=&
    Examples(iExample)%ConvergenceTestDomainSize/(NumberOfCellsInteger*(p+1))
  END DO
! -----------------------------------------------------------------------------------------------------------------------
! h-convergence
ELSEIF(TRIM(Examples(iExample)%ConvergenceTestType).EQ.'h')THEN
  SWRITE(UNIT_stdOut,'(A,E14.6,A)')&
  ' Selecting h-convergence: Expected Order of Convergence = [',Examples(iExample)%ConvergenceTestValue,']'
  SWRITE(UNIT_stdOut,'(A)')''
  ! Calc the approximate distance between the DG DOF
  DO iSubExample=1,Examples(iExample)%SubExampleNumber
    CALL str2int(ADJUSTL(TRIM(Examples(iExample)%NumberOfCellsStr(iSubExample))) &
                 ,NumberOfCellsInteger,iSTATUS) ! sanity check if the number of threads is correct
    Examples(iExample)%ConvergenceTestGridSize(iSubExample)=&
    Examples(iExample)%ConvergenceTestDomainSize/(NumberOfCellsInteger*(Examples(iExample)%ConvergenceTestValue-1.+1.))
  END DO
END IF
! -----------------------------------------------------------------------------------------------------------------------

! Calculate ConvergenceTestGridSize (average spacing between DOF)
ALLOCATE(Order(Examples(iExample)%SubExampleNumber-1,Examples(iExample)%nVar))
DO J=1,Examples(iExample)%nVar
  DO I=1,Examples(iExample)%SubExampleNumber-1
    CALL CalcOrder(2,Examples(iExample)%ConvergenceTestGridSize(I:I+1),&
                     Examples(iExample)%ConvergenceTestError(   I:I+1,J),Order(I,J))
  END DO
END DO

! Check, if the Order of Convergece is increasing with increasing polynomial degree (only important for p-convergence)
ALLOCATE(OrderIncrease(Examples(iExample)%SubExampleNumber-2,Examples(iExample)%nVar))
DO J=1,Examples(iExample)%nVar
  DO I=1,Examples(iExample)%SubExampleNumber-2
    IF(Order(I,J).LT.Order(I+1,J))THEN ! increasing order
      OrderIncrease(I,J)=1
    ELSE ! non-increasing order
      OrderIncrease(I,J)=0
    END IF
  END DO
END DO

! Calculate the averged Order of Convergence (only important for h-convergence)
ALLOCATE(OrderAveraged(Examples(iExample)%nVar))
DO J=1,Examples(iExample)%nVar
  CALL CalcOrder(Examples(iExample)%SubExampleNumber,Examples(iExample)%ConvergenceTestGridSize(:),&
                                                     Examples(iExample)%ConvergenceTestError(:,J),OrderAveraged(J))
END DO

! Check the calculated Orders of convergence
ALLOCATE(OrderReached(Examples(iExample)%nVar))
OrderReached=.FALSE. ! default
! -----------------------------------------------------------------------------------------------------------------------
! p-convergence
IF(TRIM(Examples(iExample)%ConvergenceTestType).EQ.'p')THEN
  ! 75% of the calculated values for the order of convergece must be increasing with decreasing grid spacing
  DO J=1,Examples(iExample)%nVar
    IF(REAL(SUM(OrderIncrease(:,J)))/REAL(Examples(iExample)%SubExampleNumber-2).LT.0.75)THEN
      OrderReached(J)=.FALSE.
    ELSE
      OrderReached(J)=.TRUE.
    END IF
  END DO
! -----------------------------------------------------------------------------------------------------------------------
! h-convergence
ELSEIF(TRIM(Examples(iExample)%ConvergenceTestType).EQ.'h')THEN
  ! Check Order of Convergence versus the expected value and tolerance from input
  DO J=1,Examples(iExample)%nVar
    OrderReached(J)=ALMOSTEQUALRELATIVE( OrderAveraged(J),Examples(iExample)%ConvergenceTestValue,Examples(iExample)%ConvergenceTestTolerance )
     IF((OrderReached(J).EQV..FALSE.).AND.(OrderAveraged(J).GT.0.0))THEN
       !IntegralCompare=1
       SWRITE(UNIT_stdOut,'(A)')         ' CompareConvergence does not match! Error in computation!'
       SWRITE(UNIT_stdOut,'(A,E21.14)')  ' OrderAveraged(J)                        = ',OrderAveraged(J)
       SWRITE(UNIT_stdOut,'(A,E21.14)')  ' Examples(iExample)%ConvergenceTestValue = ',Examples(iExample)%ConvergenceTestValue
       SWRITE(UNIT_stdOut,'(A,E21.14)')  ' Tolerance                               = ',Examples(iExample)%ConvergenceTestTolerance
     END IF
  END DO
END IF


IF(DoDebugOutput)THEN
  ! Write average spacing between DOF
  SWRITE(UNIT_stdOut,'(A)')' ConvergenceTestGridSize (average spacing between DOF)'
  DO I=1,Examples(iExample)%SubExampleNumber
        write(*, '(I5,E14.6)') I,Examples(iExample)%ConvergenceTestGridSize(I)
  END DO
  ! Write Order of convergence
  SWRITE(UNIT_stdOut,'(A)')''
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)')' Order of convergence for nVar=[',Examples(iExample)%nVar,&
                                 '] and SubExampleNumber-1=[',Examples(iExample)%SubExampleNumber-1,']'
  SWRITE(UNIT_stdOut,'(A5)', ADVANCE="NO") ''
  DO J=1,Examples(iExample)%nVar
    SWRITE(UNIT_stdOut, '(A10,I3,A1)',ADVANCE="NO") 'nVar=[',J,']'
  END DO
  SWRITE(UNIT_stdOut,'(A)')''
  DO I=1,Examples(iExample)%SubExampleNumber-1
    SWRITE(UNIT_stdOut,'(I5)', ADVANCE="NO") I
    DO J=1,Examples(iExample)%nVar
      SWRITE(UNIT_stdOut,'(E14.6)',ADVANCE="NO") Order(I,J)
    END DO
    SWRITE(UNIT_stdOut,'(A)')''
  END DO
  ! Write averge convergence order
  SWRITE(UNIT_stdOut,'(A5)',ADVANCE="NO")'     '
  DO J=1,Examples(iExample)%nVar
    SWRITE(UNIT_stdOut, '(A14)',ADVANCE="NO") ' -------------'
  END DO
  SWRITE(UNIT_stdOut,'(A)')''
  SWRITE(UNIT_stdOut,'(A5)',ADVANCE="NO")'mean'
  DO J=1,Examples(iExample)%nVar
    SWRITE(UNIT_stdOut,'(E14.6)',ADVANCE="NO") OrderAveraged(J)
  END DO
  SWRITE(UNIT_stdOut,'(A)')''
  SWRITE(UNIT_stdOut,'(A)')''
  !    ! Write increasing order
  !    DO I=1,Examples(iExample)%SubExampleNumber-2
  !      SWRITE(UNIT_stdOut,'(I5)', ADVANCE="NO") I
  !      DO J=1,Examples(iExample)%nVar
  !        SWRITE(UNIT_stdOut,'(I14)',ADVANCE="NO") OrderIncrease(I,J)
  !      END DO
  !      SWRITE(UNIT_stdOut,'(A)')''
  !    END DO
  !    SWRITE(UNIT_stdOut,'(A)')''

  ! Write if order of convergence was reached (h- or p-convergence)
  SWRITE(UNIT_stdOut,'(A5)',ADVANCE="NO")'Check'
  DO J=1,Examples(iExample)%nVar
    SWRITE(UNIT_stdOut,'(L14)',ADVANCE="NO") OrderReached(J)
  END DO
  SWRITE(UNIT_stdOut,'(A)')''
  SWRITE(UNIT_stdOut,'(132("-"))')
END IF

! 50% of nVar Convergence tests must succeed
DummyReal=0.
DO J=1,Examples(iExample)%nVar
  IF(OrderReached(J))DummyReal=DummyReal+1.
END DO
IF(DummyReal/REAL(Examples(iExample)%nVar).LT.0.5)THEN
  Examples(iExample)%ErrorStatus=3
ELSE
  Examples(iExample)%ErrorStatus=0
END IF

END SUBROUTINE CompareConvergence


!==================================================================================================================================
!> Compare the runtime of an example  || fixed to a specific system
!> simply extract the regressioncheck settings from the parameter_reggie.ini
!> Not yet implemented!
!==================================================================================================================================
SUBROUTINE CompareRuntime()
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================


END SUBROUTINE CompareRuntime


!==================================================================================================================================
!> Compares the L2- and LInf-Norm of an example with a reference-norm. The reference-norm is given as a constant or from a 
!> reference simulation (previous simulation.)
!> To compare the norms, the std.out file of the simulation is read-in. The last L2- and LInf-norm in the std.out file are
!> compared to the reference.
!==================================================================================================================================
SUBROUTINE CompareNorm(LNormCompare,iExample,iSubExample,ReferenceNorm)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_StringTools,           ONLY: STRICMP
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: iExample,iSubExample
REAL,INTENT(IN),OPTIONAL     :: ReferenceNorm(Examples(iExample)%nVar,2)
INTEGER,INTENT(OUT)          :: LNormCompare
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iSTATUS2,iSTATUS,iVar
INTEGER                      :: ioUnit
CHARACTER(LEN=255)           :: FileName,temp1,temp2,temp3
LOGICAL                      :: ExistFile,L2Compare,LInfCompare
REAL                         :: LNorm(Examples(iExample)%nVar),L2(Examples(iExample)%nVar),LInf(Examples(iExample)%nVar)
REAL                         :: eps
!==================================================================================================================================

! get fileID and open file
FileName=TRIM(Examples(iExample)%PATH)//'std.out'
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' CompareNorm: no File found under ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:                  ','std.out'
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:                 ',ExistFile
  ERROR STOP '-1'
END IF

! find the last L2 and LInf norm the std.out file of the example
LNorm=-1.
L2Compare=.TRUE.
LInfCompare=.TRUE.
LNormCompare=1
OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS='OLD',IOSTAT=iSTATUS,ACTION='READ') 
DO 
  READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1!,temp2,LNorm(1),LNorm(2),LNorm(3),LNorm(4),LNorm(5)
  IF(iSTATUS.EQ.-1) EXIT ! End Of File (EOF) reached: exit the loop
  
  READ(temp1,*,IOSTAT=iSTATUS2) temp2,temp3,LNorm
  IF(STRICMP(temp2,'L_2')) THEN
    L2=LNorm
  END IF
  IF(STRICMP(temp2,'L_inf')) THEN
    LInf=LNorm
  END IF
END DO
! close the file
CLOSE(ioUnit)

! when NaN is encountered set the values to HUGE
IF(ANY(ISNAN(L2)))   L2  =HUGE(1.)
IF(ANY(ISNAN(LInf))) LInf=HUGE(1.)

! Save values for ConvergenceTest
IF(Examples(iExample)%ConvergenceTest)THEN
  Examples(iExample)%ConvergenceTestError(iSubExample,1:Examples(iExample)%nVar)=L2(1:Examples(iExample)%nVar)
END IF

! compare the retrieved norms from the std.out file
IF(PRESENT(ReferenceNorm))THEN ! use user-defined norm if present, else use 0.001*SQRT(PP_RealTolerance)
  ! compare with reference file
  IF(Examples(iExample)%ReferenceTolerance.GT.0.)THEN
    eps=Examples(iExample)%ReferenceTolerance
  ELSE
    eps=0.001*SQRT(PP_RealTolerance)
  END IF
  DO iVar=1,Examples(iExample)%nVar
    IF(.NOT.ALMOSTEQUALRELATIVE(L2(iVar),ReferenceNorm(iVar,1),eps))THEN
      L2Compare=.FALSE.
      SWRITE(UNIT_stdOut,'(A)') ''
      SWRITE(UNIT_stdOut,'(A,E21.14)')  ' L2Norm                =',L2(iVar)
      SWRITE(UNIT_stdOut,'(A,E21.14)')  ' ReferenceNorm(iVar,1) =',ReferenceNorm(iVar,1)
      SWRITE(UNIT_stdOut,'(A,E21.14)')  ' eps                   =',eps
      RETURN ! fail
    END IF
  END DO ! iVar=1,Examples(iExample)%nVar
  DO iVar=1,Examples(iExample)%nVar
    IF(.NOT.ALMOSTEQUALRELATIVE(LInf(iVar),ReferenceNorm(iVar,2),eps))THEN
      LInfCompare=.FALSE.
      SWRITE(UNIT_stdOut,'(A)') ''
      SWRITE(UNIT_stdOut,'(A,E21.14)')  ' LInfNorm              =',LInf(iVar)
      SWRITE(UNIT_stdOut,'(A,E21.14)')  ' ReferenceNorm(iVar,1) =',ReferenceNorm(iVar,2)
      SWRITE(UNIT_stdOut,'(A,E21.14)')  ' eps                   =',eps
      RETURN ! fail
    END IF
  END DO ! iVar=1,Examples(iExample)%nVar
ELSE ! use user-defined norm if present, else use 100.*PP_RealTolerance
  ! compare with single value
  IF(Examples(iExample)%ReferenceTolerance.GT.0.)THEN
    eps=Examples(iExample)%ReferenceTolerance
  ELSE
    eps=1000*PP_RealTolerance ! instead of 100, use 1000 because ketchesonrk4-20 with flexi failes here
  END IF
  IF(ANY(L2.GT.eps))THEN
    L2Compare=.FALSE.
    SWRITE(UNIT_stdOut,'(A)') ''
    SWRITE(UNIT_stdOut,'(A,E21.14)')  ' L2Norm                =',MAXVAL(L2)
    SWRITE(UNIT_stdOut,'(A,E21.14)')  ' eps                   =',eps
    RETURN ! fail
  END IF
  IF(ANY(LInf.GT.eps))THEN
    LInfCompare=.FALSE.
    SWRITE(UNIT_stdOut,'(A)') ''
    SWRITE(UNIT_stdOut,'(A,E21.14)')  ' LInfNorm              =',MAXVAL(LInf)
    SWRITE(UNIT_stdOut,'(A,E21.14)')  ' eps                   =',eps
    RETURN ! fail
  END IF
END IF
IF(L2Compare.AND.LInfCompare)LNormCompare=0

END SUBROUTINE CompareNorm


!==================================================================================================================================
!> Read in the error norms (L2,Linf) from a given reference computation and reference norm file
!> The reference files contains only the L2 and Linf norm for each variable of the reference computation.
!==================================================================================================================================
SUBROUTINE ReadNorm(iExample,ReferenceNorm)
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                   :: iExample
REAL,INTENT(OUT)                     :: ReferenceNorm(Examples(iExample)%nVar,2)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iSTATUS2,iSTATUS
INTEGER                              :: ioUnit
CHARACTER(LEN=255)                   :: FileName,temp1,temp2,temp3

LOGICAL                              :: ExistFile
REAL                                 :: LNorm(Examples(iExample)%nVar)
!==================================================================================================================================
! open file and read in
FileName=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%ReferenceNormFile)
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' ReadNorm: no File found under ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:                     ',TRIM(Examples(iExample)%ReferenceNormFile)
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:                    ',ExistFile
  ERROR STOP '-1'
END IF

! read in the norms
OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS='OLD',IOSTAT=iSTATUS,ACTION='READ') 
DO 
  READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1!,temp2,LNorm(1),LNorm(2),LNorm(3),LNorm(4),LNorm(5)
  IF(iSTATUS.EQ.-1) EXIT
  
  READ(temp1,*,IOSTAT=iSTATUS2) temp2,temp3,LNorm
  IF(TRIM(temp2).EQ.'L_2') THEN
    ReferenceNorm(1:Examples(iExample)%nVar,1)=LNorm
  END IF
  IF(TRIM(temp2).EQ.'L_Inf') THEN
    ReferenceNorm(1:Examples(iExample)%nVar,2)=LNorm
  END IF
END DO
CLOSE(ioUnit)

END SUBROUTINE ReadNorm


!==================================================================================================================================
!> Compares dataset of two different h5 files
!> It uses the reference and check-state-file information as well as the dataset information from the parameter_reggie.ini
!> The two datasets in the two different files are compared by a system-call to h5diff. If h5diff finds a difference, the
!> return status of the systemcall  is >0. Additionally, a absolute tolerance is used to allow for deviation of the datasets due to
!> different compilers.
!> This routine can compare all given datasets by their name, it is not restricted to the dg_solution. Thus it can be applied to 
!> all h5-files. Attention: This subroutine requires h5diff in the path of the used shell.
!==================================================================================================================================
SUBROUTINE CompareDataSet(iExample)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: DataSet
CHARACTER(LEN=255)             :: CheckedFileName
CHARACTER(LEN=255)             :: ReferenceNormFileName
CHARACTER(LEN=500)             :: SYSCOMMAND
CHARACTER(LEN=21)              :: tmpTol
INTEGER                        :: iSTATUS
LOGICAL                        :: ExistCheckedFile,ExistReferenceNormFile
!==================================================================================================================================
CheckedFilename  =TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%CheckedStateFile)
ReferenceNormFilename=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%ReferenceStateFile)
INQUIRE(File=CheckedFilename,EXIST=ExistCheckedFile)
IF(.NOT.ExistCheckedFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)')  ' h5diff: generated state file does not exist! need ',CheckedFilename
  Examples(iExample)%ErrorStatus=3
  RETURN
END IF
INQUIRE(File=ReferenceNormFilename,EXIST=ExistReferenceNormFile)
IF(.NOT.ExistReferenceNormFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)')  ' h5diff: reference state file does not exist! need ',ReferenceNormFilename
  Examples(iExample)%ErrorStatus=3
  RETURN
END IF

DataSet=TRIM(Examples(iExample)%ReferenceDataSetName)

WRITE(tmpTol,'(E21.14)') SQRT(PP_RealTolerance)
SYSCOMMAND=H5DIFF//' --delta='//TRIM(tmpTol)//' '//TRIM(ReferenceNormFileName)//' ' &
          //TRIM(CheckedFileName)//' /'//TRIM(DataSet)//' /'//TRIM(DataSet)
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' HDF5 Datasets do not match! Error in computation!'
  Examples(iExample)%ErrorStatus=3
END IF

END SUBROUTINE CompareDataSet


!==================================================================================================================================
!> Read column number data from a file and integrates the values numerically
!==================================================================================================================================
SUBROUTINE IntegrateLine(IntegralCompare,iExample)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RegressionCheck_Vars,  ONLY: Examples
USE MOD_RegressionCheck_tools, ONLY: str2real
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
INTEGER,INTENT(OUT)            :: IntegralCompare
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=1)               :: Delimiter
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=355)             :: temp1,temp2
INTEGER                        :: iSTATUS,ioUnit,LineNumbers,I,HeaderLines,j,IndMax,CurrentColumn,IndNum,MaxColumn!,K
INTEGER                        :: IndFirstA,IndLastA,IndFirstB,IndLastB,EOL,MaxRow
LOGICAL                        :: ExistFile,IndexNotFound,IntegralValuesAreEqual
REAL,ALLOCATABLE               :: Values(:,:),Q
!==================================================================================================================================
! check if output file with data for integration over line exists
Filename=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%IntegrateLineFile)
INQUIRE(File=Filename,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)')  ' IntegrateLine: reference state file does not exist! need ',TRIM(Filename)
  Examples(iExample)%ErrorStatus=5
  RETURN
ELSE
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS='OLD',IOSTAT=iSTATUS,ACTION='READ') 
END IF
! init parameters for reading the data file
HeaderLines=Examples(iExample)%IntegrateLineHeaderLines
!HeaderLines=1
Delimiter=ADJUSTL(TRIM(Examples(iExample)%IntegrateLineDelimiter))
MaxColumn=MAXVAL(Examples(iExample)%IntegrateLineRange)
IndMax  =LEN(temp1) ! complete string length
IndFirstA=1
IndLastA =IndMax
IndFirstB=1
IndLastB =IndMax
IndexNotFound=.TRUE.
CurrentColumn=0
EOL=0
DO I=1,2 ! read the file twice in Order to determine the array size
  LineNumbers=0
  DO 
    READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1 ! get first line assuming it is something like 'nVar= 5'
    IF(iSTATUS.EQ.-1) EXIT ! end of file (EOF) reached
    IF(INDEX(temp1,'!').GT.0)temp1=temp1(1:INDEX(temp1,'!')-1) ! if temp1 contains a '!', remove it and the following characters
    LineNumbers=LineNumbers+1
    IF(I.EQ.2)THEN ! read the data on second round reading the file (in first round, collect the file length by checking each line)
      IF(LineNumbers.GT.HeaderLines)THEN ! remove header lines
        IF(IndexNotFound)THEN
            !temp2=ADJUSTL(TRIM(temp1))  ! don't use ADJUSTL because it cuts away the spaces left to the first column
            temp2=TRIM(temp1)
          ! get index range
          DO J=1,MaxColumn
            IndNum=INDEX(TRIM(temp2),Delimiter)
            IF(IndNum.EQ.1)THEN ! still is same column!!!
              DO ! while IndNum.EQ.1
                IndNum=IndNum+INDEX(TRIM(temp2(IndNum+1:IndMax)),Delimiter)
                IF(IndNum.LE.0)EXIT ! not found - exit
                IF(IndNum.GT.1)EXIT
             END DO ! while
            END IF !IndNum.EQ.1
            IF(IndNum.GT.0)THEN
              CurrentColumn=CurrentColumn+IndNum
              ! first index
              IF(J.EQ.Examples(iExample)%IntegrateLineRange(1)-1)IndFirstA=CurrentColumn+1
              IF(J.EQ.Examples(iExample)%IntegrateLineRange(2)-1)IndFirstB=CurrentColumn+1
              ! last index
              IF(J.EQ.Examples(iExample)%IntegrateLineRange(1))IndLastA=CurrentColumn-1
              IF(J.EQ.Examples(iExample)%IntegrateLineRange(2))IndLastB=CurrentColumn-1
            ELSE
              EOL=EOL+1
              IF(EOL.GT.1)THEN
                SWRITE(UNIT_stdOut,'(A)')  ' IntegrateLines failed to read data! Error in computation!'
                SWRITE(UNIT_stdOut,'(A)')  ' The chosen column for line integration is larger than the available ones!'
                Examples(iExample)%ErrorStatus=5
                stop
                RETURN
              END IF!IF(EOL.GT.1)
            END IF!IF(IndNum.GT.0)
            temp2=TRIM(temp1(CurrentColumn+1:IndMax))
          END DO!J=1,MaxColumn
        IndexNotFound=.FALSE.
        END IF ! IndexNotFound
        CALL str2real(temp1(IndFirstA:IndLastA),Values(LineNumbers-HeaderLines,1),iSTATUS) 
        CALL str2real(temp1(IndFirstB:IndLastB),Values(LineNumbers-HeaderLines,2),iSTATUS) 
      END IF!IF(LineNumbers.GT.HeaderLines)
    END IF!IF(I.EQ.2)
  END DO ! DO [WHILE]
  IF(I.EQ.1)REWIND(ioUnit)
  IF(I.EQ.2)CLOSE(ioUnit)
  IF(I.EQ.1)MaxRow=LineNumbers-HeaderLines
  IF(I.EQ.1)ALLOCATE(Values(MaxRow,MaxColumn))
  If(I.EQ.1)Values=0.
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !print*,shape(Values)
         !IF(I.EQ.2)THEN
           !DO J=1,MaxRow
             !DO K=1,2
                 !write(*,'(E21.14,A)', ADVANCE = 'NO') Values(J,K),'  '
               !IF(K.EQ.2)print*,''
             !END DO
           !END DO
         !END IF
         !read*
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END DO

! integrate the values numerically
Q=0.
DO I=1,MaxRow-1
  Q=Q+(Values(I+1,1)-Values(I,1))*(Values(I+1,2)+Values(I,2))/2.
END DO

IntegralValuesAreEqual=ALMOSTEQUALRELATIVE( Q, Examples(iExample)%IntegrateLineValue, 5.e-2 )
IF(.NOT.IntegralValuesAreEqual)THEN
  IntegralCompare=1
  SWRITE(UNIT_stdOut,'(A)')         ' IntegrateLines do not match! Error in computation!'
  SWRITE(UNIT_stdOut,'(A,E21.14)')  ' IntegrateLineValue                    = ',Q
  SWRITE(UNIT_stdOut,'(A,E21.14)')  ' Examples(iExample)%IntegrateLineValue = ',Examples(iExample)%IntegrateLineValue
  SWRITE(UNIT_stdOut,'(A,E21.14)')  ' Tolerance                             = ',1.e-2!0.1*SQRT(PP_RealTolerance)
  !SWRITE(UNIT_stdOut,'(A,E21.14)')  ' 0.1*SQRT(PP_RealTolerance)            = ',0.1*SQRT(PP_RealTolerance)
  Examples(iExample)%ErrorStatus=5
ELSE
  IntegralCompare=0
END IF

END SUBROUTINE IntegrateLine


!==================================================================================================================================
!> Read column number data from a file and integrates the values numerically
!==================================================================================================================================
SUBROUTINE CompareDatafileRow(DataCompare,iExample)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RegressionCheck_Vars,  ONLY: Examples
USE MOD_RegressionCheck_tools, ONLY: str2real
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
INTEGER,INTENT(OUT)            :: DataCompare
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=1)               :: Delimiter
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: ColumnHeaders(:)
CHARACTER(LEN=10000)           :: temp1,temp2
INTEGER                        :: iSTATUS,ioUnit,LineNumbers,HeaderLines,j,k,ColumnNumber
LOGICAL                        :: ExistFile,ReadHeaderLine,RowFound
LOGICAL,ALLOCATABLE            :: ValuesAreEqual(:)
REAL,ALLOCATABLE               :: Values(:),ValuesRef(:)
INTEGER                        :: DimValues,DimValuesRef,DimColumnHeaders
!==================================================================================================================================
RowFound=.FALSE.
DO K=1,2 ! open the data and reference file
  ! check if output file with data for integration over line exists
  SELECT CASE(K)
  CASE(1) ! reference data file
    Filename=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%CompareDatafileRowRefFile)
    ReadHeaderLine=Examples(iExample)%CompareDatafileRowReadHeader
  CASE(2) ! newly created data file
    Filename=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%CompareDatafileRowFile)
    ReadHeaderLine=.FALSE.
  END SELECT
!print*,""
!print*,""
!print*,"Filename=",Filename
!read*
  INQUIRE(File=Filename,EXIST=ExistFile)
  IF(.NOT.ExistFile) THEN
    SWRITE(UNIT_stdOut,'(A,A)')  ' CompareDatafileRow: reference state file does not exist! need ',TRIM(Filename)
    Examples(iExample)%ErrorStatus=5
    RETURN
  ELSE
    OPEN(NEWUNIT=ioUnit,FILE=TRIM(FileName),STATUS='OLD',IOSTAT=iSTATUS,ACTION='READ') 
  END IF
  ! init parameters for reading the data file
  HeaderLines=Examples(iExample)%CompareDatafileRowHeaderLines
  IF(HeaderLines.GE.Examples(iExample)%CompareDatafileRowNumber)CALL abort(&
    __STAMP__&
    ,'CompareDatafileRow: The number of header lines exceeds the number of the row for comparison!')
  Delimiter=ADJUSTL(TRIM(Examples(iExample)%CompareDatafileRowDelimiter))
  LineNumbers=0
  DO 
    READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1
    IF(iSTATUS.EQ.-1) EXIT ! end of file (EOF) reached
    temp2=ADJUSTL(temp1)
    IF(INDEX(temp2,'!').GT.0)temp2=TRIM(temp2(1:INDEX(temp2,'!')-1)) ! if temp2 contains a '!', 
                                                                     ! remove it and the following characters
    LineNumbers=LineNumbers+1
    IF((LineNumbers.EQ.1).AND.(ReadHeaderLine))THEN
      ColumnNumber=0
!print*,"K=",K,"ColumnNumber=",ColumnNumber
!read*
      CALL GetColumns(temp2,Delimiter,ColumnString=ColumnHeaders,Column=ColumnNumber)
    ELSEIF(LineNumbers.EQ.Examples(iExample)%CompareDatafileRowNumber)THEN ! remove header lines
      RowFound=.TRUE.
      EXIT
    END IF!IF(LineNumbers.GT.HeaderLines)
  END DO ! DO [WHILE]

  IF(ADJUSTL(TRIM(temp2)).NE.'')THEN ! if string is not empty
    SELECT CASE(K)
    CASE(1) ! reference data file
!print*,"K=",K,"ColumnNumber=",ColumnNumber
!read*
      CALL GetColumns(temp2,Delimiter,ColumnReal=ValuesRef,Column=ColumnNumber)
    CASE(2) ! newly created data file
!print*,"K=",K,"ColumnNumber=",ColumnNumber
!read*
      CALL GetColumns(temp2,Delimiter,ColumnReal=Values   ,Column=ColumnNumber)
    END SELECT
  END IF
  CLOSE(ioUnit)
END DO ! K=1,2


DimValues=SIZE(Values)
DimValuesRef=SIZE(ValuesRef)
IF(DimValues.NE.DimValuesRef)THEN ! dimensions of ref values and data file values is different
  SWRITE(UNIT_stdOut,'(A,A)')&
    ' CompareDatafileRow: reference and datafile vector "ValuesRef" and "Values" have different dimensions!' 
  Examples(iExample)%ErrorStatus=5
  RETURN
END IF
IF(ALLOCATED(ColumnHeaders).EQV..TRUE.)THEN
  DimColumnHeaders=SIZE(ColumnHeaders)
  IF(DimValues.NE.DimColumnHeaders)THEN
    SWRITE(UNIT_stdOut,'(A,A)')&
      ' CompareDatafileRow: Header line vector "ColumnHeaders" (1st line in ref file) and "Values" have different dimensions!' 
    Examples(iExample)%ErrorStatus=5
    RETURN
  END IF
ELSE
  ALLOCATE(ColumnHeaders(1:DimValues))
  ColumnHeaders='no header found'
END IF
!print*,"done"
!print*,"ColumnNumber=",ColumnNumber
!read*
print*,""
IF(ColumnNumber.GT.0)THEN
  ALLOCATE(ValuesAreEqual(1:ColumnNumber))
  ValuesAreEqual=.FALSE.
  DO J=1,ColumnNumber
    ValuesAreEqual(J)=ALMOSTEQUALRELATIVE( Values(J), ValuesRef(J), Examples(iExample)%CompareDatafileRowTolerance)
    IF(ValuesAreEqual(J).EQV..FALSE.)THEN
      SWRITE(UNIT_stdOut,'(A)')         ' CompareDatafileRows mismatch: '//TRIM(ColumnHeaders(J))
      SWRITE(UNIT_stdOut,'(A,E24.17)')  ' Value in Refernece            = ',ValuesRef(J)
      SWRITE(UNIT_stdOut,'(A,E24.17)')  ' Value in data file            = ',Values(J)
      SWRITE(UNIT_stdOut,'(A,E24.17)')  ' Tolerance                     = ',Examples(iExample)%CompareDatafileRowTolerance
    END IF
  END DO
END IF
IF(ANY(.NOT.ValuesAreEqual))THEN
  DataCompare=1
  SWRITE(UNIT_stdOut,'(A)')         ' CompareDatafileRows do not match! Error in computation!'
  Examples(iExample)%ErrorStatus=5
ELSE
  DataCompare=0
END IF

!print*,"DataCompare=",DataCompare
!stop "STOPPPP"

END SUBROUTINE CompareDatafileRow

!==================================================================================================================================
!> Read column data from a supplied string InputString
!==================================================================================================================================
SUBROUTINE GetColumns(InputString,Delimiter,ColumnString,ColumnReal,Column)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RegressionCheck_tools, ONLY: str2real
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL                         :: Column
CHARACTER(LEN=*),INTENT(INOUT)                      :: InputString
CHARACTER(LEN=*),ALLOCATABLE,INTENT(INOUT),OPTIONAL :: ColumnString(:)
REAL,ALLOCATABLE,INTENT(INOUT),OPTIONAL             :: ColumnReal(:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: ColumnStringLocal(:)
CHARACTER(LEN=1)               :: Delimiter
INTEGER                        :: IndNumOld,ColumnNumber
INTEGER                        :: iSTATUS,j,IndNum
LOGICAL                        :: InquireColumns
!==================================================================================================================================
!print*,"InputString=",TRIM(InputString)
!print*,"Continue?"
!read*
IndNumOld=0
IF(PRESENT(Column))THEN
  IF(Column.GT.0)THEN
    ! the number of columns in the string is pre-defined
    ColumnNumber=Column
    InquireColumns=.FALSE.
  ELSE
    InquireColumns=.TRUE.
  END IF
ELSE
  InquireColumns=.TRUE.
END IF
IF(InquireColumns)THEN
  ColumnNumber=1
  ! inquire the number of columns in the string
  IndNum=0
  DO ! while IndNum.EQ.1
    IndNum=IndNum+INDEX(TRIM(InputString(IndNum+1:LEN(InputString))),Delimiter)
    IF(IndNum.LE.0)EXIT ! not found - exit
    IF(IndNumOld.EQ.IndNum)EXIT ! EOL reached - exit
    IndNumOld=IndNum
    ColumnNumber=ColumnNumber+1
  END DO ! while
END IF
!print*,"ColumnNumber=",ColumnNumber
IF(PRESENT(Column))Column=ColumnNumber
!print*,"Continue?"
!read*
IF(ADJUSTL(TRIM(InputString)).EQ.'')ColumnNumber=0 ! if InputString is empty, no ColumnNumber information can be extracted
IF(ColumnNumber.GT.0)THEN
  ALLOCATE(ColumnStringLocal(ColumnNumber))
  ColumnStringLocal='' ! default
  IndNum=0
  DO J=1,ColumnNumber
    IndNum=INDEX(TRIM(InputString(1:LEN(InputString))),Delimiter) ! for columns 1 to ColumnNumber-1
    IF(J.EQ.ColumnNumber)IndNum=LEN(InputString)-1          ! for the last ColumnNumber
    IF(IndNum.GT.0)THEN
      ColumnStringLocal(J)=ADJUSTL(TRIM(InputString(1:IndNum-1)))
      !print*,"ColumnStringLocal(",J,")=[",TRIM(ColumnStringLocal(J)),"]"
      InputString=InputString(IndNum+1:LEN(InputString))
    END IF
  END DO
END IF
!print*,"Continue?"
!read*

IF(PRESENT(ColumnString))THEN
  ALLOCATE(ColumnString(ColumnNumber))
  ColumnString='' ! default
  DO J=1,ColumnNumber
    ColumnString(J)=ADJUSTL(TRIM(ColumnStringLocal(J)))
    !print*,"ColumnString(",J,")=",ColumnString(J)
  END DO
!print*,"GetColumns DONE"
!read*
  RETURN
END IF
IF(PRESENT(ColumnReal))THEN
  ALLOCATE(ColumnReal(ColumnNumber))
  ColumnReal=0 ! default
  DO J=1,ColumnNumber
    CALL str2real(ColumnStringLocal(J),ColumnReal(J),iSTATUS) 
    !print*,"ColumnReal(",J,")=",ColumnReal(J)
  END DO
!print*,"GetColumns DONE"
!read*
  RETURN
END IF
!print*,"GetColumns DONE"
!read*
END SUBROUTINE GetColumns

END MODULE MOD_RegressionCheck_Compare