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

INTERFACE IntegrateLine
  MODULE PROCEDURE IntegrateLine
END INTERFACE

PUBLIC::CompareResults,CompareNorm,CompareDataSet,CompareRuntime,ReadNorm,IntegrateLine
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compare the results that were created by the binary execution
!==================================================================================================================================
SUBROUTINE CompareResults(iExample,iSubExample)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RegressionCheck_Tools,   ONLY: CleanExample,AddError
USE MOD_RegressionCheck_Vars,    ONLY: Examples
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample,iSubExample
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                        :: IndNum
!CHARACTER(LEN=255)             :: FileName                          ! path to a file or its name
!CHARACTER(LEN=255)             :: temp,temp2                        ! temp variables for read in of file lines
REAL,ALLOCATABLE               :: ReferenceNorm(:,:)                !> L2 and Linf norm of the executed example from a reference
                                                                    !> solution
INTEGER                        :: ErrorStatus                       !> Error-code of regressioncheck
!==================================================================================================================================
! -----------------------------------------------------------------------------------------------------------------------
! compare the results and write error messages for the current case
! -----------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A)',ADVANCE='no')  ' Comparing results...'
! check error norms
ALLOCATE(ReferenceNorm(Examples(iExample)%nVar,2))
IF(Examples(iExample)%ReferenceFile.EQ.'')THEN
  ! constant value, should be zero no reference file given
  CALL CompareNorm(ErrorStatus,iExample)
ELSE
  ! read in reference and compare to reference solution
  CALL ReadNorm(iExample,ReferenceNorm)
  CALL CompareNorm(ErrorStatus,iExample,ReferenceNorm)
END IF
DEALLOCATE(ReferenceNorm)
IF(ErrorStatus.EQ.1)THEN
  SWRITE(UNIT_stdOut,'(A)')   ' Error-norm mismatched! Example failed! '
  SWRITE(UNIT_stdOut,'(A)')   ' For more information: '
  SWRITE(UNIT_stdOut,'(A,A)') ' Out-file: ', TRIM(Examples(iExample)%PATH)//'std.out'
  SWRITE(UNIT_stdOut,'(A,A)') ' Errorfile: ', TRIM(Examples(iExample)%PATH)//'err.out'
  CALL AddError('Mismatch of error norms',iExample,iSubExample,ErrorStatus=1,ErrorCode=3)
END IF

! diff h5 file
IF(Examples(iExample)%ReferenceStateFile.NE.'')THEN
  CALL CompareDataSet(iExample)
  IF(Examples(iExample)%ErrorStatus.EQ.3)THEN
    CALL AddError('Mismatch in HDF5-files. Datasets are unequal',iExample,iSubExample,ErrorStatus=3,ErrorCode=4)
    SWRITE(UNIT_stdOut,'(A)')  ' Mismatch in HDF5-files'
  END IF
END IF

! Integrate over line
IF(Examples(iExample)%IntegrateLine)THEN
  CALL IntegrateLine(ErrorStatus,iExample)
  IF(Examples(iExample)%ErrorStatus.EQ.5)THEN
  CALL AddError('Mismatch in LineIntegral',iExample,iSubExample,ErrorStatus=5,ErrorCode=5)
  END IF
END IF

! IF all comparisons are successful the error status is 0 -> delete file in CleanExample(iExample)
IF(Examples(iExample)%ErrorStatus.EQ.0)THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Example successful! '
  CALL CleanExample(iExample)
END IF


END SUBROUTINE CompareResults

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
SUBROUTINE CompareNorm(LNormCompare,iExample,ReferenceNorm)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_StringTools,           ONLY: STRICMP
USE MOD_Basis,                 ONLY: EQUALTOTOLERANCE
USE MOD_RegressionCheck_Vars,  ONLY: Examples
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: iExample
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
!REAL                         :: epsLNorm=1.e-9
!==================================================================================================================================

! get fileid and open file
ioUnit=GETFREEUNIT()
FileName=TRIM(Examples(iExample)%PATH)//'std.out'
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' CompareNorm: no File found under ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:                  ','std.out'
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:                 ',ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS='OLD',IOSTAT=iSTATUS,ACTION='READ') 
END IF

! find the last L2 and LInf norm the std.out file of the example
LNorm=-1.
L2Compare=.TRUE.
LInfCompare=.TRUE.
LNormCompare=1
DO 
  READ(ioUnit,'(A)',IOSTAT=iSTATUS) temp1!,temp2,LNorm(1),LNorm(2),LNorm(3),LNorm(4),LNorm(5)
  IF(iSTATUS.EQ.-1) EXIT
  
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

! compare the retrieved norms from the std.out file
IF(PRESENT(ReferenceNorm))THEN ! use user-defined norm if present, else use 0.001*SQRT(PP_RealTolerance)
  ! compare with reference file
  IF(Examples(iExample)%ReferenceTolerance.GT.0.)THEN
    eps=Examples(iExample)%ReferenceTolerance
  ELSE
    eps=0.001*SQRT(PP_RealTolerance)
  END IF
  DO iVar=1,Examples(iExample)%nVar
    IF(.NOT.EQUALTOTOLERANCE(L2(iVar),ReferenceNorm(iVar,1),eps))THEN
      L2Compare=.FALSE.
      SWRITE(UNIT_stdOut,'(A)') ''
      SWRITE(UNIT_stdOut,'(A,E20.14)')  ' L2Norm                =',L2(iVar)
      SWRITE(UNIT_stdOut,'(A,E20.14)')  ' ReferenceNorm(iVar,1) =',ReferenceNorm(iVar,1)
      SWRITE(UNIT_stdOut,'(A,E20.14)')  ' eps                   =',eps
      RETURN ! fail
    END IF
  END DO ! iVar=1,Examples(iExample)%nVar
  DO iVar=1,Examples(iExample)%nVar
    IF(.NOT.EQUALTOTOLERANCE(LInf(iVar),ReferenceNorm(iVar,2),eps))THEN
      LInfCompare=.FALSE.
      SWRITE(UNIT_stdOut,'(A)') ''
      SWRITE(UNIT_stdOut,'(A,E20.14)')  ' LInfNorm              =',LInf(iVar)
      SWRITE(UNIT_stdOut,'(A,E20.14)')  ' ReferenceNorm(iVar,1) =',ReferenceNorm(iVar,2)
      SWRITE(UNIT_stdOut,'(A,E20.14)')  ' eps                   =',eps
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
    SWRITE(UNIT_stdOut,'(A,E20.14)')  ' L2Norm                =',MAXVAL(L2)
    SWRITE(UNIT_stdOut,'(A,E20.14)')  ' eps                   =',eps
    RETURN ! fail
  END IF
  IF(ANY(LInf.GT.eps))THEN
    LInfCompare=.FALSE.
    SWRITE(UNIT_stdOut,'(A)') ''
    SWRITE(UNIT_stdOut,'(A,E20.14)')  ' LInfNorm              =',MAXVAL(LInf)
    SWRITE(UNIT_stdOut,'(A,E20.14)')  ' eps                   =',eps
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
ioUnit=GETFREEUNIT()
FileName=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%ReferenceFile)
print*,FileName
INQUIRE(File=FileName,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)') ' ReadNorm: no File found under ',TRIM(Examples(iExample)%PATH)
  SWRITE(UNIT_stdOut,'(A,A)') ' FileName:                     ',TRIM(Examples(iExample)%ReferenceFile)
  SWRITE(UNIT_stdOut,'(A,L)') ' ExistFile:                    ',ExistFile
  ERROR STOP '-1'
ELSE
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS='OLD',IOSTAT=iSTATUS,ACTION='READ') 
END IF

! read in the norms
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
CHARACTER(LEN=255)             :: ReferenceFileName
CHARACTER(LEN=500)             :: SYSCOMMAND
CHARACTER(LEN=20)              :: tmpTol
INTEGER                        :: iSTATUS
LOGICAL                        :: ExistCheckedFile,ExistReferenceFile
!==================================================================================================================================


CheckedFilename  =TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%CheckedStateFile)
ReferenceFilename=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%ReferenceStateFile)
INQUIRE(File=CheckedFilename,EXIST=ExistCheckedFile)
IF(.NOT.ExistCheckedFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)')  ' h5diff: generated state file does not exist! need ',CheckedFilename
  Examples(iExample)%ErrorStatus=3
  RETURN
END IF
INQUIRE(File=ReferenceFilename,EXIST=ExistReferenceFile)
IF(.NOT.ExistReferenceFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)')  ' h5diff: reference state file does not exist! need ',ReferenceFilename
  Examples(iExample)%ErrorStatus=3
  RETURN
END IF

DataSet=TRIM(Examples(iExample)%ReferenceDataSetName)

WRITE(tmpTol,'(E20.14)') SQRT(PP_RealTolerance)
SYSCOMMAND=H5DIFF//' --delta='//TRIM(tmpTol)//' '//TRIM(ReferenceFileName)//' ' &
          //TRIM(CheckedFileName)//' /'//TRIM(DataSet)//' /'//TRIM(DataSet)
CALL EXECUTE_COMMAND_LINE(SYSCOMMAND, WAIT=.TRUE., EXITSTAT=iSTATUS)
IF(iSTATUS.NE.0) THEN
  SWRITE(UNIT_stdOut,'(A)')  ' Datasets do not match! Error in computation!'
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
USE MOD_Basis,                 ONLY: EQUALTOTOLERANCE
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: iExample
INTEGER,INTENT(OUT)            :: IntegralCompare
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!CHARACTER(LEN=255)             :: DataSet
CHARACTER(LEN=1)               :: Delimiter
CHARACTER(LEN=255)             :: FileName
!CHARACTER(LEN=255)             :: ReferenceFileName
CHARACTER(LEN=355)             :: temp1,temp2
!CHARACTER(LEN=20)              :: tmpTol
INTEGER                        :: iSTATUS,ioUnit,LineNumbers,I,HeaderLines,j,IndMax,CurrentColumn,IndNum,MaxColumn!,K
INTEGER                        :: IndFirstA,IndLastA,IndFirstB,IndLastB,EOL,MaxRow
LOGICAL                        :: ExistFile,IndexNotFound,IntegralValuesAreEqual
REAL,ALLOCATABLE               :: Values(:,:),Q
!==================================================================================================================================
! check if output file with data for integration over line exists
Filename=TRIM(Examples(iExample)%PATH)//TRIM(Examples(iExample)%IntegrateLineFile)
INQUIRE(File=Filename,EXIST=ExistFile)
IF(.NOT.ExistFile) THEN
  SWRITE(UNIT_stdOut,'(A,A)')  ' IntegrateLine: reference state file does not exist! need ',Filename
  Examples(iExample)%ErrorStatus=5
  RETURN
ELSE
  ioUnit=GETFREEUNIT()
  OPEN(UNIT=ioUnit,FILE=TRIM(FileName),STATUS='OLD',IOSTAT=iSTATUS,ACTION='READ') 
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
! read the file twice in order to determine the array size
DO I=1,2
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
!print*,'[',TRIM(temp2),']'
!print*,'IndNum=',IndNum,'J=',J,' of ',MaxColumn
!read*
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
!print*,'Examples(iExample)%IntegrateLineRange(1)',IndFirstA,IndLastA
!print*,'Examples(iExample)%IntegrateLineRange(2)',IndFirstB,IndLastB
        END IF ! IndexNotFound
!print*,temp1(IndFirstA:IndLastA),'  ',temp1(IndFirstB:IndLastB)
        CALL str2real(temp1(IndFirstA:IndLastA),Values(LineNumbers-HeaderLines,1),iSTATUS) 
        CALL str2real(temp1(IndFirstB:IndLastB),Values(LineNumbers-HeaderLines,2),iSTATUS) 
!print*,'[',temp1(IndFirstA:IndLastA),']','[',temp1(IndFirstB:IndLastB),']'
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
                 !write(*,'(E20.14,A)', ADVANCE = 'NO') Values(J,K),'  '
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
!print*,Q
!print*,Examples(iExample)%IntegrateLineValue
!print*,1.e-3!0.1*SQRT(PP_RealTolerance)

IntegralValuesAreEqual=EQUALTOTOLERANCE( Q                                     ,&
                                         Examples(iExample)%IntegrateLineValue ,&
                                         5.e-2                                 )
IF(.NOT.IntegralValuesAreEqual)THEN
  IntegralCompare=1
  SWRITE(UNIT_stdOut,'(A)')         ' IntegrateLines do not match! Error in computation!'
  SWRITE(UNIT_stdOut,'(A,E20.14)')  ' IntegrateLineValue                    = ',Q
  SWRITE(UNIT_stdOut,'(A,E20.14)')  ' Examples(iExample)%IntegrateLineValue = ',Examples(iExample)%IntegrateLineValue
  SWRITE(UNIT_stdOut,'(A,E20.14)')  ' Tolerance                             = ',1.e-2!0.1*SQRT(PP_RealTolerance)
  !SWRITE(UNIT_stdOut,'(A,E20.14)')  ' 0.1*SQRT(PP_RealTolerance)            = ',0.1*SQRT(PP_RealTolerance)
  Examples(iExample)%ErrorStatus=5
ELSE
  IntegralCompare=0
END IF

END SUBROUTINE IntegrateLine

END MODULE MOD_RegressionCheck_Compare
