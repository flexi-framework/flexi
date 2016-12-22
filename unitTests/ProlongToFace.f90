#include "flexi.h"

!==================================================================================================================================
!> Unit test 'ProlongToFaceUnitTest'
!> Test the routine: 'ProlongToFace', from module: 'ProlongToFace'.
!> Compare against precomputed and stored values.
!==================================================================================================================================
PROGRAM ProlongToFaceUnitTest
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Unittest_Vars
USE MOD_Unittest,           ONLY: ReadInReferenceElementData
USE MOD_ProlongToFaceCons,  ONLY: ProlongToFaceCons
USE MOD_Basis,              ONLY: EQUALTOTOLERANCE
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems,FV_Elems_master,FV_Elems_slave
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Uvol(0:NRef,0:NRef,0:NRefZ,nElemsRef)
REAL                           :: Uvol_nVar(PP_nVar,0:NRef,0:NRef,0:NRefZ,nElemsRef)
REAL                           :: Uface_master(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef),Uface_master_ref(0:NRef,0:NRefZ,1:nSidesRef)
REAL                           :: Uface_slave(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef)!,Uface_slave_ref(0:9,0:9,7:6)
#if FV_ENABLED
REAL                           :: FV_Uface_master(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef),FV_Uface_master_ref(0:NRef,0:NRefZ,1:nSidesRef)
REAL                           :: FV_Uface_slave(PP_nVar,0:NRef,0:NRefZ,1:nSidesRef)
#endif
INTEGER                        :: i,j,k,l,nArgs
CHARACTER(LEN=*),PARAMETER     :: BinaryUvolString='ProlongToFaceUvol.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,equal,doGenerateUvol=.FALSE.
CHARACTER(LEN=255)             :: BinaryString,argument
!==================================================================================================================================
! Set binary file name for different node types
#if (PP_NodeType==1)
#if PP_dim == 3
BinaryString='ProlongToFace_G3D.bin'
#else
BinaryString='ProlongToFace_G2D.bin'
#endif
#elif (PP_NodeType==2)
#if PP_dim == 3
BinaryString='ProlongToFace_GL3D.bin'
#else
BinaryString='ProlongToFace_GL2D.bin'
#endif
#endif

! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) THEN
  CALL GET_COMMAND_ARGUMENT(1,argument)
  IF (argument.EQ.TRIM('--generate-reference')) THEN
    doGenerateReference = .TRUE.
  ELSE IF (argument.EQ.TRIM('--generate-uvol')) THEN
    doGenerateUvol = .TRUE.
  ELSE
    WRITE(*,*) 'ERROR - Unknown command line argument.'
    STOP -1
  END IF
END IF

IF (doGenerateUvol) THEN
  ! Generate a random volume solution
  CALL RANDOM_NUMBER(Uvol) 
  ! Save the calculated volume solution to a binary file for later input
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryUvolString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Uvol
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference Uvol to file ',BinaryUvolString
ELSE
  ! Read in the random Uvol
  OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryUvolString),FORM='unformatted')  ! open an existing file
  READ(10) Uvol
  CLOSE(10) ! close the file
  ! Build Uvol on PP_nVar as expected by ProlongToFace
  DO i=1,PP_nVar
#if PP_dim == 3
  Uvol_nVar(i,:,:,:,:) = Uvol(:,:,:,:)
#else
  Uvol_nVar(i,:,:,0,:) = Uvol(:,:,0,:)
#endif
  END DO
END IF


! Read in data from single curved element
CALL ReadInReferenceElementData()

! Call ProlongToFace
#if FV_ENABLED
ALLOCATE(FV_Elems(1:1))
ALLOCATE(FV_Elems_master(1:nSidesRef))
ALLOCATE(FV_Elems_slave (1:nSidesRef))
FV_Elems = 0
FV_Elems_master = 0
FV_Elems_slave = 0
#endif
CALL ProlongToFaceCons(NRef,Uvol_nVar,Uface_master,Uface_slave,L_Minus,L_Plus,.FALSE.)
#if FV_ENABLED
FV_Elems = 1
FV_Elems_master = 1
FV_Elems_slave = 1
CALL ProlongToFaceCons(NRef,Uvol_nVar,FV_Uface_master,FV_Uface_slave,L_Minus,L_Plus,.FALSE.)
#endif


IF (doGenerateReference) THEN
  ! Save the calculated solution to a binary file for later comparison
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) Uface_master(1,:,:,:)!,Uface_slave(1,:,:,:)
#if FV_ENABLED
  WRITE(10) FV_Uface_master(1,:,:,:)
#endif
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference to file ',BinaryString
ELSE
  ! Check if binary results file exists
  INQUIRE(FILE=TRIM(BinaryString),EXIST=binaryExists)
  
  IF (binaryExists) THEN
    ! Read the reference solution
    OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryString),FORM='unformatted')  ! open an existing file
    READ(10) Uface_master_ref!,Uface_slave_ref
#if FV_ENABLED
    READ(10) FV_Uface_master_ref
#endif
    CLOSE(10) ! close the file
    ! Check if the computed and the reference solutions are within a given tolerance
    equal =  .TRUE.
    DO i=1,PP_nVar; DO j=0,NRef; DO k=0,NRefZ; DO l=1,nSidesRef
      equal = EQUALTOTOLERANCE(Uface_master(i,j,k,l),Uface_master_ref(j,k,l),100.*PP_RealTolerance) .AND. equal
#if FV_ENABLED
      equal = EQUALTOTOLERANCE(FV_Uface_master(i,j,k,l),FV_Uface_master_ref(j,k,l),100.*PP_RealTolerance) .AND. equal
#endif      
    END DO; END DO; END DO; END DO
    ! Plus sides not needed in single element case
    !DO i=1,PP_nVar; DO j=0,9; DO k=0,9; DO l=7,6
      !equal = EQUALTOTOLERANCE(Uface_slave(i,j,k,l),Uface_slave_ref(j,k,l),100.*PP_RealTolerance) .AND. equal
    !END DO; END DO; END DO; END DO
    IF (.NOT.equal) THEN
      WRITE(*,*) 'ERROR - Calculated prolonged values deviate from reference.'
      STOP -1
    ELSE
      WRITE(*,*) 'Compared prolonged values against stored data in ',TRIM(BinaryString),' -- SUCCESSFUL.'
    END IF
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    STOP -1
  END IF
END IF

END PROGRAM ProlongToFaceUnitTest
