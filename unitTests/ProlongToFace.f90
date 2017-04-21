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
USE MOD_ProlongToFaceCons,  ONLY: ProlongToFaceCons
! Modules needed to read in reference element
USE MOD_Mesh_Vars,          ONLY: nElems,sJ,FS2M,V2S,S2V,S2V2
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR, lastMPISide_MINE, nSides
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems,FV_Elems_master,FV_Elems_slave
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Uvol(0:9,0:9,0:9,1)
REAL                           :: Uvol_nVar(PP_nVar,0:9,0:9,0:9,1)
REAL                           :: Uface_master(PP_nVar,0:9,0:9,1:6),Uface_master_ref(0:9,0:9,1:6)
REAL                           :: Uface_slave(PP_nVar,0:9,0:9,7:6)!,Uface_slave_ref(0:9,0:9,7:6)
#if FV_ENABLED
REAL                           :: FV_Uface_master(PP_nVar,0:9,0:9,1:6),FV_Uface_master_ref(0:9,0:9,1:6)
REAL                           :: FV_Uface_slave(PP_nVar,0:9,0:9,7:6)
#endif
INTEGER                        :: i,j,k,l,nArgs
CHARACTER(LEN=*),PARAMETER     :: BinaryUvolString='ProlongToFaceUvol.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,equal,doGenerateUvol=.FALSE.
CHARACTER(LEN=255)             :: BinaryString,argument
!==================================================================================================================================
! Set file name for different node types
#if (PP_NodeType==1)
BinaryString='ProlongToFace_G.bin'
#elif (PP_NodeType==2)
BinaryString='ProlongToFace_GL.bin'
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
    Uvol_nVar(i,:,:,:,:) = Uvol(:,:,:,:)
  END DO
END IF


! Read in data from single curved element
ALLOCATE(SideToElem(1:5,1:6))
ALLOCATE(L_Minus(0:9))
ALLOCATE(L_Plus(0:9))
ALLOCATE(L_HatMinus(0:9))
ALLOCATE(L_HatPlus(0:9))
ALLOCATE(sJ(0:9,0:9,0:9,0:0,1:1))
ALLOCATE(FS2M(1:2,0:9,0:9,0:4))
ALLOCATE(V2S(1:3,0:9,0:9,0:9,0:4,1:6))
ALLOCATE(S2V(1:3,0:9,0:9,0:9,0:4,1:6))
ALLOCATE(S2V2(1:2,0:9,0:9,0:4,1:6))
OPEN(UNIT = 10, STATUS='old',FILE='UnittestElementData.bin',FORM='unformatted')  ! open an existing file
READ(10) nElems,SideToElem,firstMPISide_YOUR,lastMPISide_MINE,nSides,FS2M,V2S,S2V,S2V2,L_Minus,L_Plus,L_HatPlus,L_HatMinus,sJ
CLOSE(10) ! close the file


! Call ProlongToFace
#if FV_ENABLED
ALLOCATE(FV_Elems(1:1))
ALLOCATE(FV_Elems_master(1:6))
ALLOCATE(FV_Elems_slave (1:6))
FV_Elems = 0
FV_Elems_master = 0
FV_Elems_slave = 0
#endif
CALL ProlongToFaceCons(9,Uvol_nVar,Uface_master,Uface_slave,L_Minus,L_Plus,.FALSE.)
#if FV_ENABLED
FV_Elems = 1
FV_Elems_master = 1
FV_Elems_slave = 1
CALL ProlongToFaceCons(9,Uvol_nVar,FV_Uface_master,FV_Uface_slave,L_Minus,L_Plus,.FALSE.)
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
    DO i=1,PP_nVar; DO j=0,9; DO k=0,9; DO l=1,6
      equal = ALMOSTEQUALABSORREL(Uface_master(i,j,k,l),Uface_master_ref(j,k,l),100.*PP_RealTolerance) .AND. equal
#if FV_ENABLED
      equal = ALMOSTEQUALABSORREL(FV_Uface_master(i,j,k,l),FV_Uface_master_ref(j,k,l),100.*PP_RealTolerance) .AND. equal
#endif      
    END DO; END DO; END DO; END DO
    ! Plus sides not needed in single element case
    !DO i=1,PP_nVar; DO j=0,9; DO k=0,9; DO l=7,6
      !equal = ALMOSTEQUALABSORREL(Uface_slave(i,j,k,l),Uface_slave_ref(j,k,l),100.*PP_RealTolerance) .AND. equal
    !END DO; END DO; END DO; END DO
    IF (.NOT.equal) THEN
      WRITE(*,*) 'ERROR - Calculated prolonged values deviate from reference.'
      STOP -1
    END IF
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    STOP -1
  END IF
END IF

END PROGRAM ProlongToFaceUnitTest
