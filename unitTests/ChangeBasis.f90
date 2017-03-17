#include "flexi.h"

!==================================================================================================================================
!> Unit test 'ChangeBasisUnitTest'
!> Test the routine: 'ChangeBasis', from module: 'ChangeBasis'.
!> Compare against precomputed and stored values.
!==================================================================================================================================
PROGRAM ChangeBasisUnitTest
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis,       ONLY: ChangeBasis2D,ChangeBasis3D,ChangeBasis3D_XYZ
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: nVar=3, NIn=4, NOut=5, nElems=6
INTEGER                        :: p,q,z,iVar,i,j,k,iElem
REAL                           :: Vdm(0:NOut,0:NIn)
REAL                           :: Vdm2(0:NOut,0:NIn)
REAL                           :: Vdm3(0:NOut,0:NIn)
REAL                           :: UIn(nVar,0:NIn,0:NIn,0:NIn,nElems)
REAL                           :: UOut    (nVar,0:NOut,0:NOut,0:NOut,nElems,1:4)
REAL                           :: UOut_ref(nVar,0:NOut,0:NOut,0:NOut,nElems,1:4)
REAL                           :: UIn2D(nVar,0:NIn,0:NIn,nElems)
REAL                           :: UOut2D    (nVar,0:NOut,0:NOut,nElems,1:3)
REAL                           :: UOut2D_ref(nVar,0:NOut,0:NOut,nElems,1:3)
INTEGER                        :: nArgs
CHARACTER(LEN=*),PARAMETER     :: BinaryString='ChangeBasis.bin'
LOGICAL                        :: binaryExists,doGenerateReference=.FALSE.,equal
CHARACTER(LEN=255)             :: argument
!==================================================================================================================================
! Check for command line arguments to generate the reference solution
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs.GT.0) THEN
  CALL GET_COMMAND_ARGUMENT(1,argument)
  IF (argument.EQ.TRIM('--generate-reference')) THEN
    doGenerateReference = .TRUE.
  ELSE
    WRITE(*,*) 'ERROR - Unknown command line argument.'
    STOP -1
  END IF
END IF

UOut = 0.
UOut2D = 0.
! Initialize VDMs with arbitrary data
z=1
DO q=0,NIn; DO p=0,NOut
  z=z+1
  Vdm(p,q) = 1. + 1./z
  Vdm2(p,q) = 1. - 1./z
  Vdm3(p,q) = 1. + 0.1*z
END DO; END DO ! p,q=0,NIn
z=1
! Fill Uin
DO iElem=1,nElems
  DO k=0,NIn; DO j=0,NIn; DO i=0,NIn
    DO iVar=1,nVar
      z=z+1
      UIn(iVar,i,j,k,iElem) = 1.0 + 0.1*z
    END DO
  END DO; END DO; END DO! i,j,k=0,NIn
END DO
 
z=1
! Fill Uout for mode with add to Uout 
DO iElem=1,nElems
  DO k=0,NOut; DO j=0,NOut; DO i=0,NOut
    DO iVar=1,nVar
      z=z+1
      UOut(iVar,i,j,k,iElem,2) = -1.0 - 0.2*z
    END DO
  END DO; END DO; END DO! i,j,k=0,NIn
END DO

! Check ChangeBasis3D_Mult
CALL ChangeBasis3D(nVar,nElems,NIn,NOut,Vdm,UIn,UOut(:,:,:,:,:,1),addToOutput=.FALSE.)

! Check ChangeBasis3D_Mult with add to UOut 
CALL ChangeBasis3D(nVar,nElems,NIn,NOut,Vdm,UIn,UOut(:,:,:,:,:,2),addToOutput=.TRUE.)

! Check ChangeBasis3D_Single
CALL ChangeBasis3D(nVar,NIn,NOut,Vdm,UIn(:,:,:,:,1),UOut(:,:,:,:,1,3))

! Check ChangeBasis3D_XYZ
CALL ChangeBasis3D_XYZ(nVar,NIn,NOut,Vdm,Vdm2,Vdm3,UIn(:,:,:,:,1),UOut(:,:,:,:,1,4))

UIn2D = UIn(:,:,:,0,:)
UOut2D = UOut(:,:,:,0,:,1:3)
! Check ChangeBasis2D_Single
CALL ChangeBasis2D(nVar,NIn,NOut,Vdm,UIn2D(:,:,:,1),UOut2D(:,:,:,1,1))

IF (doGenerateReference) THEN
  ! Save the calculated solution to a binary file for later comparison
  OPEN(UNIT = 10, STATUS='replace',FILE=TRIM(BinaryString),FORM='unformatted')  ! replace an existing file or create a new one
  WRITE(10) UOut,UOut2D
  CLOSE(10) ! close the file
  WRITE(*,*) 'Saved reference to file ',BinaryString
ELSE
  ! Check if binary results file exists
  INQUIRE(FILE=TRIM(BinaryString),EXIST=binaryExists)
  
  IF (binaryExists) THEN
    ! Read the reference solution
    OPEN(UNIT = 10, STATUS='old',FILE=TRIM(BinaryString),FORM='unformatted')  ! open an existing file
    READ(10) UOut_ref,UOut2D_ref
    CLOSE(10) ! close the file
    equal =  .TRUE.
    DO z=1,4
      DO iElem=1,nElems
        DO k=0,NIn; DO j=0,NIn; DO i=0,NIn
          DO iVar=1,nVar
            equal = ALMOSTEQUALABSORREL(UOut(iVar,i,j,k,iElem,z),UOut_ref(iVar,i,j,k,iElem,z),50.*PP_RealTolerance) .AND. equal
          END DO
        END DO; END DO; END DO
      END DO
    END DO
    DO z=1,3
      DO iElem=1,nElems
        DO j=0,NIn; DO i=0,NIn
          DO iVar=1,nVar
            equal = ALMOSTEQUALABSORREL(UOut2D(iVar,i,j,iElem,z),UOut2D_ref(iVar,i,j,iElem,z),50.*PP_RealTolerance) .AND. equal
          END DO
        END DO; END DO
      END DO
    END DO
    IF (.NOT.equal) THEN
      WRITE(*,*) 'ERROR - Calculated ChangeBasis deviate from reference.'
      STOP -1
    ELSE
      WRITE(*,*) 'Compared ChangeBasis against stored data -- SUCCESSFUL.'
    END IF
  ELSE
    WRITE(*,*) 'ERROR - No reference solution has been found.'
    STOP -1
  END IF
END IF

END PROGRAM ChangeBasisUnitTest
