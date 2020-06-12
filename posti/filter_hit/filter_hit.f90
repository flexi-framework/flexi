!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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

MODULE MOD_Filter_Hit
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ReadOldStateFile
  MODULE PROCEDURE ReadOldStateFile
END INTERFACE

INTERFACE WriteNewStateFile
  MODULE PROCEDURE WriteNewStateFile
END INTERFACE

INTERFACE FourierFilter
  MODULE PROCEDURE FourierFilter
END INTERFACE

PUBLIC:: ReadOldStateFile,WriteNewStateFile,FourierFilter

CONTAINS

!===================================================================================================================================
!===================================================================================================================================
SUBROUTINE InitFilterHit()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_HIT_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT BASIS...'

SWRITE(UNIT_stdOut,'(A)')' INIT BASIS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitFilterHit

!===================================================================================================================================
!> Open a state file, read the old state and store the information later needed to write a new state.
!===================================================================================================================================
SUBROUTINE ReadOldStateFile(StateFile)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_DG_Vars,         ONLY: U,Ut
USE MOD_HDF5_Input,      ONLY: OpenDataFile,CloseDataFile,ISVALIDHDF5FILE
USE MOD_HDF5_Input,      ONLY: ReadArray,ReadAttribute,GetDataProps,GetDataSize,DataSetExists
USE MOD_IO_HDF5,         ONLY: File_ID,nDims,HSize
USE MOD_Filter_Hit_Vars, ONLY: nVar_HDF5,N_HDF5,nElems_HDF5,nVarField_HDF5
USE MOD_Filter_Hit_Vars, ONLY: Time_HDF5,NodeType_HDF5,ProjectName_HDF5
USE MOD_Filter_Hit_Vars, ONLY: FieldDataExists,OverwriteMeshFile
USE MOD_ReadInTools,     ONLY: ExtractParameterFile
USE MOD_Output_Vars,     ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Output,          ONLY: insert_userblock
USE MOD_Mesh_Vars,       ONLY: MeshFile
USE MOD_StringTools,     ONLY: STRICMP,GetFileExtension
USE ISO_C_BINDING,       ONLY: C_NULL_CHAR
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)      :: StateFile !< State file to be read
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: userblockFound
CHARACTER(LEN=255)               :: prmfile=".parameter.ini"
!===================================================================================================================================
SWRITE(*,*) "READING SOLUTION FROM STATE FILE """,TRIM(StateFile), """"

! Get start index of file extension to check if it is a h5 file
IF (.NOT.STRICMP(GetFileExtension(StateFile), 'h5')) &
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid file extension of input state file!')

! Check if state file is a valid state
IF(.NOT.ISVALIDHDF5FILE(StateFile)) &
  CALL CollectiveStop(__STAMP__,'ERROR - Not a valid state file!')

! Open the data file
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Get data size of solution array
CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5)

! Allocate solution array in correct size
ALLOCATE(U( 1:nVar_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5,1:nElems_HDF5))

! Read the DG solution and store in UNew
CALL ReadArray('DG_Solution',5,&
               (/nVar_HDF5,N_HDF5+1,N_HDF5+1,N_HDF5+1,nElems_HDF5/),0,5,RealArray=U)

! Also read FieldData if present
CALL DatasetExists(File_ID,'FieldData',FieldDataExists)
IF (FieldDataExists) THEN
  CALL GetDataSize(File_ID,'FieldData',nDims,HSize)
  nVarField_HDF5 = INT(HSize(1))
  ALLOCATE(Ut(1:nVarField_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5,1:nElems_HDF5))
  CALL ReadArray('FieldData',5,&
                 (/nVarField_HDF5,N_HDF5+1,N_HDF5+1,N_HDF5+1,nElems_HDF5/),0,5,RealArray=Ut)
ENDIF

! Read the attributes from file
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_HDF5)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName_HDF5)
IF (.NOT. OverwriteMeshfile) THEN
  CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
END IF

! Extract parameter file from userblock (if found)
CALL ExtractParameterFile(StateFile,TRIM(prmfile),userblockFound)
! prepare userblock file
CALL insert_userblock(TRIM(UserBlockTmpFile)//C_NULL_CHAR,TRIM(prmfile)//C_NULL_CHAR)
INQUIRE(FILE=TRIM(UserBlockTmpFile),SIZE=userblock_total_len)

! Close the data file
CALL CloseDataFile()

SWRITE(*,*) "READING SOLUTION DONE!"
END SUBROUTINE ReadOldStateFile

!===================================================================================================================================
!> Write the new state file by calling the WriteState routine from FLEXI. All necessary variables must have been set correctly!
!===================================================================================================================================
SUBROUTINE WriteNewStateFile()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_DG_Vars,            ONLY: Ut
USE MOD_Mesh_Vars,          ONLY: MeshFile
USE MOD_IO_HDF5,            ONLY: tFieldOut,FieldOut,File_ID
USE MOD_IO_HDF5,            ONLY: OpenDataFile,CloseDataFile,AddToFieldData
USE MOD_HDF5_Output,        ONLY: WriteState,WriteAttribute
USE MOD_Interpolation_Vars, ONLY: NodeType
USE MOD_Output_Vars,        ONLY: NOut,ProjectName
USE MOD_Filter_Hit_Vars,    ONLY: N_HDF5,ProjectName_HDF5,Time_HDF5,NodeType_HDF5
USE MOD_Filter_Hit_Vars,    ONLY: FieldDataExists,nVarField_HDF5
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tFieldOut),POINTER        :: f
CHARACTER(LEN=255)             :: FileName
!===================================================================================================================================
! Set output variables according to current state file
NOut        = N_HDF5
ProjectName = TRIM(ProjectName_HDF5)//'_filtered'

! Add (filtered) FieldData to output if it exists
IF (FieldDataExists) THEN
  CALL AddToFieldData(FieldOut,(/nVarField_HDF5,N_HDF5+1,N_HDF5+1,N_HDF5+1/),'dudt',(/'dudt1','dudt2','dudt3','dudt4','dudt5'/),RealArray=Ut)
END IF

! Write new state file
CALL WriteState(TRIM(MeshFile),Time_HDF5,Time_HDF5,isErrorFile=.FALSE.)

! Annoyingly we have to overwrite the nodetype manually if it does not match the nodetype Flexi is compiled with
IF (NodeType.NE.NodeType_HDF5) THEN
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//'State',Time_HDF5))//'.h5'
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
  CALL WriteAttribute(File_ID,'NodeType',1,StrScalar=(/NodeType_HDF5/))
  CALL CloseDataFile()
END IF

! Cleanup FieldOut by deallocating last entry of list until only first element left
IF (FieldDataExists) THEN
  DO WHILE(ASSOCIATED(FieldOut%next))
    f=>FieldOut
    DO WHILE(ASSOCIATED(f%next))
      f=>f%next
    END DO
    f%next => NULL()
    DEALLOCATE(f) ! deallocate last entry of list
  END DO
  IF (ASSOCIATED(FieldOut)) DEALLOCATE(FieldOut) ! also deallocate first list entry
END IF

END SUBROUTINE WriteNewStateFile

!===================================================================================================================================
!> ?
!===================================================================================================================================
SUBROUTINE FourierFilter(nVar_In,U_In)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_Filter_Hit_Vars,    ONLY: N_FFT,Endw,Localk,N_Filter,Nc
USE MOD_FFT,                ONLY: Interpolate_DG2FFT,Interpolate_FFT2DG
USE MOD_FFT,                ONLY: ComputeFFT_R2C,ComputeFFT_C2R
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nVar_In
REAL,INTENT(INOUT)    :: U_in(1:nVar_In,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< elementwise DG solution from state file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER        :: i,j,k
REAL           :: U_Global(1:nVar_In,1:N_FFT  ,1:N_FFT  ,1:N_FFT  ) ! Real global DG solution
COMPLEX        :: U_FFT(   1:nVar_In,1:Endw(1),1:Endw(2),1:Endw(3)) ! Complex FFT solution
!===================================================================================================================================
! 1. Interpolate DG solution to equidistant points
CALL Interpolate_DG2FFT(nVar_In,U_in,U_Global)

! 2. Apply complex Fourier-Transform on solution from state file
CALL ComputeFFT_R2C(nVar_In,U_Global,U_FFT)

! 3. Fourier cutoff filter
IF (N_Filter.GT.-1) THEN
  DO k=1,Endw(3); DO j=1,Endw(2); DO i=1,Endw(1)
    IF(localk(4,i,j,k).GT.N_Filter) U_FFT(:,i,j,k) = 0.
  END DO; END DO; END DO
ELSE ! Nyquist filter
  DO k=1,Endw(3); DO j=1,Endw(2); DO i=1,Endw(1)
    IF(localk(4,i,j,k).GT.Nc) U_FFT(:,i,j,k) = 0.
  END DO; END DO; END DO
END IF

! 4. Apply inverse Fourier-Transform back into physical space
CALL ComputeFFT_C2R(nVar_In,U_FFT,U_Global)

! 5. Interpolate global solution at equidistant points back to DG points
CALL Interpolate_FFT2DG(nVar_In,U_Global,U_in)

END SUBROUTINE FourierFilter


!===================================================================================================================================
!===================================================================================================================================
SUBROUTINE FinalizeFilterHit()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Filter_Hit_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(LocalXYZ)
SDEALLOCATE(LocalK)

END SUBROUTINE FinalizeFilterHit

END MODULE MOD_Filter_Hit
