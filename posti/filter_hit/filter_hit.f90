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

PUBLIC:: ReadOldStateFile,WriteNewStateFile

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
USE MOD_HDF5_Input,      ONLY: OpenDataFile,CloseDataFile,ReadArray,ReadAttribute,GetDataProps
USE MOD_IO_HDF5,         ONLY: File_ID
USE MOD_Filter_Hit_Vars, ONLY: nVar_HDF5,N_HDF5,nElems_HDF5,U_HDF5
USE MOD_Filter_Hit_Vars, ONLY: Time_HDF5,MeshFile_HDF5,NodeType_HDF5,ProjectName_HDF5
USE MOD_ReadInTools,     ONLY: ExtractParameterFile
USE MOD_Output_Vars,     ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Output,          ONLY: insert_userblock
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
! Open the data file
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! Get data size of solution array
CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5)

! Allocate solution array in correct size
ALLOCATE(U_HDF5(1:nVar_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5,1:nElems_HDF5))

! Read the DG solution and store in UNew
CALL ReadArray('DG_Solution',5,&
               (/nVar_HDF5,N_HDF5+1,N_HDF5+1,N_HDF5+1,nElems_HDF5/),0,5,RealArray=U_HDF5)

! Read the attributes from file
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile_HDF5)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_HDF5)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName_HDF5)

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
USE MOD_HDF5_Output,        ONLY: WriteState
USE MOD_Output_Vars,        ONLY: NOut,ProjectName
USE MOD_Filter_Hit_Vars,    ONLY: N_HDF5,ProjectName_HDF5,Time_HDF5,MeshFile_HDF5
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Set output variables according to current state file
NOut        = N_HDF5
ProjectName = TRIM(ProjectName_HDF5)//'_filtered'

! Write new state file
CALL WriteState(TRIM(MeshFile_HDF5),Time_HDF5,Time_HDF5,isErrorFile=.FALSE.)

END SUBROUTINE WriteNewStateFile

!===================================================================================================================================
!===================================================================================================================================
SUBROUTINE FinalizeFilterHit()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Filter_Hit_Vars
USE MOD_DG_Vars,       ONLY: U
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(U)
SDEALLOCATE(U_j)
SDEALLOCATE(U_k)
SDEALLOCATE(U_FFT)
IF(MPIRoot) THEN
  SDEALLOCATE(Uloc)
  SDEALLOCATE(LocalXYZ)
  SDEALLOCATE(LocalK)
  SDEALLOCATE(phat)
  SDEALLOCATE(fhat)
  SDEALLOCATE(F_vv)
END IF

END SUBROUTINE FinalizeFilterHit

END MODULE MOD_Filter_Hit
