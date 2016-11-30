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
#include "flexi.h"

!==================================================================================================================================
!> Initializes HDF5 IO and sets HDF-MPI parameters, opens ans closes files.
!==================================================================================================================================
MODULE MOD_IO_HDF5
! MODULES
USE HDF5
USE MOD_Globals,ONLY: iError
IMPLICIT NONE

ABSTRACT INTERFACE
  SUBROUTINE EvalElemInt(ElemData)
  USE MOD_Mesh_Vars,ONLY:nElems
  REAL,INTENT(OUT) :: ElemData(nElems)
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE EvalFieldInt(FieldData)
  REAL,INTENT(OUT) :: FieldData(:,:,:,:,:)
  END SUBROUTINE
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                  :: gatheredWrite       !< flag whether every process should output data or data should first be gathered
                                                !< on IO processes first
INTEGER(HID_T)           :: File_ID             !< file which is currently opened
INTEGER(HSIZE_T),POINTER :: HSize(:)            !< HDF5 array size (temporary variable)
INTEGER                  :: nDims               !< 
INTEGER                  :: MPIInfo             !< hardware / storage specific / file system MPI parameters to pass to HDF5
                                                !< for optimized performance on specific systems

!> Type containing pointers to data to be written to HDF5 in an element-wise scalar fashion.
!> Alternatively a function pointer can be specified providing the desired data.
!> Only one of the pointers may be associated.
TYPE tElementOut
  CHARACTER(LEN=255)         :: VarName                !< variable name
  REAL,POINTER               :: RealArray(:) => NULL()
  REAL,POINTER               :: RealScalar   => NULL()
  INTEGER,POINTER            :: IntArray(:)  => NULL()
  INTEGER,POINTER            :: IntScalar    => NULL()
  PROCEDURE(EvalElemInt),POINTER,NOPASS :: eval  => NULL()
  TYPE(tElementOut),POINTER  :: next         => NULL() !< next list item
END TYPE

!> Type containing pointers to nodal data to be written to HDF5 in a per node fashion.
!> Alternatively a function pointer can be specified providing the desired data.
!> Only one of the pointers may be associated.
TYPE tFieldOut
  INTEGER                    :: nVal(4)                !< size of array (5th dim is nElems)
  CHARACTER(LEN=255)         :: DataSetName            !< name of dataset to be created
  CHARACTER(LEN=255),ALLOCATABLE :: VarNames(:)        !< variable names in dataset
  REAL,POINTER               :: RealArray(:,:,:,:,:) => NULL()
  PROCEDURE(EvalFieldInt),POINTER,NOPASS :: eval     => NULL()
  TYPE(tFieldOut),POINTER    :: next         => NULL() !< next list item
END TYPE

TYPE(tElementOut),POINTER    :: ElementOut   => NULL() !< linked list of output pointers
TYPE(tFieldOut),POINTER      :: FieldOut     => NULL() !< linked list of output pointers


INTERFACE InitIOHDF5
  MODULE PROCEDURE InitIOHDF5
END INTERFACE

INTERFACE OpenDataFile
  MODULE PROCEDURE OpenDataFile
END INTERFACE

INTERFACE CloseDataFile
  MODULE PROCEDURE CloseDataFile
END INTERFACE

INTERFACE AddToElemData
  MODULE PROCEDURE AddToElemData
END INTERFACE

INTERFACE AddToFieldData
  MODULE PROCEDURE AddToFieldData
END INTERFACE

!==================================================================================================================================

PUBLIC::DefineParametersIO_HDF5,InitIOHDF5,OpenDataFile,CloseDataFile
PUBLIC::AddToElemData,AddToFieldData

CONTAINS


!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersIO_HDF5()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("IO_HDF5")
CALL prms%CreateLogicalOption('gatheredWrite', "Set true to activate gathered HDF5 IO for parallel computations. "//&
                                               "Only local group masters will write data after gathering from local slaves.",&
                                               '.FALSE.')
END SUBROUTINE DefineParametersIO_HDF5

!==================================================================================================================================
!> Initialize HDF5 IO
!==================================================================================================================================
SUBROUTINE InitIOHDF5()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,ONLY:GETLOGICAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
gatheredWrite=.FALSE.
IF(nLeaderProcs.LT.nProcessors) gatheredWrite=GETLOGICAL('gatheredWrite','.FALSE.')
#if USE_MPI
CALL MPI_Info_Create(MPIInfo, iError)

!normal case:
MPIInfo=MPI_INFO_NULL

! Large block IO extremely slow on Juqeen cluster (only available on IBM clusters)
!CALL MPI_Info_set(MPIInfo, "IBM_largeblock_io", "true", ierror)
#ifdef LUSTRE
CALL MPI_Info_Create(MPIInfo, iError)
! For lustre file system:
! Disables ROMIO's data-sieving
CALL MPI_Info_set(MPIInfo, "romio_ds_read", "disable",iError)
CALL MPI_Info_set(MPIInfo, "romio_ds_write","disable",iError)
! Enable ROMIO's collective buffering
CALL MPI_Info_set(MPIInfo, "romio_cb_read", "enable", iError)
CALL MPI_Info_set(MPIInfo, "romio_cb_write","enable", iError)
#endif
#endif /*USE_MPI*/
END SUBROUTINE InitIOHDF5


!==================================================================================================================================
!> Open HDF5 file and groups
!==================================================================================================================================
SUBROUTINE OpenDataFile(FileString,create,single,readOnly,communicatorOpt,userblockSize)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileString     !< filename to be opened
LOGICAL,INTENT(IN)            :: create         !< create file if it doesn't exist. Overwrited file if already present!
LOGICAL,INTENT(IN)            :: single         !< single=T : only one processor opens file, single=F : open/create collectively
LOGICAL,INTENT(IN)            :: readOnly       !< T : file is opened in read only mode, so file system timestamp remains unchanged
                                                !< F: file is open read/write mode
INTEGER,INTENT(IN),OPTIONAL   :: communicatorOpt !< only MPI and single=F: optional communicator to be used for collective access
                                                 !< default: MPI_COMM_WORLD
INTEGER,INTENT(IN),OPTIONAL   :: userblockSize  !< size of the file to be prepended to HDF5 file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                :: Plist_ID
INTEGER                       :: comm
INTEGER(HSIZE_T)              :: userblockSize_loc, tmp, tmp2
LOGICAL                       :: fileExists
!==================================================================================================================================
LOGWRITE(*,'(A)')'  OPEN HDF5 FILE "',TRIM(FileString),'" ...'

userblockSize_loc = 0
IF (PRESENT(userblockSize)) userblockSize_loc = userblockSize

! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)

! Setup file access property list with parallel I/O access (MPI) or with default property list.
IF(create)THEN
  CALL H5PCREATE_F(H5P_FILE_CREATE_F, Plist_ID, iError)
ELSE
  CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
END IF
#if USE_MPI
comm = MERGE(communicatorOpt,MPI_COMM_WORLD,PRESENT(communicatorOpt))
IF(.NOT.single)  CALL H5PSET_FAPL_MPIO_F(Plist_ID, comm, MPIInfo, iError)
#endif /*USE_MPI*/

! Open the file collectively.
IF(create)THEN
  IF (userblockSize_loc > 0) THEN
    tmp = userblockSize_loc/512
    IF (MOD(userblockSize_loc,512).GT.0) tmp = tmp+1
    tmp2 = 512*2**CEILING(LOG(REAL(tmp))/LOG(2.))
    CALL H5PSET_USERBLOCK_F(Plist_ID, tmp2, iError)
  END IF
  CALL H5FCREATE_F(TRIM(FileString), H5F_ACC_TRUNC_F, File_ID, iError, creation_prp = Plist_ID)
ELSE
  INQUIRE(FILE=TRIM(FileString),EXIST=fileExists)
  IF(.NOT.fileExists) CALL abort(__STAMP__,&
    'ERROR: Specified file '//TRIM(FileString)//' does not exist.')
  IF (readOnly) THEN
    CALL H5FOPEN_F(  TRIM(FileString), H5F_ACC_RDONLY_F,  File_ID, iError, access_prp = Plist_ID)
  ELSE 
    CALL H5FOPEN_F(  TRIM(FileString), H5F_ACC_RDWR_F,  File_ID, iError, access_prp = Plist_ID)
  END IF
END IF
IF(iError.NE.0) CALL abort(__STAMP__,&
  'ERROR: Could not open or create file '//TRIM(FileString))

CALL H5PCLOSE_F(Plist_ID, iError)
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE OpenDataFile



!==================================================================================================================================
!> Close HDF5 file and groups
!==================================================================================================================================
SUBROUTINE CloseDataFile()
! MODULES
USE MOD_Globals,ONLY:UNIT_logOut,Logging
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
LOGWRITE(*,'(A)')'  CLOSE HDF5 FILE...'
! Close file
CALL H5FCLOSE_F(File_ID, iError)
! Close FORTRAN predefined datatypes.
CALL H5CLOSE_F(iError)
File_ID=0
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE CloseDataFile

!==================================================================================================================================
!> Set pointers to element-wise scalar arrays which will be gathered and written out
!==================================================================================================================================
SUBROUTINE AddToElemData(VarName,RealArray,RealScalar,IntArray,IntScalar,Eval)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)        :: VarName
REAL,INTENT(IN),TARGET,OPTIONAL    :: RealArray(nElems)
REAL,INTENT(IN),TARGET,OPTIONAL    :: RealScalar
INTEGER,INTENT(IN),TARGET,OPTIONAL :: IntArray(nElems)
INTEGER,INTENT(IN),TARGET,OPTIONAL :: IntScalar
PROCEDURE(EvalElemInt),POINTER,OPTIONAL :: Eval
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElementOut),POINTER          :: eout
INTEGER                            :: nOpts
!==================================================================================================================================
IF(.NOT.ASSOCIATED(ElementOut))THEN 
  ! list is empty, create first entry
  ALLOCATE(ElementOut)
  eout=>ElementOut
ELSE
  eout=>ElementOut
  ! loop until last entry
  DO WHILE(ASSOCIATED(eout%next))
    eout=>eout%next
  END DO
  ! insert new entry
  ALLOCATE(eout%next)
  eout=>eout%next
ENDIF

! set varname and data pointer
NULLIFY(eout%next)
eout%VarName=VarName
nOpts=0
IF(PRESENT(RealArray))THEN
  eout%RealArray  => RealArray
  nOpts=nOpts+1
ENDIF
IF(PRESENT(RealScalar))THEN
  eout%RealScalar => RealScalar
  nOpts=nOpts+1
ENDIF
IF(PRESENT(IntArray))THEN
  eout%IntArray   => IntArray
  nOpts=nOpts+1
ENDIF
IF(PRESENT(IntScalar))THEN
  eout%IntScalar  => IntScalar
  nOpts=nOpts+1
ENDIF
IF(PRESENT(eval))THEN
  eout%eval       => Eval
  nOpts=nOpts+1
ENDIF
IF(nOpts.NE.1) CALL Abort(__STAMP__,&
  'More then one optional argument passed to AddToElemData.')
END SUBROUTINE AddToElemData


!==================================================================================================================================
!> Set pointers to node-wise arrays for output
!==================================================================================================================================
SUBROUTINE AddToFieldData(nVal,DataSetName,VarNames,RealArray,Eval)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: nVal(4)
CHARACTER(LEN=*),INTENT(IN)        :: DataSetName
CHARACTER(LEN=*),INTENT(IN)        :: VarNames(nVal(1))
REAL,INTENT(IN),TARGET,OPTIONAL    :: RealArray(nVal(1),nVal(2),nVal(3),nVal(4),nElems)
PROCEDURE(EvalFieldInt),POINTER,OPTIONAL :: Eval
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tFieldOut),POINTER            :: nout
INTEGER                            :: nOpts
!==================================================================================================================================
IF(.NOT.ASSOCIATED(FieldOut))THEN
  ! list is empty, create first entry
  ALLOCATE(FieldOut)
  nout=>FieldOut
ELSE
  ! loop until last entry
  nout=>FieldOut
  DO WHILE(ASSOCIATED(nout%next))
    nout=>nout%next
  END DO
  ! insert new entry
  ALLOCATE(nout%next)
  nout=>nout%next
ENDIF

! set varname, array size, and data pointer
NULLIFY(nout%next)
ALLOCATE(nout%Varnames(nVal(1)))
nout%VarNames=VarNames
nout%DataSetName=DataSetName
nout%nVal=nVal
nOpts=0
IF(PRESENT(RealArray))THEN
  nout%RealArray  => RealArray
  nOpts=nOpts+1
ENDIF
IF(PRESENT(eval))THEN
  nout%eval       => Eval
  nOpts=nOpts+1
ENDIF
IF(nOpts.NE.1) CALL Abort(__STAMP__,&
  'More then one optional argument passed to AddToFieldData.')
END SUBROUTINE AddToFieldData

END MODULE MOD_IO_HDF5
