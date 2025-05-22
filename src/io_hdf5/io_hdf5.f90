!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  SUBROUTINE EvalElemInt(ElemData)
    ! MODULES
    USE MOD_Mesh_Vars,ONLY:nElems
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    REAL,INTENT(OUT) :: ElemData(nElems)
  END SUBROUTINE EvalElemInt
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE EvalFieldInt(FieldData)
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    ! INPUT / OUTPUT VARIABLES
    REAL,INTENT(OUT) :: FieldData(:,:,:,:,:)
  END SUBROUTINE EvalFieldInt
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                  :: gatheredWrite       !< flag whether every process should output data or data should first be gathered
                                                !< on IO processes first
LOGICAL                  :: output2D            !< Flag whether to use true 2D input/output or not
INTEGER(HID_T)           :: File_ID             !< file which is currently opened
INTEGER(HID_T)           :: Plist_File_ID       !< property list of file which is currently opened
INTEGER(HSIZE_T),POINTER :: HSize(:)            !< HDF5 array size (temporary variable)
INTEGER                  :: nDims               !< data size dimensions
INTEGER,PARAMETER        :: nLimit = INT(REAL(HUGE(INT(1,KIND=4)))/REAL(KIND(REAL(1))))
                                                !< Max number of entries in double arrays for HDF5 IO. Limit is computed as
                                                !< INT( Max. Bytes allowed by MPI (2GB per rank) / Size of single double entry)
#if USE_MPI
TYPE(MPI_Info)           :: MPIInfo             !< hardware / storage specific / file system MPI parameters to pass to HDF5
                                                !< for optimized performance on specific systems
#endif /*USE_MPI*/

!> Type containing pointers to data to be written to HDF5 in an element-wise scalar fashion.
!> Alternatively a function pointer can be specified providing the desired data.
!> Only one of the pointers may be associated.
TYPE tElementOut
  CHARACTER(LEN=255)                    :: VarName                    !< variable name
  REAL,POINTER                          :: RealArray(:) => NULL()
  REAL,POINTER                          :: RealScalar   => NULL()
  INTEGER,POINTER                       :: IntArray(:)  => NULL()
  INTEGER,POINTER                       :: IntScalar    => NULL()
  PROCEDURE(EvalElemInt),POINTER,NOPASS :: eval         => NULL()
  TYPE(tElementOut),POINTER             :: next         => NULL()     !< next list item
END TYPE tElementOut

!> Type containing pointers to nodal data to be written to HDF5 in a per node fashion.
!> Alternatively a function pointer can be specified providing the desired data.
!> Only one of the pointers may be associated.
TYPE tFieldOut
  INTEGER                                :: nVal(4)                        !< size of array (5th dim is nElems)
  CHARACTER(LEN=255)                     :: DataSetName                    !< name of dataset to be created
  CHARACTER(LEN=255),ALLOCATABLE         :: VarNames(:)                    !< variable names in dataset
  REAL,POINTER                           :: RealArray(:,:,:,:,:) => NULL()
  PROCEDURE(EvalFieldInt),POINTER,NOPASS :: eval                 => NULL()
  LOGICAL                                :: doSeparateOutput               !< If set, array will be written as seperate dataset,
                                                                           !< regardless of N
  TYPE(tFieldOut),POINTER                :: next                 => NULL() !< next list item
END TYPE tFieldOut

TYPE(tElementOut),POINTER    :: ElementOut   => NULL() !< linked list of output pointers
TYPE(tFieldOut),POINTER      :: FieldOut     => NULL() !< linked list of output pointers
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: DefineParametersIO_HDF5
PUBLIC:: InitIOHDF5
PUBLIC:: InitMPIInfo
PUBLIC:: OpenDataFile
PUBLIC:: CloseDataFile
PUBLIC:: AddToElemData
PUBLIC:: AddToFieldData
PUBLIC:: RemoveFromElemData
PUBLIC:: RemoveFromFieldData
PUBLIC:: FinalizeElemData
PUBLIC:: FinalizeFieldData
PUBLIC:: GetDatasetNamesInGroup
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersIO_HDF5()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("IO_HDF5")
CALL prms%CreateLogicalOption('gatheredWrite', "Set true to activate gathered HDF5 IO for parallel computations. "//&
                                               "Only local group masters will write data after gathering from local slaves.",&
                                               '.FALSE.')
#if PP_dim == 2
CALL prms%CreateLogicalOption('output2D'     , "Set true to activate hdf5 data output with flat third dimension.",'.TRUE.')
#endif /*PP_dim == 2*/
END SUBROUTINE DefineParametersIO_HDF5


!==================================================================================================================================
!> Initialize HDF5 IO
!==================================================================================================================================
SUBROUTINE InitIOHDF5()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,ONLY:GETLOGICAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
gatheredWrite=.FALSE.
IF(nLeaderProcs.LT.nProcessors) gatheredWrite=GETLOGICAL('gatheredWrite')

#if PP_dim == 3
output2D = .FALSE.
#else
output2D = GETLOGICAL('output2D')
#endif

CALL InitMPIInfo()
END SUBROUTINE InitIOHDF5


!==================================================================================================================================
!> Initialize MPIInfo variable
!==================================================================================================================================
SUBROUTINE InitMPIInfo()
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

#if USE_MPI
CALL MPI_Info_Create(MPIInfo, iError)

!normal case:
MPIInfo=MPI_INFO_NULL

! Large block IO extremely slow on Juqeen cluster (only available on IBM clusters)
!CALL MPI_Info_set(MPIInfo, "IBM_largeblock_io", "true", iError)
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
END SUBROUTINE InitMPIInfo


!==================================================================================================================================
!> Open HDF5 file and groups
!==================================================================================================================================
SUBROUTINE OpenDataFile(FileString,create,single,readOnly,userblockSize &
#if USE_MPI
                       ,communicatorOpt &
#endif /*USE_MPI*/
)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)        :: FileString      !< filename to be opened
LOGICAL,INTENT(IN)                 :: create          !< create file if it doesn't exist. Overwrited file if already present!
LOGICAL,INTENT(IN)                 :: single          !< single=T : only one processor opens file, single=F : open/create collectively
LOGICAL,INTENT(IN)                 :: readOnly        !< T : file is opened in read only mode, so file system timestamp remains unchanged
                                                      !< F: file is open read/write mode
#if USE_MPI
TYPE(MPI_Comm),INTENT(IN),OPTIONAL :: communicatorOpt !< only MPI and single=F: optional communicator to be used for collective access
#endif /*USE_MPI*/
                                                      !< default: MPI_COMM_FLEXI
INTEGER,INTENT(IN),OPTIONAL        :: userblockSize   !< size of the file to be prepended to HDF5 file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HSIZE_T)              :: userblockSize_loc, tmp, tmp2
#if USE_MPI
TYPE(MPI_Comm)                :: comm
#endif /*USE_MPI*/
!==================================================================================================================================
LOGWRITE(*,'(A)')'  OPEN HDF5 FILE "',TRIM(FileString),'" ...'

userblockSize_loc = 0
IF (PRESENT(userblockSize)) userblockSize_loc = userblockSize

! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)

! Setup file access property list with parallel I/O access (MPI) or with default property list.
IF (create) THEN; CALL H5PCREATE_F(H5P_FILE_CREATE_F, Plist_File_ID, iError)
ELSE            ; CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_File_ID, iError)
END IF

#if USE_MPI
IF (PRESENT(communicatorOpt)) THEN; comm = communicatorOpt
ELSE                              ; comm = MPI_COMM_FLEXI
END IF

#if HDF5_HAS_MPIF08
IF (.NOT.single)  CALL H5PSET_FAPL_MPIO_F(Plist_File_ID, comm        , MPIInfo        , iError)
#else
IF (.NOT.single)  CALL H5PSET_FAPL_MPIO_F(Plist_File_ID, comm%MPI_VAL, MPIInfo%MPI_VAL, iError)
#endif /*HDF5_HAS_MPIF08*/
#endif /*USE_MPI*/

! Open the file collectively.
IF(create)THEN
  IF (userblockSize_loc > 0) THEN
    tmp = userblockSize_loc/512
    IF (MOD(userblockSize_loc,INT(512,HSIZE_T)).GT.0) tmp = tmp+1
    tmp2 = 512*2**CEILING(LOG(REAL(tmp))/LOG(2.))
    CALL H5PSET_USERBLOCK_F(Plist_File_ID, tmp2, iError)
  END IF
  CALL H5FCREATE_F(TRIM(FileString), H5F_ACC_TRUNC_F, File_ID, iError, creation_prp = Plist_File_ID)
ELSE
  IF(.NOT.FILEEXISTS(FileString)) CALL Abort(__STAMP__,&
    'ERROR: Specified file '//TRIM(FileString)//' does not exist.')
  IF (readOnly) THEN
    CALL H5FOPEN_F(  TRIM(FileString), H5F_ACC_RDONLY_F,  File_ID, iError, access_prp = Plist_File_ID)
  ELSE
    CALL H5FOPEN_F(  TRIM(FileString), H5F_ACC_RDWR_F,  File_ID, iError, access_prp = Plist_File_ID)
  END IF
END IF
IF(iError.NE.0) CALL Abort(__STAMP__,&
  'ERROR: Could not open or create file '//TRIM(FileString))

LOGWRITE(*,*)'...DONE!'

END SUBROUTINE OpenDataFile


!==================================================================================================================================
!> Close HDF5 file and groups
!==================================================================================================================================
SUBROUTINE CloseDataFile()
! MODULES
USE MOD_Globals,ONLY:UNIT_logOut,Logging
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
LOGWRITE(*,'(A)')'  CLOSE HDF5 FILE...'
! Close file
CALL H5PCLOSE_F(Plist_File_ID, iError)
CALL H5FCLOSE_F(File_ID, iError)
! Close FORTRAN predefined datatypes.
CALL H5CLOSE_F(iError)
File_ID=0
LOGWRITE(*,*)'...DONE!'

END SUBROUTINE CloseDataFile


!==================================================================================================================================
!> Set pointers to element-wise arrays or scalars which will be gathered and written out. Both real or integer data types
!> are supported. It is also possible to pass a function pointer which will be evaluated to calculate the data.
!==================================================================================================================================
SUBROUTINE AddToElemData(ElementOut_In,VarName,RealArray,RealScalar,IntArray,IntScalar,Eval)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tElementOut),POINTER,INTENT(INOUT) :: ElementOut_In     !< Pointer list of element-wise data that is written to the state file
CHARACTER(LEN=*),INTENT(IN)             :: VarName           !< Name of the current array/scalar
REAL,INTENT(IN),TARGET,OPTIONAL         :: RealArray(nElems) !< Data is an array containing reals
REAL,INTENT(IN),TARGET,OPTIONAL         :: RealScalar        !< Data is a real scalar
INTEGER,INTENT(IN),TARGET,OPTIONAL      :: IntArray(nElems)  !< Data is an array containing integers
INTEGER,INTENT(IN),TARGET,OPTIONAL      :: IntScalar         !< Data is a integer scalar
PROCEDURE(EvalElemInt),POINTER,OPTIONAL :: Eval              !< Data is evaluated using a function pointer
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElementOut),POINTER          :: eout
INTEGER                            :: nOpts
!==================================================================================================================================
IF(.NOT.ASSOCIATED(ElementOut_In))THEN
  ! list is empty, create first entry
  ALLOCATE(ElementOut_In)
  eout=>ElementOut_In
ELSE
  eout=>ElementOut_In
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
!> Set pointers to node-wise arrays for output. Only real arrays or a function pointer are supported as input data.
!> Optionally, arrays can always be written to a separate dataset (even if the size is equal to the DG solution) using the
!> doSeparateOutput flag.
!==================================================================================================================================
SUBROUTINE AddToFieldData(FieldOut_In,nVal,DataSetName,VarNames,RealArray,Eval,doSeparateOutput)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tFieldOut),POINTER,INTENT(INOUT)    :: FieldOut_In       !< Pointer list of field-wise data that is written to the state file
INTEGER,INTENT(IN)                       :: nVal(4)           !< Size of array
CHARACTER(LEN=*),INTENT(IN)              :: DataSetName       !< Name of the current array (used for name of dataset)
CHARACTER(LEN=*),INTENT(IN)              :: VarNames(nVal(1)) !< Names of the variables in the array
REAL,INTENT(IN),TARGET,OPTIONAL          :: RealArray(nVal(1),nVal(2),nVal(3),nVal(4),nElems) !< Data is a real array
PROCEDURE(EvalFieldInt),POINTER,OPTIONAL :: Eval              !< Data is evaluated using a function pointer
LOGICAL,OPTIONAL,INTENT(IN)              :: doSeparateOutput  !< Flag used to always write this array to a separate dataset
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tFieldOut),POINTER            :: nout
INTEGER                            :: nOpts
INTEGER                            :: mask(3)
!==================================================================================================================================
! Mask gives the dimension of the DG solution (depending on 2D or 3D calculation)
#if PP_dim == 3
mask=(/PP_N+1,PP_N+1,PP_N+1/)
#else
mask=(/PP_N+1,PP_N+1,1/)
#endif

IF(.NOT.ASSOCIATED(FieldOut_In))THEN
  ! list is empty, create first entry
  ALLOCATE(FieldOut_In)
  nout=>FieldOut_In
ELSE
  ! loop until last entry
  nout=>FieldOut_In
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
! Optional argument that writes this array as a separate dataset even if N = PP_N. Arrays with a size different from
! the DG solution will always be written to a separate dataset.
nout%doSeparateOutput = .FALSE.
IF (PRESENT(doSeparateOutput)) nout%doSeparateOutput = doSeparateOutput
IF (ANY(nVal(2:4).NE.mask)) nout%doSeparateOutput=.TRUE.
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


!==================================================================================================================================
!>
!==================================================================================================================================
SUBROUTINE RemoveFromElemData(ElementOut_In,VarName)
! MODULES
USE MOD_Globals
USE MOD_StringTools,ONLY: STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tElementOut),POINTER,INTENT(INOUT) :: ElementOut_In     !< Pointer list of element-wise data that is written to the state file
CHARACTER(LEN=*),INTENT(IN)             :: VarName           !< Name of the current array/scalar
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElementOut),POINTER          :: current,next
!==================================================================================================================================

IF(ASSOCIATED(ElementOut_In)) THEN
  current => ElementOut_In

  ! Loop over the next entry, so we can change the current%next pointer
  IF (ASSOCIATED(current%next)) THEN
    DO WHILE (ASSOCIATED(current%next))
      next => current%next%next
      IF (STRICMP(current%next%VarName,VarName)) THEN
        IF (ASSOCIATED( current%next%RealArray  ))  NULLIFY(current%next%RealArray )
        IF (ASSOCIATED( current%next%RealScalar ))  NULLIFY(current%next%RealScalar)
        IF (ASSOCIATED( current%next%IntArray   ))  NULLIFY(current%next%IntArray  )
        IF (ASSOCIATED( current%next%IntScalar  ))  NULLIFY(current%next%IntScalar )
        IF (ASSOCIATED( current%next%eval       ))  NULLIFY(current%next%eval      )
        IF (ASSOCIATED( current%next%next       ))  NULLIFY(current%next%next      )
        DEALLOCATE(current%next)
        NULLIFY(current%next)
        RETURN
      END IF
    END DO

  ! If the list has only one entry, check if it matches and remove
  ELSE
    IF (STRICMP(current%VarName,VarName)) THEN
      IF (ASSOCIATED( current%RealArray         ))  NULLIFY(current%RealArray )
      IF (ASSOCIATED( current%RealScalar        ))  NULLIFY(current%RealScalar)
      IF (ASSOCIATED( current%IntArray          ))  NULLIFY(current%IntArray  )
      IF (ASSOCIATED( current%IntScalar         ))  NULLIFY(current%IntScalar )
      IF (ASSOCIATED( current%eval              ))  NULLIFY(current%eval      )
      IF (ASSOCIATED( current%next              ))  NULLIFY(current%next      )
      DEALLOCATE(current)
      NULLIFY(current)
      NULLIFY(ElementOut_In)
      RETURN
    END IF
  END IF
END IF

END SUBROUTINE RemoveFromElemData


!==================================================================================================================================
!>
!==================================================================================================================================
SUBROUTINE RemoveFromFieldData(FieldOut_In,DatasetName)
! MODULES
USE MOD_Globals
USE MOD_StringTools,ONLY: STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tFieldOut),POINTER,INTENT(INOUT)   :: FieldOut_In       !< Pointer list of field data that is written to the state file
CHARACTER(LEN=*),INTENT(IN)             :: DatasetName       !< Name of the dataset
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tFieldOut),POINTER                 :: current,next
!==================================================================================================================================

IF(ASSOCIATED(FieldOut_In)) THEN
  current => FieldOut_In

  ! Loop over the next entry, so we can change the current%next pointer
  IF (ASSOCIATED(current%next)) THEN
    DO WHILE (ASSOCIATED(current%next))
      next => current%next%next
      IF (STRICMP(current%next%DatasetName,DatasetName)) THEN
        IF (ASSOCIATED( current%next%RealArray  ))  NULLIFY(current%next%RealArray )
        IF (ASSOCIATED( current%next%eval       ))  NULLIFY(current%next%eval      )
        IF (ASSOCIATED( current%next%next       ))  NULLIFY(current%next%next      )
        DEALLOCATE(current%next)
        NULLIFY(current%next)
        RETURN
      END IF
    END DO

  ! If the list has only one entry, check if it matches and remove
  ELSE IF (ASSOCIATED(current)) THEN
    IF (STRICMP(current%DatasetName,DatasetName)) THEN
      IF (ASSOCIATED( current%RealArray         ))  NULLIFY(current%RealArray )
      IF (ASSOCIATED( current%eval              ))  NULLIFY(current%eval      )
      IF (ASSOCIATED( current%next              ))  NULLIFY(current%next      )
      DEALLOCATE(current)
      NULLIFY(current)
      NULLIFY(FieldOut_In)
      RETURN
    END IF
  END IF
END IF

END SUBROUTINE RemoveFromFieldData


!==================================================================================================================================
!> Takes a group and reads the names of the datasets
!==================================================================================================================================
SUBROUTINE GetDatasetNamesInGroup(group,names)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                  :: group    !< name of group
CHARACTER(LEN=255),ALLOCATABLE,INTENT(INOUT) :: names(:) !< names of datasets
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nMembers,i,type
!===================================================================================================================================
CALL H5GN_MEMBERS_F(File_ID, TRIM(group), nMembers, iError)
ALLOCATE(names(nMembers))
DO i=1,nMembers
  CALL h5gget_obj_info_idx_f(File_ID, TRIM(group), i-1, names(i), type, iError)
  IF (type.NE.H5G_DATASET_F) names(i) = ''
END DO
END SUBROUTINE GetDatasetNamesInGroup


!==================================================================================================================================
!> Finalizes HDF5 IO
!==================================================================================================================================
SUBROUTINE FinalizeIOHDF5()
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL FinalizeElemData(ElementOut)
CALL FinalizeFieldData(FieldOut)

END SUBROUTINE FinalizeIOHDF5


!==================================================================================================================================
!> Deallocate and nullify element-wise arrays or scalars which were added for writeout by means of a linked list.
!==================================================================================================================================
SUBROUTINE FinalizeElemData(ElementOut_In)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tElementOut),POINTER,INTENT(INOUT) :: ElementOut_In     !< Pointer list of element-wise data that is written to the state file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
TYPE(tElementOut),POINTER          :: current,next
!===================================================================================================================================
IF(ASSOCIATED(ElementOut_In)) THEN
  current => ElementOut_In
  NULLIFY(ElementOut_In)

  ! Finalize all entries
  DO WHILE (ASSOCIATED(current))
    next => current%next
    IF (ASSOCIATED( current%RealArray  ))  NULLIFY(current%RealArray )
    IF (ASSOCIATED( current%RealScalar ))  NULLIFY(current%RealScalar)
    IF (ASSOCIATED( current%IntArray   ))  NULLIFY(current%IntArray  )
    IF (ASSOCIATED( current%IntScalar  ))  NULLIFY(current%IntScalar )
    IF (ASSOCIATED( current%eval       ))  NULLIFY(current%eval      )
    IF (ASSOCIATED( current%next       ))  NULLIFY(current%next      )
    DEALLOCATE(current)
    NULLIFY(current)
    current => next
  END DO
END IF

END SUBROUTINE FinalizeElemData


!==================================================================================================================================
!> Deallocate and nullify additional field data arrays or scalars which were added for writeout by means of a linked list.
!==================================================================================================================================
SUBROUTINE FinalizeFieldData(FieldOut_In)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tFieldOut),POINTER,INTENT(INOUT) :: FieldOut_In     !< Pointer list of field data that is written to the state file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
TYPE(tFieldOut),POINTER          :: current,next
!===================================================================================================================================
IF(ASSOCIATED(FieldOut_In)) THEN
  current => FieldOut_In
  NULLIFY(FieldOut_In)

  ! Finalize all entries
  DO WHILE (ASSOCIATED(current))
    next => current%next
    IF (ASSOCIATED( current%RealArray ))  NULLIFY(current%RealArray)
    IF (ASSOCIATED( current%eval      ))  NULLIFY(current%eval     )
    IF (ASSOCIATED( current%next      ))  NULLIFY(current%next     )
    DEALLOCATE(current)
    NULLIFY(current)
    current => next
  END DO
END IF

END SUBROUTINE FinalizeFieldData

END MODULE MOD_IO_HDF5
