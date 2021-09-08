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
!> IO routines to gather informations from HDF5 files
!==================================================================================================================================
MODULE MOD_HDF5_Input
! MODULES
USE MOD_IO_HDF5
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ISVALIDHDF5FILE
  MODULE PROCEDURE ISVALIDHDF5FILE
END INTERFACE

INTERFACE ISVALIDMESHFILE
  MODULE PROCEDURE ISVALIDMESHFILE
END INTERFACE

INTERFACE GetNextFileName
  MODULE PROCEDURE GetNextFileName
END INTERFACE

INTERFACE DatasetExists
  MODULE PROCEDURE DatasetExists
END INTERFACE

INTERFACE GetDataSize
  MODULE PROCEDURE GetDataSize
END INTERFACE

INTERFACE GetAttributeSize
  MODULE PROCEDURE GetAttributeSize
END INTERFACE

INTERFACE GetDataProps
  MODULE PROCEDURE GetDataProps
END INTERFACE

INTERFACE ReadAttribute
  MODULE PROCEDURE ReadAttribute
END INTERFACE

INTERFACE GetVarnames
  MODULE PROCEDURE GetVarnames
END INTERFACE

PUBLIC :: Plist_File_ID,File_ID,HSize,nDims        ! Variables from MOD_IO_HDF5 that need to be public
PUBLIC :: OpenDataFile,CloseDataFile ! Subroutines from MOD_IO_HDF5 that need to be public
PUBLIC :: ISVALIDHDF5FILE,ISVALIDMESHFILE,GetDataSize,GetAttributeSize,GetDataProps,GetNextFileName
PUBLIC :: ReadArray,ReadAttribute
PUBLIC :: GetArrayAndName
PUBLIC :: GetVarnames
PUBLIC :: DatasetExists
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Subroutine to check if a file is a valid Flexi HDF5 file
!==================================================================================================================================
FUNCTION ISVALIDHDF5FILE(FileName,ProgramName,FileType,FileVersion)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName        !< name of file to be checked
CHARACTER(LEN=255),INTENT(OUT),OPTIONAL :: ProgramName !< program name
CHARACTER(LEN=255),INTENT(OUT),OPTIONAL :: FileType    !< type of the file (only if valid)
REAL,INTENT(OUT),OPTIONAL      :: FileVersion     !< desired version
LOGICAL                        :: isValidHDF5File !< result: file is valid HDF5 file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Plist_ID
LOGICAL                        :: exists
!==================================================================================================================================
isValidHDF5File=.TRUE.
iError=0

! Disable error messages
CALL H5ESET_AUTO_F(0, iError)
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Create property list
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#if USE_MPI
! Setup file access property list with parallel I/O access (MPI)
CALL H5PSET_FAPL_MPIO_F(Plist_ID,MPI_COMM_FLEXI, MPIInfo, iError)
#endif /*USE_MPI*/

! Check if file exists
exists = FILEEXISTS(FileName)
IF(.NOT.exists) THEN
  CALL abort(__STAMP__,'ERROR: HDF5 file '//TRIM(FileName)//' does not exist.')
  RETURN
END IF

! Open HDF5 file
CALL H5FOPEN_F(TRIM(FileName), H5F_ACC_RDONLY_F, File_ID, iError,access_prp = Plist_ID)
IF(iError.EQ.0) THEN
  isValidHDF5File=.TRUE.
  ! Check attributes file type and version -----------------------------
  IF(PRESENT(ProgramName))THEN
    CALL DatasetExists(File_ID,'Program',exists,attrib=.TRUE.)
    IF(exists)THEN
      CALL ReadAttribute(File_ID,'Program',1,StrScalar=ProgramName)
    ELSE
      ProgramName=''
      isValidHDF5File=.FALSE.
    END IF
  END IF
  IF(PRESENT(FileType))THEN
    CALL DatasetExists(File_ID,'File_Type',exists,attrib=.TRUE.)
    IF(exists)THEN
      CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=FileType)
    ELSE
      FileType=''
      isValidHDF5File=.FALSE.
    END IF
  END IF
  IF(PRESENT(FileVersion))THEN
    CALL DatasetExists(File_ID,'File_Version',exists,attrib=.TRUE.)
    IF(exists)THEN
      CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=FileVersion)
    ELSE
      FileVersion=-1.
      isValidHDF5File=.FALSE.
    END IF
  END IF
  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close the file.
  CALL H5FCLOSE_F(File_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  isValidHDF5File=.FALSE.
  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
END IF

END FUNCTION ISVALIDHDF5FILE

!==================================================================================================================================
!> Subroutine to check if a file is a valid mesh file
!==================================================================================================================================
FUNCTION ISVALIDMESHFILE(MeshFileName)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName    !< name of mesh file to be checked
LOGICAL                        :: isValidMeshFile !< result: file is valid mesh file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: NGeoExists
INTEGER(HID_T)                 :: Plist_ID
!==================================================================================================================================
! Disable error messages
CALL H5ESET_AUTO_F(0, iError)

! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Create property list
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#if USE_MPI
! Setup file access property list with parallel I/O access (MPI)
CALL H5PSET_FAPL_MPIO_F(Plist_ID,MPI_COMM_FLEXI, MPIInfo, iError)
#endif /*USE_MPI*/

! Check if file exists
IF(.NOT.FILEEXISTS(MeshFileName)) THEN
  CALL abort(__STAMP__,'ERROR: Mesh file '//TRIM(MeshFileName)//' does not exist.')
  isValidMeshFile = .FALSE.
  RETURN
END IF

! Open HDF5 file
CALL H5FOPEN_F(TRIM(MeshFileName), H5F_ACC_RDONLY_F, File_ID, iError,access_prp = Plist_ID)
IF(iError.EQ.0) THEN
  isValidMeshFile=.TRUE.

  ! Check NGeo attribute --------------------------------------------------------------------------------------------------------
  CALL DatasetExists(File_ID,'Ngeo',NGeoExists,attrib=.TRUE.)
  IF (.NOT.NGeoExists) isValidMeshFile = .FALSE.

  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close the file.
  CALL H5FCLOSE_F(File_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  isValidMeshFile=.FALSE.
  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
END IF
END FUNCTION ISVALIDMESHFILE

!==================================================================================================================================
!> Subroutine to determine HDF5 datasize
!==================================================================================================================================
SUBROUTINE GetDataSize(Loc_ID,DSetName,nDims,Size)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                     :: DSetName !< name if dataset to be checked
INTEGER(HID_T),INTENT(IN)            :: Loc_ID   !< ID of dataset
INTEGER,INTENT(OUT)                  :: nDims    !< found data size dimensions
INTEGER(HSIZE_T),POINTER,INTENT(OUT) :: Size(:)  !< found data size
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: DSet_ID,FileSpace
INTEGER(HSIZE_T), POINTER            :: SizeMax(:)
!==================================================================================================================================
! Open the dataset with default properties.
CALL H5DOPEN_F(Loc_ID, TRIM(DSetName) , DSet_ID, iError)
! Get the data space of the dataset.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
! Get size and max size of data space
ALLOCATE(Size(nDims),SizeMax(nDims))
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Size, SizeMax, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(DSet_ID, iError)
DEALLOCATE(SizeMax)
END SUBROUTINE GetDataSize

!==================================================================================================================================
!> Subroutine to determine HDF5 size of attribute
!==================================================================================================================================
SUBROUTINE GetAttributeSize(Loc_ID,AttribName,nDims,Size)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                     :: AttribName !< name if attribute to be checked
INTEGER(HID_T),INTENT(IN)            :: Loc_ID   !< ID of dataset
INTEGER,INTENT(OUT)                  :: nDims    !< found data size dimensions
INTEGER(HSIZE_T),POINTER,INTENT(OUT) :: Size(:)  !< found data size
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: Attr_ID,FileSpace
INTEGER(HSIZE_T), POINTER            :: SizeMax(:)
!==================================================================================================================================
! Open the dataset with default properties.
CALL H5AOPEN_F(Loc_ID, TRIM(AttribName), Attr_ID, iError)
! Get the data space of the dataset.
CALL H5AGET_SPACE_F(Attr_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
! Get size and max size of data space
ALLOCATE(Size(nDims),SizeMax(nDims))
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Size, SizeMax, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5ACLOSE_F(Attr_ID, iError)
DEALLOCATE(SizeMax)
END SUBROUTINE GetAttributeSize

!==================================================================================================================================
!> @brief Subroutine to check wheter a dataset in the hdf5 file exists
!>
!> We have no "h5dexists_f", so we use the error given by h5dopen_f.
!> this produces hdf5 error messages even if everything is ok, so we turn the error msgs off
!> during this operation.
!> auto error messages off
!==================================================================================================================================
SUBROUTINE DatasetExists(Loc_ID,DSetName,Exists,attrib)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                     :: DSetName !< name if dataset to be checked
INTEGER(HID_T),INTENT(IN)            :: Loc_ID   !< ID of dataset
LOGICAL,INTENT(IN),OPTIONAL          :: attrib   !< check dataset or attribute
LOGICAL,INTENT(OUT)                  :: Exists   !< result: dataset exists
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                              :: attrib_loc
!==================================================================================================================================

IF (PRESENT(attrib)) THEN
  attrib_loc = attrib
ELSE
  attrib_loc = .FALSE.
END IF

! Check attribute or data set. Data sets can be checked by determining the existence of the corresponding link
IF(attrib_loc)THEN
  CALL H5AEXISTS_F(Loc_ID, TRIM(DSetName), Exists, iError)
ELSE
  CALL H5LEXISTS_F(Loc_ID, TRIM(DSetName), Exists, iError)
END IF

END SUBROUTINE DatasetExists



!==================================================================================================================================
!> Subroutine to determine HDF5 dataset properties in Flexi terminology
!==================================================================================================================================
SUBROUTINE GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,ArrayName_opt)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(OUT)                     :: nVar_HDF5     !< number of variables
INTEGER,INTENT(OUT)                     :: N_HDF5        !< polynomial degree
INTEGER,INTENT(OUT)                     :: nElems_HDF5   !< inumber of elements
CHARACTER(LEN=255),OPTIONAL,INTENT(OUT) :: NodeType_HDF5 !< nodetype string
CHARACTER(LEN=*),  OPTIONAL,INTENT(IN)  :: ArrayName_opt !< array name to use in state file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                      :: ArrayName
INTEGER                                 :: Rank
INTEGER(HID_T)                          :: Dset_ID,FileSpace
INTEGER(HSIZE_T), DIMENSION(7)          :: Dims,DimsMax
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,A)')' GET SIZE OF DATA IN HDF5 FILE... '

IF(.NOT.PRESENT(ArrayName_opt)) THEN
    ArrayName = 'DG_Solution'
ELSE
    ArrayName = ArrayName_opt
END IF

! Read in attributes
! Open the dataset with default properties.
CALL H5DOPEN_F(File_ID, ArrayName, Dset_ID, iError)
! Get the data space of the dataset.
CALL H5DGET_SPACE_F(Dset_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, Rank, iError)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Rank of database',' | ',Rank,' | HDF5    | '
! Get size and max size of data space
Dims   =0
DimsMax=0
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Dims(1:Rank), DimsMax(1:Rank), iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(Dset_ID, iError)
IF(PRESENT(NodeType_HDF5)) THEN
  ! Read in NodeType
  CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType_HDF5)
END IF

! Display data
! nVar = first array index
!nVar_HDF5 = INT(Dims(1),4)
CHECKSAFEINT(Dims(1),4)
nVar_HDF5 = INT(Dims(1),4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Number of variables nVar',' | ',nVar_HDF5,' | HDF5    | '
! N = index 2-4 of array, is expected to have the same value for each direction
CHECKSAFEINT(Dims(2)-1,4)
N_HDF5 = INT(Dims(2)-1,4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Polynomial degree N',' | ',N_HDF5,' | HDF5    | '
IF(PRESENT(NodeType_HDF5)) THEN
  SWRITE(UNIT_stdOut,'(A3,A30,A3,A33,A13)')' | ','          Node type',' | ',TRIM(NodeType_HDF5),' | HDF5    | '
END IF
! nElems = index 5 of array
CHECKSAFEINT(Dims(5),4)
nElems_HDF5 = INT(Dims(5),4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','GeometricnElems',' | ',nElems_HDF5,' | HDF5    | '

SWRITE(UNIT_stdOut,'(A)')' DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE GetDataProps


SUBROUTINE GetVarNames(AttribName,VarNames,AttribExists)
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                :: AttribName
CHARACTER(LEN=255),ALLOCATABLE,INTENT(OUT) :: VarNames(:)
LOGICAL,INTENT(OUT)                        :: AttribExists
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dims, nVal
!===================================================================================================================================
SDEALLOCATE(VarNames)
CALL DatasetExists(File_ID,AttribName,AttribExists,attrib=.TRUE.)
IF (AttribExists) THEN
  ! get size of array
  CALL GetAttributeSize(File_ID,AttribName,dims,HSize)
  nVal=INT(HSize(1))
  DEALLOCATE(HSize)
  ALLOCATE(VarNames(nVal))

  ! read variable names
  CALL ReadAttribute(File_ID,TRIM(AttribName),nVal,StrArray=VarNames)
END IF
END SUBROUTINE GetVarNames

!===================================================================================================================================
!> High level wrapper to ReadArray and ReadAttrib. Check if array exists and directly
!> allocate, read array and attribs
!> Assume that the array to be read is of size (nVar,.,.,.,.,nElems) and that an associated
!> attribute containing the variable names exists
!===================================================================================================================================
SUBROUTINE GetArrayAndName(ArrayName,AttribName,nVal,Array,VarNames)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars    ,ONLY: nElems,nGlobalElems,OffsetElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)     :: ArrayName   !< name of array to be read
CHARACTER(LEN=*),INTENT(IN)     :: AttribName  !< name of varnames to be read
INTEGER,INTENT(OUT)             :: nVal(15)    !< size of array
REAL,ALLOCATABLE,INTENT(OUT)    :: Array(:)    !< array to be read
CHARACTER(LEN=255),ALLOCATABLE,INTENT(OUT) :: VarNames(:) !< variable names
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL  :: found
INTEGER  :: dims
!===================================================================================================================================
nVal=-1
SDEALLOCATE(Array)
SDEALLOCATE(VarNames)

CALL DatasetExists(File_ID, TRIM(ArrayName), found)
IF (found) THEN
  ! get size of array
  CALL GetDataSize(File_ID,TRIM(ArrayName),dims,HSize)
  nVal(1:dims)=INT(HSize)
  IF(nVal(dims).NE.nGlobalElems) STOP 'Last array dimension != nElems !'
  nVal(dims)=nElems
  DEALLOCATE(HSize)
  ALLOCATE(array(PRODUCT(nVal(1:dims))))
  ALLOCATE(VarNames(nVal(1)))

  ! read array
  CALL ReadArray(TRIM(ArrayName),dims,nVal(1:dims),OffsetElem,dims,RealArray=array)

  ! read variable names
  CALL ReadAttribute(File_ID,TRIM(AttribName),nVal(1),StrArray=VarNames)
END IF

END SUBROUTINE GetArrayAndName


!==================================================================================================================================
!> Subroutine to read arrays of rank "Rank" with dimensions "Dimsf(1:Rank)".
!==================================================================================================================================
SUBROUTINE ReadArray(ArrayName,Rank,nVal,Offset_in,Offset_dim,RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER                        :: Rank                  !< number of dimensions of the array
INTEGER                        :: offset_in             !< offset =0, start at beginning of the array
INTEGER                        :: offset_dim            !< which dimension is the offset (only one dimension possible here)
INTEGER                        :: nVal(Rank)            !< local size of array to read
CHARACTER(LEN=*),INTENT(IN)    :: ArrayName             !< name of array to be read
REAL              ,DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET :: RealArray    !< only if real array shall be read
INTEGER           ,DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET :: IntArray     !< only if integer array shall be read
CHARACTER(LEN=255),DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET :: StrArray     !< only if real string shall be read
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: DSet_ID,Type_ID,MemSpace,FileSpace,PList_ID
INTEGER(HSIZE_T)               :: Offset(Rank),Dimsf(Rank)
TYPE(C_PTR)                    :: buf
INTEGER(HID_T)                 :: driver
!==================================================================================================================================
#if USE_MPI
! HDF5 with MPI can only read max. (32 bit signed integer / size of single element) elements (2GB per MPI rank)
IF(PRODUCT(nVal).GT.nLimit)  CALL Abort(__STAMP__, &
 'Dataset "'//TRIM(ArrayName)//'" exceeds HDF5 chunksize limit of 2GB per rank! Increase number of ranks or compile without MPI!')
#endif

LOGWRITE(*,'(A,I1.1,A,A,A)')'    READ ',Rank,'D ARRAY "',TRIM(ArrayName),'"'
Dimsf=nVal
LOGWRITE(*,*)'Dimsf,Offset=',Dimsf,Offset_in
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
CALL H5DOPEN_F(File_ID, TRIM(ArrayName) , DSet_ID, iError)

IF(iError.NE.0) &
  CALL Abort(__STAMP__,'Array '//TRIM(ArrayName)//' does not exist.')

! Define and select the hyperslab to use for reading.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
Offset(:)=0
Offset(offset_dim)=Offset_in
CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
! Create property list
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#if USE_MPI
! Set property list to collective dataset read
CALL H5PGET_DRIVER_F(Plist_File_ID,driver,iError)
IF(driver.EQ.H5FD_MPIO_F) CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F, iError)
#endif

IF(PRESENT(RealArray)) Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntArray))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(StrArray))  CALL H5DGET_TYPE_F(DSet_ID, Type_ID, iError)

buf=C_NULL_PTR
IF(PRESENT(RealArray)) buf=C_LOC(RealArray)
IF(PRESENT(IntArray))  buf=C_LOC(IntArray)
IF(PRESENT(StrArray))  buf=C_LOC(StrArray(1))

! Read the data
CALL H5DREAD_F(DSet_ID,Type_ID,buf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)

! Close the datatype, property list, dataspaces and dataset.
IF(PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
CALL H5PCLOSE_F(PList_ID,iError)
CALL H5SCLOSE_F(FileSpace,iError)
CALL H5DCLOSE_F(DSet_ID, iError)
CALL H5SCLOSE_F(MemSpace,iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE ReadArray



!==================================================================================================================================
!> Subroutine to read attributes from HDF5 file.
!==================================================================================================================================
SUBROUTINE ReadAttribute(Loc_ID_in,AttribName,nVal,DatasetName,RealScalar,IntScalar,&
                                 StrScalar,LogicalScalar,RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(HID_T)    ,INTENT(IN)                  :: Loc_ID_in         !< HDF5 file id of opened file
INTEGER           ,INTENT(IN)                  :: nVal              !< number of attributes in case an array is expected
CHARACTER(LEN=*)  ,INTENT(IN)                  :: AttribName        !< name of attribute to be read
CHARACTER(LEN=*)  ,INTENT(IN) ,OPTIONAL        :: DatasetName       !< dataset name in case attribute is located in a dataset
REAL              ,INTENT(OUT),OPTIONAL,TARGET :: RealArray(nVal)   !< Array of real array attributes
INTEGER           ,INTENT(OUT),OPTIONAL,TARGET :: IntArray(nVal)    !< Array for integer array for attributes
REAL              ,INTENT(OUT),OPTIONAL,TARGET :: RealScalar        !< Scalar real attribute
INTEGER           ,INTENT(OUT),OPTIONAL,TARGET :: IntScalar         !< Scalar integer attribute
CHARACTER(LEN=255),INTENT(OUT),OPTIONAL,TARGET :: StrScalar         !< Scalar string attribute
CHARACTER(LEN=255),INTENT(OUT),OPTIONAL,TARGET :: StrArray(nVal)    !< Array for character array attributes
LOGICAL           ,INTENT(OUT),OPTIONAL        :: LogicalScalar     !< Scalar logical attribute
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Attr_ID,Type_ID,Loc_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER,TARGET                 :: IntToLog
CHARACTER(LEN=255),TARGET      :: StrTmp(1)
TYPE(C_PTR)                    :: buf
!==================================================================================================================================
LOGWRITE(*,*)' READ ATTRIBUTE "',TRIM(AttribName),'" FROM HDF5 FILE...'
Dimsf(1)=nVal
Loc_ID=Loc_ID_in
IF(PRESENT(DatasetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
END IF

! Create the attribute for group Loc_ID.
CALL H5AOPEN_F(Loc_ID, TRIM(AttribName), Attr_ID, iError)

IF(iError.NE.0) &
  CALL Abort(__STAMP__,'Attribute '//TRIM(AttribName)//' does not exist.')

IF(PRESENT(RealArray))     RealArray=0.
IF(PRESENT(RealScalar))    RealScalar=0.
IF(PRESENT(IntArray))      IntArray=0
IF(PRESENT(IntScalar))     IntScalar=0
IF(PRESENT(LogicalScalar)) LogicalScalar=.FALSE.
IF(PRESENT(StrScalar))THEN
  StrScalar=''
  StrTmp(1)=''
END IF
IF(PRESENT(StrArray))      StrArray(:)=''

IF(PRESENT(RealArray))     Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(RealScalar))    Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntArray))      Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntScalar))     Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(LogicalScalar)) Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5AGET_TYPE_F(Attr_ID, Type_ID, iError)

buf=C_NULL_PTR
IF(PRESENT(RealArray))     buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))    buf=C_LOC(RealScalar)
IF(PRESENT(IntArray))      buf=C_LOC(IntArray)
IF(PRESENT(IntScalar))     buf=C_LOC(IntScalar)
IF(PRESENT(LogicalScalar)) buf=C_LOC(IntToLog)
IF(PRESENT(StrScalar))     buf=C_LOC(StrTmp(1))
IF(PRESENT(StrArray))      buf=C_LOC(StrArray(1))

! Read the attribute data.
CALL H5AREAD_F(Attr_ID, Type_ID, buf, iError)

IF(PRESENT(LogicalScalar)) LogicalScalar=(IntToLog.EQ.1)
IF(PRESENT(StrScalar))     StrScalar=StrTmp(1)
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)

! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE ReadAttribute



!==================================================================================================================================
!> Subroutine to determine filename of next HDF5 file for FlushFiles
!==================================================================================================================================
SUBROUTINE GetNextFileName(FileName,NextFileName_HDF5,single)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName                !< filename to check
LOGICAL,INTENT(IN)             :: single                  !< switch whether file is being accessed in parallel my MPI_COMM_FLEXI
CHARACTER(LEN=255),INTENT(OUT) :: NextFileName_HDF5       !< output: follow up file according to checked file opened
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: ReadError
INTEGER(HID_T)                 :: File_ID_loc,Plist_ID
!==================================================================================================================================
LOGWRITE(*,*)' GET NEXT FILE NAME FROM HDF5 FILE ', TRIM (FileName),' ...'
ReadError=0
NextFileName_HDF5=''
! Disable error messages
CALL H5ESET_AUTO_F(0, iError)
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Setup file access property list
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#if USE_MPI
IF(.NOT.single)THEN
  ! Set property list to MPI IO
  CALL H5PSET_FAPL_MPIO_F(Plist_ID, MPI_COMM_FLEXI, MPI_INFO_NULL, iError)
END IF
#endif /*USE_MPI*/
! Open file
CALL H5FOPEN_F(TRIM(FileName), H5F_ACC_RDONLY_F, File_ID_loc, iError,access_prp = Plist_ID)
ReadError=iError
CALL H5PCLOSE_F(Plist_ID, iError)
iError=ReadError
IF (iError .EQ. 0) THEN
  ! Get Name of the mesh file, stored as third atrribute with name "NextFile"
  ! Open the attribute "NextFile" of opened file
  CALL ReadAttribute(File_ID_loc,'NextFile',1,StrScalar=NextFileName_HDF5)
  ! Close the file.
  CALL H5FCLOSE_F(File_ID_loc, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
  iError=-1
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE GetNextFileName

END MODULE MOD_HDF5_Input
