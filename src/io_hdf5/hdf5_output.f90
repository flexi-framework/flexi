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
!> Module providing IO routines for parallel output in HDF5 format: solution, time averaged files, baseflow, record points,...
!==================================================================================================================================
MODULE MOD_HDF5_Output
! MODULES
USE MOD_IO_HDF5
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE WriteState
  MODULE PROCEDURE WriteState
END INTERFACE

INTERFACE WriteTimeAverage
  MODULE PROCEDURE WriteTimeAverage
END INTERFACE

INTERFACE WriteBaseflow
  MODULE PROCEDURE WriteBaseflow
END INTERFACE

INTERFACE FlushFiles
  MODULE PROCEDURE FlushFiles
END INTERFACE

INTERFACE WriteHeader
  MODULE PROCEDURE WriteHeader
END INTERFACE

!INTERFACE WriteArray
!  MODULE PROCEDURE WriteArray
!END INTERFACE

INTERFACE WriteAttribute
  MODULE PROCEDURE WriteAttribute
END INTERFACE

INTERFACE MarkWriteSuccessfull
  MODULE PROCEDURE MarkWriteSuccessfull
END INTERFACE

INTERFACE WriteAdditionalElemData
  MODULE PROCEDURE WriteAdditionalElemData
END INTERFACE

INTERFACE
  SUBROUTINE copy_userblock(outfilename,infilename) BIND(C)
      USE ISO_C_BINDING, ONLY: C_CHAR
      CHARACTER(KIND=C_CHAR) :: outfilename(*)
      CHARACTER(KIND=C_CHAR) :: infilename(*)
  END SUBROUTINE copy_userblock
END INTERFACE

INTERFACE GenerateFileSkeleton
  MODULE PROCEDURE GenerateFileSkeleton
END INTERFACE


PUBLIC :: WriteState,FlushFiles,WriteHeader,WriteTimeAverage,WriteBaseflow,GenerateFileSkeleton
PUBLIC :: WriteArray,WriteAttribute,GatheredWriteArray,WriteAdditionalElemData,MarkWriteSuccessfull
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Subroutine to write the solution U to HDF5 format
!> Is used for postprocessing and for restart
!==================================================================================================================================
SUBROUTINE WriteState(MeshFileName,OutputTime,FutureTime,isErrorFile)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars           ,ONLY: U
USE MOD_Output_Vars       ,ONLY: ProjectName,NOut,Vdm_N_NOut,WriteStateFiles
USE MOD_Mesh_Vars         ,ONLY: offsetElem,nGlobalElems,sJ,nElems
USE MOD_ChangeBasisByDim  ,ONLY: ChangeBasisVolume
USE MOD_Equation_Vars     ,ONLY: StrVarNames
#if PP_dim == 2
USE MOD_2D                ,ONLY: ExpandArrayTo3D
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName   !< file name of mesh used for simulation
REAL,INTENT(IN)                :: OutputTime     !< simulation time when output is performed
REAL,INTENT(IN)                :: FutureTime     !< hint, when next file will be written
LOGICAL,INTENT(IN)             :: isErrorFile    !< indicate whether an error file is written in case of a crashed simulation
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName,FileType
REAL                           :: StartT,EndT
REAL,POINTER                   :: UOut(:,:,:,:,:)
#if PP_dim == 2
REAL,ALLOCATABLE               :: UOutTmp(:,:,:,:,:)
#endif
REAL                           :: Utmp(5,0:PP_N,0:PP_N,0:PP_NZ)
REAL                           :: JN(1,0:PP_N,0:PP_N,0:PP_NZ),JOut(1,0:NOut,0:NOut,0:ZDIM(NOut))
INTEGER                        :: iElem,i,j,k
INTEGER                        :: nVal(5)
!==================================================================================================================================
IF (.NOT.WriteStateFiles) RETURN
IF(MPIRoot)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE STATE TO HDF5 FILE...'
  GETTIME(StartT)
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileType=MERGE('ERROR_State','State      ',isErrorFile)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileType),OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton(TRIM(FileName),'State',PP_nVar,NOut,StrVarNames,&
                                      MeshFileName,OutputTime,FutureTime,withUserblock=.TRUE.)

! Set size of output
nVal=(/PP_nVar,NOut+1,NOut+1,ZDIM(NOut)+1,nElems/)

! build output data
IF(NOut.NE.PP_N)THEN
#if FV_ENABLED
  CALL Abort(__STAMP__, &
      "NOut not working for FV!")
#endif
  ! Project JU and J to NOut, compute U on Nout
  ALLOCATE(UOut(PP_nVar,0:NOut,0:NOut,0:ZDIM(NOut),nElems))
  DO iElem=1,nElems
    JN(1,:,:,:)=1./sJ(:,:,:,iElem,0)
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Utmp(:,i,j,k)=U(:,i,j,k,iElem)*JN(1,i,j,k)
    END DO; END DO; END DO
    CALL ChangeBasisVolume(PP_nVar,PP_N,NOut,Vdm_N_NOut,&
                           Utmp,UOut(1:PP_nVar,:,:,:,iElem))
    ! Jacobian
    CALL ChangeBasisVolume(1,PP_N,NOut,Vdm_N_NOut,JN,JOut)
    DO k=0,ZDIM(NOut); DO j=0,NOut; DO i=0,NOut
      UOut(:,i,j,k,iElem)=UOut(:,i,j,k,iElem)/JOut(1,i,j,k)
    END DO; END DO; END DO
  END DO
#if PP_dim == 2
  ! If the output should be done with a full third dimension in a two dimensional computation, we need to expand the solution
  IF (.NOT.output2D) THEN
    ALLOCATE(UOutTmp(PP_nVar,0:NOut,0:NOut,0:ZDIM(NOut),nElems))
    UOutTmp = UOut
    DEALLOCATE(UOut)
    ALLOCATE(UOut(PP_nVar,0:NOut,0:NOut,0:NOut,nElems))
    CALL ExpandArrayTo3D(5,nVal,4,Nout+1,UOutTmp,UOut)
    DEALLOCATE(UOutTmp)
    nVal=(/PP_nVar,NOut+1,NOut+1,NOut+1,nElems/)
  END IF
#endif

ELSE ! write state on same polynomial degree as the solution

#if PP_dim == 3
  UOut => U
#else
  IF (.NOT.output2D) THEN
    ! If the output should be done with a full third dimension in a two dimensional computation, we need to expand the solution
    ALLOCATE(UOut(PP_nVar,0:NOut,0:NOut,0:NOut,nElems))
    CALL ExpandArrayTo3D(5,(/PP_nVar,NOut+1,NOut+1,ZDIM(NOut)+1,nElems/),4,NOut+1,U,UOut)
    ! Correct size of the output array
    nVal=(/PP_nVar,NOut+1,NOut+1,NOut+1,nElems/)
  ELSE
    UOut => U
  END IF
#endif
END IF ! (NOut.NE.PP_N)


! Reopen file and write DG solution
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/PP_nVar,NOut+1,NOut+1,NOut+1,nGlobalElems/),&
                        nVal=nVal                                              ,&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=UOut)

! Deallocate UOut only if we did not point to U
IF((PP_N .NE. NOut).OR.((PP_dim .EQ. 2).AND.(.NOT.output2D))) DEALLOCATE(UOut)

CALL WriteAdditionalElemData(FileName,ElementOut)
CALL WriteAdditionalFieldData(FileName,FieldOut)


IF(MPIRoot)THEN
  CALL MarkWriteSuccessfull(FileName)
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF

#if USE_MPI
! Since we are going to abort directly after this wenn an error state is written, make sure that all processors are finished
! with everything or we might end up with a non-valid error state file
IF (isErrorFile) CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif
END SUBROUTINE WriteState


!==================================================================================================================================
!> This routine is a wrapper routine for WriteArray and first gathers all output arrays of an MPI sub group,
!> then only the master will write the data. Significantly reduces IO overhead for a large number of processes!
!==================================================================================================================================
SUBROUTINE GatheredWriteArray(FileName,create,DataSetName,rank,nValGlobal,nVal,offset,collective,RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName          !< Name of the file to write to
CHARACTER(LEN=*),INTENT(IN)    :: DataSetName       !< Name of the dataset to write
LOGICAL,INTENT(IN)             :: create            !< Should the file be created or not
LOGICAL,INTENT(IN)             :: collective        !< Collective write or not
INTEGER,INTENT(IN)             :: rank              !< Rank of array
INTEGER,INTENT(IN)             :: nVal(rank)        !< Local number of variables in each rank
INTEGER,INTENT(IN)             :: nValGlobal(rank)  !< Global number of variables in each rank
INTEGER,INTENT(IN)             :: offset(rank)      !< Offset in each rank
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal)) !< Real array to write
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray( PRODUCT(nVal)) !< Integer array to write
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal)) !< String array to write
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
REAL,              ALLOCATABLE :: UReal(:)
CHARACTER(LEN=255),ALLOCATABLE :: UStr(:)
INTEGER,           ALLOCATABLE :: UInt(:)
INTEGER                        :: i,nValGather(rank),nDOFLocal
INTEGER,DIMENSION(nLocalProcs) :: nDOFPerNode,offsetNode
!==================================================================================================================================
! HDF5 with MPI can only write max. (32 bit integer / 8) elements
IF(REAL(PRODUCT(nVal)).GT.((2**28-1)/8.))  CALL Abort(__STAMP__, &
    'Total size of HDF5 array "'//TRIM(DataSetName)//'" is too big! Reduce number of entries per rank or compile without MPI!')

IF(gatheredWrite)THEN
  IF(ANY(offset(1:rank-1).NE.0)) &
    CALL abort(__STAMP__,'Offset only allowed in last dimension for gathered IO.')

  ! Get last dim of each array on IO nodes
  nDOFLocal=PRODUCT(nVal)
  CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER,nDOFPerNode,1,MPI_INTEGER,0,MPI_COMM_NODE,iError)

  ! Allocate big array and compute offsets of small arrs inside big
  offsetNode=0
  IF(MPILocalRoot)THEN
    nValGather=nVal
    nValGather(rank)=SUM(nDOFPerNode)/PRODUCT(nVal(1:rank-1))
    DO i=2,nLocalProcs
      offsetNode(i)=offsetNode(i-1)+nDOFPerNode(i-1)
    END DO
    IF(PRESENT(RealArray)) ALLOCATE(UReal(PRODUCT(nValGather)))
    IF(PRESENT(IntArray))  ALLOCATE(UInt( PRODUCT(nValGather)))
    IF(PRESENT(StrArray))  ALLOCATE(UStr( PRODUCT(nValGather)))
  ELSE
    IF(PRESENT(RealArray)) ALLOCATE(UReal(1))
    IF(PRESENT(IntArray))  ALLOCATE(UInt( 1))
    IF(PRESENT(StrArray))  ALLOCATE(UStr( 1))
  ENDIF

  ! Gather small arrays on IO nodes
  IF(PRESENT(RealArray)) CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
                                          UReal,nDOFPerNode,offsetNode,MPI_DOUBLE_PRECISION,0,MPI_COMM_NODE,iError)
  IF(PRESENT(IntArray))  CALL MPI_GATHERV(IntArray, nDOFLocal,MPI_INTEGER,&
                                          UInt, nDOFPerNode,offsetNode,MPI_INTEGER,0,MPI_COMM_NODE,iError)
  !IF(PRESENT(StrArray))  CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
  !                                        UReal,nDOFPerNode, offsetNode,MPI_DOUBLE_PRECISION,0,MPI_COMM_NODE,iError)

  IF(MPILocalRoot)THEN
    ! Reopen file and write DG solution (only IO nodes)
    CALL OpenDataFile(FileName,create=create,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS)
    IF(PRESENT(RealArray)) CALL WriteArray(DataSetName,rank,nValGlobal,nValGather,&
                                                 offset,collective=collective,RealArray=UReal)
    IF(PRESENT(IntArray))  CALL WriteArray(DataSetName,rank,nValGlobal,nValGather,&
                                                 offset,collective=collective,IntArray =UInt)
    !IF(PRESENT(StrArray))  CALL WriteArray(DataSetName,rank,nValGlobal,nValGather,&
    !                                             offset,collective=collective,StrArr =UStr)
    CALL CloseDataFile()
  END IF

  SDEALLOCATE(UReal)
  SDEALLOCATE(UInt)
  SDEALLOCATE(UStr)
ELSE
#endif
  CALL OpenDataFile(FileName,create=create,single=.FALSE.,readOnly=.FALSE.)
  IF(PRESENT(RealArray)) CALL WriteArray(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,RealArray=RealArray)
  IF(PRESENT(IntArray))  CALL WriteArray(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,IntArray =IntArray)
  IF(PRESENT(StrArray))  CALL WriteArray(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,StrArray =StrArray)
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif

END SUBROUTINE GatheredWriteArray


!==================================================================================================================================
!> Write additional data for analyze purpose to HDF5.
!> The data is taken from a lists, containing either pointers to data arrays or pointers
!> to functions to generate the data, along with the respective varnames.
!>
!> Two options are available:
!>    1. WriteAdditionalElemData:
!>       Element-wise scalar data, e.g. the timestep or indicators.
!>       The data is collected in a single array and written out in one step.
!>       DO NOT MISUSE NODAL DATA FOR THIS! IT WILL DRASTICALLY INCREASE FILE SIZE AND SLOW DOWN IO!
!>    2. WriteAdditionalFieldData:
!>       Nodal data, e.g. coordinates or sgs viscosities.
!>       Each list entry is written into a separate array.
!>
!> TODO:
!>    1. Writing items separatly is slow. Maybe use multiwrite features of coming HDF5.
!>    2. Reorder dimensions, so nVar is last for all arrays.
!==================================================================================================================================
SUBROUTINE WriteAdditionalElemData(FileName,ElemList)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY: offsetElem,nGlobalElems,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)        :: FileName !< Name of the file to be written to
TYPE(tElementOut),POINTER,INTENT(IN) :: ElemList !< Linked list of arrays to write to file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: VarNames(:)
REAL,ALLOCATABLE               :: ElemData(:,:)
INTEGER                        :: nVar
TYPE(tElementOut),POINTER      :: e
!==================================================================================================================================
IF(.NOT. ASSOCIATED(ElemList)) RETURN

! Count the additional variables
nVar = 0
e=>ElemList
DO WHILE(ASSOCIATED(e))
  nVar=nVar+1
  e=>e%next
END DO

! Allocate variable names and data array
ALLOCATE(VarNames(nVar))
ALLOCATE(ElemData(nVar,nElems))

! Fill the arrays
nVar = 0
e=>ElemList
DO WHILE(ASSOCIATED(e))
  nVar=nVar+1
  VarNames(nVar)=e%VarName
  IF(ASSOCIATED(e%RealArray))  ElemData(nVar,:)=e%RealArray
  IF(ASSOCIATED(e%RealScalar)) ElemData(nVar,:)=e%RealScalar
  IF(ASSOCIATED(e%IntArray))   ElemData(nVar,:)=REAL(e%IntArray)
  IF(ASSOCIATED(e%IntScalar))  ElemData(nVar,:)=REAL(e%IntScalar)
  IF(ASSOCIATED(e%eval))       CALL e%eval(ElemData(nVar,:)) ! function fills elemdata
  e=>e%next
END DO

IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttribute(File_ID,'VarNamesAdd',nVar,StrArray=VarNames)
  CALL CloseDataFile()
END IF
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='ElemData', rank=2,  &
                        nValGlobal=(/nVar,nGlobalElems/),&
                        nVal=      (/nVar,nElems      /),&
                        offset=    (/0   ,offSetElem  /),&
                        collective=.TRUE.,RealArray=ElemData)
DEALLOCATE(ElemData,VarNames)
END SUBROUTINE WriteAdditionalElemData


!==================================================================================================================================
!> Comparable to WriteAdditionalElemData, but for field data (rank 5 arrays, where the last dimension is 1:nElems)
!> See also general comment of WriteAdditionalElemData.
!> All arrays that are of the same size as the DG solution will be written to a single dataset, since it is a lot faster than
!> writing several datasets. All arrays with a different size will be written separately. Also the optional doSeparateOutput
!> flag can be used to force the output to a separate dataset.
!==================================================================================================================================
SUBROUTINE WriteAdditionalFieldData(FileName,FieldList)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY: offsetElem,nGlobalElems,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)      :: FileName  !< Name of the file to be written to
TYPE(tFieldOut),POINTER,INTENT(IN) :: FieldList !< Linked list of arrays to write to file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: VarNames(:)
REAL,ALLOCATABLE,TARGET        :: tmp(:,:,:,:,:)
REAL,POINTER                   :: NodeData(:,:,:,:,:)
INTEGER                        :: nVar,nVarTotal
TYPE(tFieldOut),POINTER        :: f
!==================================================================================================================================
! TODO: Perform one write for each dataset.
IF(.NOT. ASSOCIATED(FieldList)) RETURN

! Count fixed size and total number of entries
nVar=0
nVarTotal=0
f=>FieldList
DO WHILE(ASSOCIATED(f))
  IF(.NOT.f%doSeparateOutput) nVar=nVar+f%nVal(1)
  nVarTotal=nVarTotal+f%nVal(1)
  f=>f%next
END DO

! --------------------------------------------------------------------------------------------- !
! First the variable size arrays or arrays that should always be written as a separate dataset
! --------------------------------------------------------------------------------------------- !
! Write the attributes
IF(MPIRoot.AND.(nVarTotal.NE.nVar))THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  f=>FieldList
  DO WHILE(ASSOCIATED(f))
  IF(f%doSeparateOutput) CALL WriteAttribute(File_ID,f%DataSetName,f%nVal(1),StrArray=f%VarNames)
    f=>f%next
  END DO
  CALL CloseDataFile()
END IF

! Write the arrays
f=>FieldList
DO WHILE(ASSOCIATED(f))
  IF(f%doSeparateOutput)THEN
    IF(ASSOCIATED(f%RealArray)) THEN ! real array
      NodeData=>f%RealArray
    ELSE IF(ASSOCIATED(f%Eval)) THEN ! eval function
      ALLOCATE(tmp(f%nVal(1),f%nVal(2),f%nVal(3),f%nVal(4),nElems))
      CALL f%eval(tmp)
      NodeData=>tmp
    END IF
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
                            DataSetName=f%DatasetName, rank=5, &
                            nValGlobal=(/f%nVal,nGlobalElems/),&
                            nVal=      (/f%nVal,nElems      /),&
                            offset=    (/0,0,0,0,  offsetElem  /),&
                            collective=.TRUE.,RealArray=NodeData)
    IF(ASSOCIATED(f%Eval)) DEALLOCATE(tmp)
  END IF
  f=>f%next
END DO


! --------------------------------------------------------------------------------------------- !
! Now process arrays with standard size PP_N
! --------------------------------------------------------------------------------------------- !
IF(nVar.LE.0) RETURN ! no standard data present

ALLOCATE(VarNames(nVar))
ALLOCATE(tmp(nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

! Write the attributes
IF(MPIRoot)THEN
  nVar=0
  f=>FieldList
  DO WHILE(ASSOCIATED(f))
    IF(.NOT.f%doSeparateOutput)THEN
      VarNames(nVar+1:nVar+f%nVal(1))=f%VarNames
      nVar=nVar+f%nVal(1)
    END IF
    f=>f%next
  END DO
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttribute(File_ID,'VarNamesAddField',nVar,StrArray=VarNames)
  CALL CloseDataFile()
END IF

! Collect all fixed size arrays in one array
nVar=0
f=>FieldList
DO WHILE(ASSOCIATED(f))
  IF(.NOT.f%doSeparateOutput)THEN
    IF(ASSOCIATED(f%RealArray))THEN ! real array
      tmp(nVar+1:nVar+f%nVal(1),:,:,:,:)=f%RealArray
    ELSEIF(ASSOCIATED(f%Eval))THEN  ! eval function
      CALL f%Eval(tmp(nVar+1:nVar+f%nVal(1),:,:,:,:))
    END IF
    nVar=nVar+f%nVal(1)
  ENDIF
  f=>f%next
END DO
! Write the arrays (fixed size)
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='FieldData', rank=5,  &
                        nValGlobal=(/nVar,PP_N+1,PP_N+1,PP_NZ+1,nGlobalElems/),&
                        nVal=      (/nVar,PP_N+1,PP_N+1,PP_NZ+1,nElems      /),&
                        offset=    (/0   ,0     ,0     ,0     ,offsetElem  /),&
                        collective=.TRUE.,RealArray=tmp)
DEALLOCATE(VarNames,tmp)

END SUBROUTINE WriteAdditionalFieldData


!==================================================================================================================================
!> Subroutine to write the baseflow to HDF5 format
!==================================================================================================================================
SUBROUTINE WriteBaseflow(MeshFileName,OutputTime)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars  ,ONLY: ProjectName
USE MOD_Mesh_Vars    ,ONLY: offsetElem,nGlobalElems,nElems
USE MOD_Sponge_Vars  ,ONLY: SpBaseFlow
USE MOD_Equation_Vars,ONLY: StrVarNames
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName       !< Name of mesh file
REAL,INTENT(IN)                :: OutputTime         !< Time of output
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
REAL                           :: StartT,EndT
REAL,POINTER                   :: UOut(:,:,:,:,:)
INTEGER                        :: NZ_loc
#if PP_dim == 2
INTEGER                        :: iElem,i,j,iVar
#endif
!==================================================================================================================================
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE BASE FLOW TO HDF5 FILE...'
  GETTIME(StartT)
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_BaseFlow',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton(TRIM(FileName),'BaseFlow',PP_nVar,PP_N,StrVarNames,MeshFileName,OutputTime)

#if PP_dim == 3
  UOut => SpBaseFlow
  NZ_loc=PP_N
#else
IF (.NOT.output2D) THEN
  ALLOCATE(UOut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
  DO iElem=1,nElems
    DO j=0,PP_N; DO i=0,PP_N
      DO iVar=1,PP_nVar
        UOut(iVar,i,j,:,iElem)=SpBaseFlow(iVar,i,j,0,iElem)
      END DO ! iVar=1,PP_nVar
    END DO; END DO
  END DO
  NZ_loc=PP_N
ELSE
  UOut => SpBaseFlow
  NZ_loc=0
END IF
#endif

! Write DG solution
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,NZ_loc+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,NZ_loc+1,nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE., RealArray=UOut)

#if PP_dim == 2
IF(.NOT.output2D) DEALLOCATE(UOut)
#endif
IF(MPIRoot)THEN
  CALL MarkWriteSuccessfull(FileName)
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
END SUBROUTINE WriteBaseflow


!==================================================================================================================================
!> Subroutine to write time averaged data and fluctuations HDF5 format
!==================================================================================================================================
SUBROUTINE WriteTimeAverage(MeshFileName,OutputTime,dtAvg,FV_Elems_In,nVal,&
                            nVarAvg,VarNamesAvg,UAvg,&
                            nVarFluc,VarNamesFluc,UFluc,&
                            FileName_In,FutureTime)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,ONLY: ProjectName
USE MOD_Mesh_Vars  ,ONLY: offsetElem,nGlobalElems,nElems
USE MOD_2D         ,ONLY: ExpandArrayTo3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVarAvg                                      !< Dimension of UAvg
INTEGER,INTENT(IN)             :: nVarFluc                                     !< Dimension of UAvg
INTEGER,INTENT(IN)             :: nVal(3)                                      !< Dimension of UAvg
INTEGER,INTENT(IN)             :: FV_Elems_In(nElems)                          !< Array with custom FV_Elem information
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName                                 !< Name of mesh file
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesAvg(nVarAvg)                         !< Average variable names
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesFluc(nVarFluc)                       !< Average variable names
REAL,INTENT(IN)                :: OutputTime                                   !< Time of output
REAL,INTENT(IN)                :: dtAvg                                        !< Timestep of averaging
REAL,INTENT(IN),TARGET         :: UAvg(nVarAvg,nVal(1),nVal(2),nVal(3),nElems) !< Averaged Solution
REAL,INTENT(IN),TARGET         :: UFluc(nVarFluc,nVal(1),nVal(2),nVal(3),nElems) !< Averaged Solution
REAL,INTENT(IN),OPTIONAL       :: FutureTime                                   !< Time of next output
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: Filename_In                            !< custom filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName,DataSet
REAL                           :: StartT,EndT
REAL,POINTER                   :: UOut(:,:,:,:,:)
#if PP_dim == 2
REAL,POINTER                   :: UOut2D(:,:,:,:,:)
#endif
TYPE(tElementOut),POINTER      :: ElementOutTimeAvg
INTEGER                        :: nVar_loc, nVal_loc(5), nVal_glob(5), i
!==================================================================================================================================
IF(ANY(nVal(1:PP_dim).EQ.0)) RETURN ! no time averaging
IF(nVarAvg.EQ.0.AND.nVarFluc.EQ.0) RETURN ! no time averaging
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE TIME AVERAGED STATE TO HDF5 FILE...'
  GETTIME(StartT)
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_TimeAvg',OutputTime))//'.h5'
IF(PRESENT(Filename_In)) Filename=TRIM(Filename_In)

! Write time averaged data --------------------------------------------------------------------------------------------------------
IF(MPIRoot)THEN
                    CALL GenerateFileSkeleton(TRIM(FileName),'TimeAvg',1 ,PP_N,(/'DUMMY_DO_NOT_VISUALIZE'/),&
                           MeshFileName,OutputTime,FutureTime,create=.TRUE.) ! dummy DG_Solution to fix Posti error, tres oegly !!!
  IF(nVarAvg .GT.0) CALL GenerateFileSkeleton(TRIM(FileName),'TimeAvg',nVarAvg ,PP_N,VarNamesAvg,&
                           MeshFileName,OutputTime,FutureTime,create=.FALSE.,Dataset='Mean')
  IF(nVarFluc.GT.0) CALL GenerateFileSkeleton(TRIM(FileName),'TimeAvg',nVarFluc,PP_N,VarNamesFluc,&
                           MeshFileName,OutputTime,FutureTime,create=.FALSE.,Dataset='MeanSquare')

  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttribute(File_ID,'AvgTime',1,RealScalar=dtAvg)
  CALL CloseDataFile()
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif

! write dummy FV array
NULLIFY(ElementOutTimeAvg)
CALL AddToElemData(ElementOutTimeAvg,'FV_Elems',IntArray=FV_Elems_In)
CALL WriteAdditionalElemData(FileName,ElementOutTimeAvg)
DEALLOCATE(ElementOutTimeAvg)

DO i=1,2
  nVar_loc =  MERGE(nVarAvg,nVarFluc,i.EQ.1)
  IF(nVar_loc.EQ.0) CYCLE
  DataSet  =  MERGE('Mean      ','MeanSquare',i.EQ.1)
  IF(i.EQ.1)THEN
    UOut   => UAvg
  ELSE
    UOut   => UFluc
  END IF
  nVal_loc =  (/nVar_loc,nVal,nElems/)
#if PP_dim == 2
  IF (.NOT.output2D) THEN
    ! If the output should be done with a full third dimension in a two dimensional computation, we need to expand the solution
    NULLIFY(UOut2D)
    ALLOCATE(UOut2D(nVal_loc(1),nVal_loc(2),nVal_loc(3),nVal(1),nVal_loc(5)))
    CALL ExpandArrayTo3D(5,nVal_loc,4,nVal(1),UOut,UOut2D)
    nVal_loc(4)=nVal(1)
    UOut=>UOut2D
  END IF
#endif
  nVal_glob=  (/nVal_loc(1:4),nGlobalElems/)

  ! Reopen file and write DG solution
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName=TRIM(DataSet), rank=5,&
                          nValGlobal=nVal_glob,&
                          nVal=      nVal_loc,&
                          offset=    (/0,0,0,0,offsetElem/),&
                          collective=.TRUE., RealArray=UOut)
#if PP_dim == 2
  ! Deallocate UOut only if we did not point to UAvg
  IF(.NOT.output2D) DEALLOCATE(UOut2D)
#endif
END DO

IF(MPIROOT) CALL MarkWriteSuccessfull(FileName)

IF(MPIROOT)THEN
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
END SUBROUTINE WriteTimeAverage


!==================================================================================================================================
!> Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!==================================================================================================================================
SUBROUTINE GenerateFileSkeleton(FileName,TypeString,nVar,NData,StrVarNames,MeshFileName,OutputTime,&
                                FutureTime,Dataset,create,withUserblock)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars        ,ONLY: ProjectName,UserBlockTmpFile,userblock_total_len
USE MOD_Mesh_Vars          ,ONLY: nGlobalElems
USE MOD_Interpolation_Vars ,ONLY: NodeType
#if FV_ENABLED
USE MOD_FV_Vars            ,ONLY: FV_X,FV_w
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName           !< Name of file to create
CHARACTER(LEN=*),INTENT(IN)    :: TypeString         !< Type of file to be created (state,timeaverage etc.)
INTEGER,INTENT(IN)             :: nVar               !< Number of variables
INTEGER,INTENT(IN)             :: NData              !< Polynomial degree of data
CHARACTER(LEN=*)               :: StrVarNames(nVar)  !< Variabel names
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName       !< Name of mesh file
REAL,INTENT(IN)                :: OutputTime         !< Time of output
REAL,INTENT(IN),OPTIONAL       :: FutureTime         !< Time of next output
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: Dataset      !< Name of the dataset
LOGICAL,INTENT(IN),OPTIONAL    :: create             !< specify whether file should be newly created
LOGICAL,INTENT(IN),OPTIONAL    :: withUserblock      !< specify whether userblock data shall be written or not
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: DSet_ID,FileSpace,HDF5DataType
INTEGER(HSIZE_T)               :: Dimsf(5)
CHARACTER(LEN=255)             :: MeshFile255
CHARACTER(LEN=255)             :: tmp255
CHARACTER(LEN=255)             :: Dataset_Str,Varname_Str
#if FV_ENABLED
REAL                           :: FV_w_array(0:PP_N)
#endif
LOGICAL                        :: withUserblock_loc,create_loc
!==================================================================================================================================
! Create file
create_loc=.TRUE.
withUserblock_loc=.FALSE.
IF(PRESENT(create))                       create_loc       =create
IF(PRESENT(withUserblock).AND.create_loc) withUserblock_loc=withUserblock
Dataset_Str='DG_Solution'
Varname_Str='VarNames'
IF(PRESENT(Dataset))THEN
  Dataset_Str=TRIM(Dataset)
  Varname_Str='VarNames_'//TRIM(DataSet)
END IF

CALL OpenDataFile(TRIM(FileName),create=create_loc,single=.TRUE.,readOnly=.FALSE.,&
                  userblockSize=MERGE(userblock_total_len,0,withUserblock_loc))

! Preallocate the data space for the dataset.
IF(output2D) THEN
  Dimsf=(/nVar,NData+1,NData+1,1,nGlobalElems/)
ELSE
  Dimsf=(/nVar,NData+1,NData+1,NData+1,nGlobalElems/)
END IF

CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_DOUBLE
CALL H5DCREATE_F(File_ID,TRIM(Dataset_Str), HDF5DataType, FileSpace, DSet_ID, iError)
! Close the filespace and the dataset
CALL H5DCLOSE_F(Dset_id, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL WriteAttribute(File_ID,TRIM(Varname_Str),nVar,StrArray=StrVarNames)

! Write default attributes only if file is created
IF(create_loc)THEN

  ! Write file header
  CALL WriteHeader(TRIM(TypeString),File_ID)

  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttribute(File_ID,'N',1,IntScalar=PP_N)
  CALL WriteAttribute(File_ID,'Dimension',1,IntScalar=PP_dim)
  CALL WriteAttribute(File_ID,'Time',1,RealScalar=OutputTime)
  tmp255=TRIM(MeshFileName)
  CALL WriteAttribute(File_ID,'MeshFile',1,StrScalar=(/tmp255/))
  IF(PRESENT(FutureTime))THEN
    MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),FutureTime))//'.h5'
    CALL WriteAttribute(File_ID,'NextFile',1,StrScalar=(/MeshFile255/))
  END IF
  tmp255=TRIM(NodeType)
  CALL WriteAttribute(File_ID,'NodeType',1,StrScalar=(/tmp255/))
#if FV_ENABLED
  CALL WriteAttribute(File_ID,'FV_Type',1,IntScalar=2)
  CALL WriteAttribute(File_ID,'FV_X',PP_N+1,RealArray=FV_X)
  FV_w_array(:)= FV_w
  CALL WriteAttribute(File_ID,'FV_w',PP_N+1,RealArray=FV_w_array)
#endif

  CALL WriteAttribute(File_ID,'NComputation',1,IntScalar=PP_N)
END IF

CALL CloseDataFile()

! Add userblock to hdf5-file (only if create)
IF(withUserblock_loc) CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)

END SUBROUTINE GenerateFileSkeleton


!==================================================================================================================================
!> Add time attribute, after all relevant data has been written to a file,
!> to indicate the writing process has been finished successfully
!==================================================================================================================================
SUBROUTINE MarkWriteSuccessfull(FileName)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName           !< Name of the file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Time(8)
!==================================================================================================================================
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
CALL DATE_AND_TIME(VALUES=time)
CALL WriteAttribute(File_ID,'TIME',8,IntArray=time)
CALL CloseDataFile()
END SUBROUTINE MarkWriteSuccessfull


!==================================================================================================================================
!> Deletes all HDF5 output files, beginning from time Flushtime. Used for cleanup at the beginning of a new simulation
!==================================================================================================================================
SUBROUTINE FlushFiles(FlushTime_In)
! MODULES
USE MOD_Globals
USE MOD_Output_Vars ,ONLY: ProjectName
USE MOD_HDF5_Input  ,ONLY: GetNextFileName
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN),OPTIONAL :: FlushTime_In     !< Time to start flush
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: stat,ioUnit
REAL                     :: FlushTime
CHARACTER(LEN=255)       :: FileName,InputFile,NextFile
!==================================================================================================================================
IF(.NOT.MPIRoot) RETURN

WRITE(UNIT_stdOut,'(a)')' DELETING OLD HDF5 FILES...'
IF (.NOT.PRESENT(FlushTime_In)) THEN
  FlushTime=0.0
ELSE
  FlushTime=FlushTime_In
END IF
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',FlushTime))//'.h5'

! Delete state files
InputFile=TRIM(FileName)
! Read calculation time from file
CALL GetNextFileName(Inputfile,NextFile,.TRUE.)
! Delete File - only root
stat=0
OPEN ( NEWUNIT= ioUnit,         &
       FILE   = InputFile,      &
       STATUS = 'OLD',          &
       ACTION = 'WRITE',        &
       ACCESS = 'SEQUENTIAL',   &
       IOSTAT = stat          )
IF(stat .EQ. 0) CLOSE ( ioUnit,STATUS = 'DELETE' )
DO
  InputFile=TRIM(NextFile)
  ! Read calculation time from file
  CALL GetNextFileName(Inputfile,NextFile,.TRUE.)
  ! Delete File - only root
  stat=0
  OPEN ( NEWUNIT= ioUnit,         &
         FILE   = InputFile,      &
         STATUS = 'OLD',          &
         ACTION = 'WRITE',        &
         ACCESS = 'SEQUENTIAL',   &
         IOSTAT = stat          )
  IF(stat .EQ. 0) CLOSE ( ioUnit,STATUS = 'DELETE' )
  IF(iError.NE.0) EXIT  ! iError is set in GetNextFileName !
END DO

WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'

END SUBROUTINE FlushFiles


!==================================================================================================================================
!> Subroutine to write a distinct file header to each HDF5 file
!==================================================================================================================================
SUBROUTINE WriteHeader(FileType_in,File_ID)
! MODULES
USE MOD_Output_Vars,ONLY:ProgramName,FileVersion,ProjectName
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)              :: FileType_in   !< Type of file (e.g. state, timeaverage)
INTEGER(HID_T),INTENT(IN)                :: File_ID       !< HDF5 file id
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255) :: tmp255
!==================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files
! Attributes are program name, file type identifier, project name and version number
tmp255=TRIM(ProgramName)
CALL WriteAttribute(File_ID,'Program'     ,1,StrScalar=(/tmp255/))
tmp255=TRIM(FileType_in)
CALL WriteAttribute(File_ID,'File_Type'   ,1,StrScalar=(/tmp255/))
tmp255=TRIM(ProjectName)
CALL WriteAttribute(File_ID,'Project_Name',1,StrScalar=(/tmp255/))
CALL WriteAttribute(File_ID,'File_Version',1,RealScalar=FileVersion)
END SUBROUTINE WriteHeader


!==================================================================================================================================
!> Low-level subroutine to actually write data to HDF5 format
!==================================================================================================================================
SUBROUTINE WriteArray(DataSetName,rank,nValGlobal,nVal,offset,&
                            collective,resizeDim,chunkSize,&
                            RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: DataSetName      !< name of the dataset to write the data into
INTEGER,INTENT(IN)            :: rank             !< number of dimensions of the array
INTEGER,INTENT(IN)            :: nValGlobal(rank) !< max size of array in offset dimension
INTEGER,INTENT(IN)            :: nVal(rank)       !< size of complete (local) array to write
INTEGER,INTENT(IN)            :: offset(rank)     !< offset =0, start at beginning of the array
LOGICAL,INTENT(IN)            :: collective       !< use collective writes from all procs
LOGICAL,INTENT(IN),OPTIONAL   :: resizeDim(rank)  !< specify dimensions which can be resized (enlarged)
INTEGER,INTENT(IN),OPTIONAL   :: chunkSize(rank)  !< specify chunksize
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal)) !< number of array entries
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray(PRODUCT(nVal))  !< number of array entries
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(PRODUCT(nVal))  !< number of array entries (length 255)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,Type_ID,dsetparams
INTEGER(HSIZE_T)               :: Dimsf(Rank),OffsetHDF(Rank),nValMax(Rank)
INTEGER(SIZE_T)                :: SizeSet=255
LOGICAL                        :: chunky
TYPE(C_PTR)                    :: buf
!==================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(DataSetName),'" TO HDF5 FILE...'

! specify chunk size if desired
nValMax=nValGlobal
chunky=.FALSE.
CALL H5PCREATE_F(H5P_DATASET_CREATE_F,dsetparams,iError)
IF(PRESENT(chunkSize))THEN
  chunky=.TRUE.
  Dimsf=chunkSize
  CALL H5PSET_CHUNK_F(dsetparams,rank,dimsf,iError)
END IF
! make array extendable in case you want to append something
IF(PRESENT(resizeDim))THEN
  IF(.NOT.PRESENT(chunkSize))&
    CALL abort(__STAMP__,&
               'Chunk size has to be specified when using resizable arrays.')
  nValMax = MERGE(H5S_UNLIMITED_F,nValMax,resizeDim)
END IF

! Create the dataset with default properties.
IF(PRESENT(RealArray)) Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntArray))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, SizeSet, iError)
END IF

Dimsf = nValGlobal ! we need the global array size
CALL H5ESET_AUTO_F(0,iError)
CALL H5DOPEN_F(File_ID, TRIM(DatasetName),DSet_ID, iError)
IF(iError.NE.0)THEN ! does not exist
  ! Create the data space for the  dataset.
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError, nValMax)
  CALL H5DCREATE_F(File_ID, TRIM(DataSetName), Type_ID, FileSpace, DSet_ID,iError,dsetparams)
  CALL H5SCLOSE_F(FileSpace, iError)
END IF
CALL H5ESET_AUTO_F(1,iError)
IF(chunky)THEN
  CALL H5DSET_EXTENT_F(DSet_ID,Dimsf,iError) ! if resizable then dataset may need to be extended
END IF

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf=nVal  ! Now we need the local array size
OffsetHDF = Offset
! Create the data space in the memory
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
ELSE
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
END IF
! Select hyperslab in the file.
CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SSELECT_NONE_F(FileSpace,iError)
ELSE
  CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, OffsetHDF, Dimsf, iError)
END IF

! Create property list for collective dataset write
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#if USE_MPI
IF(collective)THEN
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F,  iError)
ELSE
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)
END IF
#endif

!Write the dataset collectively.
IF(PRESENT(IntArray))  buf=C_LOC(IntArray)
IF(PRESENT(RealArray)) buf=C_LOC(RealArray)
IF(PRESENT(StrArray))  buf=C_LOC(StrArray(1))
CALL H5DWRITE_F(DSet_ID,Type_ID,buf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)

IF(PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close the property list, dataspaces and dataset.
CALL H5PCLOSE_F(dsetparams, iError)
CALL H5PCLOSE_F(PList_ID, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteArray


!==================================================================================================================================
!> Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,groupID. This must be opened
!> outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!==================================================================================================================================
SUBROUTINE WriteAttribute(Loc_ID_in,AttribName,nVal,DataSetname,&
                          RealScalar,IntScalar,StrScalar,LogicalScalar, &
                          RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(HID_T)    ,INTENT(IN)           :: Loc_ID_in              !< Dataset ID (only if already open)
CHARACTER(LEN=*)  ,INTENT(IN)           :: AttribName             !< name of the attribute to be written
INTEGER           ,INTENT(IN)           :: nVal                   !< number of array entries if array is written
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL  :: DatasetName            !< name of the dataset created
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealScalar       !< real scalar
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntScalar        !< integer scalar
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrScalar(1)     !< scalar string
LOGICAL           ,INTENT(IN),OPTIONAL        :: LogicalScalar    !< logical scalar
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(nVal)  !< real array of length nVal
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray(nVal)   !< integer array of length nVal
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL,TARGET :: StrArray(nVal)   !< string array of length nVal
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Rank
INTEGER(HID_T)                 :: DataSpace,Attr_ID,Loc_ID,Type_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER(SIZE_T)                :: AttrLen
INTEGER,TARGET                 :: logtoint
TYPE(C_PTR)                    :: buf
INTEGER                        :: hdferr
!==================================================================================================================================
LOGWRITE(*,*)' WRITE ATTRIBUTE "',TRIM(AttribName),'" TO HDF5 FILE...'
IF(PRESENT(DataSetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
ELSE
  Loc_ID=Loc_ID_in
END IF
! Create scalar data space for the attribute.
Rank=1
Dimsf(:)=0 !???
Dimsf(1)=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, DataSpace, iError)
! Create the attribute for group Loc_ID.
IF(PRESENT(RealScalar)) Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(RealArray))  Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntScalar))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntArray))   Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(LogicalScalar))THEN
  LogToInt=MERGE(1,0,LogicalScalar)
  Type_ID=H5T_NATIVE_INTEGER
END IF

! Create character string datatype for the attribute.
! For a attribute character, we have to build our own type with corresponding attribute length
IF(PRESENT(StrScalar))THEN
  AttrLen=LEN_TRIM(StrScalar(1))
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
END IF
IF(PRESENT(StrArray))THEN
  AttrLen=255
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
ENDIF

CALL H5ESET_AUTO_F(0, hdferr)
CALL H5AOPEN_F(    Loc_ID, TRIM(AttribName), Attr_ID, iError)
IF(iError.EQ.0)THEN
  CALL H5ACLOSE_F(Attr_ID, iError)
  CALL H5ADELETE_F(Loc_ID, TRIM(AttribName)         , iError)
END IF
CALL H5ESET_AUTO_F(1, hdferr)
CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), Type_ID, DataSpace, Attr_ID, iError)
IF(iError.NE.0) STOP 'Could not open or create attribute!'

! Write the attribute data.
buf=C_NULL_PTR
IF(PRESENT(RealArray))     buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))    buf=C_LOC(RealScalar)
IF(PRESENT(IntArray))      buf=C_LOC(IntArray)
IF(PRESENT(IntScalar))     buf=C_LOC(IntScalar)
IF(PRESENT(LogicalScalar)) buf=C_LOC(LogToInt)
IF(PRESENT(StrScalar))     buf=C_LOC(StrScalar(1))
IF(PRESENT(StrArray))      buf=C_LOC(StrArray(1))
IF(C_ASSOCIATED(buf))&
  CALL H5AWRITE_F(Attr_ID, Type_ID, buf, iError)

! Close datatype
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close dataspace
CALL H5SCLOSE_F(DataSpace, iError)
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteAttribute

END MODULE MOD_HDF5_output
