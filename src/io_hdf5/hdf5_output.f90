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
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

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

INTERFACE
  SUBROUTINE copy_userblock(outfilename,infilename) BIND(C)
      USE ISO_C_BINDING, ONLY: C_CHAR
      CHARACTER(KIND=C_CHAR) :: outfilename(*)
      CHARACTER(KIND=C_CHAR) :: infilename(*)
  END SUBROUTINE copy_userblock
END INTERFACE


PUBLIC :: WriteState,FlushFiles,WriteHeader,WriteTimeAverage,WriteBaseflow
PUBLIC :: WriteArray,WriteAttribute
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
USE MOD_Output_Vars       ,ONLY: ProjectName,NOut,Vdm_N_NOut
USE MOD_Mesh_Vars         ,ONLY: offsetElem,nGlobalElems,sJ,nElems
USE MOD_ChangeBasisByDim  ,ONLY: ChangeBasisVolume
USE MOD_Equation_Vars     ,ONLY: StrVarNames
USE MOD_2D                ,ONLY: ExpandArrayTo3D
#if FV_ENABLED && FV_RECONSTRUCT
USE MOD_FV_Vars           ,ONLY: gradUxi,gradUeta,gradUzeta,FV_dx_XI_L,FV_dx_ETA_L,FV_dx_ZETA_L
USE MOD_FV_Vars           ,ONLY: FV_dx_XI_R,FV_dx_ETA_R,FV_dx_ZETA_R
USE MOD_EOS               ,ONLY: ConsToPrim,PrimToCons
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
REAL,ALLOCATABLE               :: UOutTmp(:,:,:,:,:)
REAL                           :: Utmp(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL                           :: JN(1,0:PP_N,0:PP_N,0:PP_NZ),JOut(1,0:NOut,0:NOut,0:PP_NOutZ)
INTEGER                        :: iElem,i,j,k,iVar,nVal(5)
!#if FV_ENABLED & FV_RECONSTRUCT
!REAL                           :: UPrim(1:PP_nVarPrim)
!REAL                           :: UCons(1:PP_nVar)
!REAL                           :: gradUx(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!REAL                           :: gradUy(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!REAL                           :: gradUz(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!#endif
!==================================================================================================================================
IF(MPIRoot)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE STATE TO HDF5 FILE...'
  GETTIME(StartT)
END IF

!#if FV_ENABLED & FV_RECONSTRUCT
!! transform physical gradients of FV to reference space gradients (easier POSTI)
!DO iElem=1,nElems
  !DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    !CALL ConsToPrim(UPrim ,U(:,i,j,k,iElem))
    !CALL PrimToCons(UPrim+gradUxi(:,j,k,i,iElem)*FV_dx_XI_R(i,j,k,iElem),Ucons)
    !gradUx(:,i,j,k,iElem) = (Ucons-U(:,i,j,k,iElem))/FV_dx_XI_R  (i,j,k,iElem)* &
        !(FV_dx_XI_L  (i,j,k,iElem)+FV_dx_XI_R  (i,j,k,iElem)) * 0.5
    !CALL PrimToCons(UPrim+gradUeta(:,i,k,j,iElem)*FV_dx_ETA_R (i,j,k,iElem),Ucons)
    !gradUy(:,i,j,k,iElem) = (Ucons-U(:,i,j,k,iElem))/FV_dx_ETA_R (i,j,k,iElem)* & 
        !(FV_dx_ETA_L (i,j,k,iElem)+FV_dx_ETA_R (i,j,k,iElem)) * 0.5
    !CALL PrimToCons(UPrim+gradUzeta(:,i,j,k,iElem)*FV_dx_ZETA_R(i,j,k,iElem),Ucons)
    !gradUz(:,i,j,k,iElem) = (Ucons-U(:,i,j,k,iElem))/FV_dx_ZETA_R(i,j,k,iElem)* & 
        !(FV_dx_ZETA_L(i,j,k,iElem)+FV_dx_ZETA_R(i,j,k,iElem)) * 0.5
  !END DO; END DO; END DO! i,j,k=0,PP_N
!END DO ! iElem
!#endif

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileType=MERGE('ERROR_State','State      ',isErrorFile)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileType),OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton(TRIM(FileName),'State',PP_nVar,NOut,StrVarNames,MeshFileName,OutputTime,FutureTime)

! Set size of output 
nVal=(/PP_nVar,NOut+1,NOut+1,PP_NOutZ+1,nElems/)

! build output data
IF(NOut.NE.PP_N)THEN
#if FV_ENABLED
  CALL Abort(__STAMP__, &
      "NOut not working for FV!")
#endif
  ! Project JU and J to NOut, compute U on Nout
  ALLOCATE(UOut(PP_nVar,0:NOut,0:NOut,0:PP_NOutZ,nElems))
  DO iElem=1,nElems
    JN(1,:,:,:)=1./sJ(:,:,:,iElem,0)
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Utmp(:,i,j,k)=U(:,i,j,k,iElem)*JN(1,i,j,k)
      DO iVar=1,PP_nVar
        Utmp(iVar,i,j,:)=U(iVar,i,j,k,iElem)*JN(1,i,j,k)
      END DO ! iVar=1,PP_nVar
    END DO; END DO; END DO
    CALL ChangeBasisVolume(PP_nVar,PP_N,NOut,Vdm_N_NOut,&
                           Utmp,UOut(1:PP_nVar,:,:,:,iElem))
    ! Jacobian
    CALL ChangeBasisVolume(1,PP_N,NOut,Vdm_N_NOut,JN,JOut)
    DO k=0,PP_NOutZ; DO j=0,NOut; DO i=0,NOut
      UOut(:,i,j,k,iElem)=UOut(:,i,j,k,iElem)/JOut(1,i,j,k)
    END DO; END DO; END DO
  END DO
#if PP_dim == 2
  ! If the output should be done with a full third dimension in a two dimensional computation, we need to expand the solution
  IF (.NOT.IO_2D) THEN 
    ALLOCATE(UOutTmp(PP_nVar,0:NOut,0:NOut,0:PP_NOutZ,nElems))
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
  IF (.NOT.IO_2D) THEN
    ! If the output should be done with a full third dimension in a two dimensional computation, we need to expand the solution
    ALLOCATE(UOut(PP_nVar,0:NOut,0:NOut,0:PP_N,nElems))
    CALL ExpandArrayTo3D(5,(/PP_nVar,PP_N+1,PP_N+1,PP_NZ+1,nElems/),4,Nout+1,U,UOut)
    ! Correct size of the output array
    nVal=(/PP_nVar,NOut+1,NOut+1,NOut+1,nElems/)
  ELSE
    UOut => U
  END IF
#endif
END IF ! (NOut.NE.PP_N)


! Reopen file and write DG solution
#if MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/PP_nVar,NOut+1,NOut+1,NOut+1,nGlobalElems/),&
                        nVal=nVal                                              ,&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=UOut)

!#if FV_ENABLED & FV_RECONSTRUCT
!CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        !DataSetName='gradUxi', rank=5,&
                        !nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        !nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nElems/),   &
                        !offset=    (/0,      0,     0,     0,     offsetElem/),  &
                        !collective=.TRUE., RealArray=gradUx)

!CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        !DataSetName='gradUeta', rank=5,&
                        !nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        !nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nElems/),   &
                        !offset=    (/0,      0,     0,     0,     offsetElem/),  &
                        !collective=.TRUE., RealArray=gradUy)

!CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        !DataSetName='gradUzeta', rank=5,&
                        !nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        !nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nElems/),   &
                        !offset=    (/0,      0,     0,     0,     offsetElem/),  &
                        !collective=.TRUE., RealArray=gradUz)
!#endif
                    
! Deallocate UOut only if we did not point to U
IF((PP_N .NE. NOut).OR.((PP_dim .EQ. 2).AND.(.NOT.IO_2D))) DEALLOCATE(UOut)

CALL WriteAdditionalElemData(FileName,ElementOut)
CALL WriteAdditionalFieldData(FileName,FieldOut)

IF(MPIRoot)THEN
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
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
#if MPI
REAL,              ALLOCATABLE :: UReal(:)
CHARACTER(LEN=255),ALLOCATABLE :: UStr(:)
INTEGER,           ALLOCATABLE :: UInt(:)
INTEGER                        :: i,nValGather(rank),nDOFLocal
INTEGER,DIMENSION(nLocalProcs) :: nDOFPerNode,offsetNode
!==================================================================================================================================
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
#if MPI
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
CHARACTER(LEN=255),INTENT(IN)  :: FileName           !< Name of the file to be written to
TYPE(tElementOut),POINTER,INTENT(IN) :: ElemList     !< Linked list of arrays to write to file
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
!==================================================================================================================================
SUBROUTINE WriteAdditionalFieldData(FileName,FieldList)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY: offsetElem,nGlobalElems,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName          !< Name of the file to be written to
TYPE(tFieldOut),POINTER,INTENT(IN) :: FieldList     !< Linked list of arrays to write to file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: VarNames(:)
REAL,ALLOCATABLE,TARGET        :: tmp(:,:,:,:,:)
REAL,POINTER                   :: NodeData(:,:,:,:,:)
INTEGER                        :: nVar,nVarTotal
TYPE(tFieldOut),POINTER        :: f
INTEGER                        :: mask(3)
!==================================================================================================================================
! TODO: Perform one write for each dataset.
IF(.NOT. ASSOCIATED(FieldList)) RETURN

#if PP_dim == 3
mask=(/PP_N+1,PP_N+1,PP_N+1/)
#else
mask=(/PP_N+1,PP_N+1,1/)
#endif

! Count fixed size and total number of entries 
nVar=0
nVarTotal=0
f=>FieldList
DO WHILE(ASSOCIATED(f))
  IF(ALL(f%nVal(2:4).EQ.mask)) nVar=nVar+f%nVal(1)
  nVarTotal=nVarTotal+f%nVal(1)
  f=>f%next
END DO

! First the variable size arrays 
! Write the attributes (variable size)
IF(MPIRoot.AND.(nVarTotal.NE.nVar))THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  f=>FieldList
  DO WHILE(ASSOCIATED(f))
    IF(ANY(f%nVal(2:4).NE.PP_N+1))&
      CALL WriteAttribute(File_ID,f%DataSetName,f%nVal(1),StrArray=f%VarNames)
    f=>f%next
  END DO
  CALL CloseDataFile()
END IF

! Write the arrays (variable size)
f=>FieldList
DO WHILE(ASSOCIATED(f))
  IF(.NOT. ALL(f%nVal(2:4).EQ.mask))THEN ! not fixed size
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


! Now process arrays with standard size PP_N
IF(nVar.LE.0) RETURN ! no standard data present

ALLOCATE(VarNames(nVar))
ALLOCATE(tmp(nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

! Write the attributes (fixed size)
IF(MPIRoot)THEN
  nVar=0
  f=>FieldList
  DO WHILE(ASSOCIATED(f))
    IF(ALL(f%nVal(2:4).EQ.mask))THEN
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
  IF(ALL(f%nVal(2:4).EQ.mask))THEN
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
!==================================================================================================================================
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE BASE FLOW TO HDF5 FILE...'
  GETTIME(StartT)
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_BaseFlow',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton(TRIM(FileName),'BaseFlow',PP_nVar,PP_N,StrVarNames,MeshFileName,OutputTime)

! Write DG solution
#if MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE., RealArray=SpBaseFlow)

IF(MPIRoot)THEN
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
END SUBROUTINE WriteBaseflow


!==================================================================================================================================
!> Subroutine to write time averaged data and fluctuations HDF5 format
!==================================================================================================================================
SUBROUTINE WriteTimeAverage(MeshFileName,OutputTime,FutureTime,VarNamesAvg,VarNamesFluc,UAvg,UFluc,dtAvg,nVar_Avg,nVar_Fluc)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,ONLY: ProjectName
USE MOD_Mesh_Vars  ,ONLY: offsetElem,nGlobalElems,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName                                 !< Name of mesh file
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesAvg(nVar_Avg)                        !< Average variable names
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesFluc(nVar_Fluc)                      !< Fluctuations variable names
REAL,INTENT(IN)                :: OutputTime                                   !< Time of output
REAL,INTENT(IN),OPTIONAL       :: FutureTime                                   !< Time of next output
REAL,INTENT(IN),TARGET         :: UAvg(nVar_Avg,0:PP_N,0:PP_N,0:PP_N,nElems)   !< Averaged Solution
REAL,INTENT(IN),TARGET         :: UFluc(nVar_Fluc,0:PP_N,0:PP_N,0:PP_N,nElems) !< Fluctuations
REAL,INTENT(IN)                :: dtAvg                                        !< Timestep of averaging
INTEGER,INTENT(IN)             :: nVar_Avg                                     !< Number of averaged variables
INTEGER,INTENT(IN)             :: nVar_Fluc                                    !< Number of fluctuations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
REAL                           :: StartT,EndT
!==================================================================================================================================
IF((nVar_Avg.EQ.0).AND.(nVar_Fluc.EQ.0)) RETURN ! no time averaging
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE TIME AVERAGED STATE AND FLUCTUATIONS TO HDF5 FILE...'
  GETTIME(StartT)
END IF

! Write timeaverages ---------------------------------------------------------------------------------------------------------------
IF(nVar_Avg.GT.0)THEN
  ! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_TimeAvg',OutputTime))//'.h5'
  IF(MPIRoot)THEN
    CALL GenerateFileSkeleton(TRIM(FileName),'TimeAvg',nVar_Avg,PP_N,VarNamesAvg,MeshFileName,OutputTime,FutureTime)
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttribute(File_ID,'AvgTime',1,RealScalar=dtAvg)
    CALL CloseDataFile()
  END IF
#if MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

  ! Reopen file and write DG solution
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution', rank=5,&
                          nValGlobal=(/nVar_Avg,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                          nVal=      (/nVar_Avg,PP_N+1,PP_N+1,PP_N+1,nElems/),&
                          offset=    (/0,       0,     0,     0,     offsetElem/),&
                          collective=.TRUE., RealArray=UAvg)
END IF

! Write fluctuations ---------------------------------------------------------------------------------------------------------------
IF(nVar_Fluc.GT.0)THEN
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Fluc',OutputTime))//'.h5'
  IF(MPIRoot)THEN
    CALL GenerateFileSkeleton(TRIM(FileName),'Fluc',nVar_Fluc,PP_N,VarNamesFluc,MeshFileName,OutputTime,FutureTime)
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttribute(File_ID,'AvgTime',1,RealScalar=dtAvg)
    CALL CloseDataFile()
  END IF
#if MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

  ! Reopen file and write DG solution
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution', rank=5,&
                          nValGlobal=(/nVar_Fluc,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                          nVal=      (/nVar_Fluc,PP_N+1,PP_N+1,PP_N+1,nElems/),&
                          offset=    (/0,        0,     0,     0,     offsetElem/),&
                          collective=.TRUE., RealArray=UFluc)
END IF

IF(MPIROOT)THEN
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
END SUBROUTINE WriteTimeAverage


!==================================================================================================================================
!> Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!==================================================================================================================================
SUBROUTINE GenerateFileSkeleton(FileName,TypeString,nVar,NData,StrVarNames,MeshFileName,OutputTime,FutureTime)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars  ,ONLY: ProjectName,UserBlockTmpFile,userblock_total_len
USE MOD_Mesh_Vars  ,ONLY: nGlobalElems
USE MOD_Interpolation_Vars,ONLY: NodeType
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_X,FV_w
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
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: DSet_ID,FileSpace,HDF5DataType
INTEGER(HSIZE_T)               :: Dimsf(5)
CHARACTER(LEN=255)             :: MeshFile255
#if FV_ENABLED
REAL                           :: FV_w_array(0:PP_N)
#endif
!==================================================================================================================================
! Create file
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)

! Write file header
CALL WriteHeader(TRIM(TypeString),File_ID)

! Preallocate the data space for the dataset.
IF(IO_2D) THEN
  Dimsf=(/nVar,NData+1,NData+1,1,nGlobalElems/)
ELSE
  Dimsf=(/nVar,NData+1,NData+1,NData+1,nGlobalElems/)
END IF

CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_DOUBLE
CALL H5DCREATE_F(File_ID,'DG_Solution', HDF5DataType, FileSpace, DSet_ID, iError)
! Close the filespace and the dataset
CALL H5DCLOSE_F(Dset_id, iError)
CALL H5SCLOSE_F(FileSpace, iError)

!#if FV_ENABLED
!Dimsf=(/nVar,NData+1,NData+1,NData+1,nGlobalElems/)
!CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
!! Create the dataset with default properties.
!HDF5DataType=H5T_NATIVE_DOUBLE
!CALL H5DCREATE_F(File_ID,'gradUxi', HDF5DataType, FileSpace, DSet_ID, iError)
!! Close the filespace and the dataset
!CALL H5DCLOSE_F(Dset_id, iError)
!CALL H5SCLOSE_F(FileSpace, iError)
!Dimsf=(/nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/)
!CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
!! Create the dataset with default properties.
!HDF5DataType=H5T_NATIVE_DOUBLE
!CALL H5DCREATE_F(File_ID,'gradUeta', HDF5DataType, FileSpace, DSet_ID, iError)
!! Close the filespace and the dataset
!CALL H5DCLOSE_F(Dset_id, iError)
!CALL H5SCLOSE_F(FileSpace, iError)
!Dimsf=(/nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/)
!CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
!! Create the dataset with default properties.
!HDF5DataType=H5T_NATIVE_DOUBLE
!CALL H5DCREATE_F(File_ID,'gradUzeta', HDF5DataType, FileSpace, DSet_ID, iError)
!! Close the filespace and the dataset
!CALL H5DCLOSE_F(Dset_id, iError)
!CALL H5SCLOSE_F(FileSpace, iError)
!#endif

! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
CALL WriteAttribute(File_ID,'N',1,IntScalar=PP_N)
CALL WriteAttribute(File_ID,'Dimension',1,IntScalar=PP_dim)
CALL WriteAttribute(File_ID,'Time',1,RealScalar=OutputTime)
CALL WriteAttribute(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
IF(PRESENT(FutureTime))THEN
  MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),FutureTime))//'.h5'
  CALL WriteAttribute(File_ID,'NextFile',1,StrScalar=(/MeshFile255/))
END IF
CALL WriteAttribute(File_ID,'NodeType',1,StrScalar=(/NodeType/))
CALL WriteAttribute(File_ID,'VarNames',nVar,StrArray=StrVarNames)
#if FV_ENABLED
CALL WriteAttribute(File_ID,'FV_Type',1,IntScalar=2)
CALL WriteAttribute(File_ID,'FV_X',PP_N+1,RealArray=FV_X)
FV_w_array(:)= FV_w
CALL WriteAttribute(File_ID,'FV_w',PP_N+1,RealArray=FV_w_array)
#endif

CALL WriteAttribute(File_ID,'NComputation',1,IntScalar=PP_N)

CALL CloseDataFile()

! Add userblock to hdf5-file
CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)

END SUBROUTINE GenerateFileSkeleton


!==================================================================================================================================
!> Deletes all HDF5 output files, beginning from time Flushtime. Used for cleanup at the beginning of a new simulation
!==================================================================================================================================
SUBROUTINE FlushFiles(FlushTime_In)
! MODULES
!USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,ONLY:ProjectName
USE MOD_HDF5_Input,ONLY:GetNextFileName
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
ioUnit=GETFREEUNIT()
OPEN ( UNIT   = ioUnit,            &
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
  ioUnit=GETFREEUNIT()
  OPEN ( UNIT   = ioUnit,            &
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
!==================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files
! Attributes are program name, file type identifier, project name and version number
CALL WriteAttribute(File_ID,'Program'     ,1,StrScalar=(/TRIM(ProgramName)/))
CALL WriteAttribute(File_ID,'File_Type'   ,1,StrScalar=(/TRIM(FileType_in)/))
CALL WriteAttribute(File_ID,'Project_Name',1,StrScalar=(/TRIM(ProjectName)/))
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
#if MPI
IF(collective)THEN
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F,  iError)
ELSE
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)
END IF
#endif

!Write the dataset collectively.
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(IntArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,IntArray, Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(RealArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,RealArray,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,StrArray, Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
#else
IF(PRESENT(IntArray))  buf=C_LOC(IntArray)
IF(PRESENT(RealArray)) buf=C_LOC(RealArray)
IF(PRESENT(StrArray))  buf=C_LOC(StrArray(1))
CALL H5DWRITE_F(DSet_ID,Type_ID,buf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
#endif /* HDF5_F90 */

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
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL,TARGET :: StrScalar(1)     !< scalar string
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
#ifndef HDF5_F90
TYPE(C_PTR)                    :: buf
#endif
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
  AttrLen=LEN(StrScalar(1))
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
END IF
IF(PRESENT(StrArray))THEN
  AttrLen=LEN(StrArray(1))
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
ENDIF

CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), Type_ID, DataSpace, Attr_ID, iError)
! Write the attribute data.
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(RealArray))     CALL H5AWRITE_F(Attr_ID, Type_ID, RealArray,  Dimsf, iError)
IF(PRESENT(RealScalar))    CALL H5AWRITE_F(Attr_ID, Type_ID, RealScalar, Dimsf, iError)
IF(PRESENT(IntArray))      CALL H5AWRITE_F(Attr_ID, Type_ID, IntArray,   Dimsf, iError)
IF(PRESENT(IntScalar))     CALL H5AWRITE_F(Attr_ID, Type_ID, IntScalar,  Dimsf, iError)
IF(PRESENT(LogicalScalar)) CALL H5AWRITE_F(Attr_ID, Type_ID, LogToInt,   Dimsf, iError)
IF(PRESENT(StrScalar))     CALL H5AWRITE_F(Attr_ID, Type_ID, StrScalar,  Dimsf, iError)
IF(PRESENT(StrArray))      CALL H5AWRITE_F(Attr_ID, Type_ID, StrArray,   Dimsf, iError)
#else /* HDF5_F90 */
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
#endif /* HDF5_F90 */

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
