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

!===================================================================================================================================
!> Module for generic data output in vtk xml fromat
!> WARNING: WriteDataToVTK works only for POSTPROCESSING or for debug output during runtime
!===================================================================================================================================
MODULE MOD_VTK
USE ISO_C_BINDING
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE, BIND(C) :: CARRAY
  INTEGER (C_INT) :: len
  INTEGER (C_INT) :: dim
  TYPE (C_PTR)    :: data
END TYPE CARRAY

INTERFACE WriteDataToVTK
  MODULE PROCEDURE WriteDataToVTK
END INTERFACE

INTERFACE WriteVTKMultiBlockDataSet
  MODULE PROCEDURE WriteVTKMultiBlockDataSet
END INTERFACE

INTERFACE WriteCoordsToVTK_array
  MODULE PROCEDURE WriteCoordsToVTK_array
END INTERFACE

INTERFACE WriteDataToVTK_array
  MODULE PROCEDURE WriteDataToVTK_array
END INTERFACE

INTERFACE WriteVarnamesToVTK_array
  MODULE PROCEDURE WriteVarnamesToVTK_array
END INTERFACE

#if USE_MPI
INTERFACE WriteParallelVTK
  MODULE PROCEDURE WriteParallelVTK
END INTERFACE
#endif /*USE_MPI*/

PUBLIC::WriteDataToVTK
PUBLIC::WriteVTKMultiBlockDataSet
PUBLIC::WriteCoordsToVTK_array
PUBLIC::WriteDataToVTK_array
PUBLIC::WriteVarnamesToVTK_array
#if USE_MPI
PUBLIC::WriteParallelVTK
#endif /*USE_MPI*/
PUBLIC::CARRAY
!===================================================================================================================================

CONTAINS

SUBROUTINE CreateConnectivity(NVisu,nElems,nodeids,dim,DGFV,HighOrder)
! MODULES
USE ISO_C_BINDING
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                       :: NVisu
INTEGER,INTENT(IN)                       :: nElems
INTEGER,ALLOCATABLE,TARGET,INTENT(INOUT) :: nodeids(:)        !< stores the connectivity
INTEGER,INTENT(IN)                       :: dim               !< 3 = 3d connectivity, 2 = 2d connectivity
INTEGER,INTENT(IN)                       :: DGFV              !< flag indicating DG = 0 or FV = 1 data
INTEGER,INTENT(IN)                       :: HighOrder         !< flag indicating high order
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
INTEGER             :: NodeID,NodeIDElem
INTEGER             :: NVisu_k, NVisu_j, NVisu_elem, NVisu_p1_2
INTEGER             :: nVTKCells
! High-order elements
INTEGER,ALLOCATABLE :: NodeIndizes(:,:,:)
!===================================================================================================================================

SELECT CASE(dim)
  CASE(3)
    NVisu_k = NVisu
    NVisu_j = NVisu
  CASE(2)
    NVisu_k = 1
    NVisu_j = NVisu
  CASE(1)
    NVisu_k = 1
    NVisu_j = 1
  CASE DEFAULT
    ! Dummy variables to make GCC happy
    NVisu_k = -1
    NVisu_j = -1
    CALL CollectiveStop(__STAMP__, "Only 2D and 3D connectivity can be created. dim must be 2 or 3.")
END SELECT

NVisu_elem = (NVisu+1)**dim
NVisu_p1_2 = (NVisu+1)**2

IF (HighOrder.EQ.1 .AND. DGFV.EQ.0) THEN
  IF (.NOT.ALLOCATED(NodeIndizes))   ALLOCATE(NodeIndizes  ( 1:NVisu+1  ,1:NVisu_j+1 , 1:NVisu_k+1))

  nVTKCells  = nElems
  NodeID     = 0
  DO k=1,NVisu_k+1; DO j=1,NVisu_j+1; DO i=1,NVisu+1
    NodeIndizes(i,j,k) = NodeID
    NodeID             = NodeID+1
  END DO; END DO; END DO

  SDEALLOCATE(nodeids)
  ALLOCATE(nodeids(((NVisu+1)**dim)*nVTKCells))

  ! create connectivity
  NodeID     = 0
  NodeIDElem = 0
  nodeids    = 0

  ! Create cell for each element
  DO iElem = 1,nElems
    ! Vertices
    nodeids(NodeID+1) = NodeIDElem + NodeIndizes(1      ,1,1)
    nodeids(NodeID+2) = NodeIDElem + NodeIndizes(NVisu+1,1,1)
    IF (dim.GE.2) THEN
      nodeids(NodeID+3) = NodeIDElem +NodeIndizes(NVisu+1,NVisu_j+1,1)
      nodeids(NodeID+4) = NodeIDElem +NodeIndizes(1      ,NVisu_j+1,1)
    END IF
    IF (dim.EQ.3) THEN
      nodeids(NodeID+5) = NodeIDElem +NodeIndizes(1      ,1        ,NVisu_k+1)
      nodeids(NodeID+6) = NodeIDElem +NodeIndizes(NVisu+1,1        ,NVisu_k+1)
      nodeids(NodeID+7) = NodeIDElem +NodeIndizes(NVisu+1,NVisu_j+1,NVisu_k+1)
      nodeids(NodeID+8) = NodeIDElem +NodeIndizes(1      ,NVisu_j+1,NVisu_k+1)
    END IF
    NodeID = NodeID + 2**dim

    ! Edges
    ! Note from trixi: This order doesn't make any sense. This is completely different
    ! from what is shown in
    ! https://blog.kitware.com/wp-content/uploads/2018/09/Source_Issue_43.pdf
    ! but this is the way it works.
    IF (nVisu.GT.1) THEN
      nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(2:NVisu,1,1)
      NodeID = NodeID + NVisu-1
      IF (dim.GE.2) THEN
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(NVisu+1,2:NVisu_j,1)
        NodeID = NodeID + NVisu_j-1
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(2:NVisu,NVisu_j+1,1)
        NodeID = NodeID + NVisu  -1
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(1      ,2:NVisu_j,1)
        NodeID = NodeID + NVisu_j-1
      END IF
      IF (dim.EQ.3) THEN
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(2:NVisu,1        ,NVisu_k+1)
        NodeID = NodeID + NVisu  -1
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(NVisu+1,2:NVisu_j,NVisu_k+1)
        NodeID = NodeID + NVisu_j-1
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(2:NVisu,NVisu_j+1,NVisu_k+1)
        NodeID = NodeID + NVisu  -1
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(1      ,2:NVisu_j,NVisu_k+1)
        NodeID = NodeID + NVisu_j-1
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(1      ,1        ,2:NVisu_k)
        NodeID = NodeID + NVisu_k-1
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(NVisu+1,1        ,2:NVisu_k)
        NodeID = NodeID + NVisu_k-1
        ! The following two are switched compared to trixi because ParaView changed the ordering from VTK8 to VTK9 convention
        ! https://gitlab.kitware.com/paraview/paraview/-/issues/20728
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(NVisu+1,NVisu_j+1,2:NVisu_k)
        NodeID = NodeID + NVisu_k-1
        nodeids(NodeID+1:NodeID+NVisu-1) = NodeIDElem + NodeIndizes(1      ,NVisu_j+1,2:NVisu_k)
        NodeID = NodeID + NVisu_k-1
      END IF

      ! Faces
      ! Note from trixi: See above
      IF (dim.EQ.3) THEN
        nodeids(NodeID+1:NodeID+((NVisu_j-1)*(NVisu_k-1))) = NodeIDElem &
                                                           + RESHAPE(NodeIndizes(1      ,2:NVisu_j,2:NVisu_k),(/(NVisu_j-1)*(NVisu_k-1)/))
        NodeID = NodeID + (NVisu_j-1)*(NVisu_k-1)
        nodeids(NodeID+1:NodeID+((NVisu_j-1)*(NVisu_k-1))) = NodeIDElem &
                                                           + RESHAPE(NodeIndizes(NVisu+1,2:NVisu_j,2:NVisu_k),(/(NVisu_j-1)*(NVisu_k-1)/))
        NodeID = NodeID + (NVisu_j-1)*(NVisu_k-1)
        nodeids(NodeID+1:NodeID+((NVisu  -1)*(NVisu_k-1))) = NodeIDElem &
                                                           + RESHAPE(NodeIndizes(2:NVisu,1        ,2:NVisu_k),(/(NVisu  -1)*(NVisu_k-1)/))
        NodeID = NodeID + (NVisu  -1)*(NVisu_k-1)
        nodeids(NodeID+1:NodeID+((NVisu_j-1)*(NVisu_k-1))) = NodeIDElem &
                                                           + RESHAPE(NodeIndizes(2:NVisu,NVisu_j+1,2:NVisu_k),(/(NVisu_j-1)*(NVisu_k-1)/))
        NodeID = NodeID + (NVisu  -1)*(NVisu_k-1)
      END IF
      IF (dim.GE.2) THEN
        nodeids(NodeID+1:NodeID+((NVisu  -1)*(NVisu_j-1))) = NodeIDElem &
                                                           + RESHAPE(NodeIndizes(2:NVisu,2:NVisu_j,1        ),(/(NVisu  -1)*(NVisu_j-1)/))
        NodeID = NodeID + (NVisu  -1)*(NVisu_j-1)
      END IF
      IF (dim.EQ.3) THEN
        nodeids(NodeID+1:NodeID+((NVisu  -1)*(NVisu_j-1))) = NodeIDElem &
                                                           + RESHAPE(NodeIndizes(2:NVisu,2:NVisu_j,NVisu_k+1),(/(NVisu  -1)*(NVisu_j-1)/))
        NodeID = NodeID + (NVisu  -1)*(NVisu_j-1)
      END IF

      ! Volume
      IF (dim.EQ.3) THEN
        nodeids(NodeID+1:NodeID+((NVisu  -1)*(NVisu_j-1)*(NVisu_k-1))) = NodeIDElem &
                                                           + RESHAPE(NodeIndizes(2:NVisu,2:NVisu_j,2:NVisu_k),(/(NVisu-1)*(NVisu_j-1)*(NVisu_k-1)/))
        NodeID = NodeID + (NVisu  -1)*(NVisu_j-1)*(NVisu_k-1)
      END IF
    END IF

    NodeIDElem = NodeIDElem + NVisu_elem
  END DO

  DEALLOCATE(NodeIndizes)

ELSE
  nVTKCells  = ((NVisu+DGFV)/(1+DGFV))**dim*nElems
  SDEALLOCATE(nodeids)
  ALLOCATE(nodeids((2**dim)*nVTKCells))

  ! create connectivity
  NodeID     = 0
  NodeIDElem = 0
  DO iElem=1,nElems
    DO k=1,NVisu_k,(DGFV+1)
      DO j=1,NVisu_j,(DGFV+1)
        DO i=1,NVisu,(DGFV+1)
          IF (dim.GE.2) THEN
            NodeID=NodeID+1
            nodeids(NodeID) = NodeIDElem+i+   j   *(NVisu+1)+(k-1)*NVisu_p1_2-1 !P4(CGNS=tecVisu standard)
          END IF
          NodeID=NodeID+1
          nodeids(NodeID) = NodeIDElem+i+  (j-1)*(NVisu+1)+(k-1)*NVisu_p1_2-1 !P1
          NodeID=NodeID+1
          nodeids(NodeID) = NodeIDElem+i+1+(j-1)*(NVisu+1)+(k-1)*NVisu_p1_2-1 !P2
          IF (dim.GE.2) THEN
            NodeID=NodeID+1
            nodeids(NodeID) = NodeIDElem+i+1+ j   *(NVisu+1)+(k-1)*NVisu_p1_2-1 !P3
          END IF
          IF (dim.EQ.3) THEN
            NodeID=NodeID+1
            nodeids(NodeID)=NodeIDElem+i+   j   *(NVisu+1)+ k   *NVisu_p1_2-1 !P8
            NodeID=NodeID+1
            nodeids(NodeID)=NodeIDElem+i+  (j-1)*(NVisu+1)+ k   *NVisu_p1_2-1 !P5
            NodeID=NodeID+1
            nodeids(NodeID)=NodeIDElem+i+1+(j-1)*(NVisu+1)+ k   *NVisu_p1_2-1 !P6
            NodeID=NodeID+1
            nodeids(NodeID)=NodeIDElem+i+1+ j   *(NVisu+1)+ k   *NVisu_p1_2-1 !P7
          END IF
        END DO
      END DO
    END DO
    NodeIDElem=NodeIDElem+NVisu_elem
  END DO
END IF

END SUBROUTINE CreateConnectivity

!===================================================================================================================================
!> Subroutine to write 2D or 3D point data to VTK format
!===================================================================================================================================
SUBROUTINE WriteDataToVTK(nVal,NVisu,nElems,VarNames,Coord,Value,FileString,dim,DGFV,nValAtLastDimension,HighOrder,PostiParallel,&
                          OutputDirectory)
! MODULES
USE MOD_Globals
USE MOD_Restart_Vars   ,ONLY: RestartTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nVal                 !< Number of nodal output variables
INTEGER,INTENT(IN)          :: NVisu                !< Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)          :: nElems               !< Number of output elements
INTEGER,INTENT(IN)          :: dim                  !< dimension: 2 or 3
REAL,INTENT(IN)             :: Coord(1:3,0:NVisu,0:NVisu*(MERGE(1,0,dim.GT.1)),0:NVisu*(MERGE(1,0,dim.GT.2)),nElems)     !< CoordsVector
CHARACTER(LEN=*),INTENT(IN) :: VarNames(nVal)       !< Names of all variables that will be written out
REAL,INTENT(IN)             :: Value(:,:,:,:,:)     !< Statevector
CHARACTER(LEN=*),INTENT(IN) :: FileString           !< Output file name
INTEGER,OPTIONAL,INTENT(IN) :: DGFV                 !< flag indicating DG = 0 or FV =1 data
LOGICAL,OPTIONAL,INTENT(IN) :: nValAtLastDimension  !< if TRUE, nVal is stored in the last index of value
LOGICAL,OPTIONAL,INTENT(IN) :: HighOrder            !< if TRUE, posti uses high-order element representation
LOGICAL,OPTIONAL,INTENT(IN) :: PostiParallel        !< if TRUE, posti runs parallel
CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: OutputDirectory      !< Custom output directory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iVal,ivtk
INTEGER                     :: nElems_glob(0:nProcessors-1)
INTEGER                     :: NVisu_elem,nVTKPoints,nVTKCells
INTEGER                     :: nTotalElems
INTEGER                     :: nBytes,Offset
INTEGER                     :: INTdummy
REAL(KIND=4)                :: FLOATdummy
CHARACTER(LEN=35)           :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=1)            :: lf
INTEGER                     :: ElemType,iElem
INTEGER,ALLOCATABLE,TARGET  :: nodeids(:)
INTEGER                     :: NVisu_k,NVisu_j,PointsPerVTKCell
#if USE_MPI
INTEGER                     :: iProc,nElems_proc,nElemsMax
REAL,ALLOCATABLE            :: buf(:,:,:,:), buf2(:,:,:,:,:)
#endif /*USE_MPI*/
INTEGER                     :: DGFV_loc
LOGICAL                     :: nValAtLastDimension_loc
INTEGER                     :: HighOrder_loc                                          ! INTEGER to be consistent with visu_Cwrapper
LOGICAL                     :: PostiParallel_loc
CHARACTER(LEN=255)          :: FileString_loc
!===================================================================================================================================
IF (PRESENT(DGFV)) THEN
  DGFV_loc = DGFV
ELSE
  DGFV_loc = 0
END IF

IF (PRESENT(nValAtLastDimension)) THEN
  nValAtLastDimension_loc = nValAtLastDimension
ELSE
  nValAtLastDimension_loc = .FALSE.
END IF

IF (PRESENT(HighOrder)) THEN
  HighOrder_loc = MERGE(1,0,HighOrder.AND. DGFV_loc.EQ.0)
ELSE
  HighOrder_loc = 0
END IF

IF (PRESENT(PostiParallel)) THEN
  IF(nProcessors.GT.1)THEN
    PostiParallel_loc=.TRUE.
    FileString_loc='visu/'//TRIM(INTSTAMP(TRIM(FileString),myRank))//'.vtu'
  ELSE
    PostiParallel_loc = .FALSE.
    FileString_loc=TRIM(FileString)//'.vtu'
  END IF
ELSE
  PostiParallel_loc = .FALSE.
  FileString_loc=TRIM(FileString)//'.vtu'
END IF

! Prepend output directory
IF (PRESENT(OutputDirectory)) THEN
  IF (TRIM(OutputDirectory).NE.'') THEN
    FileString_loc=TRIM(OutputDirectory)//'/'//TRIM(FileString_loc)
    ! create visu dir, where all vtu files are placed
    IF (PostiParallel_loc) THEN
      IF (MPIRoot) CALL SYSTEM('mkdir -p '//TRIM(OutputDirectory)//'/visu')
    END IF
  ELSE
    ! create visu dir, where all vtu files are placed
    IF (PostiParallel_loc) THEN
      IF (MPIRoot) CALL SYSTEM('mkdir -p visu')
    END IF
  END IF ! TRIM(OutputDirectory).NE.''
ELSE
  ! create visu dir, where all vtu files are placed
  IF (PostiParallel_loc) THEN
    IF (MPIRoot) CALL SYSTEM('mkdir -p visu')
  END IF
END IF ! PRESENT(OutputDirectory)

#if USE_MPI
IF (PostiParallel_loc) CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

IF (dim.EQ.3) THEN
  NVisu_k = NVisu
  NVisu_j = NVisu
ELSE IF (dim.EQ.2) THEN
  NVisu_k = 0
  NVisu_j = NVisu
ELSE IF (dim.EQ.1) THEN
  NVisu_k = 0
  NVisu_j = 0
ELSE
  CALL Abort(__STAMP__, &
      "Only 2D and 3D connectivity can be created. dim must be 1, 2 or 3.")
END IF
PointsPerVTKCell = MERGE((NVisu+1)**dim,2**dim,HighOrder_loc.EQ.1)

SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO') " WRITE ",dim,"D DATA TO VTX XML BINARY (VTU) FILE..."

! get total number of elements on all processors
#if USE_MPI
CALL MPI_GATHER(nElems,1,MPI_INTEGER,nElems_glob,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
#else
nElems_glob(0) = nElems
#endif /*USE_MPI*/
IF(.NOT.PostiParallel_loc)THEN
  nTotalElems = SUM(nElems_glob)
ELSE
  nTotalElems = nElems
END IF

NVisu_elem = (NVisu+1)**dim
nVTKPoints = NVisu_elem * nTotalElems
IF (HighOrder_loc.EQ.1) THEN
  nVTKCells  = nTotalElems
ELSE
  nVTKCells  = ((NVisu+DGFV_loc)/(1+DGFV_loc))**dim*nTotalElems
END IF

! write header of VTK file
IF((.NOT.PostiParallel_loc.AND.MPIRoot).OR.PostiParallel_loc)THEN
  ! Line feed character
  lf = char(10)

  ! Write file
  OPEN(NEWUNIT=ivtk,FILE=TRIM(FileString_loc),ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! This version is important for high-order elements because ParaView automatically converts the node numbering when switching from VTK8 to VTK9 convention
  ! https://gitlab.kitware.com/paraview/paraview/-/issues/20728
  Buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify output time
  Buffer='    <FieldData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;WRITE(ivtk) TRIM(Buffer)
  WRITE(TempStr1,'(F17.9)')RestartTime
  Buffer='        '//TRIM(ADJUSTL(TempStr1))//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </DataArray>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </FieldData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify field pieces
  WRITE(TempStr1,'(I16)')nVTKPoints
  WRITE(TempStr2,'(I16)')nVTKCells
  Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
         &NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  DO iVal=1,nVal
    Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVal))//'" '// &
                     'format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    Offset=Offset+SIZEOF_F(INTdummy)+nVTKPoints*SIZEOF_F(FLOATdummy)
    WRITE(StrOffset,'(I16)')Offset
  END DO
  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+3*nVTKPoints*SIZEOF_F(FLOATdummy)
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify necessary cell data
  Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Connectivity
  Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+PointsPerVTKCell*nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Offsets
  Buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Elem types
  Buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
                   'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)
END IF


#if USE_MPI
IF(.NOT.PostiParallel_loc)THEN
  IF(MPIRoot)THEN
    !ALLOCATE buffer for Root
    nElemsMax=MAXVAL(nElems_glob)
    ALLOCATE(buf(   0:NVisu,0:NVisu_j,0:NVisu_k,nElemsMax))
  END IF
END IF
#endif /*USE_MPI*/

! Write binary raw data into append section
! Solution data
IF(.NOT.PostiParallel_loc)THEN
  DO iVal=1,nVal
    IF(MPIRoot)THEN
      nBytes = nVTKPoints*SIZEOF_F(FLOATdummy)
      IF (nValAtLastDimension_loc) THEN
        WRITE(ivtk) nBytes,REAL(Value(:,:,:,:,iVal),4)
      ELSE
        WRITE(ivtk) nBytes,REAL(Value(iVal,:,:,:,:),4)
      END IF
#if USE_MPI
      DO iProc=1,nProcessors-1
        nElems_proc=nElems_glob(iProc)
        IF (nElems_proc.GT.0) THEN
          CALL MPI_RECV(buf(:,:,:,1:nElems_proc),nElems_proc*NVisu_elem,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_FLEXI,MPIstatus,iError)
          WRITE(ivtk) REAL(buf(:,:,:,1:nElems_proc),4)
        END IF
      END DO !iProc
    ELSE
      IF (nElems.GT.0) THEN
        IF (nValAtLastDimension_loc) THEN
          CALL MPI_SEND(Value(:,:,:,:,iVal),nElems*NVisu_elem,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_FLEXI,iError)
        ELSE
          CALL MPI_SEND(Value(iVal,:,:,:,:),nElems*NVisu_elem,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_FLEXI,iError)
        END IF
      END IF
#endif /*USE_MPI*/
    END IF !MPIRoot
  END DO       ! iVar
ELSE
  DO iVal=1,nVal
    nBytes = nVTKPoints*SIZEOF_F(FLOATdummy)
    IF (nValAtLastDimension_loc) THEN
      WRITE(ivtk) nBytes,REAL(Value(:,:,:,:,iVal),4)
    ELSE
      WRITE(ivtk) nBytes,REAL(Value(iVal,:,:,:,:),4)
    END IF
  END DO
END IF

#if USE_MPI
IF(.NOT.PostiParallel_loc)THEN
  IF(MPIRoot)THEN
    SDEALLOCATE(buf)
    ALLOCATE(buf2(3,0:NVisu,0:NVisu_j,0:NVisu_k,nElemsMax))
  END IF
END IF
#endif /*USE_MPI*/

! Coordinates
IF(.NOT.PostiParallel_loc)THEN
  IF(MPIRoot)THEN
    nBytes = nVTKPoints*SIZEOF_F(FLOATdummy) * 3
    WRITE(ivtk) nBytes
    WRITE(ivtk) REAL(Coord(:,:,:,:,:),4)
#if USE_MPI
    DO iProc=1,nProcessors-1
      nElems_proc=nElems_glob(iProc)
      IF (nElems_proc.GT.0) THEN
        CALL MPI_RECV(buf2(:,:,:,:,1:nElems_proc),nElems_proc*NVisu_elem*3,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_FLEXI,MPIstatus,iError)
        WRITE(ivtk) REAL(buf2(:,:,:,:,1:nElems_proc),4)
      END IF
    END DO !iProc
  ELSE
    IF (nElems.GT.0) THEN
      CALL MPI_SEND(Coord(:,:,:,:,:),nElems*NVisu_elem*3,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_FLEXI,iError)
    END IF
#endif /*USE_MPI*/
  END IF !MPIRoot
ELSE
  nBytes = nVTKPoints*SIZEOF_F(FLOATdummy) * 3
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL(Coord(:,:,:,:,:),4)
END IF

#if USE_MPI
IF(MPIRoot.AND..NOT.PostiParallel_loc)THEN
  SDEALLOCATE(buf2)
END IF
#endif /*USE_MPI*/

! Connectivity and footer
IF((.NOT.PostiParallel_loc.AND.MPIRoot).OR.PostiParallel_loc)THEN
  CALL CreateConnectivity(NVisu=NVisu,nElems=nTotalElems,nodeids=nodeids,dim=dim,DGFV=DGFV_loc,HighOrder=HighOrder_loc)

  nBytes = PointsPerVTKCell*nVTKCells*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) nodeids
  ! Offset
  nBytes = nVTKCells*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) (Offset,Offset=PointsPerVTKCell,PointsPerVTKCell*nVTKCells,PointsPerVTKCell)
  ! Elem type
  IF (HighOrder_loc.EQ.1) THEN
    IF (dim.EQ.3) THEN
      ElemType = 72 ! VTK_LAGRANGE_HEXAHEDRON
    ELSE IF (dim.EQ.2) THEN
      ElemType = 70 ! VTK_LAGRANGE_QUADRILATERAL
    ELSE IF (dim.EQ.1) THEN
      ElemType = 68 ! VTK_LAGRANGE_CURVE
    END IF
  ELSE
    IF (dim.EQ.3) THEN
      ElemType = 12 ! VTK_HEXAHEDRON
    ELSE IF (dim.EQ.2) THEN
      ElemType = 9  ! VTK_QUAD
    ELSE IF (dim.EQ.1) THEN
      ElemType = 3  ! VTK_LINE
    END IF
  END IF
  WRITE(ivtk) nBytes
  WRITE(ivtk) (ElemType,iElem=1,nVTKCells)

  DEALLOCATE(nodeids)

  ! Footer
  lf = char(10)
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
ENDIF

#if USE_MPI
IF(PostiParallel_loc.AND.MPIRoot)THEN
  IF(PRESENT(OutputDirectory))THEN
    CALL WriteParallelVTK(FileString,nVal,VarNames,OutputDirectory)
  ELSE
    CALL WriteParallelVTK(FileString,nVal,VarNames)
  END IF
ENDIF
#endif


SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK

!===================================================================================================================================
!> Links DG and FV VTK files together
!===================================================================================================================================
SUBROUTINE WriteVTKMultiBlockDataSet(FileString,FileString_DG,FileString_FV,OutputDirectory)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileString     !< Output file name
CHARACTER(LEN=*),INTENT(IN) :: FileString_DG  !< Filename of DG VTU file
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FileString_FV    !< Filename of FV VTU file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: OutputDirectory  !< Custom output directory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ivtk
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
!===================================================================================================================================
IF (MPIRoot) THEN
  IF(nProcessors.GT.1)THEN
    ! write '.vtu">'//multiblock file
    IF (PRESENT(OutputDirectory)) THEN
      IF (TRIM(OutputDirectory).NE.'') THEN
        OPEN(NEWUNIT=ivtk,FILE=TRIM(OutputDirectory)//'/'//TRIM(FileString)//'.pvd',STATUS='REPLACE',ACCESS='STREAM')
      ELSE
        OPEN(NEWUNIT=ivtk,FILE=                            TRIM(FileString)//'.pvd',STATUS='REPLACE',ACCESS='STREAM')
      END IF
    ELSE
      OPEN(NEWUNIT=ivtk,FILE=                            TRIM(FileString)//'.pvd',STATUS='REPLACE',ACCESS='STREAM')
    END IF
    ! Line feed character
    lf = char(10)
    Buffer='<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf
    WRITE(ivtk) TRIM(BUFFER)
    Buffer='  <Collection>'//lf;WRITE(ivtk) TRIM(BUFFER)
    Buffer='    <DataSet part="0" name="DG" file="'//TRIM(FileString_DG)//'.pvtu"/>'//lf;WRITE(ivtk) TRIM(BUFFER)
    Buffer='    <DataSet part="1" name="FV" file="'//TRIM(FileString_FV)//'.pvtu"/>'//lf;WRITE(ivtk) TRIM(BUFFER)
    Buffer='  </Collection>'//lf;WRITE(ivtk) TRIM(BUFFER)
    Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(BUFFER)
    CLOSE(ivtk)
  ELSE
    ! write '.vtu">'//multiblock file
    OPEN(NEWUNIT=ivtk,FILE=TRIM(FileString)//'.vtm',STATUS='REPLACE',ACCESS='STREAM')
    ! Line feed character
    lf = char(10)
    Buffer='<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf
    WRITE(ivtk) TRIM(BUFFER)
    Buffer='  <vtkMultiBlockDataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
    Buffer='    <DataSet index="0" name="DG" file="'//TRIM(FileString_DG)//'.vtu"/>'//lf;WRITE(ivtk) TRIM(BUFFER)
    Buffer='    <DataSet index="1" name="FV" file="'//TRIM(FileString_FV)//'.vtu"/>'//lf;WRITE(ivtk) TRIM(BUFFER)
    Buffer='  </vtkMultiBlockDataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
    Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(BUFFER)
    CLOSE(ivtk)
  ENDIF
ENDIF
END SUBROUTINE WriteVTKMultiBlockDataSet

#if USE_MPI
!===================================================================================================================================
!> Writes PVTU files
!===================================================================================================================================
SUBROUTINE WriteParallelVTK(FileString,nVal,VarNames,OutputDirectory)
! MODULES
USE MOD_Globals
USE MOD_Restart_Vars   ,ONLY: RestartTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileString     !< Output file name
INTEGER,INTENT(IN)          :: nVal                 !< Number of nodal output variables
CHARACTER(LEN=*),INTENT(IN) :: VarNames(nVal)       !< Names of all variables that will be written out
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: OutputDirectory  !< Custom output directory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ivtk, iVal, iProc
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
CHARACTER(LEN=35)  :: StrProc,TempStr1
!===================================================================================================================================
IF (MPIRoot) THEN
  ! write multiblock file
  IF (PRESENT(OutputDirectory)) THEN
    IF (TRIM(OutputDirectory).NE.'') THEN
      OPEN(NEWUNIT=ivtk,FILE=TRIM(OutputDirectory)//'/'//TRIM(FileString)//'.pvtu',STATUS='REPLACE',ACCESS='STREAM')
    ELSE
      OPEN(NEWUNIT=ivtk,FILE=                            TRIM(FileString)//'.pvtu',STATUS='REPLACE',ACCESS='STREAM')
    END IF
  ELSE
    OPEN(NEWUNIT=ivtk,FILE=                            TRIM(FileString)//'.pvtu',STATUS='REPLACE',ACCESS='STREAM')
  END IF
  ! Line feed character
  lf = char(10)
  Buffer='<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf
  WRITE(ivtk) TRIM(BUFFER)
  Buffer='  <PUnstructuredGrid GhostLevel="0">'//lf;WRITE(ivtk) TRIM(BUFFER)
  ! Specify output time
  Buffer='    <FieldData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;WRITE(ivtk) TRIM(Buffer)
  WRITE(TempStr1,'(F17.9)')RestartTime
  Buffer='        '//TRIM(ADJUSTL(TempStr1))//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </DataArray>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </FieldData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify point data
  Buffer='    <PPointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  DO iVal=1,nVal
    Buffer='      <PDataArray type="Float32" Name="'//TRIM(VarNames(iVal))//'" format="appended"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  END DO
  Buffer='    </PPointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='    <PCellData> </PCellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='    <PPoints>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      <PDataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </PPoints>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify necessary cell data
  Buffer='    <PCells>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Connectivity
  Buffer='      <PDataArray type="Int32" Name="connectivity" format="appended"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Offsets
  Buffer='      <PDataArray type="Int32" Name="offsets" format="appended"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Elem types
  Buffer='      <PDataArray type="Int32" Name="types" format="appended"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </PCells>'//lf;WRITE(ivtk) TRIM(Buffer)
  DO iProc=0,nProcessors-1
    WRITE(StrProc,'(I16)')iProc
    Buffer='    <Piece Source="visu/'//TRIM(INTSTAMP(TRIM(FileString),iProc))//'.vtu"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  END DO
  Buffer='  </PUnstructuredGrid>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(BUFFER)
  CLOSE(ivtk)
ENDIF
END SUBROUTINE WriteParallelVTK
#endif /*USE_MPI*/

!===================================================================================================================================
!> Subroutine to write 2D or 3D coordinates to VTK format
!===================================================================================================================================
SUBROUTINE WriteCoordsToVTK_array(NVisu,nElems,coords_out,nodeids_out,coords,nodeids,dim,DGFV,HighOrder)
USE ISO_C_BINDING
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: NVisu                        !< Polynomial degree for visualization
INTEGER,INTENT(IN)                   :: nElems                       !< Number of elements
INTEGER,INTENT(IN)                   :: dim                          !< Spacial dimension (2D or 3D)
INTEGER,INTENT(IN)                   :: DGFV                         !< flag indicating DG = 0 or FV =1 data
REAL,ALLOCATABLE,TARGET,INTENT(IN)   :: coords(:,:,:,:,:)            !< Array containing coordinates
INTEGER,INTENT(IN)                   :: HighOrder                    !< flag indicating high order
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,ALLOCATABLE,TARGET,INTENT(INOUT) :: nodeids(:)
TYPE (CARRAY), INTENT(INOUT)         :: coords_out
TYPE (CARRAY), INTENT(INOUT)         :: nodeids_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
coords_out%dim  = dim

IF (nElems.EQ.0) THEN
  coords_out%len  = 0
  nodeids_out%len = 0
  RETURN
END IF

SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO') " WRITE ",dim,"D COORDS TO VTX XML BINARY (VTU) ARRAY..."
! values and coords are already in the correct structure of VTK/Paraview

! create connectivity
CALL CreateConnectivity(NVisu=NVisu,nElems=nElems,nodeids=nodeids,dim=dim,DGFV=DGFV,HighOrder=HighOrder)

IF (HighOrder.EQ.1 .AND.DGFV.EQ.0) THEN
  ! set the sizes of the arrays
  coords_out%len  = 3*(NVisu+1)**dim*nElems
  nodeids_out%len =   (NVisu+1)**dim*nElems
ELSE
  ! set the sizes of the arrays
  coords_out%len = 3*(NVisu+1)**dim*nElems
  nodeids_out%len = (2**dim)*((NVisu+DGFV)/(1+DGFV))**dim*nElems
END IF

! assign data to the arrays (no copy!!!)
coords_out%data  = C_LOC(coords(1,0,0,0,1))
nodeids_out%data = C_LOC(nodeids(1))

SWRITE(UNIT_stdOut,'(A)')" Done!"

END SUBROUTINE WriteCoordsToVTK_array

!===================================================================================================================================
!> Subroutine to write actual 2D or 3D point data to VTK format
!===================================================================================================================================
SUBROUTINE WriteDataToVTK_array(nVal,NVisu,nElems,Values_out,values,dim)
USE ISO_C_BINDING
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: nVal                         !> Number of nodal output variables
INTEGER,INTENT(IN)                :: NVisu                        !> Polynomial degree for visualization
INTEGER,INTENT(IN)                :: nElems                       !> Number of elements
INTEGER,INTENT(IN)                :: dim                          !> Spacial dimension (2D or 3D)
REAL(C_DOUBLE),ALLOCATABLE,TARGET,INTENT(IN) :: values(:,:,:,:,:) !> Array containing the points values
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE (CARRAY), INTENT(INOUT)      :: Values_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
values_out%dim  = dim
IF (nElems.EQ.0) THEN
  values_out%len  = 0
  RETURN
END IF
SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO')" WRITE ",dim,"D DATA   TO VTX XML BINARY (VTU) ARRAY..."

! values and coords are already in the correct structure of VTK/Paraview
! set the sizes of the arrays
values_out%len = nVal*(NVisu+1)**dim*nElems

IF (nVal*(NVisu+1)**dim*nElems.GT.0) THEN
  ! assign data to the arrays (no copy!!!)
  values_out%data = C_LOC(values(0,0,0,1,1))
END IF

SWRITE(UNIT_stdOut,'(A)')" Done!"
END SUBROUTINE WriteDataToVTK_array

!===================================================================================================================================
!> Subroutine to write variable names to VTK format
!===================================================================================================================================
SUBROUTINE WriteVarnamesToVTK_array(nVarTotal,mapVisu,varnames_out,VarNamesTotal,nVarVisu)
USE ISO_C_BINDING
! MODULES
USE MOD_Globals
USE MOD_StringTools    ,ONLY: STRICMP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVarTotal
INTEGER,INTENT(IN)             :: mapVisu(nVarTotal)
TYPE (CARRAY), INTENT(INOUT)   :: varnames_out
CHARACTER(LEN=255),INTENT(IN)  :: VarNamesTotal(nVarTotal)
INTEGER,INTENT(IN)             :: nVarVisu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(C_CHAR),POINTER    :: VarNames_loc(:,:)
INTEGER                      :: i,iVar
!===================================================================================================================================
! copy varnames
ALLOCATE(VarNames_loc(255,nVarVisu))
varnames_out%len  = nVarVisu*255
IF (nVarVisu.GT.0) THEN
  varnames_out%data = C_LOC(VarNames_loc(1,1))

  DO iVar=1,nVarTotal
    IF (mapVisu(iVar).GT.0) THEN
      DO i=1,255
        VarNames_loc(i,mapVisu(iVar)) = VarNamesTotal(iVar)(i:i)
      END DO
    END IF
  END DO
END IF

END SUBROUTINE WriteVarnamesToVTK_array

END MODULE MOD_VTK
