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

!===================================================================================================================================
!> Module for generic data output in vtk xml fromat
!> WARNING: WriteDataToVTK works only for POSTPROCESSING
!===================================================================================================================================
MODULE MOD_VTK
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDataToVTK3D
  MODULE PROCEDURE WriteDataToVTK3D
END INTERFACE

INTERFACE WriteVTKMultiBlockDataSet
  MODULE PROCEDURE WriteVTKMultiBlockDataSet
END INTERFACE

PUBLIC::WriteDataToVTK3D
PUBLIC::WriteVTKMultiBlockDataSet
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write 3D point data to VTK format
!===================================================================================================================================
SUBROUTINE WriteDataToVTK3D(NPlot,nElems,nVal,VarNames,Coord,Value,FileString)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nVal                                        !< Number of nodal output variables
INTEGER,INTENT(IN)          :: NPlot                                       !< Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)          :: nElems                                      !< Number of output elements
REAL,INTENT(IN)             :: Coord(3,0:NPlot,0:NPlot,0:NPlot,nElems)     !< CoordsVector
CHARACTER(LEN=*),INTENT(IN) :: VarNames(nVal)                              !< Names of all variables that will be written out
REAL,INTENT(IN)             :: Value(nVal,0:NPlot,0:NPlot,0:NPlot,nElems)  !< Statevector
CHARACTER(LEN=*),INTENT(IN) :: FileString                                  !< Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iVal,iElem,Offset,nBytes,nVTKPoints,nVTKCells,ivtk=44
INTEGER            :: nGlobalElems_loc
INTEGER            :: INTdummy
INTEGER            :: NPlot_p1_3,NPlot_p1_2,PointID,CellID,ElemType
INTEGER,ALLOCATABLE:: Vertex(:,:)
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
REAL(KIND=4)       :: FLOATdummy
#if MPI
INTEGER            :: iProc,nElems_proc,nElemsMax
REAL,ALLOCATABLE   :: buf(:,:,:,:), buf2(:,:,:,:,:)
#endif /*MPI*/
INTEGER            :: nElems_glob(0:nProcessors-1)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE 3D DATA TO VTX XML BINARY (VTU) FILE..."

NPlot_p1_3=(NPlot+1)**3
NPlot_p1_2=(NPlot+1)**2

#if MPI
CALL MPI_GATHER(nElems,1,MPI_INTEGER,nElems_glob,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
#else
nElems_glob(0) = nElems
#endif

IF(MPIROOT)THEN
  ! here comes the MPI stuff
  nGlobalElems_loc=nElems
#if MPI
  !ALLOCATE buffer for Root
  nElemsMax=MAXVAL(nElems_glob)
  ALLOCATE(buf(   0:Nplot,0:Nplot,0:Nplot,nElemsMax))
  ALLOCATE(buf2(3,0:Nplot,0:Nplot,0:Nplot,nElemsMax))
  nGlobalElems_loc=SUM(nElems_glob)
#endif /*MPI*/

  ! Line feed character
  lf = char(10)

  ! Write file
  OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  nVTKPoints=NPlot_p1_3*nGlobalElems_loc
  nVTKCells =NPlot**3  *nGlobalElems_loc
  Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  WRITE(TempStr1,'(I16)')nVTKPoints
  WRITE(TempStr2,'(I16)')nVTKCells
  Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
         &NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  DO iVal=1,nVal
    Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVal))//'" &
           &format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    Offset=Offset+SIZEOF_F(INTdummy)+nVTKPoints*SIZEOF_F(FLOATdummy)
    WRITE(StrOffset,'(I16)')Offset
  END DO
  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+3*nVTKPoints*SIZEOF_F(FLOATdummy)
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify necessary cell data
  Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Connectivity
  Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+8*nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Offsets
  Buffer='        <DataArray type="Int32" Name="offsets" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Elem types
  Buffer='        <DataArray type="Int32" Name="types" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)

END IF

! Write binary raw data into append section
! Solution data
DO iVal=1,nVal
  IF(MPIroot)THEN
    nBytes = nVTKPoints*SIZEOF_F(FLOATdummy)
    WRITE(ivtk) nBytes,REAL(Value(iVal,:,:,:,:),4)
#if MPI
    DO iProc=1,nProcessors-1
      nElems_proc=nElems_glob(iProc)
      CALL MPI_RECV(buf(:,:,:,1:nElems_proc),nElems_proc*NPlot_p1_3,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
      WRITE(ivtk) REAL(buf(:,:,:,1:nElems_proc),4)
    END DO !iProc
  ELSE
    CALL MPI_SEND(Value(iVal,:,:,:,:),nElems*NPlot_p1_3,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_WORLD,iError)
#endif /*MPI*/
  END IF !MPIroot
END DO       ! iVar

! Coordinates
IF(MPIRoot)THEN
  nBytes = nVTKPoints*SIZEOF_F(FLOATdummy) * 3
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL(Coord(:,:,:,:,:),4)
#if MPI
  DO iProc=1,nProcessors-1
    nElems_proc=nElems_glob(iProc)
    CALL MPI_RECV(buf2(:,:,:,:,1:nElems_proc),nElems_proc*NPlot_p1_3*3,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
    WRITE(ivtk) REAL(buf2(:,:,:,:,1:nElems_proc),4)
  END DO !iProc
ELSE
  CALL MPI_SEND(Coord(:,:,:,:,:),nElems*NPlot_p1_3*3,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_WORLD,iError)
#endif /*MPI*/
END IF !MPIroot


! Connectivity
IF(MPIROOT)THEN
  PointID=0
  CellID=0
  ALLOCATE(Vertex(8,nVTKCells))
  DO iElem=1,nGlobalElems_loc
    DO k=1,NPlot; DO j=1,NPlot; DO i=1,NPlot
      CellID=CellID+1
      !
      Vertex(:,CellID)=(/                                       &
        PointID+i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P4(CGNS=tecplot standard)
        PointID+i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P1
        PointID+i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P2
        PointID+i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P3
        PointID+i+   j   *(NPlot+1)+ k   *NPlot_p1_2-1,      & !P8
        PointID+i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2-1,      & !P5
        PointID+i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2-1,      & !P6
        PointID+i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2-1      /) !P7
    END DO; END DO; END DO
    PointID=PointID+NPlot_p1_3
  END DO
  nBytes = 8*nVTKCells*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) Vertex(:,:)
  ! Offset
  nBytes = nVTKCells*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) (Offset,Offset=8,8*nVTKCells,8)
  ! Elem type
  ElemType = 12 ! VTK_HEXAHEDRON
  WRITE(ivtk) nBytes
  WRITE(ivtk) (ElemType,iElem=1,nVTKCells)
  ! Write footer
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
  SDEALLOCATE(Vertex)
#if MPI
  SDEALLOCATE(buf)
  SDEALLOCATE(buf2)
#endif /*MPI*/
ENDIF
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK3D



SUBROUTINE WriteVTKMultiBlockDataSet(FileString,FileString_DG,FileString_FV)
!===================================================================================================================================
! Linkes VTK data- und mesh-files together
!===================================================================================================================================
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileString     !< Output file name
CHARACTER(LEN=*),INTENT(IN) :: FileString_DG  !< Filename of DG VTU file 
CHARACTER(LEN=*),INTENT(IN) :: FileString_FV  !< Filename of FV VTU file 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ivtk=44
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
!===================================================================================================================================
IF (MPIRoot) THEN                   
  ! write multiblock file
  OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
  ! Line feed character
  lf = char(10)
  Buffer='<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf
  WRITE(ivtk) TRIM(BUFFER)
  Buffer='  <vtkMultiBlockDataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    <DataSet index="0" name="DG" file="'//TRIM(FileString_DG)//'">'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    </DataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    <DataSet index="1" name="FV" file="'//TRIM(FileString_FV)//'">'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='    </DataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='  </vtkMultiBlockDataSet>'//lf;WRITE(ivtk) TRIM(BUFFER)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(BUFFER)
  CLOSE(ivtk)
ENDIF
END SUBROUTINE WriteVTKMultiBlockDataSet

END MODULE MOD_VTK
