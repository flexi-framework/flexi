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

!===================================================================================================================================
!> Module containing the main procedures for the visu tool: visu_requestInformation is called by ParaView to create a
!> list of available variables and visu is the main routine which is either called by ParaView to get the data it visualizes
!> or by the standalone tool.
!===================================================================================================================================
MODULE MOD_Visu_HDF5_Output
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE visu_WriteHDF5
  MODULE PROCEDURE visu_WriteHDF5
END INTERFACE

PUBLIC:: visu_WriteHDF5
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Subroutine to write the solution U to HDF5 format
!> Is used for postprocessing
!==================================================================================================================================
SUBROUTINE visu_WriteHDF5(nVarVisu,NVisu,nElems_loc,FileString,MeshFileName,VarNames_loc  &
                         ,Coords_DG,Coords_DG2D            ,dim                           & ! ,Coords_DG1D
                         ,UVisu_DG ,UVisu_DG2D             )                                ! ,UVisu_DG1D
! MODULES
USE MOD_Globals               !,ONLY: ABORT,TIMESTAMP,MPIROOT,MPI_COMM_FLEXI,UNIT_stdOut
USE MOD_PreProc
USE MOD_2D                    ,ONLY: ExpandArrayTo3D
USE MOD_HDF5_Output           ,ONLY: GenerateFileSkeleton,WriteAttribute,GatheredWriteArray,MarkWriteSuccessful
USE MOD_IO_HDF5               ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Mesh_Vars             ,ONLY: nElems,nGlobalElems,offsetElem
!USE MOD_Output_Vars           ,ONLY: ProjectName
USE MOD_Visu_Vars             ,ONLY: OutputTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVarVisu
INTEGER,INTENT(IN)             :: NVisu
INTEGER,INTENT(IN)             :: nElems_loc
CHARACTER(LEN=*),INTENT(IN)    :: FileString
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
CHARACTER(LEN=*),INTENT(IN)    :: VarNames_loc(nVarVisu)
INTEGER,INTENT(IN)             :: dim
REAL,INTENT(IN),TARGET,OPTIONAL:: Coords_DG  (1:3,0:nVisu,0:nVisu,0:MERGE(nVisu,0,dim.GT.2),1:nElems)
REAL,INTENT(IN),TARGET,OPTIONAL:: Coords_DG2D(1:3,0:nVisu,0:ZDIM(NVisu),1:nElems_loc)
! REAL,INTENT(IN),TARGET,OPTIONAL:: Coords_DG1D(1:3)
REAL,INTENT(IN),TARGET,OPTIONAL:: UVisu_DG   (0:nVisu,0:nVisu,0:MERGE(nVisu,0,dim.GT.2),1:nElems,1:nVarVisu)
REAL,INTENT(IN),TARGET,OPTIONAL:: UVisu_DG2D (0:NVisu,0:ZDIM(NVisu),1:nElems_loc,1:nVarVisu)
! REAL,INTENT(IN),TARGET,OPTIONAL:: UVisu_DG1D (1:nVarVisu)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName!,FileType
INTEGER                        :: iElem,iVar!,k
REAL,POINTER                   :: UOut(:,:,:,:,:),UOut2D(:,:,:,:)
INTEGER,ALLOCATABLE            :: nVal(:),nValGlobal(:),offset(:)
#if USE_MPI
INTEGER                        :: nGlobalElems_loc,offsetElem_loc
INTEGER                        :: recvbuf,sendbuf
#endif /*USE_MPI*/
!==================================================================================================================================

SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO')" WRITE ",dim,"D DATA TO HDF5 FILE..."

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
!FileType = 'Solution'
!FileName = TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileType),OutputTime))//'.h5'
FileName = FileString

SELECT CASE(dim)
  CASE(3)
    IF(MPIRoot) CALL GenerateFileSkeleton( TRIM(FileName)                                   &
                                         , 'State'                                          &
                                         , nVarVisu                                         &
                                         , NVisu                                            &
                                         , VarNames_loc                                     &
                                         , MeshFileName                                     &
                                         , OutputTime                                       &
                                         , OutputTime                                       &
                                         , create        = .TRUE.                           &
                                         , Dataset       = 'Visu')
    ALLOCATE( nVal      (dim+2) &
            , nValGlobal(dim+2) &
            , offset    (dim+2))

  CASE(2)
    IF(MPIRoot) THEN
      CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteAttribute(File_ID,'VarNames_Visu2D',nVarVisu,StrArray=VarNames_loc)
      CALL CloseDataFile()
    END IF
    ALLOCATE( nVal      (dim+2) &
            , nValGlobal(dim+2) &
            , offset    (dim+2))

  CASE(1)
    IF(MPIRoot) THEN
      CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteAttribute(File_ID,'VarNames_Visu1D',nVarVisu,StrArray=VarNames_loc)
      CALL CloseDataFile()
    END IF
    ALLOCATE( nVal      (1) &
            , nValGlobal(1) &
            , offset    (1))

  CASE DEFAULT
    CALL ABORT(__STAMP__,'Invalid dimension for hdf5 output!')

END SELECT

! Set size of output
SELECT CASE(dim)
  CASE(3)
    nVal       = (/nVarVisu,NVisu+1,NVisu+1,NVisu+1,nElems      /)
    nValGlobal = (/nVarVisu,NVisu+1,NVisu+1,NVisu+1,nGlobalElems/)
    offset     = (/0       ,0      ,0      ,0      ,offsetElem  /)
    ! UVisu_DG is sorted in the old style, re-sort for HDF5
    ALLOCATE(UOut(1:nVarVisu,0:NVisu,0:NVisu,0:NVisu,1:nElems))
    DO iElem = 1,nElems
  !    DO k = 0,NVisu; DO j = 0,NVisu; DO i = 0,NVisu;
        DO iVar = 1,nVarVisu
          UOut(iVar,:,:,:,iElem) = UVisu_DG(:,:,:,iElem,iVar)
        END DO
  !    END DO; END DO; END DO
    END DO

  CASE (2)
#if USE_MPI
    ASSOCIATE( nElems       => nElems_loc       &
             , nGlobalElems => nGlobalElems_loc &
             , offsetElem   => offsetElem_loc)

    ! Calculate local offset
    sendbuf = nElems
    recvbuf = 0
    CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_FLEXI,iError)
    offsetElem = recvbuf
    sendbuf    = recvbuf + nElems
    CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nProcessors-1,MPI_COMM_FLEXI,iError)
    nGlobalElems = sendbuf
#else
    ASSOCIATE( nElems       => nElems_loc       &
             , nGlobalElems => nElems_loc       &
             , offsetElem   => 0             )
#endif /*USE_MPI*/

    ! USurfVisu_DG is sorted in the old style, re-sort for HDF5
    nVal       = (/nVarVisu,NVisu+1,NVisu+1,nElems      /)
    nValGlobal = (/nVarVisu,NVisu+1,NVisu+1,nGlobalElems/)
    offset     = (/0       ,0      ,0      ,offsetElem  /)
  !  ! The output should be done with a full third dimension in a two dimensional computation, we need to expand the solution
  !  ALLOCATE(UOut(nVarVisu,0:NVisu,0:NVisu,0:NVisu,nElems))
  !  CALL ExpandArrayTo3D(5,(/nVarVisu,NVisu+1,NVisu+1,1,nElems/),4,NVisu+1,UVisu_DG(:,:,:,1,:),UOut)
    ALLOCATE(UOut2D(1:nVarVisu,0:NVisu,0:NVisu,1:nElems))
    DO iElem = 1,nElems
      DO iVar = 1,nVarVisu
        UOut2D(iVar,:,:,iElem) = UVisu_DG2D(:,:,iElem,iVar)
      END DO
    END DO

    END ASSOCIATE

  ! CASE (1)
  !   nVal       = (/nVarVisu/)
  !   nValGlobal = (/nVarVisu/)
  !   offset     = (/0       /)
  CASE(1)
    CALL Abort(__STAMP__,'1D HDF5 output currently not implemented')

END SELECT

! Reopen file and write DG solution
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

SELECT CASE(dim)
  CASE(3)
    ASSOCIATE(nValGlobal => (/3       ,NVisu+1,NVisu+1,NVisu+1,nGlobalElems/), &
              nVal       => (/3       ,NVisu+1,NVisu+1,NVisu+1,nElems      /))

    CALL GatheredWriteArray( TRIM(FileName)            &
                           , create      = .FALSE.     &
                           , DataSetName = 'Coords'    &
                           , rank        = dim+2       &
                           , nValGlobal  = nValGlobal  &
                           , nVal        = nVal        &
                           , offset      = offset      &
                           , collective  =.TRUE.       &
                           , RealArray   = Coords_DG)
    END ASSOCIATE

    CALL GatheredWriteArray( TRIM(FileName)            &
                           , create      = .FALSE.     &
                           , DataSetName = 'Visu'      &
                           , rank        = dim+2       &
                           , nValGlobal  = nValGlobal  &
                           , nVal        = nVal        &
                           , offset      = offset      &
                           , collective  =.TRUE.       &
                           , RealArray   = UOut)
    DEALLOCATE(UOut)

  CASE(2)
#if USE_MPI
    ASSOCIATE(nValGlobal => (/3       ,NVisu+1,NVisu+1,nGlobalElems_loc/), &
              nVal       => (/3       ,NVisu+1,NVisu+1,nElems_loc      /))
#else
    ASSOCIATE(nValGlobal => (/3       ,NVisu+1,NVisu+1,nElems_loc      /), &
              nVal       => (/3       ,NVisu+1,NVisu+1,nElems_loc      /))
#endif /*USE_MPI*/

    CALL GatheredWriteArray( TRIM(FileName)            &
                           , create      = .FALSE.     &
                           , DataSetName = 'Coords2D'  &
                           , rank        = dim+2       &
                           , nValGlobal  = nValGlobal  &
                           , nVal        = nVal        &
                           , offset      = offset      &
                           , collective  =.TRUE.       &
                           , RealArray   = Coords_DG2D)
    END ASSOCIATE

    CALL GatheredWriteArray( TRIM(FileName)            &
                           , create      = .FALSE.     &
                           , DataSetName = 'Visu2D'    &
                           , rank        = dim+2       &
                           , nValGlobal  = nValGlobal  &
                           , nVal        = nVal        &
                           , offset      = offset      &
                           , collective  =.TRUE.       &
                           , RealArray   = UOut2D)
    DEALLOCATE(UOut2D)

  ! CASE(1)
  !   ASSOCIATE(nValGlobal => (/3/), &
  !             nVal       => (/3/))

  !   CALL GatheredWriteArray( TRIM(FileName)            &
  !                          , create      = .FALSE.     &
  !                          , DataSetName = 'Coords1D'  &
  !                          , rank        = 1           &
  !                          , nValGlobal  = nValGlobal  &
  !                          , nVal        = nVal        &
  !                          , offset      = offset      &
  !                          , collective  =.TRUE.       &
  !                          , RealArray   = Coords_DG1D)
  !   END ASSOCIATE

  !   CALL GatheredWriteArray( TRIM(FileName)            &
  !                          , create      = .FALSE.     &
  !                          , DataSetName = 'Visu1D'    &
  !                          , rank        = 1           &
  !                          , nValGlobal  = nValGlobal  &
  !                          , nVal        = nVal        &
  !                          , offset      = offset      &
  !                          , collective  =.TRUE.       &
  !                          , RealArray   = UVisu_DG1D)


END SELECT

DEALLOCATE(nVal,nValGlobal,offset)

IF(MPIRoot)THEN
  CALL MarkWriteSuccessful(FileName)
  WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END IF

END SUBROUTINE visu_WriteHDF5

END MODULE MOD_Visu_HDF5_Output
