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
!> Module to handle the input of the recordpoints
!===================================================================================================================================
MODULE MOD_RPData
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ReadRPData
  MODULE PROCEDURE ReadRPData
END INTERFACE

INTERFACE AssembleRPData
  MODULE PROCEDURE AssembleRPData
END INTERFACE

INTERFACE FinalizeRPData
  MODULE PROCEDURE FinalizeRPData
END INTERFACE

PUBLIC :: ReadRPData
PUBLIC :: AssembleRPData
PUBLIC :: FinalizeRPData
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Read in the RP data from a single .h5 file
!===================================================================================================================================
SUBROUTINE ReadRPData(FileString,firstFile)
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input
USE MOD_RPData_Vars
USE MOD_RPSetVisu          ,ONLY: InitRPSet,ChangeRPSet
USE MOD_RPSetVisuVisu_Vars ,ONLY: nRP_HDF5,RPOutMap
USE MOD_ParametersVisu     ,ONLY: ProjectName
USE MOD_ParametersVisu     ,ONLY: skip,RP_DefFile,RP_SET_defined
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),INTENT(IN) :: FileString
LOGICAL,INTENT(IN),OPTIONAL   :: firstFile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: firstFile_loc
CHARACTER(LEN=255)            :: FileType
CHARACTER(LEN=255)            :: RP_DefFile_loc
INTEGER                       :: nSamples_loc,nSamples_skip
INTEGER                       :: iSample
REAL,ALLOCATABLE              :: temparray(:,:,:)
#if USE_MPI
INTEGER                       :: nChunks,iChunk,nSamples_chunk,offset
REAL                          :: nTotal,limit
#endif
!===================================================================================================================================

firstFile_loc=.FALSE.
IF (PRESENT(firstFile)) THEN
  IF (firstFile.EQV..TRUE.) firstFile_loc=.TRUE.
END IF

IF (firstFile_loc) THEN
  WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' Read recordpoint data from data file "'//TRIM(FileString)//'" ...'
ELSE
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO' )' Read recordpoint data from data file "'//TRIM(FileString)//'" ...'
END IF

! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
IF(firstFile_loc.AND.(TRIM(ProjectName).EQ.TRIM(''))) THEN
  CALL ReadAttribute(File_ID,'ProjectName',1,StrScalar=ProjectName)
END IF
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=FileType)
IF(TRIM(FileType) .NE. 'RecordPoints_Data') THEN
  WRITE(UNIT_stdOut,'(A)')' No valid RP data file, skipping!'
  RETURN
END IF
! Check the RP definition file path in the dataset file
CALL ReadAttribute(File_ID,'RPDefFile',1,StrScalar=RP_DefFile_loc)
CALL CloseDataFile()
IF(firstFile_loc) THEN
  IF(.NOT.RP_SET_defined) THEN
    RP_DefFile=RP_DefFile_loc
  END IF
  WRITE(UNIT_stdOut,'(A,A)')' Using RP set file ',TRIM(RP_DefFile)
  CALL InitRPSet(RP_DefFile)
END IF
IF((RP_DefFile_loc.NE.RP_DefFile).AND.(.NOT.RP_SET_defined)) THEN
  WRITE(UNIT_stdOut,'(A,A,A,A)')' RP set file differs from previous one: ',&
                                 TRIM(RP_DefFile_loc),' -- ',TRIM(RP_DefFile)
  RP_DefFile=RP_DefFile_loc
  WRITE(UNIT_stdOut,'(A)')' Initializing new RP set.'
  CALL ChangeRPSet(RP_DefFile)
END IF

CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
! Readin RP Data
CALL GetDataSize(File_ID,'RP_Data',nDims,HSize)
IF(nRP_HDF5 .NE. HSize(2)) THEN
  WRITE(UNIT_stdOut,'(A)')' Number of RPs do not match RP definition file, skipping!'
  RETURN
END IF
IF(nDims.EQ.2) THEN
  nSamples_loc=1
ELSE
  nSamples_loc=INT(HSize(3) )
END IF
! first file: get VarNames and create first dataset pointer.
IF(firstFile_loc.EQV..TRUE.) THEN
  nVar_HDF5 = INT(HSize(1) -1)
  ALLOCATE(VarNames_HDF5(nVar_HDF5))
  CALL ReadAttribute(File_ID,'VarNames',nVar_HDF5,StrArray=VarNames_HDF5)
  nSamples_global=1
  IF(MOD(nSamples_loc,skip).EQ.0) THEN
    nSamples_skip=nSamples_loc/skip
  ELSE
    nSamples_skip=nSamples_loc/skip+1
  END IF
  ! in this case, assign first data set pointer
  CALL getNewRPDataSet(firstset,nSamples_skip)
  actualset=>firstset
ELSE
  !create next dataset
  IF(MOD(nSamples_loc,skip).EQ.0) THEN
    nSamples_skip=nSamples_loc/skip
  ELSE
    nSamples_skip=nSamples_loc/skip+1
  END IF
  CALL getNewRPDataSet(actualset%nextset,nSamples_skip)
  actualset=>actualset%nextset
END IF
DEALLOCATE(HSize)

ALLOCATE(temparray(0:nVar_HDF5,1:nRP_HDF5,1:nSamples_loc)) ! storing complete sample set
#if USE_MPI
! check array data size.
nTotal=REAL((nVar_HDF5+1)*nRP_HDF5*nSamples_loc)
!limit=(2**31-1)/8.
limit=2**28-1/8. ! max. 32 bit integer / 8
IF(nTotal.GT.limit)THEN
  ! read in chunkwise
  nChunks=2*CEILING(nTotal/limit)
  nSamples_chunk=CEILING(REAL(nSamples_loc)/REAL(nChunks))
  offset=0
  DO iChunk=1,nChunks
    IF(iChunk.EQ.nChunks) nSamples_chunk=nSamples_loc-offset
    CALL ReadArray('RP_Data',3,(/nVar_HDF5+1,nRP_HDF5,nSamples_chunk/),offset,3, &
                   RealArray=temparray(:,:,offset+1:offset+nSamples_chunk))
    offset=offset+nSamples_chunk
  END DO
ELSE
#endif
CALL ReadArray('RP_Data',3,(/nVar_HDF5+1,nRP_HDF5,nSamples_loc/),0,3,RealArray=temparray)
#if USE_MPI
END IF
#endif
DO iSample=1,nSamples_skip
  actualset%data(:,:,iSample)=temparray(:,RPOutMap(:),(iSample-1)*skip+1) !RPOutMap filters out RPs which are to be visualized
END DO
DEALLOCATE(temparray)
nSamples_global=nSamples_global+nSamples_skip-1


CALL CloseDataFile()

IF (firstFile_loc) THEN
  WRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' Read recordpoint data from data file "'//TRIM(FileString)//'" ... done'
ELSE
  WRITE(UNIT_stdOut,'(A)',ADVANCE='YES') ' done'
END IF

END SUBROUTINE ReadRPData


!===================================================================================================================================
!> Assemble the data from the seperate .h5 recordpoint files into the global RPData array
!===================================================================================================================================
SUBROUTINE AssembleRPData()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars        ,ONLY: firstset,actualset,nSamples_global, RPData, RPTime,nVar_HDF5
USE MOD_RPSetVisuVisu_Vars ,ONLY: nRP_global
USE MOD_ParametersVisu     ,ONLY: skip
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: nSamples_loc
INTEGER                  :: iStart,iEnd
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')   ' Assemble recordpoint data...'
WRITE(UNIT_stdOut,'(A,I8)')' | total number of samples : ',nSamples_global

ALLOCATE(RPData(1:nVar_HDF5,nRP_global,nSamples_global))
ALLOCATE(RPTime(nSamples_global))
actualset=>firstset
iStart=1
RPTime(1)=firstset%data(0,1,1)
IF(skip.EQ.1) THEN
  DO WHILE(ASSOCIATED(actualset))
    IF(RPTime(iStart).EQ.actualset%data(0,1,1)) THEN          ! dataset only valid if its first sample has same time as the last
      nSamples_loc=actualset%nSamples                         ! sample of its precursor
      iEnd=MIN(iStart-1+nSamples_loc,nSamples_global)
      RPData(1:nVar_HDF5,:,iStart:iEnd)=actualset%data(1:nVar_HDF5,:,1:nSamples_loc)
      RPTime(iStart:iEnd)=actualset%data(0,1,1:nSamples_loc)  ! each RP has same time...
      iStart=iEnd
      actualset=>actualset%nextset
    ELSE
      CALL Abort(__STAMP__,'ERROR - Time History of RP Data Files does not match!')
    END IF
  END DO
ELSE
  DO WHILE(ASSOCIATED(actualset))
    nSamples_loc=actualset%nSamples                         ! sample of its precursor
    iEnd=MIN(iStart-1+nSamples_loc,nSamples_global)
    RPData(1:nVar_HDF5,:,iStart:iEnd)=actualset%data(1:nVar_HDF5,:,1:nSamples_loc)
    RPTime(iStart:iEnd)=actualset%data(0,1,1:nSamples_loc)  ! each RP has same time...
    iStart=iEnd
    actualset=>actualset%nextset
  END DO
END IF

WRITE(UNIT_stdOut,'(A)')  ' Assemble recordpoint data done'
END SUBROUTINE AssembleRPData


!===================================================================================================================================
!> Deallocate global variable for Recordpoints data module
!===================================================================================================================================
SUBROUTINE FinalizeRPData()
! MODULES
USE MOD_RPData_Vars
IMPLICIT NONE
!===================================================================================================================================
SDEALLOCATE(VarNames_HDF5)
SDEALLOCATE(RPData)
SDEALLOCATE(RPTime)
END SUBROUTINE FinalizeRPData

END MODULE MOD_RPData
