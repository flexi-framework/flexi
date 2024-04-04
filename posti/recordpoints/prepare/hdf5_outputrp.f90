!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
!> Module to handle the output of the RPs
!===================================================================================================================================
MODULE MOD_HDF5_OutputRP
! MODULES
USE MOD_IO_HDF5, ONLY: FILE_ID,OpenDataFile,CloseDataFile
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE WriteRecordPointstoHDF5
  MODULE PROCEDURE WriteRecordPointstoHDF5
END INTERFACE

PUBLIC :: WriteRecordPointstoHDF5
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Subroutine to write the recordpoints to HDF5 format
!===================================================================================================================================
SUBROUTINE WriteRecordPointstoHDF5(ProjectName,MeshFileName)
! MODULES
USE MOD_Globals
USE MOD_HDF5_Output
USE MOD_Mesh_Vars       ,ONLY:NGeo,nGlobalElems
USE MOD_Interpolation_Vars,ONLY: NodeType
USE MOD_RPSet_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: ProjectName
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iGr,iL,iRP,iPl,iBx
INTEGER                        :: i,j,k
CHARACTER(LEN=255)             :: FileName,FileString,tmp255
CHARACTER(LEN=5)               :: PlaneType
CHARACTER(LEN=3)               :: BoxType
TYPE(tGroup),POINTER           :: Group
TYPE(tLine),POINTER            :: Line
TYPE(tPlane),POINTER           :: Plane
TYPE(tBox),POINTER             :: Box
INTEGER,ALLOCATABLE            :: RPset(:),RPset2D(:,:),RPset3D(:,:,:)
REAL,ALLOCATABLE               :: x_dummy(:,:)
INTEGER                        :: s1(1),s2(2),s3(3)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE RECORDPOINTS TO HDF5 FILE...'
FileName=TRIM(ProjectName)//'_RPSet'
FileString=TRIM(FileName)//'.h5'
CALL OpenDataFile(TRIM(Filestring),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

CALL WriteAttribute(File_ID,'File_Type',1,StrScalar=(/CHARACTER(LEN=255)::'RecordPoints'/))
tmp255=TRIM(MeshFileName)
CALL WriteAttribute(File_ID,'MeshFile',1,StrScalar=(/tmp255/))
CALL WriteAttribute(File_ID,'NGeo',1,IntScalar=NGeo)
tmp255=TRIM(NodeType)
CALL WriteAttribute(File_ID,'NodeType',1,StrScalar=(/tmp255/))

! Write RP Set to File -------------------------------------------------------------------------------------------------------------
! 1. Groups
s1=nGroups
CALL WriteArray('GroupNames',1,s1,s1,(/0/),.FALSE.,StrArray=Groups(:)%Name)
DO iGr=1,nGroups
  Group=>Groups(iGr)
  IF(Group%nRP.GT.0) THEN
    ALLOCATE(RPset(1:Group%nRP))
    DO iRP=1,Group%nRP
      RPset(iRP)=Group%RP_ptr(iRP)%RP%ID
    END DO
    s1=Group%nRP
    CALL WriteArray(TRIM(Group%Name),1,s1,s1,(/0/),.FALSE.,IntArray=RPset)
    DEALLOCATE(RPset)
  END IF !Group%nRP.GT.0
END DO !iGr

! 2. Lines
IF(nLines.GT.0) THEN
  s1=nLines
  CALL WriteArray('LineNames',1,s1,s1,(/0/),.FALSE.,StrArray=Lines(:)%Name)
  DO iL=1,nLines
    Line=>Lines(iL)
    ALLOCATE(RPset(1:Line%nRP))
    DO iRP=1,Line%nRP
      RPset(iRP)=Line%RP_ptr(iRP)%RP%ID
    END DO
    s1=Line%nRP
    CALL WriteArray(TRIM(Line%Name),1,s1,s1,(/0/),.FALSE.,IntArray=RPset)
    DEALLOCATE(RPset)
    ! groupID
    CALL WriteAttribute(FILE_ID,'GroupID',1,DataSetname=Line%Name,IntScalar=Line%GroupID)
  END DO !iL
END IF

! 3. Points
IF(nPoints.GT.0) THEN
  ALLOCATE(RPset(1:nPoints))
  DO iRP=1,nPoints
    RPset(iRP)=Points(iRP)%RP%ID
  END DO
  s1=nPoints
  CALL WriteArray('Points_IDlist'     ,1,s1,s1,(/0/),.FALSE.,IntArray=RPSet)
  CALL WriteArray('Points_GroupIDlist',1,s1,s1,(/0/),.FALSE.,IntArray=Points(:)%GroupID)
END IF

! 4. Planes
IF(nPlanes.GT.0) THEN
  s1=nPlanes
  CALL WriteArray('PlaneNames',1,s1,s1,(/0/),.FALSE.,StrArray=Planes(:)%Name)
  DO iPl=1,nPlanes
    Plane=>Planes(iPl)
    ALLOCATE(RPset2D(1:Plane%nRP(1),1:Plane%nRP(2)))
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        RPset2D(i,j)=Plane%RP_ptr(i,j)%RP%ID
      END DO !i
    END DO !j
    CALL WriteArray(TRIM(Plane%Name),2,Plane%nRP,Plane%nRP,(/0,0/),.FALSE.,IntArray=RPset2D)
    DEALLOCATE(RPset2D)
    ! groupID
    CALL WriteAttribute(FILE_ID,'GroupID',1,DataSetname=Plane%Name,IntScalar=Plane%GroupID)

    ! write normal and tangent vectors in case of the BLPlane
    PlaneType=TRIM(Plane%Name(1:5))
    IF(PlaneType.EQ.TRIM("BLPla")) THEN
      s2=(/3,Plane%nRP(1)/)
      WRITE(tmp255,'(A,A)')TRIM(Plane%Name),'_NormVec'
      CALL WriteArray(tmp255,2,s2,s2,(/0,0/),.FALSE.,RealArray=Plane%NormVec)
      WRITE(tmp255,'(A,A)')TRIM(Plane%Name),'_TangVec'
      CALL WriteArray(tmp255,2,s2,s2,(/0,0/),.FALSE.,RealArray=Plane%TangVec)
    END IF
  END DO !iPl
END IF

! 5. Boxes
IF(nBoxes.GT.0) THEN
  s1=nBoxes
  CALL WriteArray('BoxNames',1,s1,s1,(/0/),.FALSE.,StrArray=Boxes(:)%Name)
  DO iBx=1,nBoxes
    Box=>Boxes(iBx)
    ALLOCATE(RPset3D(1:Box%nRP(1),1:Box%nRP(2),1:Box%nRP(3)))
    DO k=1,Box%nRP(3)
      DO j=1,Box%nRP(2)
        DO i=1,Box%nRP(1)
          RPset3D(i,j,k)=Box%RP_ptr(i,j,k)%RP%ID
        END DO !i
      END DO !j
    END DO !k
    CALL WriteArray(TRIM(Box%Name),3,Box%nRP,Box%nRP,(/0,0,0/),.FALSE.,IntArray=RPset3D)
    DEALLOCATE(RPset3D)
    ! groupID
    CALL WriteAttribute(FILE_ID,'GroupID',1,DataSetname=Box%Name,IntScalar=Box%GroupID)

    ! write normal and tangent vectors in case of the BLBox
    BoxType=TRIM(Box%Name(1:3))
    IF(BoxType.EQ.TRIM("BLB")) THEN
      s3=(/3,Box%nRP(1),Box%nRP(3)/)
      WRITE(tmp255,'(A,A)')TRIM(Box%Name),'_NormVec'
      CALL WriteArray(tmp255,3,s3,s3,(/0,0,0/),.FALSE.,RealArray=Box%NormVec)
      WRITE(tmp255,'(A,A)')TRIM(Box%Name),'_TangVec'
      CALL WriteArray(tmp255,3,s3,s3,(/0,0,0/),.FALSE.,RealArray=Box%TangVec)
    END IF
  END DO !iBx
END IF

! Offset Array
s2=(/2,nGlobalElems/)
CALL WriteArray('OffsetRP',2,s2,s2,(/0,0/),.FALSE.,IntArray=OffsetRP)

! Coordinates
ALLOCATE(x_dummy(1:3,1:nRP_global))
s2=(/3,nRP_global/)
DO iRP=1,nRP_global
  x_dummy(:,iRP)=RPlist(iRP)%RP%xF
END DO
CALL WriteArray('xF_RP',2,s2,s2,(/0,0/),.FALSE.,RealArray=x_dummy)
DO iRP=1,nRP_global
  x_dummy(:,iRP)=RPlist(iRP)%RP%x
END DO
CALL WriteArray('x_RP' ,2,s2,s2,(/0,0/),.FALSE.,RealArray=x_dummy)
DO iRP=1,nRP_global
  x_dummy(:,iRP)=RPlist(iRP)%RP%xi
END DO
CALL WriteArray('xi_RP',2,s2,s2,(/0,0/),.FALSE.,RealArray=x_dummy)
DEALLOCATE(x_dummy)
! Close the file.
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteRecordPointstoHDF5

END MODULE MOD_HDF5_OutputRP
