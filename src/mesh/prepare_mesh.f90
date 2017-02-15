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
!> Contains subroutines to collect and condense mesh topology information (master,slave,flip)
!> and define communication neighbours.
!==================================================================================================================================
MODULE MOD_Prepare_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE setLocalSideIDs
  MODULE PROCEDURE setLocalSideIDs
END INTERFACE

INTERFACE fillMeshInfo
  MODULE PROCEDURE fillMeshInfo
END INTERFACE

PUBLIC::setLocalSideIDs,fillMeshInfo

#if MPI
INTERFACE exchangeFlip
  MODULE PROCEDURE exchangeFlip
END INTERFACE

PUBLIC::exchangeFlip
#endif
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine sorts sides into three groups containing BC sides, inner sides and MPI sides in the following manner:
!> 
!> * BCSides         : sides with boundary conditions (excluding periodic BCs)
!> * InnerMortars    : "virtual" sides introduced for collecting the data of the small sides at a non-conforming interface
!> * InnerSides      : normal inner sides
!> * MPI_MINE sides  : MPI sides to be processed by the current processor (e.g. flux computation)
!> * MPI_YOUR sides  : MPI sides to be processed by the neighbour processor
!> * MPIMortars      : mortar interfaces to be comunicated
!> 
!> Each side can be accessed through its SideID defining its position in the processor local side list.
!> The routine furthermore sets the MPI masters and slave sides.
!==================================================================================================================================
SUBROUTINE setLocalSideIDs()
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,  ONLY:tElem,tSide
USE MOD_Mesh_Vars,  ONLY: nElems,nInnerSides,nSides,nBCSides,offsetElem
USE MOD_Mesh_Vars,  ONLY: Elems,nMPISides_MINE,nMPISides_YOUR,BoundaryType,nBCs
USE MOD_Mesh_Vars,  ONLY: nMortarSides,nMortarInnerSides,nMortarMPISides
#if MPI
USE MOD_ReadInTools,ONLY: GETLOGICAL
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
INTEGER   :: iElem,FirstElemInd,LastElemInd
INTEGER   :: iLocSide,iSide,iInnerSide,iBCSide
INTEGER   :: iMortar,iMortarInnerSide,iMortarMPISide,nMortars
INTEGER   :: i,j
INTEGER   :: PeriodicBCMap(nBCs)       !connected periodic BCs
#if MPI
INTEGER   :: nSmallMortarSides
INTEGER   :: nSmallMortarInnerSides
INTEGER   :: nSmallMortarMPISides_MINE
INTEGER   :: nSmallMortarMPISides_YOUR
INTEGER               :: iNbProc,ioUnit,addToInnerMortars
INTEGER               :: ProcInfo(9),nNBmax      !for output only
INTEGER,ALLOCATABLE   :: SideIDMap(:)
INTEGER,ALLOCATABLE   :: NBinfo(:,:),NBinfo_glob(:,:,:),nNBProcs_glob(:),Procinfo_glob(:,:),tmparray(:,:)  !for output only
REAL,ALLOCATABLE      :: tmpreal(:,:)
CHARACTER(LEN=10)     :: formatstr
CHARACTER(LEN=255)    :: PartitionInfoFileName
LOGICAL               :: writePartitionInfo
#endif
!==================================================================================================================================
FirstElemInd= offsetElem+1
LastElemInd = offsetElem+nElems
! ----------------------------------------
! Set side IDs to arrange sides:
! 1. BC sides
! 2. inner sides
! 3. MPI sides
! MPI Sides are not included here!
! ----------------------------------------
! Get connection between periodic BCs
PeriodicBCMap=-2
DO i=1,nBCs
  IF((BoundaryType(i,BC_TYPE).NE.1)) PeriodicBCMap(i)=-1 ! not periodic
  IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).GT.0)) PeriodicBCMap(i)=-1 ! slave
  IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).LT.0))THEN
    DO j=1,nBCs
      IF(BoundaryType(j,BC_TYPE).NE.1) CYCLE
      IF(BoundaryType(j,BC_ALPHA).EQ.(-BoundaryType(i,BC_ALPHA))) PeriodicBCMap(i)=j
    END DO
  END IF
END DO
IF(ANY(PeriodicBCMap.EQ.-2))&
  CALL abort(__STAMP__,'Periodic connection not found.')

DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp

      aSide%sideID=-1
      ! periodics have two bcs: set to (positive) master bc (e.g. from -1 to 1)
      IF(aSide%BCIndex.GE.1)THEN
        IF(PeriodicBCMap(aSide%BCIndex).NE.-1)&
          aSide%BCIndex=PeriodicBCMap(aSide%BCIndex)
      END IF
    END DO !iMortar
  END DO
END DO

nMortarInnerSides=0
nMortarMPISides=0
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    aSide%tmp=0
    IF(aSide%nMortars.GT.0)THEN
      DO iMortar=1,aSide%nMortars
        IF(aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp%nbProc.NE.-1)THEN
          aSide%tmp=-1
          EXIT
        END IF
      END DO !iMortar
      IF(aSide%tmp.EQ.-1) THEN
        nMortarMPISides=nMortarMPISides+1
      ELSE
        nMortarInnerSides=nMortarInnerSides+1
      END IF
    END IF !nMortars>0
  END DO
END DO
IF((nMortarInnerSides+nMortarMPISides).NE.nMortarSides) &
   CALL abort(__STAMP__,'nInner+nMPI mortars <> nMortars.')

iSide=0
iBCSide=0
iMortarInnerSide=nBCSides
iInnerSide=nBCSides+nMortarInnerSides
iMortarMPISide=nSides-nMortarMPISides
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp

      IF(aSide%sideID.EQ.-1)THEN
        IF(aSide%NbProc.EQ.-1)THEN ! no MPI Sides
          IF(ASSOCIATED(aSide%connection))THEN
            iInnerSide=iInnerSide+1
            iSide=iSide+1
            aSide%SideID=iInnerSide
            aSide%connection%SideID=iInnerSide
          ELSE
            IF(aSide%MortarType.GT.0) THEN
              IF(aSide%tmp.EQ.-1)THEN !MPI mortar side
                iMortarMPISide=iMortarMPISide+1
                aSide%SideID=iMortarMPISide
              ELSE
                iMortarInnerSide=iMortarInnerSide+1
                iSide=iSide+1
                aSide%SideID=iMortarInnerSide
              END IF !mpi mortar
            ELSE !this is now a BC side, really!
              iBCSide=iBCSide+1
              iSide=iSide+1
              aSide%SideID=iBCSide
            END IF !mortar
          END IF !associated connection
        END IF ! .NOT. MPISide
      END IF !sideID NE -1
    END DO !iMortar
  END DO ! iLocSide=1,6
END DO !iElem
IF(iSide.NE.nInnerSides+nBCSides+nMortarInnerSides) STOP 'not all SideIDs are set!'
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A22,I8)')'nMortarSides:',nMortarSides
LOGWRITE(*,'(A22,I8)')'nMortarInnerSides:',nMortarInnerSides
LOGWRITE(*,'(A22,I8)')'nMortarMPISides:',nMortarMPISides
LOGWRITE(*,*)'-------------------------------------------------------'

nMPISides_MINE=0
nMPISides_YOUR=0
#if MPI
! SPLITTING MPISides in MINE and YOURS
ALLOCATE(nMPISides_MINE_Proc(1:nNbProcs),nMPISides_YOUR_Proc(1:nNbProcs))
nMPISides_MINE_Proc=0
nMPISides_YOUR_Proc=0
DO iNbProc=1,nNbProcs
  IF(myRank.LT.NbProc(iNbProc)) THEN
    nMPISides_MINE_Proc(iNbProc)=nMPISides_Proc(iNbProc)/2
  ELSE
    nMPISides_MINE_Proc(iNbProc)=nMPISides_Proc(iNbProc)-nMPISides_Proc(iNbProc)/2
  END IF
  nMPISides_YOUR_Proc(iNbProc)=nMPISides_Proc(iNbProc)-nMPISides_MINE_Proc(iNbProc)
END DO
nMPISides_MINE=SUM(nMPISides_MINE_Proc)
nMPISides_YOUR=SUM(nMPISides_YOUR_Proc)

ALLOCATE(offsetMPISides_YOUR(0:nNbProcs),offsetMPISides_MINE(0:nNbProcs))
offsetMPISides_MINE=0
offsetMPISides_YOUR=0
! compute offset, first all MINE , then all YOUR MPISides
offsetMPISides_MINE(0)=nInnerSides+nBCSides+nMortarInnerSides
DO iNbProc=1,nNbProcs
  offsetMPISides_MINE(iNbProc)=offsetMPISides_MINE(iNbProc-1)+nMPISides_MINE_Proc(iNbProc)
END DO
offsetMPISides_YOUR(0)=offsetMPISides_MINE(nNbProcs)
DO iNbProc=1,nNbProcs
  offsetMPISides_YOUR(iNbProc)=offsetMPISides_YOUR(iNbProc-1)+nMPISides_YOUR_Proc(iNbProc)
END DO

DO iNbProc=1,nNbProcs
  ALLOCATE(SideIDMap(nMPISides_Proc(iNbProc)))
  iSide=0
  DO iElem=FirstElemInd,LastElemInd
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      nMortars=aSide%nMortars
      DO iMortar=0,nMortars
        IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
        IF(aSide%NbProc.NE.NbProc(iNbProc))CYCLE
        iSide=iSide+1
        !trick: put non-mortars first to optimize addtoInnerMortars (also small mortar sides are marked with MortarType<0)
        IF((iMortar.EQ.0).AND.(aSide%MortarType.EQ.0)) aSide%ind=-aSide%ind
        !
        SideIDMap(iSide)=aSide%ind !global Side Index
      END DO !iMortar
    END DO !iLocSide
  END DO !iElem
  ! TODO: Check if this is OK with mortars
  IF(iSide.GT.1) CALL MergeSort(SideIDMap(1:iSide),iSide) !sort by global side index
  DO iElem=FirstElemInd,LastElemInd
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      nMortars=aSide%nMortars
      DO iMortar=0,nMortars
        IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
        IF(aSide%NbProc.NE.NbProc(iNbProc))CYCLE
        aSide%SideID=INVMAP(aSide%ind,nMPISides_Proc(iNbProc),SideIDMap) ! get sorted iSide
        IF(myRank.LT.aSide%NbProc)THEN
          IF(aSide%SideID.LE.nMPISides_MINE_Proc(iNbProc))THEN !MINE
            aSide%SideID=aSide%SideID +offsetMPISides_MINE(iNbProc-1)
          ELSE !YOUR
            aSide%SideID=(aSide%SideID-nMPISides_MINE_Proc(iNbProc))+offsetMPISides_YOUR(iNbProc-1)
          END IF
        ELSE
          IF(aSide%SideID.LE.nMPISides_YOUR_Proc(iNbProc))THEN !MINE
            aSide%SideID=aSide%SideID +offsetMPISides_YOUR(iNbProc-1)
          ELSE !YOUR
            aSide%SideID=(aSide%SideID-nMPISides_YOUR_Proc(iNbProc))+offsetMPISides_MINE(iNbProc-1)
          END IF
        END IF !myrank<NbProc
      END DO !iMortar
    END DO !iLocSide
  END DO !iElem
  DEALLOCATE(SideIDMap)
END DO !nbProc(i)
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
      aSide%ind=ABS(aSide%ind) ! set back trick
    END DO !iMortar
  END DO !iLocSide
END DO !iElem

! optimize mortars: search for mortars being fully MPI_MINE and add them to innerMortars
IF(nMortarSides.GT.0)THEN
  addToInnerMortars=0
  DO iElem=FirstElemInd,LastElemInd
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      aSide%tmp=0
      DO iMortar=1,aSide%nMortars
        aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp%tmp=0
      END DO !iMortar
    END DO !iLocSide
  END DO !iElem
  DO iElem=FirstElemInd,LastElemInd
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      IF(aSide%nMortars.GT.0)THEN
        aSide%tmp=-1 !mortar side
        DO iMortar=1,aSide%nMortars
          IF(aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp%SideID.GT.offsetMPISides_YOUR(0))THEN
            aSide%tmp=-2  !mortar side with side used in MPI_YOUR
            EXIT
          END IF
        END DO !iMortar
        IF(aSide%tmp.EQ.-1) THEN
          addToInnerMortars=addToInnerMortars+1
        END IF
      END IF !nMortars>0
    END DO !iLocSide
  END DO !iElem
  addToInnerMortars=addToInnerMortars-nMortarInnerSides
  IF(addToInnerMortars.GT.0)THEN
    iMortarInnerSide=nBCSides+nMortarInnerSides
    DO iElem=FirstElemInd,LastElemInd
      aElem=>Elems(iElem)%ep
      DO iLocSide=1,6
        aSide=>aElem%Side(iLocSide)%sp
        IF((aSide%tmp.EQ.0).AND.(aSide%SideID.GT.iMortarInnerSide))THEN
          !shift SideID
          aSide%SideID=aSide%SideID+addToInnerMortars
          aSide%tmp=1
        END IF
        nMortars=aSide%nMortars
        DO iMortar=1,nMortars
          aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
          IF((aSide%tmp.EQ.0).AND.(aSide%SideID.GT.iMortarInnerSide))THEN
            aSide%SideID=aSide%SideID+addToInnerMortars
            aSide%tmp=1
          END IF
        END DO !iMortar
      END DO !iLocSide
    END DO !iElem
    offsetMPISides_MINE=offsetMPISides_MINE+addToInnerMortars
    offsetMPISides_YOUR=offsetMPISides_YOUR+addToInnerMortars

    nMortarInnerSides=nMortarInnerSides+addToInnerMortars
    nMortarMPISides  =nMortarMPISides-addToInnerMortars
    iMortarMPISide=nSides-nMortarMPISides
    iMortarInnerSide=nBCSides

    DO iElem=FirstElemInd,LastElemInd
      aElem=>Elems(iElem)%ep
      DO iLocSide=1,6
        aSide=>aElem%Side(iLocSide)%sp
        IF(aSide%tmp.EQ.-2)THEN !MPI mortars, renumber SideIDs
          iMortarMPISide=iMortarMPISide+1
          aSide%SideID=iMortarMPISide
        ELSEIF(aSide%tmp.EQ.-1)THEN !innermortars mortars, renumber SideIDs
          iMortarInnerSide=iMortarInnerSide+1
          aSide%SideID=iMortarInnerSide
        END IF !aSide%tmp==-1
      END DO !iLocSide
    END DO !iElem
  END IF !addToInnerMortars>0
  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A22,I8)')'addToInnerMortars:',addToInnerMortars
  LOGWRITE(*,'(A22,I8)')'new nMortarSides:',nMortarSides
  LOGWRITE(*,'(A22,I8)')'new nMortarInnerSides:',nMortarInnerSides
  LOGWRITE(*,'(A22,I8)')'new nMortarMPISides:',nMortarMPISides
  LOGWRITE(*,*)'-------------------------------------------------------'
END IF !nMortarSides>0

nSmallMortarSides=0
nSmallMortarMPIsides_MINE=0
nSmallMortarMPIsides_YOUR=0
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%nMortars.GT.0)THEN !mortar side
      nSmallMortarSides=nSmallMortarSides+aSide%nMortars
      DO iMortar=1,aSide%nMortars
        IF (aSide%MortarSide(iMortar)%sp%SideID.GT.offsetMPISides_YOUR(0))THEN
          nSmallMortarMPIsides_YOUR=nSmallMortarMPISides_YOUR+1
        ELSE
          IF(aSide%MortarSide(iMortar)%sp%SideID.GT.offsetMPISides_MINE(0))THEN
            nSmallMortarMPISides_MINE=nSmallMortarMPISides_MINE+1
          END IF
        END IF
      END DO !iMortar
    END IF !mortarSide
  END DO ! LocSideID
END DO ! iElem
nSmallMortarInnerSides=nSmallMortarSides-nSmallMortarMPISides_MINE-nSmallMortarMPISides_YOUR

!------------------------------------------------------
! Copy data into some MPI arrays
!------------------------------------------------------

ALLOCATE(nMPISides_send(       nNbProcs,2))
ALLOCATE(OffsetMPISides_send(0:nNbProcs,2))
ALLOCATE(nMPISides_rec(        nNbProcs,2))
ALLOCATE(OffsetMPISides_rec( 0:nNbProcs,2))
! Set number of sides and offset for SEND MINE - RECEIVE YOUR case
nMPISides_send(:,1)     =nMPISides_MINE_Proc
OffsetMPISides_send(:,1)=OffsetMPISides_MINE
nMPISides_rec(:,1)      =nMPISides_YOUR_Proc
OffsetMPISides_rec(:,1) =OffsetMPISides_YOUR
! Set number of sides and offset for SEND YOUR - RECEIVE MINE case
nMPISides_send(:,2)     =nMPISides_YOUR_Proc
OffsetMPISides_send(:,2)=OffsetMPISides_YOUR
nMPISides_rec(:,2)      =nMPISides_MINE_Proc
OffsetMPISides_rec(:,2) =OffsetMPISides_MINE

!------------------------------------------------------
! From this point on only debug output is performed
!------------------------------------------------------

WRITE(formatstr,'(a5,I2,a3)')'(A22,',nNBProcs,'I8)'
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A22,I8)')'nNbProcs:',nNbProcs
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'NbProc:'   ,NbProc
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'nMPISides_Proc:',nMPISides_Proc
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'nMPISides_MINE_Proc:',nMPISides_MINE_Proc
LOGWRITE(*,formatstr)'nMPISides_YOUR_Proc:',nMPISides_YOUR_Proc
WRITE(formatstr,'(a5,I2,a3)')'(A22,',nNBProcs+1,'I8)'
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'offsetMPISides_MINE:',offsetMPISides_MINE
LOGWRITE(*,formatstr)'offsetMPISides_YOUR:',offsetMPISides_YOUR
LOGWRITE(*,*)'-------------------------------------------------------'

writePartitionInfo = GETLOGICAL('writePartitionInfo','.FALSE.')
IF(.NOT.writePartitionInfo) RETURN
!output partitioning info
ProcInfo(1)=nElems
ProcInfo(2)=nSides
ProcInfo(3)=nInnerSides
ProcInfo(4)=nBCSides
ProcInfo(5)=nMortarInnerSides
ProcInfo(6)=nMortarMPISides
ProcInfo(7)=nSmallMortarInnerSides
ProcInfo(8)=nSmallMortarMPISides_MINE
ProcInfo(9)=nSmallMortarMPISides_YOUR
IF(MPIroot)THEN
  ALLOCATE(nNBProcs_glob(0:nProcessors-1))
  ALLOCATE(ProcInfo_glob(9,0:nProcessors-1))
  nNBProcs_glob=-99999
  Procinfo_glob=-88888
ELSE
  ALLOCATE(nNBProcs_glob(1)) !dummy for debug
  ALLOCATE(ProcInfo_glob(1,1)) !dummy for debug
END IF !MPIroot
CALL MPI_GATHER(nNBProcs,1,MPI_INTEGER,nNBProcs_glob,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
CALL MPI_GATHER(ProcInfo,9,MPI_INTEGER,ProcInfo_glob,9,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
IF(MPIroot)THEN
  nNBmax=MAXVAL(nNBProcs_glob) !count, total number of columns in table
  ALLOCATE(NBinfo_glob(6,nNBmax,0:nProcessors))
  NBinfo_glob=-77777
ELSE
  ALLOCATE(NBinfo_glob(1,1,1)) !dummy for debug
END IF
CALL MPI_BCAST(nNBmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
ALLOCATE(NBinfo(6,nNbmax))
NBinfo=0
NBinfo(1,1:nNBProcs)=NBProc
NBinfo(2,1:nNBProcs)=nMPISides_Proc
NBinfo(3,1:nNBProcs)=nMPISides_MINE_Proc
NBinfo(4,1:nNBProcs)=nMPISides_YOUR_Proc
NBinfo(5,1:nNBProcs)=offsetMPISides_MINE(0:nNBProcs-1)
NBinfo(6,1:nNBProcs)=offsetMPISides_YOUR(0:nNBProcs-1)
CALL MPI_GATHER(NBinfo,6*nNBmax,MPI_INTEGER,NBinfo_glob,6*nNBmax,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
DEALLOCATE(NBinfo)
IF(MPIroot)THEN
  ioUnit=GETFREEUNIT()
  WRITE(PartitionInfoFileName,'(A21,I6.6,A4)')'partitionInfo_nRanks_',nProcessors,'.out'
  OPEN(UNIT=ioUnit,FILE=TRIM(PartitionInfoFileName),STATUS='REPLACE')
  WRITE(ioUnit,*)'Partition Information:'
  WRITE(ioUnit,*)'total number of Procs,',nProcessors
  WRITE(ioUnit,*)'total number of Elems,',SUM(Procinfo_glob(1,:))

  WRITE(ioUnit,'(13(A23))')'Rank','nElems','nSides','nInnerSides','nBCSides','nMPISides','nMPISides_MINE','nNBProcs', &
              'nMortarInnerSides', 'nMortarMPISides', 'nSmallMortInnerSides', 'nSmallMortMPISidesMINE', 'nSmallMortMPISidesYOUR'
  WRITE(ioUnit,'(300("="))')
  !statistics
  ALLOCATE(tmparray(12,0:3),tmpreal(12,2))
  tmparray(:,0)=0      !tmp
  tmparray(:,1)=0      !mean
  tmparray(:,2)=-HUGE(1) !max
  tmparray(:,3)= HUGE(1) !min
  DO i=0,nProcessors-1
    !actual proc
    tmparray( 1,0)=Procinfo_glob(1,i)
    tmparray( 2,0)=Procinfo_glob(2,i)
    tmparray( 3,0)=Procinfo_glob(3,i)
    tmparray( 4,0)=Procinfo_glob(4,i)
    tmparray( 5,0)=SUM(NBinfo_glob(2,:,i))
    tmparray( 6,0)=SUM(NBinfo_glob(3,:,i))
    tmparray( 7,0)=nNBProcs_glob(i)
    tmparray( 8,0)=Procinfo_glob(5,i)
    tmparray( 9,0)=Procinfo_glob(6,i)
    tmparray(10,0)=Procinfo_glob(7,i)
    tmparray(11,0)=Procinfo_glob(8,i)
    tmparray(12,0)=Procinfo_glob(9,i)
    DO j=1,12
      !mean
      tmparray(j,1)=tmparray(j,1)+tmparray(j,0)
      !max
      tmparray(j,2)=MAX(tmparray(j,2),tmparray(j,0))
      tmparray(j,3)=MIN(tmparray(j,3),tmparray(j,0))
    END DO
  END DO
  tmpreal(:,1)=REAL(tmparray(:,1))/REAL(nProcessors) !mean in REAL
  tmpreal(:,2)=0.   !RMS
  DO i=0,nProcessors-1
    !actual proc
    tmparray( 1,0)=Procinfo_glob(1,i)
    tmparray( 2,0)=Procinfo_glob(2,i)
    tmparray( 3,0)=Procinfo_glob(3,i)
    tmparray( 4,0)=Procinfo_glob(4,i)
    tmparray( 5,0)=SUM(NBinfo_glob(2,:,i))
    tmparray( 6,0)=SUM(NBinfo_glob(3,:,i))
    tmparray( 7,0)=nNBProcs_glob(i)
    tmparray( 8,0)=Procinfo_glob(5,i)
    tmparray( 9,0)=Procinfo_glob(6,i)
    tmparray(10,0)=Procinfo_glob(7,i)
    tmparray(11,0)=Procinfo_glob(8,i)
    tmparray(12,0)=Procinfo_glob(9,i)
    DO j=1,12
      tmpreal(j,2)=tmpreal(j,2)+(tmparray(j,0)-tmpreal(j,1))**2
    END DO
  END DO
  tmpreal(:,2)=SQRT(tmpreal(:,2)/REAL(nProcessors))
  WRITE(ioUnit,'(A23,12(13X,F10.2))')'   MEAN        ',tmpreal(:,1)
  WRITE(ioUnit,'(300("-"))')
  WRITE(ioUnit,'(A23,12(13X,F10.2))')'   RMS         ',tmpreal(:,2)
  WRITE(ioUnit,'(300("-"))')
  WRITE(ioUnit,'(A23,12(13X,I10))')'   MIN         ',tmparray(:,3)
  WRITE(ioUnit,'(300("-"))')
  WRITE(ioUnit,'(A23,12(13X,I10))')'   MAX         ',tmparray(:,2)
  WRITE(ioUnit,'(300("="))')
  DO i=0,nProcessors-1
    WRITE(ioUnit,'(13(13X,I10))')i, &
          Procinfo_glob(1:4,i),SUM(NBinfo_glob(2,:,i)),SUM(NBinfo_glob(3,:,i)),nNBProcs_glob(i),Procinfo_glob(5:9,i)
    WRITE(ioUnit,'(300("-"))')
  END DO
  WRITE(ioUnit,*)' '
  WRITE(ioUnit,*)'Information per neighbor processor'
  WRITE(ioUnit,*)' '
  WRITE(ioUnit,'(7(A15))')'Rank','NBProc','nMPISides_Proc','nMPISides_MINE','nMPISides_YOUR','offset_MINE','offset_YOUR'
  WRITE(ioUnit,'(105("="))')
  DO i=0,nProcessors-1
    WRITE(ioUnit,'(7(5X,I10))')i,NBinfo_glob(:,1,i)
    DO j=2,nNBProcs_glob(i)
      WRITE(ioUnit,'(A15,6(5X,I10))')' ',NBinfo_glob(:,j,i)
    END DO
    WRITE(ioUnit,'(105("-"))')
  END DO
  DEALLOCATE(tmparray,tmpreal)
  CLOSE(ioUnit)
END IF !MPIroot
DEALLOCATE(NBinfo_glob,nNBProcs_glob,ProcInfo_glob)
#endif /*MPI*/
END SUBROUTINE setLocalSideIDs



!==================================================================================================================================
!> This routine condenses the mesh topology from a pointer-based structure into arrays.
!> The array ElemToSide contains for each elements local side the global SideID and its
!> flip with regard to the neighbour side.
!> The SideToElem array contains for each side the neighbour elements (master and slave)
!> as well as the local side IDs of the side within those elements.
!> The last entry is the flip of the slave with regard to the master element.
!==================================================================================================================================
SUBROUTINE fillMeshInfo()
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY: tElem,tSide,Elems
USE MOD_Mesh_Vars,ONLY: nElems,offsetElem,nBCSides
USE MOD_Mesh_Vars,ONLY: firstMortarInnerSide,lastMortarInnerSide,nMortarInnerSides,firstMortarMPISide
USE MOD_Mesh_Vars,ONLY: ElemToSide,SideToElem,BC,AnalyzeSide
USE MOD_Mesh_Vars,ONLY: MortarType,MortarInfo
#if MPI
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide,mSide
INTEGER             :: iElem,LocSideID,nSides_flip(0:4),SideID
INTEGER             :: nSides_MortarType(1:3),iMortar
!==================================================================================================================================
! ELement to Side mapping
nSides_flip=0
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    ElemToSide(E2S_SIDE_ID,LocSideID,iElem)=aSide%SideID
    ElemToSide(E2S_FLIP,LocSideID,iElem)   =aSide%Flip
    nSides_flip(aSide%flip)=nSides_flip(aSide%flip)+1
  END DO ! LocSideID
END DO ! iElem

! Side to Element mapping, sorted by SideID
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%Flip.EQ.0)THEN !root side
      SideToElem(S2E_ELEM_ID,aSide%SideID)         = iElem !root Element
      SideToElem(S2E_LOC_SIDE_ID,aSide%SideID)     = LocSideID
      AnalyzeSide(aSide%sideID)                    = aSide%BCIndex
    ELSE
      SideToElem(S2E_NB_ELEM_ID,aSide%SideID)      = iElem ! element with flipped side
      SideToElem(S2E_NB_LOC_SIDE_ID,aSide%SideID)  = LocSideID
      SideToElem(S2E_FLIP,aSide%SideID)            = aSide%Flip
    END IF
    IF(aSide%sideID .LE. nBCSides) BC(aSide%sideID)=aSide%BCIndex
  END DO ! LocSideID
END DO ! iElem

! Mapping of Mortar Master Side to Mortar Slave Side
nSides_MortarType=0

DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%nMortars.GT.0)THEN !mortar side
      ! compute index of big mortar in MortarInfo = [1:nMortarSides]
      SideID=aSide%SideID+1-MERGE(firstMortarInnerSide,firstMortarMPISide-nMortarInnerSides,&
                                  aSide%SideID.LE.lastMortarInnerSide)
      MortarType(1,aSide%SideID)=aSide%MortarType
      MortarType(2,aSide%SideID)=SideID
      DO iMortar=1,aSide%nMortars
        mSide=>aSide%MortarSide(iMortar)%sp
        MortarInfo(MI_SIDEID,iMortar,SideID)=mSide%SideID
        MortarInfo(MI_FLIP,iMortar,SideID)=mSide%Flip
      END DO !iMortar
      nSides_MortarType(aSide%MortarType)=nSides_MortarType(aSide%MortarType)+1
    END IF !mortarSide
  END DO ! LocSideID
END DO ! iElem

#if MPI
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nSides_flip,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(nSides_flip, nSides_flip,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE     ,nSides_MortarType,3,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(nSides_MortarType,nSides_MortarType,3,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=0     | ',nSides_flip(0)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=1     | ',nSides_flip(1)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=2     | ',nSides_flip(2)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=3     | ',nSides_flip(3)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=4     | ',nSides_flip(4)
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(*,'(A,A34,I0)')' |','nSides of MortarType=1 | ',nSides_MortarType(1)
SWRITE(*,'(A,A34,I0)')' |','nSides of MortarType=2 | ',nSides_MortarType(2)
SWRITE(*,'(A,A34,I0)')' |','nSides of MortarType=3 | ',nSides_MortarType(3)
SWRITE(UNIT_StdOut,'(132("."))')

LOGWRITE(*,*)'============================= START SIDE CHECKER ==================='
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  LOGWRITE(*,*)'=============== iElem= ',iElem, '==================='
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    LOGWRITE(*,'(5(A,I4))')'globSideID= ',aSide%ind, &
                 ', flip= ',aSide%flip ,&
                 ', SideID= ', aSide%SideID,', nMortars= ',aSide%nMortars,', nbProc= ',aSide%nbProc
    IF(aSide%nMortars.GT.0)THEN !mortar side
      LOGWRITE(*,*)'   --- Mortars ---'
      DO iMortar=1,aSide%nMortars
        LOGWRITE(*,'(I4,4(A,I4))') iMortar,', globSideID= ',aSide%MortarSide(iMortar)%sp%ind, &
                     ', flip= ',aSide%MortarSide(iMortar)%sp%Flip, &
                     ', SideID= ',aSide%MortarSide(iMortar)%sp%SideID, &
                     ', nbProc= ',aSide%MortarSide(iMortar)%sp%nbProc

      END DO !iMortar
    END IF !mortarSide
  END DO ! LocSideID
END DO ! iElem
LOGWRITE(*,*)'============================= END SIDE CHECKER ==================='
END SUBROUTINE fillMeshInfo


#if MPI
!==================================================================================================================================
!> This routine communicates the flip between MPI sides, as the flip determines wheter
!> a side is a master or a slave side. The flip of MINE sides is set to zero, therefore
!> send MINE flip to other processor, so YOUR sides get their corresponding flip>0.
!==================================================================================================================================
SUBROUTINE exchangeFlip()
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY: nElems,offsetElem
USE MOD_Mesh_Vars,ONLY: tElem,tSide,Elems
USE MOD_MPI_Vars
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
INTEGER             :: iElem,LocSideID,iNbProc
INTEGER             :: iMortar,nMortars
INTEGER             :: Flip_MINE(offsetMPISides_MINE(0)+1:offsetMPISides_MINE(nNBProcs))
INTEGER             :: Flip_YOUR(offsetMPISides_YOUR(0)+1:offsetMPISides_YOUR(nNBProcs))
INTEGER             :: SendRequest(nNbProcs),RecRequest(nNbProcs)
!==================================================================================================================================
IF(nProcessors.EQ.1) RETURN
!fill MINE flip info
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(LocSideID)%sp%mortarSide(iMortar)%sp
      IF((aSide%SideID.GT.offsetMPISides_MINE(0)       ).AND.&
         (aSide%SideID.LE.offsetMPISides_MINE(nNBProcs)))THEN
        Flip_MINE(aSide%sideID)=aSide%flip
      END IF
    END DO ! iMortar
  END DO ! LocSideID
END DO ! iElem
DO iNbProc=1,nNbProcs
  ! Start send flip from MINE
  IF(nMPISides_MINE_Proc(iNbProc).GT.0)THEN
    nSendVal    =nMPISides_MINE_Proc(iNbProc)
    SideID_start=OffsetMPISides_MINE(iNbProc-1)+1
    SideID_end  =OffsetMPISides_MINE(iNbProc)
    CALL MPI_ISEND(Flip_MINE(SideID_start:SideID_end),nSendVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,SendRequest(iNbProc),iError)
  END IF
  ! Start receive flip to YOUR
  IF(nMPISides_YOUR_Proc(iNbProc).GT.0)THEN
    nRecVal     =nMPISides_YOUR_Proc(iNbProc)
    SideID_start=OffsetMPISides_YOUR(iNbProc-1)+1
    SideID_end  =OffsetMPISides_YOUR(iNbProc)
    CALL MPI_IRECV(Flip_YOUR(SideID_start:SideID_end),nRecVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,RecRequest(iNbProc),iError)
  END IF
END DO !iProc=1,nNBProcs
DO iNbProc=1,nNbProcs
  IF(nMPISides_YOUR_Proc(iNbProc).GT.0)CALL MPI_WAIT(RecRequest(iNbProc) ,MPIStatus,iError)
  IF(nMPISides_MINE_Proc(iNBProc).GT.0)CALL MPI_WAIT(SendRequest(iNbProc),MPIStatus,iError)
END DO !iProc=1,nNBProcs
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(LocSideID)%sp%mortarSide(iMortar)%sp
      IF(aSide%NbProc.EQ.-1) CYCLE !no MPISide
      IF(aSide%SideID.GT.offsetMPISides_YOUR(0))THEN
        IF(aSide%flip.EQ.0)THEN
          IF(Flip_YOUR(aSide%SideID).EQ.0) STOP 'problem in exchangeflip'
          aSide%flip=Flip_YOUR(aSide%sideID)
        END IF
      ELSE
        aSide%flip=0 !MINE MPISides flip=0
      END IF
    END DO ! iMortar
  END DO ! LocSideID
END DO ! iElem

END SUBROUTINE exchangeFlip
#endif

!==================================================================================================================================
!> Fast recursive sorting algorithm for integer arrays
!==================================================================================================================================
RECURSIVE SUBROUTINE MergeSort(A,nTotal)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal    !< size of array to be sorted
INTEGER,INTENT(INOUT) :: A(nTotal) !< array to be sorted
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: nA,nB,tmp
!==================================================================================================================================
IF(nTotal.LT.2) RETURN
IF(nTotal.EQ.2)THEN
  IF(A(1).GT.A(2))THEN
    tmp  = A(1)
    A(1) = A(2)
    A(2) = tmp
  ENDIF
  RETURN
ENDIF
nA=(nTotal+1)/2
CALL MergeSort(A,nA)
nB=nTotal-nA
CALL MergeSort(A(nA+1:nTotal),nB)
! Performed first on lowest level
IF(A(nA).GT.A(nA+1)) CALL DoMerge(A,nA,nB)
END SUBROUTINE MergeSort


!==================================================================================================================================
!> Merge subarrays (part of mergesort)
!==================================================================================================================================
SUBROUTINE DoMerge(A,nA,nB)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nA        !< number of items in A
INTEGER,INTENT(IN)    :: nB        !< number of items in B
INTEGER,INTENT(INOUT) :: A(nA+nB)  !< subarray to be merged
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k
INTEGER :: part1(nA),part2(nB)
!==================================================================================================================================
part1(1:nA)=A(1:nA)
part2(1:nB)=A(nA+1:nA+nB)
i=1; j=1; k=1;
DO WHILE((i.LE.nA).AND.(j.LE.nB))
  IF(part1(i).LE.part2(j))THEN
    A(k)=part1(i)
    i=i+1
  ELSE
    A(k)=part2(j)
    j=j+1
  ENDIF
  k=k+1
END DO
j=nA-i
A(k:k+nA-i)=part1(i:nA)
END SUBROUTINE DoMerge


!==================================================================================================================================
!> Find the inverse Mapping p.e. NodeID-> entry in NodeMap (a sorted array of unique NodeIDs), using bisection.
!> If Index is not in the range, -1 will be returned. If it is in the range, but is not found, 0 will be returned!!
!==================================================================================================================================
FUNCTION INVMAP(ID,nIDs,ArrID)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)                :: ID            !< ID to search for
INTEGER, INTENT(IN)                :: nIDs          !< size of ArrID
INTEGER, INTENT(IN)                :: ArrID(nIDs)   !< 1D array of IDs
INTEGER                            :: INVMAP        !< index of ID in NodeMap array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,maxSteps,low,up,mid
!==================================================================================================================================
INVMAP=0
maxSteps=INT(LOG(REAL(nIDs))*1.4426950408889634556)+1    !1/LOG(2.)=1.4426950408889634556
low=1
up=nIDs
IF((ID.LT.ArrID(low)).OR.(ID.GT.ArrID(up))) THEN
  !WRITE(*,*)'WARNING, Node Index Not in local range -> set to -1'
  INVMAP=-1  ! not in the range!
  RETURN
END IF
IF(ID.EQ.ArrID(low))THEN
  INVMAP=low
ELSEIF(ID.EQ.ArrID(up))THEN
  INVMAP=up
ELSE
  !bisection
  DO i=1,maxSteps
    mid=(up-low)/2+low
    IF(ID .EQ. ArrID(mid))THEN
      INVMAP=mid                     !index found!
      EXIT
    ELSEIF(ID .GT. ArrID(mid))THEN ! seek in upper half
      low=mid
    ELSE
      up=mid
    END IF
  END DO
END IF
END FUNCTION INVMAP

END MODULE MOD_Prepare_Mesh
