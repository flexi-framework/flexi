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
!> Contains subroutines to build mappings for easier handling of 3D data-structures and their connectivity
!==================================================================================================================================
MODULE MOD_Mappings
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE BuildMappings
  MODULE PROCEDURE BuildMappings
END INTERFACE

INTERFACE Flip_S2M
  MODULE PROCEDURE Flip_S2M
END INTERFACE

INTERFACE Flip_M2S
  MODULE PROCEDURE Flip_M2S
END INTERFACE

INTERFACE CGNS_SideToVol
  MODULE PROCEDURE CGNS_SideToVol
END INTERFACE

INTERFACE CGNS_SideToVol2
  MODULE PROCEDURE CGNS_SideToVol2
END INTERFACE

INTERFACE CGNS_VolToSide
  MODULE PROCEDURE CGNS_VolToSide
END INTERFACE

INTERFACE SideToVol
  MODULE PROCEDURE SideToVol
END INTERFACE

INTERFACE SideToVol2
  MODULE PROCEDURE SideToVol2
END INTERFACE

INTERFACE VolToSide
  MODULE PROCEDURE VolToSide
END INTERFACE

INTERFACE VolToSide2
  MODULE PROCEDURE VolToSide2
END INTERFACE

INTERFACE VolToVol
  MODULE PROCEDURE VolToVol
END INTERFACE

INTERFACE ElemToNbElem
  MODULE PROCEDURE ElemToNbElem
END INTERFACE

PUBLIC::BuildMappings
PUBLIC::Flip_S2M
PUBLIC::Flip_M2S
PUBLIC::CGNS_SideToVol
PUBLIC::CGNS_SideToVol2
PUBLIC::CGNS_VolToSide
PUBLIC::SideToVol
PUBLIC::SideToVol2
PUBLIC::VolToSide
PUBLIC::VolToSide2
PUBLIC::VolToVol
PUBLIC::ElemToNbElem
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Routine which prebuilds mappings for a specific polynomial degree and allocates and stores them in given mapping arrays.
!==================================================================================================================================
SUBROUTINE buildMappings(Nloc,V2S,V2S2,CV2S,S2V,S2V2,CS2V2,FS2M)
! MODULES
USE MOD_Globals,           ONLY:CollectiveStop
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                       :: Nloc              !< Polynomial degree to build mappings on
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: V2S(:,:,:,:,:,:)  !< VolumeToSide mapping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: V2S2(:,:,:,:,:)   !< VolumeToSide2 mapping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: S2V(:,:,:,:,:,:)  !< SideToVolume mapping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: CV2S(:,:,:,:,:)   !< CGNS VolumeToSide mappping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: S2V2(:,:,:,:,:)   !< SideToVolume2 mapping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: CS2V2(:,:,:,:)    !< CGNS SideToVolume2 mapping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: FS2M(:,:,:,:)     !< FlipSlaveToMaster mapping
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,p,q,f,s,ijk(3),pq(3)
INTEGER,ALLOCATABLE :: V2S_check(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: S2V_check(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: S2V2_check(:,:,:,:,:)
LOGICAL             :: correct
!==================================================================================================================================

! VolToSide
ALLOCATE(V2S_check(3,0:Nloc,0:Nloc,0:Nloc,0:4,1:6)) ! used for sanity check
DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
  DO f=0,4
    DO s=1,6
      V2S_check(:,i,j,k,f,s) = VolToSide(Nloc,i,j,k,f,s)
    END DO
  END DO
END DO; END DO; END DO

IF(PRESENT(V2S))THEN
  ALLOCATE(V2S(3,0:Nloc,0:Nloc,0:Nloc,0:4,1:6))
  V2S = V2S_check
END IF

! VolToSide2
IF(PRESENT(V2S2))THEN
  ALLOCATE(V2S2(2,0:Nloc,0:Nloc,0:4,1:6))
  DO j=0,Nloc; DO i=0,Nloc
    DO f=0,4
      DO s=1,6
        V2S2(:,i,j,f,s) = VolToSide2(Nloc,i,j,f,s)
      END DO
    END DO
  END DO; END DO
END IF

! CGNS_VolToSide
IF(PRESENT(CV2S))THEN
  ALLOCATE(CV2S(3,0:Nloc,0:Nloc,0:Nloc,1:6))
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    DO s=1,6
      CV2S(:,i,j,k,s) = CGNS_VolToSide(Nloc,i,j,k,s)
    END DO
  END DO; END DO; END DO
END IF

! SideToVol
ALLOCATE(S2V_check(3,0:Nloc,0:Nloc,0:Nloc,0:4,1:6)) ! used for sanity check
DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
  DO f=0,4
    DO s=1,6
      S2V_check(:,i,j,k,f,s) = SideToVol(Nloc,i,j,k,f,s)
    END DO
  END DO
END DO; END DO; END DO

IF(PRESENT(S2V))THEN
  ALLOCATE(S2V(3,0:Nloc,0:Nloc,0:Nloc,0:4,1:6))
  S2V = S2V_check
END IF

! SideToVol2
ALLOCATE(S2V2_check(2,0:Nloc,0:Nloc,0:4,1:6)) ! used for sanity check
DO j=0,Nloc; DO i=0,Nloc
  DO f=0,4
    DO s=1,6
      S2V2_check(:,i,j,f,s) = SideToVol2(Nloc,i,j,f,s)
    END DO
  END DO
END DO; END DO

IF(PRESENT(S2V2))THEN
  ALLOCATE(S2V2(2,0:Nloc,0:Nloc,0:4,1:6))
  S2V2 = S2V2_check
END IF

! CGNS_SideToVol2
IF(PRESENT(CS2V2))THEN
  ALLOCATE(CS2V2(2,0:Nloc,0:Nloc,1:6))
  DO j=0,Nloc; DO i=0,Nloc
    DO s=1,6
      CS2V2(:,i,j,s) = CGNS_SideToVol2(Nloc,i,j,s)
    END DO
  END DO; END DO
END IF

! Flip_S2M
IF(PRESENT(FS2M))THEN
  ALLOCATE(FS2M(2,0:Nloc,0:Nloc,0:4))
  DO j=0,Nloc; DO i=0,Nloc
    DO f=0,4
      FS2M(:,i,j,f) = Flip_S2M(Nloc,i,j,f)
    END DO
  END DO; END DO
END IF

! Perform sanity checks
#if PP_dim == 3
DO f = 0, 4
  DO s = 1, 6
#else    
DO f = 0, 1
  DO s = 2, 5
#endif    
    DO q = 0,PP_NlocZ; DO p = 0,Nloc
      ijk = S2V_check(:,0,p,q,f,s)
      pq = V2S_check(:,ijk(1),ijk(2),ijk(3),f,s)
      IF ((pq(1).NE.p).OR.(pq(2).NE.q)) THEN
        CALL CollectiveStop(__STAMP__,&
          'SideToVol does not fit to VolToSide')
      END IF
    END DO; END DO
  END DO ! s = 1, 6
END DO ! f = 0, 4

#if PP_dim == 3
DO f = 0, 4
  DO s = 1, 6
#else    
DO f = 0, 1
  DO s = 2, 5
#endif    
    DO k=0,PP_NlocZ; DO j=0,Nloc; DO i=0,Nloc
      pq = V2S_check(:,i,j,k,f,s)
      ijk(1:2) = S2V2_check(:,pq(1),pq(2),f,s)
      correct=.TRUE.
      SELECT CASE(s)
      CASE(XI_MINUS,XI_PLUS)
        IF ((ijk(1).NE.j).OR.(ijk(2).NE.k)) correct=.FALSE.
      CASE(ETA_MINUS,ETA_PLUS)
        IF ((ijk(1).NE.i).OR.(ijk(2).NE.k)) correct=.FALSE.
      CASE(ZETA_MINUS,ZETA_PLUS)
        IF ((ijk(1).NE.i).OR.(ijk(2).NE.j)) correct=.FALSE.
      END SELECT
      IF(.NOT.correct)THEN
        CALL CollectiveStop(__STAMP__,&
          'SideToVol2 does not fit to VolToSide')
      END IF
    END DO; END DO; END DO! i,j,k=0,Nloc
  END DO ! s = 1, 6
END DO ! f = 0, 4

! Deallocate arrays used for sanity check
DEALLOCATE(V2S_check,S2V_check,S2V2_check)
END SUBROUTINE BuildMappings


!==================================================================================================================================
!> Transforms Coordinates from RHS of Slave to RHS of Master
!>    input: p,q in Slave-RHS, flip;
!>   output: indices in Master-RHS
!==================================================================================================================================
FUNCTION Flip_S2M(Nloc, p, q, flip)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: p,q,flip,Nloc
INTEGER,DIMENSION(2) :: Flip_S2M
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SELECT CASE(flip)
#if PP_dim == 3
  CASE(0)
    Flip_S2M = (/     p,     q/)
  CASE(1)
    Flip_S2M = (/     q,     p/)
  CASE(2)
    Flip_S2M = (/Nloc-p,     q/)
  CASE(3)
    Flip_S2M = (/Nloc-q,Nloc-p/)
  CASE(4)
    Flip_S2M = (/     p,Nloc-q/)
#else
  CASE(0)
    Flip_S2M = (/     p,     0/)
  CASE(1)
    Flip_S2M = (/Nloc-p,     0/)
#endif    
END SELECT
END FUNCTION Flip_S2M


!==================================================================================================================================
!> Transforms Coordinates from RHS of Master to RHS of Slave
!>   actualy this is the same function as Flip_S2M
!==================================================================================================================================
FUNCTION Flip_M2S(Nloc, p, q, flip)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: p,q,flip,Nloc
INTEGER,DIMENSION(2) :: Flip_M2S
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Flip_M2S=Flip_S2M(Nloc,p,q,flip)
END FUNCTION Flip_M2S


!==================================================================================================================================
!> Transforms Volume-Coordinates into RHS of the Side (uses CGNS-Notation for side orientation)
!> input: i,j,k, locSideID
!>   where: i,j,k = volume-indices
!> output: indices in Master-RHS  +  volume-index which is not used (depending on locSideID)
!==================================================================================================================================
FUNCTION CGNS_VolToSide(Nloc, i, j, k, locSideID)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: i,j,k,locSideID,Nloc
INTEGER,DIMENSION(3) :: CGNS_VolToSide
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SELECT CASE(locSideID)
#if PP_dim == 3
  CASE(XI_MINUS)
    CGNS_VolToSide = (/k,j,i/)
  CASE(XI_PLUS)
    CGNS_VolToSide = (/j,k,Nloc-i/)
  CASE(ETA_MINUS)
    CGNS_VolToSide = (/i,k,j/)
  CASE(ETA_PLUS)
    CGNS_VolToSide = (/Nloc-i,k,Nloc-j/)
  CASE(ZETA_MINUS)
    CGNS_VolToSide = (/j,i,k/)
  CASE(ZETA_PLUS)
    CGNS_VolToSide = (/i,j,Nloc-k/)
#else
  CASE(XI_MINUS)
    CGNS_VolToSide = (/Nloc-j,i,0/)
  CASE(XI_PLUS)
    CGNS_VolToSide = (/j,Nloc-i,0/)
  CASE(ETA_MINUS)
    CGNS_VolToSide = (/i,j,0/)
  CASE(ETA_PLUS)
    CGNS_VolToSide = (/Nloc-i,Nloc-j,0/)
#endif  
END SELECT
END FUNCTION CGNS_VolToSide


!==================================================================================================================================
!> Transforms RHS-Coordinates of Side (CGNS-Notation for side orientation) into Volume-Coordinates
!> input: l, p,q, locSideID
!>   where: p,q are in Master-RHS;
!>          l is the xi-,eta- or zeta-index in 0:Nloc corresponding to locSideID
!> output: volume-indices
!==================================================================================================================================
FUNCTION CGNS_SideToVol(Nloc, l, p, q, locSideID)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: l,p,q,locSideID,Nloc
INTEGER,DIMENSION(3) :: CGNS_SideToVol
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SELECT CASE(locSideID)
#if PP_dim == 3
  CASE(XI_MINUS)
    CGNS_SideToVol = (/l,q,p/)
  CASE(XI_PLUS)
    CGNS_SideToVol = (/Nloc-l,p,q/)
  CASE(ETA_MINUS)
    CGNS_SideToVol = (/p,l,q/)
  CASE(ETA_PLUS)
    CGNS_SideToVol = (/Nloc-p,Nloc-l,q/)
  CASE(ZETA_MINUS)
    CGNS_SideToVol = (/q,p,l/)
  CASE(ZETA_PLUS)
    CGNS_SideToVol = (/p,q,Nloc-l/)
#else
  CASE(XI_MINUS)
    CGNS_SideToVol = (/l,Nloc-p,0/)
  CASE(XI_PLUS)
    CGNS_SideToVol = (/Nloc-l,p,0/)
  CASE(ETA_MINUS)
    CGNS_SideToVol = (/p,l,0/)
  CASE(ETA_PLUS)
    CGNS_SideToVol = (/Nloc-p,Nloc-l,0/)
#endif    
END SELECT
END FUNCTION CGNS_SideToVol


!==================================================================================================================================
!> Transforms RHS-Coordinates of Side (CGNS-Notation for side orientation) into side-local tensor product Volume-Coordinates
!> input: l, p,q, locSideID
!>   where: p,q are in Master-RHS;
!>          l is the xi-,eta- or zeta-index in 0:Nloc corresponding to locSideID
!> output: Surface coordinates in volume frame
!==================================================================================================================================
FUNCTION CGNS_SideToVol2(Nloc, p, q, locSideID)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: p,q,locSideID,Nloc
INTEGER,DIMENSION(2) :: CGNS_SideToVol2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SELECT CASE(locSideID)
#if PP_dim == 3
  CASE(XI_MINUS)
    CGNS_SideToVol2 = (/q,p/)
  CASE(XI_PLUS)
    CGNS_SideToVol2 = (/p,q/)
  CASE(ETA_MINUS)
    CGNS_SideToVol2 = (/p,q/)
  CASE(ETA_PLUS)
    CGNS_SideToVol2 = (/Nloc-p,q/)
  CASE(ZETA_MINUS)
    CGNS_SideToVol2 = (/q,p/)
  CASE(ZETA_PLUS)
    CGNS_SideToVol2 = (/p,q/)
#else
  CASE(XI_MINUS)
    CGNS_SideToVol2 = (/Nloc-p,0/)
  CASE(XI_PLUS)
    CGNS_SideToVol2 = (/p,0/)
  CASE(ETA_MINUS)
    CGNS_SideToVol2 = (/p,0/)
  CASE(ETA_PLUS)
    CGNS_SideToVol2 = (/Nloc-p,0/)
#endif    
END SELECT
END FUNCTION CGNS_SideToVol2


!==================================================================================================================================
!> Transform Volume-Coordinates to RHS-Coordinates of Master. This is: VolToSide = Flip_S2M(CGNS_VolToSide(...))
!> input: i,j,k, flip, locSideID
!>   where: i,j,k = volume-indices
!> output: indices in Master-RHS
!==================================================================================================================================
FUNCTION VolToSide(Nloc, i, j, k, flip, locSideID)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: i,j,k,flip,locSideID,Nloc
INTEGER,DIMENSION(3) :: VolToSide
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(3) :: pq
!==================================================================================================================================
pq = CGNS_VolToSide(Nloc,i,j,k,locSideID)
VolToSide(1:2) = Flip_S2M(Nloc,pq(1),pq(2),flip)
VolToSide(3) = pq(3)
END FUNCTION VolToSide


!==================================================================================================================================
!> Transform Surface Coordinates in Volume frame to RHS-Coordinates of Slave.
!> This is: VolToSide2 = Flip_M2S(CGNS_VolToSide2(...))
!> input: i,j, flip, locSideID
!>   where: i,j = volume-indices
!> output: indices in Master-RHS
!==================================================================================================================================
FUNCTION VolToSide2(Nloc, i, j, flip, locSideID)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: i,j,flip,locSideID,Nloc
INTEGER,DIMENSION(2) :: VolToSide2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!==================================================================================================================================
pq = Flip_M2S(Nloc,i,j,flip)
VolToSide2 = CGNS_SideToVol2(Nloc,pq(1),pq(2),locSideID)
END FUNCTION VolToSide2


!==================================================================================================================================
!> Transform RHS-Coordinates of Master to Volume-Coordinates. This is: SideToVol = CGNS_SideToVol(Flip_M2S(...))
!> input: l, p,q, flip, locSideID
!>     where: p,q are in Master-RHS;
!>            l is the xi-,eta- or zeta-index in 0:Nloc corresponding to locSideID
!> output: volume-indices
!==================================================================================================================================
FUNCTION SideToVol(Nloc, l, p, q, flip, locSideID)
! MODULES
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)   :: l,p,q,flip,locSideID,Nloc
INTEGER,DIMENSION(3) :: SideToVol
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!==================================================================================================================================
pq = Flip_M2S(Nloc,p,q,flip)
SideToVol = CGNS_SideToVol(Nloc,l,pq(1),pq(2),locSideID)
END FUNCTION SideToVol

!==================================================================================================================================
!> Transform RHS-Coordinates of Master to Volume-Coordinates. This is: SideToVol2 = CGNS_SideToVol2(Flip_M2S(...))
!> input:  p,q, flip, locSideID
!>     where: p,q are in Master-RHS;
!> output: volume-indicies
!==================================================================================================================================
FUNCTION SideToVol2(Nloc, p, q, flip, locSideID)
! MODULES
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)   :: p,q,flip,locSideID,Nloc
INTEGER,DIMENSION(2) :: SideToVol2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!==================================================================================================================================
pq = Flip_M2S(Nloc,p,q,flip)
SideToVol2 = CGNS_SideToVol2(Nloc,pq(1),pq(2),locSideID)
END FUNCTION SideToVol2

!==================================================================================================================================
!> Get the index of neighbor element, return -1 if none exists
!==================================================================================================================================
FUNCTION ElemToNbElem(locSideID,iElem)
! MODULES
USE MOD_Mesh_Vars,ONLY:ElemToSide,SideToElem,firstInnerSide,lastInnerSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: locSideID,iElem
INTEGER              :: ElemToNBElem
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: flip, SideID
!==================================================================================================================================
SideID=ElemToSide(E2S_SIDE_ID,locSideID,iElem)
IF ((SideID.GE.firstInnerSide).AND.(SideID.LE.lastInnerSide)) THEN
  flip  =ElemToSide(E2S_FLIP,locSideID,iElem)
  IF (flip.EQ.0) THEN
    ElemToNbElem = SideToElem(S2E_NB_ELEM_ID,SideID)
  ELSE
    ElemToNbElem = SideToElem(S2E_ELEM_ID,SideID)
  END IF
ELSE
  ElemToNbElem = -1
END IF
END FUNCTION ElemToNbElem


!==================================================================================================================================
!> Transform Volume-Coordinates to neighboring Volume-Coordinates.  This is: VolToVol = SideToVol(VolToSide(...))
!> input: i,j,k, iElem, locSideID
!>     where: i,j,k  are Volume-Indizices of element iElem;
!>            locSideID  side to the neighboring element
!> output: volume-indices of neighboring element, that are next to ijk in direction of locSideID
!==================================================================================================================================
FUNCTION VolToVol(Nloc,i,j,k,locSideID,iElem)
! MODULES
USE MOD_Mesh_Vars,ONLY:ElemToSide,SideToElem
IMPLICIT NONE
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)   :: i,j,k,locSideID,iElem,Nloc
INTEGER,DIMENSION(3) :: VolToVol
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(3) :: pq
INTEGER              :: flip, l, SideID, neighbor_flip, neighbor_locSideID
!==================================================================================================================================
SideID=ElemToSide(E2S_SIDE_ID,locSideID,iElem)
flip  =ElemToSide(E2S_FLIP,locSideID,iElem)
pq    =VolToSide(Nloc, i,j,k, flip, locSideID)
l = pq(3)
IF (flip.EQ.0) THEN
  neighbor_locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  neighbor_flip      = SideToElem(S2E_FLIP,SideID)
ELSE
  neighbor_locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  neighbor_flip      = 0
END IF
VolToVol = SideToVol(Nloc, l, pq(1), pq(2), neighbor_flip, neighbor_locSideID)
END FUNCTION VolToVol

END MODULE MOD_Mappings
