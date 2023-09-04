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

INTERFACE SideToVol
  MODULE PROCEDURE SideToVol
END INTERFACE

INTERFACE SideToVol2
  MODULE PROCEDURE SideToVol2
END INTERFACE

INTERFACE VolToSide
  MODULE PROCEDURE VolToSide
END INTERFACE

!INTERFACE ElemToNbElem
!  MODULE PROCEDURE ElemToNbElem
!END INTERFACE

INTERFACE FinalizeMappings
  MODULE PROCEDURE FinalizeMappings
END INTERFACE


PUBLIC::BuildMappings
PUBLIC::Flip_S2M
PUBLIC::Flip_M2S
PUBLIC::SideToVol
PUBLIC::SideToVol2
PUBLIC::VolToSide
!PUBLIC::ElemToNbElem
PUBLIC::FinalizeMappings
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Routine which prebuilds mappings for a specific polynomial degree and allocates and stores them in given mapping arrays.
!==================================================================================================================================
SUBROUTINE buildMappings(Nloc,V2S,S2V,S2V2,FS2M,dim)
! MODULES
USE MOD_Globals,           ONLY:CollectiveStop
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                       :: Nloc              !< Polynomial degree to build mappings on
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: V2S(:,:,:,:,:,:)  !< VolumeToSide mapping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: S2V(:,:,:,:,:,:)  !< SideToVolume mapping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: S2V2(:,:,:,:,:)   !< SideToVolume2 mapping
INTEGER,ALLOCATABLE,INTENT(OUT),OPTIONAL :: FS2M(:,:,:,:)     !< FlipSlaveToMaster mapping
INTEGER,INTENT(IN),OPTIONAL              :: dim               !< dimension (2 or 3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,p,q,l,f,s,ijk(3),pq(3),NlocZ
INTEGER,ALLOCATABLE :: V2S_check(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: S2V_check(:,:,:,:,:,:)
INTEGER,ALLOCATABLE :: S2V2_check(:,:,:,:,:)
LOGICAL             :: correct
INTEGER             :: Flip_lower,Flip_upper,locSide_lower,locSide_upper,dim_loc
!==================================================================================================================================
IF (PRESENT(dim)) THEN
  dim_loc = dim
ELSE
  dim_loc = 3
END IF

IF (dim_loc.EQ.2) THEN
  NlocZ = 0
  Flip_lower = 0
  Flip_upper = 1
  locSide_lower = 2
  locSide_upper = 5
ELSE
  NlocZ = Nloc
  Flip_lower = 0
  Flip_upper = 4
  locSide_lower = 1
  locSide_upper = 6
END IF

! VolToSide
ALLOCATE(V2S_check(3,0:Nloc,0:Nloc,0:NlocZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper)) ! used for sanity check
DO k=0,NlocZ; DO j=0,Nloc; DO i=0,Nloc
  DO f = Flip_lower,Flip_upper
    DO s = locSide_lower,locSide_upper
      V2S_check(:,i,j,k,f,s) = VolToSide(Nloc,i,j,k,f,s,dim_loc)
    END DO
  END DO
END DO; END DO; END DO

IF(PRESENT(V2S))THEN
  SDEALLOCATE(V2S)
  ALLOCATE(V2S(3,0:Nloc,0:Nloc,0:NlocZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper))
  V2S = V2S_check
END IF

! SideToVol
ALLOCATE(S2V_check(3,0:Nloc,0:Nloc,0:NlocZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper)) ! used for sanity check
DO q=0,NlocZ; DO p=0,Nloc; DO l=0,Nloc
  DO f = Flip_lower,Flip_upper
    DO s = locSide_lower,locSide_upper
      S2V_check(:,l,p,q,f,s) = SideToVol(Nloc,l,p,q,f,s, dim_loc)
    END DO
  END DO
END DO; END DO; END DO

IF(PRESENT(S2V))THEN
  SDEALLOCATE(S2V)
  ALLOCATE(S2V(3,0:Nloc,0:Nloc,0:NlocZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper))
  S2V = S2V_check
END IF

! SideToVol2
ALLOCATE(S2V2_check(2,0:Nloc,0:NlocZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper)) ! used for sanity check
DO j=0,NlocZ; DO i=0,Nloc
  DO f = Flip_lower,Flip_upper
    DO s = locSide_lower,locSide_upper
      S2V2_check(:,i,j,f,s) = SideToVol2(Nloc,i,j,f,s,dim_loc)
    END DO
  END DO
END DO; END DO

IF(PRESENT(S2V2))THEN
  SDEALLOCATE(S2V2)
  ALLOCATE(S2V2(2,0:Nloc,0:NlocZ,Flip_lower:Flip_upper,locSide_lower:locSide_upper))
  S2V2 = S2V2_check
END IF

! Flip_S2M
IF(PRESENT(FS2M))THEN
  SDEALLOCATE(FS2M)
  ALLOCATE(FS2M(2,0:Nloc,0:NlocZ,Flip_lower:Flip_upper))
  DO j=0,NlocZ; DO i=0,Nloc
    DO f = Flip_lower,Flip_upper
      FS2M(:,i,j,f) = Flip_S2M(Nloc,i,j,f,dim_loc)
    END DO
  END DO; END DO
END IF

! Perform sanity checks
DO f = Flip_lower,Flip_upper
  DO s = locSide_lower,locSide_upper
    DO q = 0,NlocZ; DO p = 0,Nloc
      ijk = S2V_check(:,0,p,q,f,s)
      pq = V2S_check(:,ijk(1),ijk(2),ijk(3),f,s)
      IF ((pq(1).NE.p).OR.(pq(2).NE.q)) THEN
        CALL CollectiveStop(__STAMP__,&
          'SideToVol does not fit to VolToSide')
      END IF
    END DO; END DO
  END DO ! s = locSide_lower,locSide_upper
END DO ! f = Flip_lower,Flip_upper

DO f = Flip_lower,Flip_upper
  DO s = locSide_lower,locSide_upper
    DO k=0,NlocZ; DO j=0,Nloc; DO i=0,Nloc
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
  END DO ! s = locSide_lower,locSide_upper
END DO ! f = Flip_lower,Flip_upper

! Deallocate arrays used for sanity check
DEALLOCATE(V2S_check,S2V_check,S2V2_check)
END SUBROUTINE BuildMappings


!==================================================================================================================================
!> Transforms Coordinates from RHS of Slave to RHS of Master
!>    input: p,q in Slave-RHS, flip;
!>   output: indices in Master-RHS
!==================================================================================================================================
FUNCTION Flip_S2M(Nloc, p, q, flip, dim)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: p,q,flip,Nloc,dim
INTEGER,DIMENSION(2) :: Flip_S2M
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF (dim.EQ.3) THEN
  SELECT CASE(flip)
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
  END SELECT
ELSE
  SELECT CASE(flip)
  CASE(0)
    Flip_S2M = (/     p,     0/)
  CASE(1)
    Flip_S2M = (/Nloc-p,     0/)
  END SELECT
END IF
END FUNCTION Flip_S2M


!==================================================================================================================================
!> Transforms Coordinates from RHS of Master to RHS of Slave
!>   actualy this is the same function as Flip_S2M
!==================================================================================================================================
FUNCTION Flip_M2S(Nloc, p, q, flip, dim)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: p,q,flip,Nloc,dim
INTEGER,DIMENSION(2) :: Flip_M2S
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Flip_M2S=Flip_S2M(Nloc,p,q,flip,dim)
END FUNCTION Flip_M2S


!==================================================================================================================================
!> Transforms Volume-Coordinates into RHS of the Side (uses CGNS-Notation for side orientation)
!> input: i,j,k, locSideID
!>   where: i,j,k = volume-indices
!> output: indices in Master-RHS  +  volume-index which is not used (depending on locSideID)
!==================================================================================================================================
FUNCTION CGNS_VolToSide(Nloc, i, j, k, locSideID, dim)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: i,j,k,locSideID,Nloc,dim
INTEGER,DIMENSION(3) :: CGNS_VolToSide
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF (dim.EQ.3) THEN
  SELECT CASE(locSideID)
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
  END SELECT
ELSE
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    CGNS_VolToSide = (/Nloc-j,i,0/)
  CASE(XI_PLUS)
    CGNS_VolToSide = (/j,Nloc-i,0/)
  CASE(ETA_MINUS)
    CGNS_VolToSide = (/i,j,0/)
  CASE(ETA_PLUS)
    CGNS_VolToSide = (/Nloc-i,Nloc-j,0/)
  END SELECT
END IF
END FUNCTION CGNS_VolToSide


!==================================================================================================================================
!> Transforms RHS-Coordinates of Side (CGNS-Notation for side orientation) into Volume-Coordinates
!> input: l, p,q, locSideID
!>   where: p,q are in Master-RHS;
!>          l is the xi-,eta- or zeta-index in 0:Nloc corresponding to locSideID
!> output: volume-indices
!==================================================================================================================================
FUNCTION CGNS_SideToVol(Nloc, l, p, q, locSideID, dim)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: l,p,q,locSideID,Nloc,dim
INTEGER,DIMENSION(3) :: CGNS_SideToVol
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF (dim.EQ.3) THEN
  SELECT CASE(locSideID)
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
  END SELECT
ELSE
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    CGNS_SideToVol = (/l,Nloc-p,0/)
  CASE(XI_PLUS)
    CGNS_SideToVol = (/Nloc-l,p,0/)
  CASE(ETA_MINUS)
    CGNS_SideToVol = (/p,l,0/)
  CASE(ETA_PLUS)
    CGNS_SideToVol = (/Nloc-p,Nloc-l,0/)
  END SELECT
END IF
END FUNCTION CGNS_SideToVol


!==================================================================================================================================
!> Transforms RHS-Coordinates of Side (CGNS-Notation for side orientation) into side-local tensor product Volume-Coordinates
!> input: p,q, locSideID
!>   where: p,q are in Master-RHS;
!> output: Surface coordinates in volume frame
!==================================================================================================================================
FUNCTION CGNS_SideToVol2(Nloc, p, q, locSideID, dim)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: p,q,locSideID,Nloc,dim
INTEGER,DIMENSION(2) :: CGNS_SideToVol2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF (dim.EQ.3) THEN
  SELECT CASE(locSideID)
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
  END SELECT
ELSE
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    CGNS_SideToVol2 = (/Nloc-p,0/)
  CASE(XI_PLUS)
    CGNS_SideToVol2 = (/p,0/)
  CASE(ETA_MINUS)
    CGNS_SideToVol2 = (/p,0/)
  CASE(ETA_PLUS)
    CGNS_SideToVol2 = (/Nloc-p,0/)
  END SELECT
END IF
END FUNCTION CGNS_SideToVol2


!==================================================================================================================================
!> Transform Volume-Coordinates to RHS-Coordinates of Master. This is: VolToSide = Flip_S2M(CGNS_VolToSide(...))
!> input: i,j,k, flip, locSideID
!>   where: i,j,k = volume-indices
!> output: indices in Master-RHS
!==================================================================================================================================
FUNCTION VolToSide(Nloc, i, j, k, flip, locSideID, dim)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: i,j,k,flip,locSideID,Nloc,dim
INTEGER,DIMENSION(3) :: VolToSide
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(3) :: pq
!==================================================================================================================================
pq = CGNS_VolToSide(Nloc,i,j,k,locSideID, dim)
VolToSide(1:2) = Flip_S2M(Nloc,pq(1),pq(2),flip, dim)
VolToSide(3) = pq(3)
END FUNCTION VolToSide

!==================================================================================================================================
!> Transform RHS-Coordinates of Master to Volume-Coordinates. This is: SideToVol = CGNS_SideToVol(Flip_M2S(...))
!> input: l, p,q, flip, locSideID
!>     where: p,q are in Master-RHS;
!>            l is the xi-,eta- or zeta-index in 0:Nloc corresponding to locSideID
!> output: volume-indices
!==================================================================================================================================
FUNCTION SideToVol(Nloc, l, p, q, flip, locSideID,dim)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: l,p,q,flip,locSideID,Nloc,dim
INTEGER,DIMENSION(3) :: SideToVol
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!==================================================================================================================================
pq = Flip_M2S(Nloc,p,q,flip,dim)
SideToVol = CGNS_SideToVol(Nloc,l,pq(1),pq(2),locSideID,dim)
END FUNCTION SideToVol

!==================================================================================================================================
!> Transform RHS-Coordinates of Master to Volume-Coordinates. This is: SideToVol2 = CGNS_SideToVol2(Flip_M2S(...))
!> input:  p,q, flip, locSideID
!>     where: p,q are in Master-RHS;
!> output: volume-indicies (two indices on the specific local side)
!==================================================================================================================================
FUNCTION SideToVol2(Nloc, p, q, flip, locSideID,dim)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: p,q,flip,locSideID,Nloc,dim
INTEGER,DIMENSION(2) :: SideToVol2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!==================================================================================================================================
pq = Flip_M2S(Nloc,p,q,flip,dim)
SideToVol2 = CGNS_SideToVol2(Nloc,pq(1),pq(2),locSideID,dim)
END FUNCTION SideToVol2

!!!==================================================================================================================================
!!> Get the index of neighbor element, return -1 if none exists
!!==================================================================================================================================
!FUNCTION ElemToNbElem(locSideID,iElem)
!! MODULES
!USE MOD_Mesh_Vars,ONLY:ElemToSide
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!INTEGER,INTENT(IN)   :: locSideID,iElem
!INTEGER              :: ElemToNBElem
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER              :: SideID
!!==================================================================================================================================
!SideID=ElemToSide(E2S_SIDE_ID,locSideID,iElem)
!END FUNCTION ElemToNbElem

SUBROUTINE FinalizeMappings()
! MODULES
USE MOD_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!===================================================================================================================================
SDEALLOCATE(FS2M)
SDEALLOCATE(V2S)
SDEALLOCATE(S2V)
SDEALLOCATE(S2V2)
SDEALLOCATE(FS2M)
END SUBROUTINE FinalizeMappings

END MODULE MOD_Mappings
