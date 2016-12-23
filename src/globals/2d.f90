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
!==================================================================================================================================
MODULE MOD_2D
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE to2D_rank4
  MODULE PROCEDURE to2D_rank4
END INTERFACE

INTERFACE to2D_rank5
  MODULE PROCEDURE to2D_rank5
END INTERFACE

INTERFACE to2D_rank6
  MODULE PROCEDURE to2D_rank6
END INTERFACE

PUBLIC::to2D_rank4
PUBLIC::to2D_rank5
PUBLIC::to2D_rank6

CONTAINS

SUBROUTINE to2D_rank4(lbound_in,ubound_in,index3D,array) 
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)             :: lbound_in(4)
INTEGER,INTENT(IN)             :: ubound_in(4)
INTEGER,INTENT(IN)             :: index3D
REAL,INTENT(INOUT),ALLOCATABLE :: array(:,:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: array_loc(:,:,:,:)
INTEGER                        :: ubound_loc(4)
!===================================================================================================================================
ubound_loc = ubound_in
ubound_loc(index3D) = lbound_in(index3D)
ALLOCATE(array_loc(lbound_in(1):ubound_loc(1),&
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4)))
array_loc = array( lbound_in(1):ubound_loc(1),&             
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4))
DEALLOCATE(array)               
ALLOCATE(array    (lbound_in(1):ubound_loc(1),&
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4)))
array = array_loc
DEALLOCATE(array_loc)

END SUBROUTINE to2D_rank4

SUBROUTINE to2D_rank5(lbound_in,ubound_in,index3D,array) 
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)             :: lbound_in(5)
INTEGER,INTENT(IN)             :: ubound_in(5)
INTEGER,INTENT(IN)             :: index3D
REAL,INTENT(INOUT),ALLOCATABLE :: array(:,:,:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: array_loc(:,:,:,:,:)
INTEGER                        :: ubound_loc(5)
!===================================================================================================================================
ubound_loc = ubound_in
ubound_loc(index3D) = lbound_in(index3D)
ALLOCATE(array_loc(lbound_in(1):ubound_loc(1),&
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4),&
                   lbound_in(5):ubound_loc(5)))
array_loc = array( lbound_in(1):ubound_loc(1),&             
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4),&
                   lbound_in(5):ubound_loc(5))
DEALLOCATE(array)               
ALLOCATE(array    (lbound_in(1):ubound_loc(1),&
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4),&
                   lbound_in(5):ubound_loc(5)))
array = array_loc
DEALLOCATE(array_loc)

END SUBROUTINE to2D_rank5

SUBROUTINE to2D_rank6(lbound_in,ubound_in,index3D,array) 
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)             :: lbound_in(6)
INTEGER,INTENT(IN)             :: ubound_in(6)
INTEGER,INTENT(IN)             :: index3D
REAL,INTENT(INOUT),ALLOCATABLE :: array(:,:,:,:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: array_loc(:,:,:,:,:,:)
INTEGER                        :: ubound_loc(6)
!===================================================================================================================================
ubound_loc = ubound_in
ubound_loc(index3D) = lbound_in(index3D)
ALLOCATE(array_loc(lbound_in(1):ubound_loc(1),&
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4),&
                   lbound_in(5):ubound_loc(5),&
                   lbound_in(6):ubound_loc(6)))
array_loc = array( lbound_in(1):ubound_loc(1),&             
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4),&
                   lbound_in(5):ubound_loc(5),&
                   lbound_in(6):ubound_loc(6))
DEALLOCATE(array)               
ALLOCATE(array    (lbound_in(1):ubound_loc(1),&
                   lbound_in(2):ubound_loc(2),&
                   lbound_in(3):ubound_loc(3),&
                   lbound_in(4):ubound_loc(4),&
                   lbound_in(5):ubound_loc(5),&
                   lbound_in(6):ubound_loc(6)))
array = array_loc
DEALLOCATE(array_loc)

END SUBROUTINE to2D_rank6

END MODULE MOD_2D
