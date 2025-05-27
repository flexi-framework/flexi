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

!==================================================================================================================================
! Posti related variables
!==================================================================================================================================
MODULE MOD_EOS_Posti_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

INTEGER,PARAMETER :: nVarDepEOS=3

! r u v w
INTEGER,DIMENSION(1:nVarDepEOS,0:nVarDepEOS),PARAMETER :: DepTableEOS = TRANSPOSE(RESHAPE(&
(/&
  0,1,0,0 ,& !1  u
  0,0,1,0 ,& !2  v
  0,0,0,1  & !3  w
/),(/nVarDepEOS+1,nVarDepEOS/)))

CHARACTER(LEN=255),DIMENSION(nVarDepEOS),PARAMETER :: DepNames = &
(/ CHARACTER(LEN=255) :: &
"u"       ,& !1
"v"       ,& !2
"w"        & !3
/)

! Only dummy variables
INTEGER,DIMENSION(1:nVarDepEOS),PARAMETER :: DepSurfaceOnlyEOS = &
(/ 0,0,0 /)
INTEGER,DIMENSION(1:nVarDepEOS),PARAMETER :: DepVolumeOnlyEOS = &
(/ 0,0,0 /)

END MODULE MOD_EOS_Posti_Vars
