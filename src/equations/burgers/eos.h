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
! Define variables for normal and extended state vector
! Normal   U(1:3)  with conservative variables
! Extended U(1:6) with conservative and primitive variables

#define CONS 1:PP_nVar          /* all cons variables */
#define PRIM 1:PP_nVarPrim      /* all prim variables */

#define PP_2Var PP_nVar+PP_nVarPrim
! Lifting
#define PP_nVarLifting               3
#define LIFT_VARS                    (/1,2,3/)
#define PRIM_LIFT                    (/1,2,3/)
