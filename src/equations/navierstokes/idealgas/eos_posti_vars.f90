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
#include "eos.h"

!==================================================================================================================================
!==================================================================================================================================
MODULE MOD_EOS_Posti_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

INTEGER,PARAMETER :: nVarTotalEOS=27
! ATTENTION: The first     5 variables must be the conservative ones
!            The following 5 variables must be the primitive ones
!           E
!           n
!           e
!           r
!           g
!           y                    E              V
!           S            V       n              o
!           t            e     E t   T          r
!           a            l     n h   o          t
!           g            o     e a   t          i
!           n            c V   r l   a T        c
!           a            i e   g p   l o        i
!           t         T  t l   y y   T t        t
!           i         e  y o   S S   e a  V V V y     D Q
! G   M M M o V V V   m  M c   t t   m l  o o o M     i C S
! r   o o o n e e e P p  a i   a a   p P  r r r a H   l r c
! a D m m m D l l l r e  g t   g g E e r  t t t g e L a i h
! d e e e e e o o o e r  n y   n n n r e  i i i n l a t t l
! i n n n n n c c c s a  i S   a a t a s  c c c i i m a e i
! e s t t t s i i i s t  t o M t t r t s  i i i t c b t r e
! n i u u u i t t t u u  u u a i i o u u  t t t u i d i i r
! t t m m m t y y y r r  d n c o o p r r  y y y d t a o o e
! s y X Y Z y X Y Z e e  e d h n n y e e  X Y Z e y 2 n n n
INTEGER,DIMENSION(1:nVarTotalEOS,0:nVarTotalEOS),PARAMETER :: DepTableEOS = TRANSPOSE(RESHAPE(&
(/&
  0,1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !1  Density
  0,0,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !2  MomentumX
  0,0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !3  MomentumY
  0,0,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !4  MomentumZ
  0,0,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !5  EnergyStagnationDensity
  0,1,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !6  VelocityX
  0,1,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !7  VelocityY
  0,1,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !8  VelocityZ
  0,1,1,1,1,1,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !9  Pressure
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !10 Temperature
  0,0,0,0,0,0,1,1,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !11 VelocityMagnitude
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !12 VelocitySound
  0,0,0,0,0,0,0,0,0,0,0, 1,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !13 Mach
  0,1,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !14 EnergyStagnation
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,1,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !15 EnthalpyStagnation
  0,1,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !16 Entropy
  0,0,0,0,0,0,0,0,0,0,1, 1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !17 TotalTemperature
  0,1,0,0,0,0,0,0,0,1,0, 1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !18 TotalPressure
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !19 VorticityX
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !20 VorticityY
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !21 VorticityZ
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 1,1,1,0,0,0,0,0,0 ,& !22 VorticityMagnitude
  1,0,0,0,0,0,1,1,1,0,0, 1,0,0,0,0,0,0,0, 1,1,1,1,0,0,0,0,0 ,& !23 Helicity
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !24 Lambda2
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !25 Dilatation
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !26 QCriterion
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0  & !27 Schlieren
/),(/nVarTotalEOS+1,nVarTotalEOS/)))

#if FV_ENABLED && FV_RECONSTRUCT
!           E
!           n
!           e
!           r
!           g
!           y
!           S
!           t
!           a
!           g
!           n
!           a
!           t         T
!           i         e
! G   M M M o V V V   m
! r   o o o n e e e P p
! a D m m m D l l l r e
! d e e e e e o o o e r
! i n n n n n c c c s a
! e s t t t s i i i s t
! n i u u u i t t t u u
! t t m m m t y y y r r
! s y X Y Z y X Y Z e e
INTEGER,DIMENSION(1:PP_nVar,0:nVarTotalEOS),PARAMETER :: DepTablePrimToCons =TRANSPOSE(RESHAPE(&
(/&
  0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !1 Density
  0,1,0,0,0,0,1,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !2 MomentumX
  0,1,0,0,0,0,0,1,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !3 MomentumY
  0,1,0,0,0,0,0,0,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0 ,& !4 MomentumZ
  0,1,0,0,0,0,1,1,1,1,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0  & !5 EnergyStagnationDensity
/),(/nVarTotalEOS+1,5/)))
#endif

CHARACTER(LEN=255),DIMENSION(nVarTotalEOS),PARAMETER :: DepNames = &
(/ CHARACTER(LEN=255) ::    &
"Density"                  ,& !1
"MomentumX"                ,& !2
"MomentumY"                ,& !3
"MomentumZ"                ,& !4
"EnergyStagnationDensity"  ,& !5
"VelocityX"                ,& !6
"VelocityY"                ,& !7
"VelocityZ"                ,& !8
"Pressure"                 ,& !9
"Temperature"              ,& !10
"VelocityMagnitude"        ,& !11
"VelocitySound"            ,& !12
"Mach"                     ,& !13
"EnergyStagnation"         ,& !14
"EnthalpyStagnation"       ,& !15
"Entropy"                  ,& !16
"TotalTemperature"         ,& !17
"TotalPressure"            ,& !18
"VorticityX"               ,& !19
"VorticityY"               ,& !20
"VorticityZ"               ,& !21
"VorticityMagnitude"       ,& !22
"Helicity"                 ,& !23
"Lambda2"                  ,& !24
"Dilatation"               ,& !25
"QCriterion"               ,& !26
"Schlieren"                 & !27
/)

END MODULE MOD_EOS_Posti_Vars
