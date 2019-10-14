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
!> Module containing the routines needed for the inexact sparse ILU(0) inversion of the preconditioner.
!===================================================================================================================================
MODULE MOD_SparseILU
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitSparseILU
  MODULE PROCEDURE InitSparseILU
END INTERFACE

INTERFACE FinalizeSparseILU
  MODULE PROCEDURE FinalizeSparseILU
END INTERFACE

INTERFACE BuildILU0
  MODULE PROCEDURE BuildILU0
END INTERFACE

INTERFACE ApplyILU
  MODULE PROCEDURE ApplyILU
END INTERFACE

PUBLIC:: InitSparseILU,FinalizeSparseILU,BuildILU0,ApplyILU
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Init Sparse ILU
!===================================================================================================================================
SUBROUTINE InitSparseILU()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_SparseILU_Vars
USE MOD_Mesh_Vars              ,ONLY:nElems
USE MOD_Implicit_Vars          ,ONLY:nDOFVarElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ILU0...'
nMTriangle=nDOFVarElem-1
ALLOCATE(Dinv(nDOFVarElem,nElems))
ALLOCATE(nUNonZeros(nElems) &
        ,nLNonZeros(nElems) )

ALLOCATE(IU(nElems)         &
        ,IL(nElems)         )

DO iElem=1,nElems
  ALLOCATE(IU(iElem)%IEntry(nDOFVarElem))
  ALLOCATE(IL(iElem)%IEntry(nDOFVarElem))
END DO 

! machine accuracy
epsZero=EPSILON(0.0d0)

SWRITE(UNIT_stdOut,'(A)')' INIT ILU0 DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitSparseILU

!===================================================================================================================================
!> Build the ILU0 per Block in the csr format
!===================================================================================================================================
SUBROUTINE BuildILU0(ILU0,iElem)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SparseILU_Vars
USE MOD_Implicit_Vars              ,ONLY:nDOFVarElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)          :: ILU0(nDOFVarElem,nDOFVarElem)
INTEGER,INTENT(IN)          :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: ii,kk,jj,iEntry
LOGICAL                     :: first
REAL                        :: Sparsity
INTEGER                     :: lastLine,lineNonZero
LOGICAL                     :: singleValue
!===================================================================================================================================
!----------------------------------------------------------------------------------------------------------------------------------
! ILU(0) // drop all zero values // Saad 'Iterative Methods for sparse linear systems' p.307 Algorithm 10.4 // p.327 in pdf
!----------------------------------------------------------------------------------------------------------------------------------
! Saving of L and U within ILU0, L unit lower triangular, U upper triangular
!----------------------------------------------------------------------------------------------------------------------------------
! Saad

DO ii=2,nDOFVarElem
  DO kk=1,ii-1
    IF(ABS(ILU0(ii,kk)).GT.epsZero)THEN
      ILU0(ii,kk) = ILU0(ii,kk)/ILU0(kk,kk)
      DO jj=kk+1,nDOFVarElem
        IF(ABS(ILU0(ii,jj)).GT.epsZero)THEN
          ILU0(ii,jj) = ILU0(ii,jj) - ILU0(ii,kk)*ILU0(kk,jj)
        END IF ! ii,jj element NZ
      END DO ! jj
    END IF ! ii,kk element NZ
  END DO ! ii
END DO ! kk

! Diagonal entries
! It is not possible to store them before the ILU decomposition because they will be changed
DO ii=1,nDOFVarElem
  Dinv(ii,iElem)=1./ILU0(ii,ii)
END DO

!----------------------------------------------------------------------------------------------------------------------------------
! Extended CSR FORMAT
!----------------------------------------------------------------------------------------------------------------------------------

! get number of non-zero entries
nUNonZeros(iElem)=0
nLNonZeros(iElem)=0

DO ii=1,nDOFVarElem
  DO kk=ii+1,nDOFVarElem
    ! upper
    IF(ABS(ILU0(ii,kk)).GT.epsZero)THEN
      nUNonZeros(iElem)=nUNonZeros(iElem)+1
    END IF
  END DO
  ! lower
  DO kk=1,ii-1
    IF(ABS(ILU0(ii,kk)).GT.epsZero)THEN
      nLNonZeros(iElem)=nLNonZeros(iElem)+1
    END IF
  END DO
END DO ! ii


Sparsity= REAL(nUNonZeros(iElem))+REAL(nLNonZeros(iElem)) + REAL(nDOFVarElem)
Sparsity=Sparsity/REAL(nDOFVarElem)/REAL(nDOFVarElem)

SDEALLOCATE(IU(iElem)%Entry)
SDEALLOCATE(IU(iElem)%JEntry)
SDEALLOCATE(IL(iElem)%Entry)
SDEALLOCATE(IL(iElem)%JEntry)
ALLOCATE( IU(iElem)%Entry(nUNonZeros(iElem))  &
        , IU(iELEM)%JEntry(nUNonZeros(iElem)) &
        , IL(iELEM)%Entry(nLNonZeros(iElem))  &
        , IL(iELEM)%JEntry(nLNonZeros(iElem)))

! nullify     
IL(iELem)%Entry=0.
IU(iELem)%Entry=0.

! U part of ILU0
! simple version
jj=0
singleValue=.FALSE.
DO ii=1,nMTriangle,1
  first=.TRUE.
  lineNonZero=0
  DO kk=1,nDOFVarElem
    IF(kk.GT.ii)THEN
      IF(ABS(ILU0(ii,kk)).GT.epsZero)THEN
        jj=jj+1
        IU(iElem)%Entry(jj)=ILU0(ii,kk)
        IU(iELem)%JENTRY(jj)=kk
        lineNonZero=lineNonZero+1
        IF(first)THEN
          IU(iElem)%IENTRY(ii)=jj
          first=.FALSE.
        END IF !first
      END IF ! zero
    END IF ! upper part
  END DO ! kk
  ! modification for zero line or single value lines
  IF(lineNonZero.EQ.0)THEN
    IF(ii.EQ.1)THEN
      IU(iElem)%iEntry(1)=1
    ELSE
      IU(iElem)%iEntry(ii)=IU(iElem)%iEntry(ii-1)+Lastline
    END IF
  END IF
  ! last line
  LastLine=LineNonZero
END DO ! ii
 IU(iElem)%iEntry(nDOFVarElem)= IU(iElem)%iEntry(1)+nUNonZeros(iElem)


! L part of ILU0
! simple version
jj=0
singleValue=.FALSE.
DO ii=1,nMTriangle,1
  first=.TRUE.
  lineNonZero=0
  iEntry=ii+1
  DO kk=1,nDOFVarElem
    IF(kk.LT.iEntry)THEN
      IF(ABS(ILU0(iEntry,kk)).GT.epsZero)THEN
        jj=jj+1
        IL(iElem)%Entry(jj) =ILU0(iEntry,kk)
        IL(iElem)%jEntry(jj)=kk
        lineNonZero=lineNonZero+1
        IF(first)THEN
          IL(iElem)%iEntry(ii)=jj
          first=.FALSE.
        END IF ! first
      END IF ! zero
    ELSE ! kk GE iEntry
      CYCLE
    END IF ! lower
  END DO ! kk
  ! modification for zero line or single value lines
  IF(lineNonZero.EQ.0)THEN
    IF(ii.EQ.1)THEN
      IL(iElem)%iEntry(ii)=1
    ELSE
      IL(iElem)%iEntry(ii)=IL(iElem)%iEntry(ii-1)+Lastline
    END IF
  END IF
  ! last line
  LastLine=LineNonZero
END DO ! ii
IL(iElem)%iEntry(nDOFVarElem)=IL(iElem)%iEntry(1)+nLNonZeros(iElem)

END SUBROUTINE BuildILU0

!==================================================================================================================================
!> Application of block ILU0 preconditioner
!==================================================================================================================================
SUBROUTINE ApplyILU(Vin,Vout)
! MODULES
USE MOD_PreProc
USE MOD_Implicit_Vars     ,ONLY: nDOFVarElem
USE MOD_SparseILU_Vars    ,ONLY: Dinv,IL,IU,nMTriangle
USE MOD_Mesh_Vars         ,ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                                    :: Vin(1:nDOFVarElem,1:nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                   :: Vout(1:nDOFVarElem,1:nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                               :: Vcalc(1:nDOFVarElem,1:nElems)
INTEGER                                            :: ii,k1,k2,jj,iEntry,iElem
!==================================================================================================================================
! Forward eliminiation
Vcalc=Vin
DO iElem=1,nElems
  DO ii=1,nMTriangle,1
    iEntry=ii+1
    k1=IL(iElem)%iEntry(ii)
    k2=IL(iElem)%iEntry(ii+1)-1
    DO jj=k1,k2
      Vcalc(iEntry,iElem)=Vcalc(iEntry,iElem)-IL(iElem)%Entry(jj)*Vcalc(IL(iELEM)%jEntry(jj),iElem)
    END DO ! jj
  END DO ! ii
  ! backward elimination
  ! init backward Gauss
  Vout(nDOFVarElem,iElem) = Vcalc(nDOFVarElem,iElem)*Dinv(nDOFVarElem,iElem)
  DO ii=nMTriangle,1,-1
    k1=IU(iElem)%iEntry(ii)
    k2=IU(iElem)%iEntry(ii+1)-1
    iEntry=ii
    DO jj=k1,k2
      Vcalc(iEntry,iElem)=Vcalc(iEntry,iElem)-IU(iElem)%Entry(jj)*Vout(IU(iElem)%jEntry(jj),iElem)
    END DO ! jj
    Vout(iEntry,iElem)=Vcalc(iEntry,iElem)*Dinv(iEntry,iElem)
  END DO ! ii
END DO

END SUBROUTINE ApplyILU

!===================================================================================================================================
!> Finalize Sparse ILU
!===================================================================================================================================
SUBROUTINE FinalizeSparseILU()
! MODULES
USE MOD_SparseILU_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Dinv)
SDEALLOCATE(nUNonZeros)
SDEALLOCATE(nLNonZeros)
SDEALLOCATE(IU)
SDEALLOCATE(IL)
END SUBROUTINE FinalizeSparseILU

END MODULE MOD_SparseILU
