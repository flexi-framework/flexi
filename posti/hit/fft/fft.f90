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

MODULE MOD_HIT_FFT
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitFFT
  MODULE PROCEDURE InitFFT
END INTERFACE

INTERFACE ComputeFFT_R2C
  MODULE PROCEDURE ComputeFFT_R2C
END INTERFACE

INTERFACE ComputeFFT_C2R
  MODULE PROCEDURE ComputeFFT_C2R
END INTERFACE

INTERFACE Interpolate_DG2FFT
  MODULE PROCEDURE Interpolate_DG2FFT
END INTERFACE

INTERFACE Interpolate_FFT2DG
  MODULE PROCEDURE Interpolate_FFT2DG
END INTERFACE

INTERFACE EvalFourierAtDGCoords
  MODULE PROCEDURE EvalFourierAtDGCoords
END INTERFACE

INTERFACE FinalizeFFT
  MODULE PROCEDURE FinalizeFFT
END INTERFACE

PUBLIC:: InitFFT,FinalizeFFT
PUBLIC:: ComputeFFT_R2C,ComputeFFT_C2R
PUBLIC:: Interpolate_DG2FFT,Interpolate_FFT2DG
PUBLIC:: EvalFourierAtDGCoords
!==================================================================================================================================



CONTAINS

!===================================================================================================================================
!> Initialize FFT. Define necessary parameters, allocate arrays for solution and auxiliary FFT variables.
!===================================================================================================================================
SUBROUTINE InitFFT()
! MODULES
USE MOD_Globals
USE MOD_HIT_FFT_Vars
#if USE_OPENMP
USE OMP_Lib
USE FFTW3
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: i,j,k
#if USE_OPENMP
INTEGER  :: void
#endif
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT FFT...'

#if USE_OPENMP
! Initialize FFTW with maximum number of OpenMP threads if enabled
void = FFTW_INIT_THREADS()
CALL FFTW_PLAN_WITH_NTHREADS(OMP_GET_MAX_THREADS())
#endif

Nc = FLOOR(REAL(N_FFT)/2.)
Endw(1)   = Nc+1
Endw(2:3) = N_FFT

kmax=NINT(SQRT(REAL((3*N_FFT**2))))+1

! Allocate array for wave numbers
ALLOCATE(Localk(1:4,1:Endw(1),1:Endw(2),1:Endw(3)))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
!$OMP DO
! fill the wave space
DO i=1,Endw(1)
  DO j=1,Endw(2)
    DO k=1,Endw(3)
      ! due to symmetry, only nx/2+1 kx are stored, ranging from 0..1..2... to Nx/2 (Nx/2+1 waves) e.g. Nx=7: k1=:0,1,2,3; Nx=4:0,1,2
      localk(1,i,j,k)=i-1
      ! ny waves are stored, order: 0,1,2,..,Nyquist,..,..,-2,-1, e.g Ny=4: 0,1,2,-1, Ny=7: 0,1,2,3,-3,-2,-1
      IF (j.LE.Nc+1) THEN
        localk(2,i,j,k)=j-1
      ELSE
        localk(2,i,j,k)=-(N_FFT+1-j)
      END IF
      IF (k.LE.Nc+1) THEN
        localk(3,i,j,k)=k-1
      ELSE
        localk(3,i,j,k)=-(N_FFT+1-k)
      END IF
      localk(4,i,j,k)=NINT(SQRT(REAL(localk(1,i,j,k)*localk(1,i,j,k) &
                                    +localk(2,i,j,k)*localk(2,i,j,k) &
                                    +localk(3,i,j,k)*localk(3,i,j,k))))
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

SWRITE(UNIT_stdOut,'(A)')' INIT FFT DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitFFT


!===================================================================================================================================
!> Wrapper routine for performing a 3D FFT (physical to Fourier) by means of the FFTW library on a array on global coordinates and
!> nVar_In separate variables. The computation time can be printed to the standard output by setting the flag doPrintTime.
!===================================================================================================================================
SUBROUTINE ComputeFFT_R2C(nVar_In,U_Global,U_FFT,doPrintTime)
! MODULES
USE MOD_Globals
USE MOD_HIT_FFT_Vars       ,ONLY: N_FFT,Endw
USE FFTW3
#if USE_OPENMP
USE OMP_Lib
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nVar_In         !< Number of independent variables in first dimension of input array
REAL,INTENT(IN)       :: U_Global(nVar_In,1:N_FFT  ,1:N_FFT  ,1:N_FFT  ) !< Global DG solution interpolated to equidistant points
COMPLEX,INTENT(OUT)   :: U_FFT(   nVar_In,1:Endw(1),1:Endw(2),1:Endw(3)) !< Complex Fourier representation of input array
LOGICAL,INTENT(IN),OPTIONAL :: doPrintTime                               !< Print execution time to standard output
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)  :: plan
INTEGER          :: iVar
REAL             :: Time
REAL             :: U_r(1:N_FFT  ,1:N_FFT  ,1:N_FFT  ) ! Real global DG solution per variable
COMPLEX          :: U_c(1:Endw(1),1:Endw(2),1:Endw(3)) ! Complex FFT solution per variable
LOGICAL          :: doPrintTimeLoc
!===================================================================================================================================
doPrintTimeLoc = .FALSE.
IF (PRESENT(doPrintTime)) doPrintTimeLoc = doPrintTime

IF(doPrintTimeLoc) WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' COMPUTE FFT FROM PHYSICAL TO FOURIER...'
Time = OMP_FLEXITIME()

! Create plan once for local arrays and reuse it for each variable
CALL DFFTW_PLAN_DFT_R2C_3D(plan,N_FFT,N_FFT,N_FFT,U_r,U_c,FFTW_ESTIMATE)

! Compute FFT for each variable. Local arrays ensure data to be contiguous in memory.
DO iVar=1,nVar_In
  U_r = U_Global(iVar,:,:,:)
  CALL DFFTW_Execute(plan)   ! Input/Output Arrays are already specified in plan generation
  U_FFT(iVar,:,:,:) = U_c
END DO

! Normalize data (FFTW does not normalize FFT results)
U_FFT=U_FFT/REAL(N_FFT**3)

! Release resources allocated with plan
CALL DFFTW_DESTROY_PLAN(plan)

IF(doPrintTimeLoc) WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',OMP_FLEXITIME()-Time,'s]'
END SUBROUTINE ComputeFFT_R2C


!===================================================================================================================================
!> Wrapper routine for performing an inverse 3D FFT (Fourier to physical) by means of the FFTW library on a array with a complex
!> Fourier representation of the data with nVar_In separate variables. The computation time can be printed to the standard output by
!> setting the flag doPrintTime.
!===================================================================================================================================
SUBROUTINE ComputeFFT_C2R(nVar_In,U_FFT,U_Global,doPrintTime)
! MODULES
USE MOD_Globals
USE MOD_HIT_FFT_Vars       ,ONLY: N_FFT,Endw
USE FFTW3
#if USE_OPENMP
USE OMP_Lib
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nVar_In         !< Number of independent variables in first dimension of input array
COMPLEX,INTENT(IN)    :: U_FFT(   nVar_In,1:Endw(1),1:Endw(2),1:Endw(3)) !< Global DG solution interpolated to equidistant points
REAL,INTENT(OUT)      :: U_Global(nVar_In,1:N_FFT  ,1:N_FFT  ,1:N_FFT  ) !< Complex Fourier representation of input array
LOGICAL,INTENT(IN),OPTIONAL :: doPrintTime                               !< Print execution time to standard output
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8) :: plan
INTEGER         :: iVar
REAL            :: Time
REAL            :: U_r(1:N_FFT  ,1:N_FFT  ,1:N_FFT  ) ! Real global DG solution per variable
COMPLEX         :: U_c(1:Endw(1),1:Endw(2),1:Endw(3)) ! Complex FFT solution per variable
LOGICAL         :: doPrintTimeLoc
!===================================================================================================================================
doPrintTimeLoc = .FALSE.
IF (PRESENT(doPrintTime)) doPrintTimeLoc = doPrintTime

IF(doPrintTimeLoc) WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' COMPUTE FFT FROM FOURIER TO PHYSICAL...'
Time = OMP_FLEXITIME()

! Create plan once for local arrays and reuse it for each variable
CALL DFFTW_PLAN_DFT_C2R_3D(plan,N_FFT,N_FFT,N_FFT,U_c,U_r,FFTW_ESTIMATE)

! Compute FFT for each variable. Local arrays ensure data to be contiguous in memory.
DO iVar=1,nVar_In
  U_c = U_FFT(iVar,:,:,:)
  CALL DFFTW_Execute(plan)   ! Input/Output Arrays are already specified in plan generation
  U_Global(iVar,:,:,:) = U_r
END DO

! Release resources allocated with plan
CALL DFFTW_DESTROY_PLAN(plan)

IF(doPrintTimeLoc) WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',OMP_FLEXITIME()-Time,'s]'
END SUBROUTINE ComputeFFT_C2R


!===================================================================================================================================
!> Interpolate via ChangeBasis() state from DG to visu (equidistant) coordinates and assemble the elementwise DG representation to
!> a global representation. This requires the IJK-sorting of the Mesh, i.e. the array Elem_IJK, to be known.
!===================================================================================================================================
SUBROUTINE Interpolate_DG2FFT(NodeType_In,nVar_In,U_DG,U_Global)
! MODULES
USE MOD_Globals
USE MOD_PreProc               ,ONLY: PP_N
USE MOD_HIT_FFT_Vars          ,ONLY: N_FFT,NCalc
USE MOD_Mesh_Vars             ,ONLY: Elem_IJK,nElems
USE MOD_Interpolation         ,ONLY: GetVandermonde
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: nVar_In          !< Number of independent variables in first dimension of input array
CHARACTER(LEN=255),INTENT(IN)    :: NodeType_In      !< NodeType of input array U_DG
REAL,INTENT(INOUT)               :: U_DG(    nVar_In,0:PP_N ,0:PP_N ,0:PP_N ,nElems) !< Elementwise input array
REAL,INTENT(OUT)                 :: U_Global(nVar_In,1:N_FFT,1:N_FFT,1:N_FFT       ) !< Global solution representation
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: iElem
INTEGER       :: i,ii,iGlob
INTEGER       :: j,jj,jGlob
INTEGER       :: k,kk,kGlob
REAL          :: U_aux(nVar_In,0:NCalc,0:NCalc,0:NCalc)
REAL          :: VdmGaussEqui(0:NCalc,0:PP_N)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' INTERPOLATE DG SOLUTION TO FFT COORDS...'

! Vandermonde to interpolate from HDF5_Nodetype to equidistant points
CALL GetVandermonde(PP_N,NodeType_In,NCalc,'VISU_INNER',VdmGaussEqui)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iElem,i,ii,iGlob,j,jj,jGlob,k,kk,kGlob,U_aux)
!$OMP DO
! Loop to get the nodal U_DG solution into a global solution in ijk form
DO iElem=1,nElems
  ! Get solution in each element on equidistant points
  CALL ChangeBasis3D(nVar_In,PP_N,NCalc,VdmGaussEqui,U_DG(:,:,:,:,iElem),U_aux(:,:,:,:))

  ! Fill the global solution array using the ijk sorting
  ii=Elem_IJK(1,iElem)
  jj=Elem_IJK(2,iElem)
  kk=Elem_IJK(3,iElem)
  DO k=0,NCalc
    kGlob=(kk-1)*(NCalc+1)+k+1
    DO j=0,NCalc
      jGlob=(jj-1)*(NCalc+1)+j+1
      DO i=0,NCalc
        iGlob=(ii-1)*(NCalc+1)+i+1
        U_Global(:,iGlob,jGlob,kGlob) = U_aux(:,i,j,k)
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')'DONE'
END SUBROUTINE Interpolate_DG2FFT


!===================================================================================================================================
!> Sort a global solution representation back into the elementwise DG representation and interpolate via ChangeBasis() state from
!> visu (equidistant) to DG coordinates. This requires the IJK-sorting of the Mesh, i.e. the array Elem_IJK, to be known.
!===================================================================================================================================
SUBROUTINE Interpolate_FFT2DG(NodeType_In,nVar_In,U_Global,U_DG)
! MODULES
USE MOD_Globals
USE MOD_PreProc               ,ONLY: PP_N
USE MOD_HIT_FFT_Vars          ,ONLY: N_FFT,NCalc
USE MOD_Mesh_Vars             ,ONLY: Elem_IJK,nElems
USE MOD_Interpolation         ,ONLY: GetVandermonde
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: nVar_In         !< Number of independent variables in first dimension of input array
CHARACTER(LEN=255),INTENT(IN)    :: NodeType_In     !< NodeType of input array U_DG
REAL,INTENT(INOUT)               :: U_Global(nVar_In,1:N_FFT,1:N_FFT,1:N_FFT       ) !< Elementwise input array
REAL,INTENT(OUT)                 :: U_DG(    nVar_In,0:PP_N ,0:PP_N ,0:PP_N ,nElems) !< Global solution representation
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: iElem
INTEGER       :: i,ii,iGlob
INTEGER       :: j,jj,jGlob
INTEGER       :: k,kk,kGlob
REAL          :: U_aux(nVar_In,0:NCalc,0:NCalc,0:NCalc)
REAL          :: VdmEquiGauss(0:NCalc,0:PP_N)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' INTERPOLATE FFT SOLUTION TO DG COORDS...'

! Vandermonde to interpolate from equidistant points to HDF5_Nodetype
CALL GetVandermonde(NCalc,'VISU_INNER',PP_N,NodeType_In,VdmEquiGauss)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iElem,i,ii,iGlob,j,jj,jGlob,k,kk,kGlob,U_aux)
!$OMP DO
! Loop to get the global solution in ijk form back to a nodal U_DG solution
DO iElem=1,nElems
  ! Fill the global solution array using the ijk sorting
  ii=Elem_IJK(1,iElem)
  jj=Elem_IJK(2,iElem)
  kk=Elem_IJK(3,iElem)
  DO k=0,NCalc
    kGlob=(kk-1)*(NCalc+1)+k+1
    DO j=0,NCalc
      jGlob=(jj-1)*(NCalc+1)+j+1
      DO i=0,NCalc
        iGlob=(ii-1)*(NCalc+1)+i+1
        U_aux(:,i,j,k) = U_Global(:,iGlob,jGlob,kGlob)
      END DO
    END DO
  END DO

  ! Get solution in each element on equidistant points
  CALL ChangeBasis3D(nVar_In,NCalc,PP_N,VdmEquiGauss,U_aux(:,:,:,:),U_DG(:,:,:,:,iElem))
END DO
!$OMP END DO
!$OMP END PARALLEL

SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')'DONE'
END SUBROUTINE Interpolate_FFT2DG


!===================================================================================================================================
!> Transform a Fourier representation of the solution back to DG coords by evaluating the global Fourier basis/solution at the DG
!> interpolation points.
!===================================================================================================================================
SUBROUTINE EvalFourierAtDGCoords(nVar_In,U_FFT,U_DG,doPrintTime)
! MODULES
USE MOD_Globals
USE MOD_PreProc               ,ONLY: PP_N
USE MOD_Mesh_Vars             ,ONLY: Elem_xGP,nElems
USE MOD_HIT_FFT_Vars          ,ONLY: II,Nc,endw
#if USE_OPENMP
USE OMP_Lib
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nVar_In      !< Number of independent variables in first dimension of input array
COMPLEX,INTENT(IN)    :: U_FFT(nVar_In,1:Endw(1),1:Endw(2),1:Endw(3)       ) !< Global DG solution interpolated to equidist. points
REAL,   INTENT(OUT)   :: U_DG( nVar_In,0:PP_N   ,0:PP_N   ,0:PP_N   ,nElems) !< Complex Fourier representation of input array
LOGICAL,INTENT(IN),OPTIONAL :: doPrintTime                                   !< Print execution time to standard output
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL         :: Time
INTEGER      :: i,j,k,iw,jw,kw,iElem
INTEGER      :: wavenumber
COMPLEX      :: basis
COMPLEX      :: U_k(1:nVar_In,1:endw(1),1:endw(2),0:PP_N)
COMPLEX      :: U_j(1:nVar_In,1:endw(1),0:PP_N   ,0:PP_N)
COMPLEX      :: U_i(1:nVar_In,0:PP_N   ,0:PP_N   ,0:PP_N)
LOGICAL      :: doPrintTimeLoc
!===================================================================================================================================
doPrintTimeLoc = .FALSE.
IF (PRESENT(doPrintTime)) doPrintTimeLoc = doPrintTime

IF(doPrintTimeLoc) WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' EVALUATE FOURIER BASIS AT DG COORDINATES...'
Time = OMP_FLEXITIME()

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iElem,U_i,i,iw,U_j,j,jw,U_k,k,kw,wavenumber,basis)
!$OMP DO
DO iElem=1,nElems
  ! Evaluate basis in z
  U_k=0.
  DO k=0,PP_N
    DO kw=1,endw(3)
      wavenumber=kw-1
      IF (kw.GE.Nc+2) wavenumber=-(2*Nc+1-kw)
      basis=EXP(II*(wavenumber)*Elem_xGP(3,0,0,k,iElem))
      DO jw=1,endw(2)
        DO iw=1,endw(1)
          U_k(:,iw,jw,k) = U_k(:,iw,jw,k)+U_FFT(:,iw,jw,kw)*basis
        END DO !k
      END DO !kw
    END DO !jw
  END DO !iw

  ! Evaluate basis in y
  U_j=0.
  DO j=0,PP_N
    DO jw=1,endw(2)
      wavenumber=jw-1
      IF (jw.GE.Nc+2) wavenumber=-(2*Nc+1-jw)
      basis=EXP(II*(wavenumber)*Elem_xGP(2,0,j,0,iElem))
      DO k=0,PP_N
        DO iw=1,endw(1)
          U_j(:,iw,j,k) = U_j(:,iw,j,k)+U_k(:,iw,jw,k)*basis
        END DO !k
      END DO !j
    END DO !jw
  END DO !iw

  ! Evaluate basis in x
  U_i=0.
  U_j(:,1,:,:)=U_j(:,1,:,:)/2.
  DO i=0,PP_N
    DO iw=1,endw(1)
      basis=EXP(II*(iw-1)*Elem_xGP(1,i,0,0,iElem))
      DO k=0,PP_N
        DO j=0,PP_N
          U_i(:,i,j,k) = U_i(:,i,j,k)+U_j(:,iw,j,k)*basis
        END DO !k
      END DO !j
    END DO !i
  END DO !iw

  U_DG(:,:,:,:,iElem) = 2*REAL(U_i)
END DO !iElem
!$OMP END DO
!$OMP END PARALLEL

IF(doPrintTimeLoc) WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',OMP_FLEXITIME()-Time,'s]'

END SUBROUTINE EvalFourierAtDGCoords


!===================================================================================================================================
!> Finalize FFT
!===================================================================================================================================
SUBROUTINE FinalizeFFT()
! MODULES                                                                                                                          !
USE MOD_HIT_FFT_Vars
USE FFTW3
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_OPENMP
CALL FFTW_CLEANUP_THREADS()
#else
CALL FFTW_CLEANUP()
#endif
SDEALLOCATE(LocalK)

END SUBROUTINE FinalizeFFT

END MODULE MOD_HIT_FFT
