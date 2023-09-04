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
!> Module containing evaluation routines regarding spectral analysis. This includes the central routine that calculates the fast
!> Fourier transform or the power density spectra. Also included are the routines to apply windowing or calculate a
!> fourt derivative.
!===================================================================================================================================
MODULE MOD_Spec
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitSpec
  MODULE PROCEDURE InitSpec
END INTERFACE

INTERFACE Spec
  MODULE PROCEDURE Spec
END INTERFACE

INTERFACE FinalizeSpec
  MODULE PROCEDURE FinalizeSpec
END INTERFACE
PUBLIC :: InitSpec,Spec,FinalizeSpec
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize all necessary information for interpolation and FFT
!===================================================================================================================================
SUBROUTINE InitSpec()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars          ,ONLY: nSamples_global,RPTime
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out
USE MOD_RPInterpolation_Vars ,ONLY: dt_out,TEnd
USE MOD_ParametersVisu       ,ONLY: nBlocks,samplingFreq,BlockSize,fourthDeriv
USE MOD_Spec_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL         :: block_fac
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')' INIT SPECTRAL ANALYSIS...'
! fourth Derivative truncation error allows only relatively large dt's.
! interpolate to equidistant time grid with dt~0.001
IF (fourthDeriv) THEN
  dt_out = (RPTime(nSamples_global)-RPTime(1))/REAL(nSamples_out-1)
  IF (dt_out .LE. 1E-3) THEN
   nSamples_out = INT(nSamples_out * dt_out/1E-3)
  END IF
END IF

! Make BlockSize even
IF(MOD(Blocksize,2).NE.0) THEN
  WRITE(UNIT_stdOut,'(A)')'   WARNING: Blocksize should be even!!! Equalizing...'
  BlockSize=BlockSize+1
END IF

! if sampling frequency is prescribed, we calculate the number of blocks and
! cut off the signal to fit the number of blocks
IF(samplingFreq.GT.0) THEN
  dt_out=1./samplingFreq
  nSamples_out=INT((RPTime(nSamples_global)-RPTime(1))/dt_out)+1
  nBlocks=INT(2*nSamples_out/BlockSize)-1
  nBlocks=MAX(1,nBlocks)
  nSamples_out=BlockSize*(nBlocks+1)/2
  TEnd   = RPTime(1)+REAL(nSamples_out-1)*dt_out
  IF(TEnd.GT.RPTime(nSamples_global)) &
   STOP 'Error: Desired blockinterval is greater than available time data! Increase sampling frequency or decrease blocksize!'
!if number of blocks is prescribed, we calculate the BlockSize and the dt
ELSE
  ! Create nSamples such as to be multiple of block_fac
  block_fac = (nBlocks+1)/2.
  nSamples_out = NINT(nSamples_out/(nBlocks+1.))*(nBlocks+1)+1
  TEnd =RPTime(nSamples_global)
END IF
dt_out = (TEnd-RPTime(1))/REAL(nSamples_out-1)

WRITE(UNIT_stdOut,'(A)')'  Equidistant time grid:'
WRITE(UNIT_stdOut,'(A18,ES16.7)')'             dt = ',dt_out
WRITE(UNIT_stdOut,'(A18,I8    )')'      Blocksize = ',Blocksize
WRITE(UNIT_stdOut,'(A18,I8    )')'        nBlocks = ',nBlocks
WRITE(UNIT_stdOut,'(A18,I8    )')'       nSamples = ',nSamples_out
WRITE(UNIT_stdOut,'(A18,ES16.7)')'       new TEnd = ',TEnd
WRITE(UNIT_stdOut,'(A18,ES16.7)')'TEnd next block = ',TEnd+REAL(BlockSize)*0.5*dt_out
WRITE(UNIT_stdOut,'(A)')' DONE.'
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitSpec

!===================================================================================================================================
!> Main routine used in spectral analysis. Will perform the fast Fourier transform using the FFTW library in blocks that will then
!> be averaged.
!> Additional calculations include the application of a window function, a conversion to PSD values or the calculation of
!> the 1/3 octave spectrum. The fourth derivatives may also be calculated.
!===================================================================================================================================
SUBROUTINE Spec()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_RPSetVisuVisu_Vars   ,ONLY: nRP_global
USE MOD_RPInterpolation_Vars
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out,RPData_out
USE MOD_ParametersVisu       ,ONLY: doPSD,doFFT,nVarVisu,nBlocks,cutoffFreq,doHanning,fourthDeriv,thirdOct
USE MOD_ParametersVisu       ,ONLY: u_infPhys,chordPhys
USE FFTW3
USE MOD_Spec_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample,iVar,iRP,iStart,iEnd
REAL                            :: df
COMPLEX,ALLOCATABLE             :: in(:),out(:)
INTEGER(KIND=8)                 :: plan
REAL                            :: Time
REAL                            :: Time_Block
REAL,ALLOCATABLE                :: RPData_tmp(:)
REAL                            :: M_t,RMS_t,RMS_PSD
INTEGER                         :: nSamples_block,iBlock
!1/3 Oct Stuff
REAL                            :: iOct, octaveLimUpper, octaveLimLower, counter
REAL,ALLOCATABLE                :: PhysFreq(:)
INTEGER                         :: iSample_Oct
REAL                            :: maxdev
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')

IF (doPSD) THEN
  WRITE(UNIT_stdOut,'(A)')' CALCULATING PSD ...'
ELSE
  WRITE(UNIT_stdOut,'(A)')' CALCULATING FFT ...'
END IF

! Time interval block
nSamples_block = INT(2*(nSamples_out)/(nBlocks+1)) ! Blocks with 50% overlap!!
nSamples_spec  = INT((nSamples_block)/2)+1
Time_Block     = dt_out*REAL(nSamples_block)
df = 1./Time_Block

ALLOCATE(in(nSamples_block))
ALLOCATE(out(nSamples_block))
ALLOCATE(RPData_spec(1:nVarVisu,nRP_global,nSamples_spec))
ALLOCATE(RPData_tmp(nSamples_out))
ALLOCATE(RPData_freq(nSamples_spec))
DO iSample=1,nSamples_spec
  RPData_freq(iSample)=(iSample-1)*df
END DO

WRITE(UNIT_stdOut,'(A,I8)')   '       Number of averaging blocks:',nBlocks
WRITE(UNIT_stdOut,'(A,I8)')   '      Number of samples per block:',nSamples_block
WRITE(UNIT_stdOut,'(A,F16.7)')'              Time interval block:',Time_Block
WRITE(UNIT_stdOut,'(A,F16.7)')'             Frequency resolution:',df
IF(cutoffFreq.NE.-999)THEN
  WRITE(UNIT_stdOut,'(A,F16.7)')'    User defined Cutoff Frequency:',cutoffFreq
  nSamples_spec=MIN(NINT(cutoffFreq/dF)+1,nSamples_spec)
END IF
WRITE(UNIT_stdOut,'(A,I8)')   '  No. spectrum output samples FFT:',nSamples_spec
WRITE(UNIT_stdOut,'(A,F16.7)')'                Nyquist frequency:',0.5*df*REAL(nSamples_block-1)
WRITE(UNIT_stdOut,'(A,F16.7)')'          Max. resolved frequency:',RPData_freq(nSamples_spec)

CALL DFFTW_PLAN_DFT_1D(plan,nSamples_block,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
GETTIME(StartTime)

IF(doPSD .OR. doFFT) THEN
  maxdev=0.
  RPData_spec=0.
  DO iRP=1,nRP_global
    IF(MOD(iRP,10).EQ.0) THEN
      WRITE(UNIT_stdOut,*)'   Processing RP ',iRP,' of ',nRP_global
    END IF
    DO iVar=1,nVarVisu
      RPData_tmp=RPData_out(iVar,iRP,:)
      IF(fourthDeriv) CALL Deriv(RPData_tmp)
      DO iBlock=1,nBlocks
        iStart=INT(0.5*REAL((iBlock-1)*nSamples_block))+1
        iEnd  =INT(REAL(0.5*(iBlock+1)*nSamples_block))
        DO iSample=1,nSamples_block
          in(iSample)= RPData_tmp(iStart-1+iSample)
        END DO
        ! Hann Window
        IF(doHanning) CALL hanning(nSamples_block,in)
        CALL DFFTW_EXECUTE_DFT(plan, in, out)
        out(1:nSamples_spec)=2./REAL(nSamples_block)*ABS(out(1:nSamples_spec))
        out(1)=0.5*out(1) !mean value
        IF (doPSD) THEN
          RPData_spec(iVar,iRP,:)=REAL(RPData_spec(iVar,iRP,:)+out(1:nSamples_spec)**2)
        ELSE
          RPData_spec(iVar,iRP,:)=REAL(RPData_spec(iVar,iRP,:)+out(1:nSamples_spec))
        END IF
      END DO ! nBlocks
      IF (doPSD) THEN
        RPData_spec(iVar,iRP,2:nSamples_spec)=0.5*RPData_spec(iVar,iRP,2:nSamples_spec)/nBlocks/df
        IF (fourthDeriv) THEN
          RPData_spec(iVar,iRP,:)= RPData_spec(iVar,iRP,:)/(RPData_freq(:)*2*PP_Pi)**8
        END IF
        ! MS value on the first index
        RPData_spec(iVar,iRP,1)=SQRT(SUM(RPData_spec(iVar,iRP,2:nSamples_spec))*df)
        RMS_PSD=SQRT(SUM(RPData_spec(iVar,iRP,2:nSamples_spec))*df)
      ELSE
        RPData_spec(iVar,iRP,:)=RPData_spec(iVar,iRP,:)/nBlocks
        IF (fourthDeriv) THEN
          RPData_spec(iVar,iRP,:)= RPData_spec(iVar,iRP,:)/(RPData_freq(:)*2*PP_Pi)**4
        END IF
        ! MS value on the first index
      !  RPData_spec(iVar,iRP,1)=SUM(RPData_spec(iVar,iRP,2:nSamples_spec))**2/nSamples_spec
        RPData_spec(iVar,iRP,1)=SQRT(SUM(RPData_spec(iVar,iRP,2:nSamples_spec)**2))
      END IF
    !  !calculate RMS with time signal
      M_t=SUM(RPData_tmp(1:nSamples_out))/nSamples_out
      RMS_t=SQRT(SUM((RPData_tmp(1:nSamples_out)-M_t)*(RPData_tmp(1:nSamples_out)-M_t))/nSamples_out)
    !  WRITE(UNIT_stdOut,*)'------------------------------'
    !  WRITE(UNIT_stdOut,*)' mean timedata',SUM(RPData_tmp(1:nSamples_out))/nSamples_out
    !  WRITE(UNIT_stdOut,*)' mean PSD     ',SQRT(RPData_spec(iVar,iRP,1))
    !  WRITE(UNIT_stdOut,*)' RMS timedata ',RMS_t
    !  WRITE(UNIT_stdOut,*)' RMS from PSD ',RMS_PSD
      maxdev=MAX(maxdev,(RMS_t-RMS_PSD)**2/RMS_t**2)
    END DO ! iVar
  END DO   ! iRP
  GETTIME(Time)
  WRITE(UNIT_stdOut,*)'------------------------------'
  WRITE(UNIT_stdOut,*)' Max. relative error in RMS: ', SQRT(maxdev)
  WRITE(UNIT_stdOut,'(A,F8.2,A)')' DONE! [',Time-StartTime,' sec ]'
END IF!(doPSD .OR. doFFT)

!1/3 octave average
IF(ThirdOct) THEN
  ALLOCATE(PhysFreq(nSamples_spec))
  PhysFreq = RPData_freq*u_infPhys/chordPhys
  nSamples_Oct = 30

  ALLOCATE(RPData_Oct(1:nVarVisu,nRP_global,nSamples_Oct))
  ALLOCATE(RPData_FreqOct(nSamples_Oct))
  RPData_Oct=0.

  iOct = -18.
  DO iSample_Oct=1,nSamples_Oct
    RPData_FreqOct(iSample_Oct) = 1000.*2**(iOct/3.)
    iOct = iOct+1
  END DO

  DO iRP=1,nRP_global
    DO iVar=1,nVarVisu
      iOct = -18.
      DO iSample_Oct=1,nSamples_Oct
        octaveLimUpper=1000.*2**(iOct/3.+1./6.)
        octaveLimLower=1000.*2**(iOct/3.-1./6.)
        iOct = iOct+1
        counter=0.
        IF(doPSD) THEN ! PSD
          DO iSample=1,nSamples_Spec
            IF((PhysFreq(iSample).GE.octaveLimLower).AND.(PhysFreq(iSample).LT.octaveLimUpper))THEN
              RPData_Oct(iVar,iRP,iSample_Oct) = RPData_Oct(iVar,iRP,iSample_Oct) + RPData_spec(iVar,iRP,iSample)*df
              counter=counter+1
            END IF
          END DO
        ELSE ! FFT
          DO iSample=1,nSamples_Spec
            IF((PhysFreq(iSample).GE.octaveLimLower).AND.(PhysFreq(iSample).LT.octaveLimUpper))THEN
              RPData_Oct(iVar,iRP,iSample_Oct) = RPData_Oct(iVar,iRP,iSample_Oct) + RPData_spec(iVar,iRP,iSample)**2
              counter=counter+1
            END IF
          END DO
          RPData_Oct(iVar,iRP,iSample_Oct)=RPData_Oct(iVar,iRP,iSample_Oct)*0.5 ! FFT case: add factor 0.5
        END IF
        IF(counter.EQ.0.)THEN
!          write(*,*)'NO frequence found in range',octaveLimUpper,'--',octaveLimLower
!          write(*,*)'SET Variable to 99999999'
          RPData_Oct(iVar,iRP,iSample_Oct)=99999999.
        ELSE
!          write(*,*)'For 1/3Octave frequency',RPData_FreqOct(iSample_Oct),counter,'samples found'
        END IF
      END DO
    END DO ! iVar
  END DO   ! iRP
END IF!1/3 Oct
DEALLOCATE(in,out,RPData_tmp)
CALL DFFTW_DESTROY_PLAN(plan)
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE Spec

!===================================================================================================================================
!> Apply a window function for spectra. Uses the window function named after Julius von Hann, often referred to as Henning function.
!> Normalized such that SUM(windowCoeff)=1. For power spectra SUM(windowCoeff)²=1.
!> For details, see e.g. Harris, F. J. (1978). "On the use of windows for harmonic analysis with the discrete Fourier transform",
!> Proceedings of the IEEE. doi:10.1109/PROC.1978.10837.
!===================================================================================================================================
SUBROUTINE Hanning(nSamples,RPData)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)          :: nSamples
COMPLEX,INTENT(INOUT)       :: RPData(nSamples)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSample
REAL               :: windowCoeff,windowNorm
!===================================================================================================================================
windowNorm=0.
DO iSample=1,nSamples
  windowCoeff=0.5*(1.-COS(2.*PP_Pi*(iSample)/(nSamples)))
  RPData(iSample)= RPData(iSample)*windowCoeff
  windowCoeff=windowCoeff**2
  windowNorm=windowNorm+windowCoeff
END DO
RPData(:)=RPData(:)/SQRT(windowNorm/nSamples)
END SUBROUTINE Hanning

!===================================================================================================================================
!> Compute the fourth derivative w.r.t time of the RP signal using a central difference of second order accuray.
!===================================================================================================================================
SUBROUTINE Deriv(RPData)
! MODULES
USE MOD_OutputRPVisu_Vars    ,ONLY: nSamples_out
USE MOD_RPInterpolation_Vars ,ONLY: dt_out
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(INOUT) :: RPData(nSamples_out)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSample
REAL               :: RPData_out(nSamples_out)
REAL               :: dh,km1,km2,k0,k1,k2
!===================================================================================================================================
dh = dt_out**4
DO iSample=4,nSamples_out-3
  km2 = RPData(iSample-2)
  km1 = RPData(iSample-1)
  k0  = RPData(iSample)
  k1  = RPData(iSample+1)
  k2  = RPData(iSample+2)
  RPData_out(iSample)= dh*(km2-4*km1+6*k0-4*k1+k2)
END DO
RPData=RPData_out
END SUBROUTINE Deriv

!===================================================================================================================================
!> Deallocate global variables
!===================================================================================================================================
SUBROUTINE FinalizeSpec()
! MODULES
USE MOD_Globals
USE MOD_Spec_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(RPData_spec)
SDEALLOCATE(RPData_freq)
WRITE(UNIT_stdOut,'(A)') ' SPECTRAL ANALYSIS FINALIZED'

END SUBROUTINE FinalizeSpec

END MODULE MOD_Spec
