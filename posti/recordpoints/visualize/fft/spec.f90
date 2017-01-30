#include "flexi.h"
MODULE MOD_spec
!===================================================================================================================================
! Module to handle the Recordpoints
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitSpec
  MODULE PROCEDURE InitSpec
END INTERFACE

INTERFACE spec
  MODULE PROCEDURE spec
END INTERFACE

INTERFACE FinalizeSpec
  MODULE PROCEDURE FinalizeSpec
END INTERFACE
PUBLIC :: InitSpec,spec,FinalizeSpec
!===================================================================================================================================

CONTAINS

SUBROUTINE InitSpec()
!===================================================================================================================================
! Initialize all necessary information for interpolation and FFT
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars          ,ONLY: nSamples_global,RPTime
USE MOD_RPSet_Vars           ,ONLY: nRP_global
USE MOD_OutputRPVisu_Vars          ,ONLY: nSamples_out
USE MOD_RPInterpolation_Vars ,ONLY: dt_out,TEnd
USE MOD_Parameters       ,ONLY: nBlocks,samplingFreq,BlockSize,cutoffFreq,doPSD,doHanning,fourthDeriv
USE MOD_spec_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL         :: block_fac
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')' INIT SPECTRAL ANALYSIS...'
! fourth derivative truncation error allows only relatively large dt's.
! interpolate to equidistant time grid with dt~0.001
IF (fourthDeriv) THEN
  dt_out = (RPTime(nSamples_global)-RPTime(1))/REAL(nSamples_out-1)
  IF (dt_out .LE. 1E-3) THEN
   nSamples_out = INT(nSamples_out * dt_out/1E-3)
  END IF
END IF
    
! Make BlockSize even
IF(MOD(Blocksize,2).NE.0) THEN
  WRITE(UNIT_StdOut,'(A)')'   WARNING: Blocksize should be even!!! Equalizing...'
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
  nSamples_out = NINT(nSamples_out/(block_fac*2))*block_fac*2+1
  TEnd =RPTime(nSamples_global) 
END IF
dt_out = (TEnd-RPTime(1))/REAL(nSamples_out-1)

WRITE(UNIT_stdOut,'(A)')'  Equidistant time grid:'
WRITE(UNIT_StdOut,'(A18,ES16.7)')'             dt = ',dt_out
WRITE(UNIT_StdOut,'(A18,I8    )')'      Blocksize = ',Blocksize
WRITE(UNIT_StdOut,'(A18,I8    )')'        nBlocks = ',nBlocks
WRITE(UNIT_StdOut,'(A18,I8    )')'       nSamples = ',nSamples_out
WRITE(UNIT_StdOut,'(A18,ES16.7)')'       new TEnd = ',TEnd
WRITE(UNIT_StdOut,'(A18,ES16.7)')'TEnd next block = ',TEnd+REAL(BlockSize)*0.5*dt_out
WRITE(UNIT_stdOut,'(A)')' DONE.'
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitSpec



SUBROUTINE spec()
!===================================================================================================================================
! Initialize all necessary information to perform interpolation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars          ,ONLY: RPTime
USE MOD_RPSet_Vars           ,ONLY: nRP_global
USE MOD_RPInterpolation_Vars
USE MOD_OutputRPVisu_Vars          ,ONLY: nSamples_out,RPData_out
USE MOD_Parameters       ,ONLY: doPSD,doFFT,nVarVisu,nBlocks,cutoffFreq,doHanning,fourthDeriv,thirdOct
USE MOD_Parameters       ,ONLY: u_infPhys,chordPhys 
USE MOD_RPSet_Vars           ,ONLY: nLines,Lines
USE FFTW3
USE MOD_spec_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample,iVar,iRP,iStart,iEnd
REAL                            :: df
COMPLEX,ALLOCATABLE             :: in(:),out(:)
INTEGER(KIND=8)                 :: plan
REAL                            :: Time
REAL                            :: Time_Block
REAL,ALLOCATABLE                :: RPData_tmp(:)
REAL                            :: pi
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

pi=acos(-1.)

!4th deriv difference according to ulli
!Test Beispiel
!DO iSample=1,nSamples_out
!!  RPData_out(:,:,iSample) = cos(10.*pi*2.*REAL(iSample-1)/REAL(nSamples_out-1))+cos(2*pi*100*RPTime(iSample))*1E-10
!!  RPData_out(:,:,iSample) = cos(10.*pi*2.*REAL(iSample-1)/REAL(nSamples_block))
!END DO


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
      IF(fourthDeriv) CALL deriv(RPData_tmp)
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
          RPData_spec(iVar,iRP,:)=RPData_spec(iVar,iRP,:)+out(1:nSamples_spec)**2
        ELSE
          RPData_spec(iVar,iRP,:)=RPData_spec(iVar,iRP,:)+out(1:nSamples_spec)
        END IF
      END DO ! nBlocks
      IF (doPSD) THEN
        RPData_spec(iVar,iRP,2:nSamples_spec)=0.5*RPData_spec(iVar,iRP,2:nSamples_spec)/nBlocks/df
        IF (fourthDeriv) THEN
          RPData_spec(iVar,iRP,:)= RPData_spec(iVar,iRP,:)/(RPData_freq(:)*2*pi)**8 
        END IF
        ! MS value on the first index
        RPData_spec(iVar,iRP,1)=SQRT(SUM(RPData_spec(iVar,iRP,2:nSamples_spec))*df)
        RMS_PSD=SQRT(SUM(RPData_spec(iVar,iRP,2:nSamples_spec))*df)
      ELSE
        RPData_spec(iVar,iRP,:)=RPData_spec(iVar,iRP,:)/nBlocks
        IF (fourthDeriv) THEN
          RPData_spec(iVar,iRP,:)= RPData_spec(iVar,iRP,:)/(RPData_freq(:)*2*pi)**4 
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
  nSamples_Oct = 30.

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
END SUBROUTINE spec



SUBROUTINE hanning(nSamples,RPData)
!===================================================================================================================================
! window function for spectra. normalized such that SUM(windowCoeff)=1. For power spectra SUM(windowCoeff)²=1. 
!===================================================================================================================================
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)          :: nSamples
COMPLEX,INTENT(INOUT)       :: RPData(nSamples)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSample
REAL               :: pi,windowCoeff,windowNorm
!===================================================================================================================================
pi=ACOS(-1.)
windowNorm=0.
DO iSample=1,nSamples
  windowCoeff=0.5*(1.-COS(2.*pi*(iSample)/(nSamples)))
  RPData(iSample)= RPData(iSample)*windowCoeff
  windowCoeff=windowCoeff**2
  windowNorm=windowNorm+windowCoeff
END DO
RPData(:)=RPData(:)/SQRT(windowNorm/nSamples)
END SUBROUTINE hanning



SUBROUTINE deriv(RPData)
!===================================================================================================================================
! Initialize all necessary information to perform interpolation
!===================================================================================================================================
! MODULES
USE MOD_OutputRPVisu_Vars          ,ONLY: nSamples_out
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
END SUBROUTINE deriv



SUBROUTINE FinalizeSpec()
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Spec_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)') '  SPECTRAL ANALYSIS FINALIZED'
SDEALLOCATE(RPData_spec)
SDEALLOCATE(RPData_freq)
END SUBROUTINE FinalizeSpec


END MODULE MOD_spec

