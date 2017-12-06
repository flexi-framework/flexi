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
!> Module to handle filter operations
!==================================================================================================================================
MODULE MOD_Filter
! MODULES
IMPLICIT NONE
PRIVATE

ABSTRACT INTERFACE
  SUBROUTINE FilterInt(U_in,FilterMat)
    USE MOD_PreProc
    USE MOD_Mesh_Vars,ONLY: nElems
    REAL,INTENT(INOUT) :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)
    REAL,INTENT(IN)    :: FilterMat(   0:PP_N,0:PP_N)
  END SUBROUTINE
END INTERFACE

PROCEDURE(FilterInt),POINTER :: Filter_pointer     !< Point to the filter routine to be used

!----------------------------------------------------------------------------------------------------------------------------------

INTEGER,PARAMETER      :: FILTERTYPE_NONE   = 0
INTEGER,PARAMETER      :: FILTERTYPE_CUTOFF = 1
INTEGER,PARAMETER      :: FILTERTYPE_MODAL  = 2
INTEGER,PARAMETER      :: FILTERTYPE_LAF    = 3

INTERFACE InitFilter
  MODULE PROCEDURE InitFilter
END INTERFACE

INTERFACE FinalizeFilter
  MODULE PROCEDURE FinalizeFilter
END INTERFACE

INTERFACE Filter_Selective
  MODULE PROCEDURE Filter_Selective
END INTERFACE

PUBLIC :: InitFilter
PUBLIC :: Filter_pointer
PUBLIC :: Filter_Selective
PUBLIC :: FinalizeFilter
PUBLIC :: DefineParametersFilter
!==================================================================================================================================

#ifdef DEBUG
#if EQNSYSNR==2
! Add dummy interfaces to unused subroutines to suppress compiler warnings.
INTERFACE DUMMY_Filter
  MODULE PROCEDURE Filter
END INTERFACE
INTERFACE DUMMY_Filter_LAF
  MODULE PROCEDURE Filter_LAF
END INTERFACE
#endif
#endif
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersFilter()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Filter")
CALL prms%CreateIntFromStringOption(   'FilterType',        "Type of filter to be applied. None, CutOff, Modal, LAF",&
                                                            'None')
CALL addStrListEntry('FilterType','none',  FILTERTYPE_NONE)
CALL addStrListEntry('FilterType','cutoff',FILTERTYPE_CUTOFF)
CALL addStrListEntry('FilterType','modal', FILTERTYPE_MODAL)
CALL addStrListEntry('FilterType','laf',   FILTERTYPE_LAF)
CALL prms%CreateIntOption(             'NFilter',           "Cut-off mode (FilterType==CutOff or LAF)")
CALL prms%CreateRealOption(            'LAF_alpha',         "Relaxation factor for LAF, see Flad et al. JCP 2016")
CALL prms%CreateRealArrayOption(       'HestFilterParam',   "Parameters for Hesthaven filter (FilterType==Modal)")
END SUBROUTINE DefineParametersFilter


!==================================================================================================================================
!> Initialize all necessary information to perform filtering and filter dealiasing
!==================================================================================================================================
SUBROUTINE InitFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone,Vdm_Leg,sVdm_Leg,NodeType
USE MOD_ChangeBasis       ,ONLY:ChangeBasis3D
USE MOD_ReadInTools       ,ONLY:GETINT,GETREAL,GETREALARRAY,GETLOGICAL,GETINTFROMSTR
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_IO_HDF5           ,ONLY:AddToElemData,ElementOut
#if EQNSYSNR==2
USE MOD_Interpolation_Vars,ONLY:wGP
USE MOD_Mesh_Vars         ,ONLY:nElems,sJ
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iDeg
#if EQNSYSNR==2
INTEGER             :: iElem,i,j,k
#endif
!==================================================================================================================================
IF(FilterInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
   CALL CollectiveStop(__STAMP__,'InitFilter not ready to be called or already called.')
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT FILTER...'

FilterType = GETINTFROMSTR('FilterType')

! set the filterfunction pointer to the cut-off filter, even if the filter itself is not active
! necessary for the overintegration, which also uses the filter_pointer
NFilter = PP_N
Filter_Pointer=>Filter
IF(FilterType.GT.0) THEN
  ALLOCATE(FilterMat(0:PP_N,0:PP_N))
  FilterMat = 0.
  SELECT CASE (FilterType)
  CASE (FILTERTYPE_CUTOFF) ! Modal Filter cut-off 
    NFilter = GETINT('NFilter')
    DO iDeg=0,NFilter
      FilterMat(iDeg,iDeg) = 1.
    END DO
  CASE (FILTERTYPE_MODAL) ! Modal Filter (e.g. Hesthaven book)
    ! Read in modal filter parameter
    HestFilterParam = GETREALARRAY('HestFilterParam',3,'(/36.,12.,1./)')
    CALL HestFilter()
#if EQNSYSNR==2
  CASE (FILTERTYPE_LAF) ! Modal Filter cut-off, adaptive (LAF), only Euler/Navier-Stokes
    NFilter = GETINT('NFilter')
    LAF_alpha= GETREAL('LAF_alpha','1.0')
    ! LAF uses a special filter routine
    Filter_Pointer=>Filter_LAF
    DO iDeg=0,NFilter
      FilterMat(iDeg,iDeg) = 1.
    END DO
    ALLOCATE(J_N(0:PP_N,0:PP_N,0:PP_NZ))
    ALLOCATE(lim(nElems))
    ALLOCATE(eRatio(nElems))
    ALLOCATE(r(nElems))
    ALLOCATE(ekin_avg_old(nElems))
    ALLOCATE(ekin_fluc_avg_old(nElems))
    ALLOCATE(Vol(nElems))
    ALLOCATE(Integrationweight(0:PP_N,0:PP_N,0:PP_NZ,nElems))
    CALL AddToElemData(ElementOut,'LAF_eRatio',RealArray=eRatio)
    CALL AddToElemData(ElementOut,'LAF_lim'   ,RealArray=lim)
    CALL AddToElemData(ElementOut,'LAF_r'     ,RealArray=r)
    Vol = 0.
    DO iElem=1,nElems
      J_N(0:PP_N,0:PP_N,0:PP_NZ)=1./sJ(:,:,:,iElem,0)
      DO k=0,PP_NZ
        DO j=0,PP_N
          DO i=0,PP_N
#if PP_dim == 3
            IntegrationWeight(i,j,k,iElem) = wGP(i)*wGP(j)*wGP(k)*J_N(i,j,k)
#else
            IntegrationWeight(i,j,k,iElem) = wGP(i)*wGP(j)*J_N(i,j,k)
#endif
            Vol(iElem) = Vol(iElem) + IntegrationWeight(i,j,k,iElem)
          END DO ! i
        END DO ! j
      END DO ! k
    END DO !iElem
    DEALLOCATE(J_N)

  ! Compute normalization for LAF 
    normMod=((REAL(NFilter)+1)**(-2./3.)-2.**(-2./3.))/(REAL(PP_N+1)**(-2./3.)-REAL(NFilter+1)**(-2./3.))
    lim=1E-8
    eRatio=0.
    r=0.
    ekin_avg_old=1.E-16
    ekin_fluc_avg_old=1.E-16
#endif /*EQNSYSNR==2*/

  CASE DEFAULT 
    CALL CollectiveStop(__STAMP__,&
                        "FilterType unknown!")
  END SELECT

  !INFO
  SWRITE(*,'(A)',ADVANCE='NO')'FILTER DIAGONAL: '
  DO iDeg=0,PP_N-1
    SWRITE(*,'(F7.3)',ADVANCE='NO')FilterMat(iDeg,iDeg)
  END DO
  SWRITE(*,'(F7.3)')FilterMat(PP_N,PP_N)

  ! Assemble filter matrix in nodal space
  FilterMat=MATMUL(MATMUL(Vdm_Leg,FilterMat),sVdm_Leg)
END IF !FilterType=0

FilterInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT FILTER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitFilter


!==================================================================================================================================
!> Builds filter matrix for modal filter function. Could be used to controll aliasing instabilities.
!> For details see NODALBOOK. The magnitude of the modal DOF are reduced, where DOF belonging to higher order are reduced more.
!> The first DOF (=mean value) is NOT reduced to keep conservation. It is also possible to combine the Filter with a modal-based
!> indicator (Resolution/Persson indicator), to keep accuracy in resolved regions.
!==================================================================================================================================
SUBROUTINE HestFilter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars,ONLY:HestFilterParam,FilterMat
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: alpha,s,eta,etac ! etac is the modal cutoff. Example:
                                        !   filtering all DOF from x ORDER to highest ORDER: alpha<0,s arbitrary, param(3)=x-1
                                        ! good start set of parameters: alpha=36, s=12, etac=1
                                        !   HestFilterParam=(/36.,12.,1./) for ini file
                                        !   Default is HestFilterParam(1)=0.
INTEGER             :: iDeg
!==================================================================================================================================
alpha = HestFilterParam(1)
s     = HestFilterParam(2)
etac  = HestFilterParam(3)/REAL(PP_N+1)

FilterMat = 0.
DO iDeg=0,MIN(INT(HestFilterParam(3))-1,PP_N)
  FilterMat(iDeg,iDeg) = 1.
END DO
IF(alpha.GE.0.) THEN
  DO iDeg=INT(HestFilterParam(3)),PP_N
    eta = REAL(iDeg+1)/REAL(PP_N+1)
    FilterMat(iDeg,iDeg) = EXP(-alpha*((eta-etac)/(1.-etac))**s)
  END DO
END IF
END SUBROUTINE HestFilter



!==================================================================================================================================
!> interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
!> to another 3D tensor product node positions (number of nodes N_out+1)
!> defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!>  xi is defined in the 1DrefElem xi=[-1,1]
!==================================================================================================================================
PURE SUBROUTINE Filter(U_in,FilterMat)
! MODULES
USE MOD_PreProc
USE MOD_ChangeBasisByDim,  ONLY: ChangeBasisVolume
USE MOD_Mesh_Vars,         ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars,           ONLY: FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution vector to be filtered
REAL,INTENT(IN)     :: FilterMat(0:PP_N,0:PP_N)                  !< filter matrix to be used
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem
!==================================================================================================================================
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) CYCLE ! Do only, when DG element
#endif
  CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FilterMat,U_in(:,:,:,:,iElem))
END DO ! iElem
END SUBROUTINE Filter




#if EQNSYSNR==2
!===============================================================================================================================
!> LAF implementation via filter (only for Euler/Navier-Stokes)
!===============================================================================================================================
SUBROUTINE Filter_LAF(U_in,FilterMat)
! MODULES
USE MOD_PreProc
USE MOD_Filter_Vars,ONLY: lim,eRatio,r,ekin_avg_old,Vol,normMod,IntegrationWeight,ekin_fluc_avg_old,LAF_alpha
USE MOD_Mesh_Vars,  ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars,    ONLY: FV_Elems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: FilterMat(0:PP_N,0:PP_N)                  !< filter matrix to be used
!-------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution vector to be filtered
!-------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                      :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ):: U_Xi, U_filtered
#if PP_dim == 3
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ):: U_Eta
#endif
REAL,DIMENSION(1:3,0:PP_N,0:PP_N,0:PP_NZ)    :: U_fluc
REAL                                         :: ekin_avg,U_avg(5)
REAL                                         :: ekin_fluc_avg,dedt_lim
!===============================================================================================================================
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) CYCLE ! Do only, when DG element
#endif
  U_Xi = 0.
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_Xi(:,i,j,k)       = U_Xi(:,i,j,k)       + FilterMat(i,l)*U_in(:,l,j,k,iElem)
        END DO 
      END DO !i
    END DO !j
  END DO !k
#if PP_dim == 2
  U_filtered= 0.
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_filtered(:,i,j,k)      = U_filtered(:,i,j,k)      + FilterMat(j,l)*U_Xi(:,i,l,k) 
        END DO !l
      END DO !i
    END DO !j
  END DO !k
#else
  U_Eta= 0.
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_Eta(:,i,j,k)      = U_Eta(:,i,j,k)      + FilterMat(j,l)*U_Xi(:,i,l,k) 
        END DO !l
      END DO !i
    END DO !j
  END DO !k
  ! We need the filtered (low mode solution) and the unfiltered (full mode solution)
  U_filtered=0.
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_filtered(:,i,j,k) = U_filtered(:,i,j,k) + FilterMat(k,l)*U_Eta(:,i,j,l)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
#endif
  
  ! Compute the small scale contribution: u_small=u_full - u_large
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
#if PP_dim==3
        u_fluc(1:3,i,j,k) = (U_in(2:4,i,j,k,iElem)-U_filtered(2:4,i,j,k))
#else
        u_fluc(1:2,i,j,k) = (U_in(2:3,i,j,k,iElem)-U_filtered(2:3,i,j,k))
#endif
      END DO !i
    END DO !j
  END DO !k

  ! Compute the average velocities 
  U_Avg=0.
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
#if PP_dim==3
        U_Avg(2:4)=U_Avg(2:4) + U_filtered(2:4,i,j,k)/U_in(1,i,j,k,iElem)*IntegrationWeight(i,j,k,iElem)
#else
        U_Avg(2:3)=U_Avg(2:3) + U_filtered(2:3,i,j,k)/U_in(1,i,j,k,iElem)*IntegrationWeight(i,j,k,iElem)
#endif
      END DO ! i
    END DO ! j
  END DO ! k
  U_Avg(:) = U_Avg(:)/Vol(iElem) !CELL average
  ekin_avg=0.
  ekin_fluc_avg=0.
  ! compute average and fluc average kinetic energy
  DO k=0,PP_NZ
    DO j=0,PP_N
      DO i=0,PP_N
        ekin_fluc_avg = ekin_fluc_avg + (u_fluc(1,i,j,k)*u_fluc(1,i,j,k)/U_in(1,i,j,k,iElem)**2+& 
#if PP_dim==3
                                         u_fluc(3,i,j,k)*u_fluc(3,i,j,k)/U_in(1,i,j,k,iElem)**2+&
#endif
                                         u_fluc(2,i,j,k)*u_fluc(2,i,j,k)/U_in(1,i,j,k,iElem)**2)*IntegrationWeight(i,j,k,iElem)
        ekin_avg=ekin_avg + ((U_filtered(2,i,j,k)-U_Avg(2))*(U_filtered(2,i,j,k)-U_Avg(2))/U_in(1,i,j,k,iElem)**2 & 
                            +(U_filtered(3,i,j,k)-U_Avg(3))*(U_filtered(3,i,j,k)-U_Avg(3))/U_in(1,i,j,k,iElem)**2 & 
#if PP_dim==3
                            +(U_filtered(4,i,j,k)-U_Avg(4))*(U_filtered(4,i,j,k)-U_Avg(4))/U_in(1,i,j,k,iElem)**2 &
#endif
                             )*IntegrationWeight(i,j,k,iElem)
          END DO ! i
        END DO ! j
      END DO ! k
      !Calculated normalized kinetic energy ratio 
      eRatio(iElem) = (ekin_fluc_avg/ekin_avg)*normMod
      !define threshhold value for eRatio by linear weighting over low mode energy change (low mode trigger LMT)
!      lim(iElem) = lim(iElem)*ekin_avg_old(iElem)/ekin_avg
      !limit may not become smaller than 1.0E-8 --> tends to 0 for from scratch computations with LMT
      lim(iElem) = MAX(lim(iElem),1.0E-8)
      !first order de/dt approx for flucs and large scales as ratio for filter cond.
      dedt_lim = ABS(ekin_avg - ekin_avg_old(iElem))/((ekin_fluc_avg - ekin_fluc_avg_old(iElem)))!*(12.-8.)/(8.-2.)


      IF ((dedt_lim .GT. 1.0) .AND. ((ekin_fluc_avg - ekin_fluc_avg_old(iElem)).GT.0.)) THEN
        lim(iElem) = eRatio(iElem)
!      ! ETB (energy transfer balance): Case rise of energy in fluctuations --> allowed if rise in low modes is
!      ! equal or higher, or drop in low modes is equal or higher (transfer)
!     ! By setting lim=eRatio, limit follows eRatio for valid states, also after filtering therefor a margin of the
!      ! limit is beneficial, see below
!      ! Cases with energy loss in fluctuations need not be considered --> limit LE Eratio
      ELSEIF ((ekin_fluc_avg - ekin_fluc_avg_old(iElem)).LT.0.) THEN
        lim(iElem) = eRatio(iElem)
!        ! energy drain in high modes is not dangerous 
      ELSEIF ((eRatio(iElem) .GT. lim(iElem)*LAF_alpha)) THEN 
!      IF ((eRatio(iElem) .GT. 1.0)) THEN ! For LAF_Fixed 
!      IF (MOD(r(iElem)+1,10000.).EQ.0) THEN ! For LAF_Clocked 
        ! Filter to NUnder whenever eRatio is to high, we allow for a small margin --> less margin > more filtering
        U_in(:,:,:,:,iElem) = U_filtered(:,:,:,:)
        r(iElem) = r(iElem)+1. 

      END IF
      ekin_avg_old(iElem)=ekin_avg
      ekin_fluc_avg_old(iElem)=ekin_fluc_avg
    END DO !iElem
END SUBROUTINE Filter_LAF
#endif /*EQNSYSNR==2*/



SUBROUTINE Filter_Selective(NVar,FilterMat,U_in,filter_ind)
! MODULES
USE MOD_PreProc
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(NVar,0:PP_N,0:PP_N,0:PP_N) ! < solution vector to be filtered
REAL,INTENT(IN)     :: FilterMat(0:PP_N,0:PP_N)        ! < filter matrix to be used
INTEGER,INTENT(IN)  :: NVar
LOGICAL, INTENT(IN) :: filter_ind(:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                   :: i,j,k,l
REAL,DIMENSION(NVar,0:PP_N,0:PP_N,0:PP_N) :: U_Xi,U_Eta
!==================================================================================================================================
! Perform filtering
#if FV_ENABLED
stop
#endif
IF(filter_ind(1)) THEN
  U_Xi = 0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_Xi(:,i,j,k) = U_Xi(:,i,j,k) + FilterMat(i,l)*U_in(:,l,j,k)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
ELSE
  U_Xi = U_in
END IF
IF(filter_ind(2)) THEN
  U_Eta= 0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_Eta(:,i,j,k) = U_Eta(:,i,j,k) + FilterMat(j,l)*U_Xi(:,i,l,k)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
ELSE
  U_Eta = U_Xi
END IF
IF(filter_ind(3)) THEN
  U_in(:,:,:,:)=0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO l=0,PP_N
          U_in(:,i,j,k) = U_in(:,i,j,k) + FilterMat(k,l)*U_Eta(:,i,j,l)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
ELSE
  U_in = U_Eta
END IF
END SUBROUTINE Filter_Selective

!==================================================================================================================================
!> Deallocate filter arrays
!==================================================================================================================================
SUBROUTINE FinalizeFilter()
! MODULES
USE MOD_Filter_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(FilterMat)
#if EQNSYSNR==2
SDEALLOCATE(lim)
SDEALLOCATE(eRatio)
SDEALLOCATE(r)
SDEALLOCATE(ekin_avg_old)
SDEALLOCATE(Integrationweight)
SDEALLOCATE(ekin_fluc_avg_old)
SDEALLOCATE(Vol)
#endif /*EQNSYSNR==2*/
FilterInitIsDone = .FALSE.
END SUBROUTINE FinalizeFilter

END MODULE MOD_Filter
