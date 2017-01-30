#include "flexi.h"

MODULE MOD_dynsmag
#ifdef EDDYVISCOSITY
!===================================================================================================================================
! Soubroutines necessary for calculating dynsmag Eddy-Viscosity
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Initdynsmag
  MODULE PROCEDURE Initdynsmag
END INTERFACE

INTERFACE dynsmag
  MODULE PROCEDURE dynsmag
END INTERFACE

INTERFACE dynsmag_surf
  MODULE PROCEDURE dynsmag_surf
END INTERFACE

INTERFACE compute_cd 
  MODULE PROCEDURE compute_cd 
END INTERFACE

INTERFACE Finalizedynsmag
  MODULE PROCEDURE Finalizedynsmag
END INTERFACE

PUBLIC::Initdynsmag,dynsmag,dynsmag_surf,compute_cd,Finalizedynsmag
INTEGER::N_testfilter
!===================================================================================================================================

CONTAINS

SUBROUTINE Initdynsmag()
!===================================================================================================================================
! Get some parameters needed by dynsmag modules and initialize dynsmags
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc                
USE MOD_EddyVisc_Vars              
USE MOD_ReadInTools       ,   ONLY:GETREAL,GETLOGICAL,GETINT
USE MOD_Interpolation_Vars,   ONLY:InterpolationInitIsDone,Vdm_Leg,sVdm_Leg,NodeType,wGP
USE MOD_Mesh_Vars         ,   ONLY:MeshInitIsDone
USE MOD_Interpolation_Vars,   ONLY:wGP
USE MOD_Mesh_Vars,            ONLY:sJ,nSides,nElems
USE MOD_Mesh_Vars,            ONLY:ElemToSide
USE MOD_Testcase_Vars,        ONLY:testcase
#if MPI
USE MOD_MPI,                  ONLY:StartReceiveMPIData,FinishExchangeMPIData,StartSendMPIData
USE MOD_MPI_Vars,             ONLY:MPIRequest_DeltaS,nNbProcs
#endif /*MPI*/ 
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iElem,j,k
REAL    :: CellVol
INTEGER :: iLocSide,SideID,FlipID
!===================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.dynsmagInitIsDone)THEN
   SWRITE(UNIT_StdOut,'(A)') "Initdynsmag not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT dynsmag...'

! Read the variables used for LES model
! Smagorinsky model
CS     = GETREAL('CS')
PrSGS  = GETREAL('PrSGS','0.7')
IF(testcase.EQ."channel") THEN
  ! Do Van Driest style damping or not
  VanDriest = GETLOGICAL('VanDriest','.FALSE.')
END IF

! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )

DO iElem=1,nElems                                        
  CellVol = 0.
  DO i=0,PP_N
    DO j=0,PP_N
      DO k=0,PP_N
        CellVol = CellVol +wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
      END DO
    END DO
  END DO
  DeltaS(iElem) = ( CellVol)**(1./3.)  / (REAL(PP_N)+1.)
  DO iLocSide=1,6
     SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
     FlipID=ElemToSide(E2S_FLIP,iLocSide,iElem) 
     IF(FlipID.EQ.0) THEN
       DeltaS_master(SideID)=DeltaS(iElem)
     ELSE
       DeltaS_slave(SideID)=DeltaS(iElem)
     END IF
  END DO
END DO
#if MPI
! Send YOUR - receive MINE
CALL StartReceiveMPIData(DeltaS_slave, 1, 1,nSides,MPIRequest_DeltaS( :,SEND),SendID=1)
CALL StartSendMPIData(   DeltaS_slave, 1, 1,nSides,MPIRequest_DeltaS( :,RECV),SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_DeltaS ) !Send MINE -receive YOUR
#endif /*MPI*/
N_Testfilter = GETINT('N_Testfilter')
DO i=0,N_Testfilter
  FilterMat_Testfilter(i,i) = 1.
END DO
!FilterMat_Testfilter(N_testfilter+1,N_testfilter+1) = 0.5
SWRITE(*,'(A)',ADVANCE='NO')'TEST FILTER, FILTER DIAGONAL: '
DO i=0,PP_N
  SWRITE(*,'(F7.3)',ADVANCE='NO')FilterMat_Testfilter(i,i)
END DO
FilterMat_Testfilter=MATMUL(MATMUL(Vdm_Leg,FilterMat_Testfilter),sVdm_Leg)
dynsmagInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT dynsmag DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE Initdynsmag

SUBROUTINE dynsmag(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,iElem,i,j,k,muSGS)
!===================================================================================================================================
!Compute dynsmag Eddy-Visosity at a given point
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EddyVisc_Vars,     ONLY:deltaS,CS,VanDriest,SGS_Ind,muSGSmax
USE MOD_Mesh_Vars,         ONLY:Elem_xGP

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem,i,j,k
REAL,INTENT(IN)                           :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho
REAL,INTENT(INOUT)                          :: muSGS 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: S_eN
REAL                :: divv
!===================================================================================================================================
! die zwei aus der Wurzel gleich hier oben verarbeitet, spart eine Operation
divv    = grad11+grad22+grad33
S_eN = 2*((grad11-2./3.*divv)**2. + (grad22-2./3.*divv)**2. + (grad33-2./3.*divv)**2.)
S_eN = S_eN + ( grad12 + grad21 )**2.
S_eN = S_eN + ( grad13 + grad31 )**2.
S_eN = S_eN + ( grad23 + grad32 )**2.
S_eN = sqrt(S_eN)
! dynsmag model
muSGS = S_eN*rho*SGS_Ind(1,i,j,k,iElem) 
SGS_Ind(2,i,j,k,iElem) = muSGS 
muSGSmax(iElem) = MAX(muSGS,muSGSmax(iElem))
END SUBROUTINE dynsmag

SUBROUTINE dynsmag_surf(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,DeltaSS,cd,muSGS,Face_xGP)
!===================================================================================================================================
!Compute dynsmag Eddy-Visosity at a given point
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EddyVisc_Vars,     ONLY:CS

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                           :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho
REAL,INTENT(IN)                           :: DeltaSS 
REAL,INTENT(IN)                           :: cd,Face_xGP
REAL,INTENT(OUT)                          :: muSGS 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: S_eN
REAL                :: divv
!===================================================================================================================================
divv    = grad11+grad22+grad33
S_eN = 2*((grad11-2./3.*divv)**2. + (grad22-2./3.*divv)**2. + (grad33-2./3.*divv)**2.)
S_eN = S_eN + ( grad12 + grad21 )**2.
S_eN = S_eN + ( grad13 + grad31 )**2.
S_eN = S_eN + ( grad23 + grad32 )**2.
S_eN = sqrt(S_eN)
! dynsmag model
muSGS= cd  * S_eN*rho
END SUBROUTINE dynsmag_surf

SUBROUTINE compute_cd(U_in)
!===============================================================================================================================
!compute TKE indicator needed for the modification of Smagorinskys viscosity.
!===============================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_EddyVisc_Vars,          ONLY:SGS_Ind,FilterMat_testfilter
USE MOD_EddyVisc_Vars,          ONLY:SGS_Ind!,SGS_Ind_Slave,SGS_Ind_Master
USE MOD_EddyVisc_Vars,          ONLY:MM_Avg, ML_Avg
!USE MOD_ProlongToFace1,         ONLY: ProlongToFace1
USE MOD_Filter,                 ONLY:Filter_General
USE MOD_Lifting_Vars,           ONLY:gradUx,gradUy,gradUz
USE MOD_Interpolation_Vars,     ONLY:wGP
USE MOD_Mesh_Vars,              ONLY:sJ, nSides, nElems
USE MOD_Timedisc_Vars,          ONLY:dt
!USE MOD_Interpolation_Vars,     ONLY: L_Minus,L_Plus
!#if MPI
!USE MOD_MPI_Vars
!USE MOD_MPI,                ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
!#endif

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!-------------------------------------------------------------------------------------------------------------------------------
! aOCAL VARIABLES 
INTEGER                                      :: i,j,k,l,m,iElem
REAL,DIMENSION(3,0:PP_N,0:PP_N,0:PP_N)       :: V_filtered
REAL                                         :: srho
REAL                :: v1,v2,v3
REAL,DIMENSION(3,3,0:PP_N,0:PP_N,0:PP_N) :: S_ik,M_ik,L_ik
REAL                                     :: D_Ratio
REAL, DIMENSION(0:PP_N,0:PP_N,0:PP_N)    :: C_d,MM,ML,S_eN_filtered,S_eN
REAL, DIMENSION(0:PP_N)                  :: dummy1D
REAL                                     :: Vol, Integrationweight,dummy
REAL,DIMENSION(1:3,1:3,0:PP_N,0:PP_N,0:PP_N)         :: gradv_all
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N)             :: gradv11_elem,gradv12_elem,gradv13_elem 
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N)             :: gradv21_elem,gradv22_elem,gradv23_elem
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N)             :: gradv31_elem,gradv32_elem,gradv33_elem
!===============================================================================================================================
DO iElem=1,nElems
    DO k=0,PP_N;  DO j=0,PP_N; DO i=0,PP_N
      srho = 1. / U_in(1,i,j,k,iElem) ! 1/rho
      v1   = U_in(2,i,j,k,iElem)*srho
      v2   = U_in(3,i,j,k,iElem)*srho
      v3   = U_in(4,i,j,k,iElem)*srho
      gradv11_elem(i,j,k) = srho*(gradUx(2,i,j,k,iElem) - v1*gradUx(1,i,j,k,iElem))
      gradv21_elem(i,j,k) = srho*(gradUx(3,i,j,k,iElem) - v2*gradUx(1,i,j,k,iElem))
      gradv31_elem(i,j,k) = srho*(gradUx(4,i,j,k,iElem) - v3*gradUx(1,i,j,k,iElem))
      gradv12_elem(i,j,k) = srho*(gradUy(2,i,j,k,iElem) - v1*gradUy(1,i,j,k,iElem))
      gradv22_elem(i,j,k) = srho*(gradUy(3,i,j,k,iElem) - v2*gradUy(1,i,j,k,iElem))
      gradv32_elem(i,j,k) = srho*(gradUy(4,i,j,k,iElem) - v3*gradUy(1,i,j,k,iElem))
      gradv13_elem(i,j,k) = srho*(gradUz(2,i,j,k,iElem) - v1*gradUz(1,i,j,k,iElem))
      gradv23_elem(i,j,k) = srho*(gradUz(3,i,j,k,iElem) - v2*gradUz(1,i,j,k,iElem))
      gradv33_elem(i,j,k) = srho*(gradUz(4,i,j,k,iElem) - v3*gradUz(1,i,j,k,iElem))
  END DO; END DO; END DO ! i,j,k
  DO k=0,PP_N;  DO j=0,PP_N; DO i=0,PP_N
      ! die zwei aus der Wurzel gleich hier oben verarbeitet, spart eine Operation
      S_eN(i,j,k)= 2*(gradv11_elem(i,j,k)**2. + gradv22_elem(i,j,k)**2. + gradv33_elem(i,j,k)**2.)
      S_eN(i,j,k)= S_eN(i,j,k) + ( gradv12_elem(i,j,k) + gradv21_elem(i,j,k) )**2.
      S_eN(i,j,k)= S_eN(i,j,k) + ( gradv13_elem(i,j,k) + gradv31_elem(i,j,k) )**2.
      S_eN(i,j,k)= S_eN(i,j,k) + ( gradv23_elem(i,j,k) + gradv32_elem(i,j,k) )**2.
      S_eN(i,j,k)= sqrt( S_eN(i,j,k) )
  END DO; END DO; END DO ! i,j,k
  ! dynamic Smagorinsky model
  DO i=1,3
    V_Filtered(i,:,:,:) = U_in(i+1,:,:,:,iElem)/U_in(1,:,:,:,iElem)
  END DO !i

  !Filter velocities
  CALL Filter_General(3,FilterMat_testfilter,V_filtered) 
  !           _ _   __
  !Compute L=-u u + uu  !!TENSOR ik
  DO i=1,3
    DO k=1,3
      L_ik(i,k,:,:,:) = U_in(i+1,:,:,:,iElem)*U_in(k+1,:,:,:,iElem)/U_in(1,:,:,:,iElem)**2
    END DO ! k
    CALL Filter_General(3,FilterMat_testfilter,L_ik(i,1:3,:,:,:)) 
  END DO ! i
  DO i=1,3
    DO k=1,3
      L_ik(i,k,:,:,:) = L_ik(i,k,:,:,:) - V_filtered(i,:,:,:)*V_filtered(k,:,:,:)
    END DO ! k
  END DO ! i
  !           _____            _ _   ____
  !Compute M=(Delta/Delta)**2*|S|S - |S|S  !!TENSOR ik
  S_ik(1,1,:,:,:) = gradv11_elem(:,:,:)
  S_ik(2,2,:,:,:) = gradv22_elem(:,:,:)
  S_ik(3,3,:,:,:) = gradv33_elem(:,:,:)
  S_ik(1,2,:,:,:) = (gradv21_elem(:,:,:)+gradv12_elem(:,:,:))*0.5
  S_ik(1,3,:,:,:) = (gradv31_elem(:,:,:)+gradv13_elem(:,:,:))*0.5
  S_ik(2,3,:,:,:) = (gradv23_elem(:,:,:)+gradv32_elem(:,:,:))*0.5
  S_ik(2,1,:,:,:) = S_ik(1,2,:,:,:)
  S_ik(3,1,:,:,:) = S_ik(1,3,:,:,:)
  S_ik(3,2,:,:,:) = S_ik(2,3,:,:,:)

  DO i=1,3
    DO k=1,3
      M_ik(i,k,:,:,:) = S_eN(:,:,:)*S_ik(i,k,:,:,:)
    END DO ! k
    CALL Filter_General(3,FilterMat_testfilter,M_ik(i,1:3,:,:,:)) 
  END DO ! i

  !filtered gradients

  gradv_all(1,1,:,:,:) = gradv11_elem(:,:,:) 
  gradv_all(1,2,:,:,:) = gradv12_elem(:,:,:)
  gradv_all(1,3,:,:,:) = gradv13_elem(:,:,:)
  gradv_all(2,1,:,:,:) = gradv21_elem(:,:,:)
  gradv_all(2,2,:,:,:) = gradv22_elem(:,:,:)
  gradv_all(2,3,:,:,:) = gradv23_elem(:,:,:)
  gradv_all(3,1,:,:,:) = gradv31_elem(:,:,:)
  gradv_all(3,2,:,:,:) = gradv32_elem(:,:,:)
  gradv_all(3,3,:,:,:) = gradv33_elem(:,:,:)
  CALL Filter_General(3,FilterMat_Testfilter,gradv_all(:,1,:,:,:))
  CALL Filter_General(3,FilterMat_Testfilter,gradv_all(:,2,:,:,:))
  CALL Filter_General(3,FilterMat_Testfilter,gradv_all(:,3,:,:,:))
  DO i=1,3
    DO j=1,3
      S_ik(i,j,:,:,:) = (gradv_all(i,j,:,:,:)+gradv_all(j,i,:,:,:))*0.5
    END DO
  END DO

  S_eN_filtered(:,:,:)= (S_ik(1,1,:,:,:)**2 + S_ik(2,2,:,:,:)**2. + S_ik(3,3,:,:,:)**2.)
  S_eN_filtered(:,:,:)= S_eN_filtered(:,:,:) +  2*(S_ik(2,1,:,:,:) )**2.
  S_eN_filtered(:,:,:)= S_eN_filtered(:,:,:) +  2*(S_ik(3,1,:,:,:) )**2.
  S_eN_filtered(:,:,:)= S_eN_filtered(:,:,:) +  2*(S_ik(3,2,:,:,:) )**2.
  S_eN_filtered(:,:,:)= sqrt(2* S_eN_filtered(:,:,:) )
  !D_Ratio=(REAL(N_testfilter+1)/REAL(PP_N+1))**2
  D_Ratio=((1./REAL(N_testfilter+1)) / (1./REAL(PP_N+1)))**2
  DO i=1,3
    DO k=1,3
      M_ik(i,k,:,:,:) = M_ik(i,k,:,:,:) - D_Ratio * S_eN_filtered(:,:,:)*S_ik(i,k,:,:,:)
    END DO ! k
  END DO ! i
  !
  !Compute C_d**2 =1/2* M_ik*L_ik / M_ik*M_ik  !!SUM OVER ik !!!for muSGS=C_d**2|S|
  MM=0.
  ML=0.
  DO i=0,PP_N
    DO j=0,PP_N
      DO k=0,PP_N
        DO l=1,3
          DO m=1,3
            MM(i,j,k)  = MM(i,j,k) + (M_ik(l,m,i,j,k)*M_ik(l,m,i,j,k))
            ML(i,j,k)  = ML(i,j,k) + (M_ik(l,m,i,j,k)*L_ik(l,m,i,j,k))
          END DO ! k
        END DO ! j
      END DO ! k
    END DO ! j
  END DO ! i


  !Average C_d**2 over homogeneous directions, or cell, or in time, or whatever
!  !-----CELL AVERAGE
!  dummy=0.
!  Vol = 0.
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        IntegrationWeight = wGP(i)*wGP(j)*wGP(k)*1/sJ(i,j,k,iElem)
!        dummy=dummy + C_d(i,j,k)*IntegrationWeight
!        Vol = Vol + Integrationweight
!      END DO ! i
!    END DO ! j
!  END DO ! k
!  SGS_Ind(1,:,:,:,iElem) =dummy/Vol !CELL average
!  SGS_Ind(1,:,:,:,iElem) =0.001
!print*,SGS_Ind(1,:,:,:,iElem)

!  !-----CELL LOCAL 
!  dummy=0.
!  Vol = 0.
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        IntegrationWeight = wGP(i)*wGP(j)
!        dummy = dummy + ML(i,j,k)*IntegrationWeight
!        Vol = Vol + Integrationweight
!      END DO ! i
!    END DO ! j
!  END DO ! j
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        ML(i,j,k) = dummy/Vol !CELL average
!      END DO ! i
!    END DO ! j
!  END DO ! j
!  dummy=0.
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        IntegrationWeight = wGP(i)*wGP(j)
!        dummy = dummy + MM(i,j,k)*IntegrationWeight
!      END DO ! i
!    END DO ! j
!  END DO ! j
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        MM(i,j,k) = dummy/Vol !CELL average
!      END DO ! i
!    END DO ! j
!  END DO ! j
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        C_d(i,k,j) = -0.5*ML(i,j,k)/MM(i,j,k) !CELL average
!        IF (ABS(MM(i,j,k)) .LE. 1e-15) C_d(i,j,k)=0. 
!      END DO ! i
!    END DO ! j
!  END DO ! j
!!  SGS_Ind(1,:,:,:,iElem) = MIN(C_d(:,:,:),CS***2)
!  SGS_Ind(1,:,:,:,iElem) = C_d(:,:,:)
!
  !-----CELL LOCAL IN HOMOGENEOUS DIRECTIONS
  !ACHTUNG NUR KARTESISCH UND IN/FUER Y=J!!!
!  dummy1D=0.
!  Vol = 0.
!  DO j=0,PP_N
!    DO i=0,PP_N
!      IntegrationWeight = wGP(i)*wGP(j)
!      dummy1D(:) = dummy1D(:) + ML(i,:,j)*IntegrationWeight
!      Vol = Vol + Integrationweight
!    END DO ! i
!  END DO ! j
!  DO j=0,PP_N
!    DO i=0,PP_N
!      ML(i,:,j) = dummy1D(:)/Vol !CELL average
!    END DO ! i
!  END DO ! j
!  dummy1D=0.
!  DO j=0,PP_N
!    DO i=0,PP_N
!      IntegrationWeight = wGP(i)*wGP(j)
!      dummy1D(:) = dummy1D(:) + MM(i,:,j)*IntegrationWeight
!    END DO ! i
!  END DO ! j
!  DO j=0,PP_N
!    DO i=0,PP_N
!      MM(i,:,j) = dummy1D(:)/Vol !CELL average
!    END DO ! i
!  END DO ! j
!  DO j=0,PP_N
!    DO i=0,PP_N
!      C_d(i,:,j) = - 0.5*ML(i,:,j)/MM(i,:,j) !CELL average
!      !C_d(i,:,j) = + 0.5*ML(i,:,j)/MM(i,:,j) !CELL average
!    END DO ! i
!  END DO ! j
!  SGS_Ind(1,:,:,:,iElem) = C_d(:,:,:)
!
!------------------------------------------
!Time Average ala pruett
!Im Nenner steht r=FilterWidth/dt(RK)
!Simple shot is taking the filter width equal 1, assuming meaningful dimensionless timescale
ML_Avg(:,:,:,iElem) = ML_Avg(:,:,:,iElem) + (ML(:,:,:)-ML_Avg(:,:,:,iElem))*(dt/1.)
MM_Avg(:,:,:,iElem) = MM_Avg(:,:,:,iElem) + (MM(:,:,:)-MM_Avg(:,:,:,iElem))*(dt/1.)
!C_d=0.
DO i=0,PP_N
  DO j=0,PP_N
    DO k=0,PP_N
!      IF (ABS(MM_Avg(i,j,k,iElem)) .LE. 1e-15) CYCLE 
      C_d(i,j,k) =  0.5*ML_Avg(i,j,k,iElem)/MM_Avg(i,j,k,iElem)
    END DO ! i
  END DO ! j
END DO ! k
SGS_Ind(1,:,:,:,iElem) = SGS_Ind(1,:,:,:,iElem) + (C_d(:,:,:)-SGS_Ind(1,:,:,:,iElem))/1.*dt
!print*,SGS_Ind(1,:,:,:,iElem)
!  SGS_Ind(1,:,:,:,iElem) = C_d(:,:,:)
!------------------------------------------
END DO
!#if MPI
!! 4.2)
!CALL StartReceiveMPIData(SGS_Ind_master,DataSizeSideScalar,1,nSides,MPIRequest_SGS_Ind(:,SEND),SendID=1)
!                                                         ! Receive YOUR / FV_surf_gradU_slave: master -> slave
!CALL ProlongToFace1(PP_N,SGS_Ind(1:1,:,:,:,:),SGS_Ind_master(:,:,:,:),SGS_Ind_Slave(:,:,:,:),L_Minus,L_Plus,.TRUE.)
!CALL StartSendMPIData(SGS_Ind_master,DataSizeSideScalar,1,nSides,MPIRequest_SGS_Ind(:,RECV),SendID=1)
!#endif
!! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
!CALL ProlongToFace1(PP_N,SGS_Ind(1:1,:,:,:,:),SGS_Ind_master(:,:,:,:),SGS_Ind_Slave(:,:,:,:),L_Minus,L_Plus,.FALSE.)
!#if MPI  
!CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_SGS_Ind)  ! U_slave: slave -> master 
!#endif
END SUBROUTINE compute_cd 

SUBROUTINE Finalizedynsmag()
!===============================================================================================================================
! Get the constant advection velocity vector from the ini file
!===============================================================================================================================
! MODULES
USE MOD_EddyVisc_Vars,ONLY:dynsmagInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===============================================================================================================================
dynsmagInitIsDone = .FALSE.
END SUBROUTINE Finalizedynsmag
#endif
END MODULE MOD_dynsmag
