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

INTERFACE compute_cd 
  MODULE PROCEDURE compute_cd 
END INTERFACE

INTERFACE Finalizedynsmag
  MODULE PROCEDURE Finalizedynsmag
END INTERFACE

PUBLIC::Initdynsmag,dynsmag,compute_cd,Finalizedynsmag
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
USE MOD_ReadInTools       ,   ONLY:GETREAL,GETLOGICAL,GETINT, GETSTR
USE MOD_Mesh_Vars         ,   ONLY:MeshInitIsDone, Elem_xgp, offsetElem
USE MOD_Interpolation_Vars,   ONLY:InterpolationInitIsDone,Vdm_Leg,sVdm_Leg,NodeType,wGP,NodeTypeCL
USE MOD_Interpolation,        ONLY:GetVandermonde
USE MOD_ChangeBasis,          ONLY:changeBasis3D
USE MOD_Interpolation_Vars,   ONLY:wGP
USE MOD_Mesh_Vars,            ONLY:sJ,nElems
USE MOD_Mesh_Vars,            ONLY:dXCL_N
USE MOD_Mesh_Vars,            ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Testcase_Vars,        ONLY:testcase
USE MOD_HDF5_input,           ONLY: OpenDataFile,CloseDataFile,ReadArray
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: isn
INTEGER :: i,j,k,iElem
REAL    :: WallDistTrsh, distNorm
REAL    :: vec(3)
REAL    :: EwallDist(3,nElems)
REAL    :: Vdm_CLN_N(         0:PP_N   ,0:PP_N)
REAL    :: dX_N(3, 0:PP_N,  0:PP_N   ,0:PP_N)
!===================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.dynsmagInitIsDone)THEN
   SWRITE(UNIT_StdOut,'(A)') "Initdynsmag not ready to be called or already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT dynsmag...'

! Read the variables used for LES model
! Smagorinsky model
PrSGS  = GETREAL('PrSGS','0.7')
!get wall distance for zonal averaging approach if needed
WallDistFile = GETSTR('WallDistFile','')
WallDistTrsh = GETREAL('WallDistanceTreshold','0.0')
IF (WallDistFile .NE. '') THEN
  CALL OpenDataFile(WallDistFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadArray('Walldistance',2,(/3,nElems/),OffsetElem,2,RealArray=EwallDist)
  CALL CloseDataFile()
ENDIF


!Find in which x, y, z direction are the i, j ,k index pointing, and
!then decide which index to filter
!MATTEO: DEBUG
filtdir_out = 0.0
DO iElem=1,nElems
  distNorm = SQRT(SUM(EwallDist(:,iElem)**2)) 
  !MATTEO: debug
  walldist_out(iElem) = distNorm
  walldist_x(iElem) = EwallDist(1,iElem)
  walldist_y(iElem) = EwallDist(2,iElem)
  walldist_z(iElem) = EwallDist(3,iElem)
  IF (distNorm .GE. WallDistTrsh) THEN !out of the boundary layer
    ! Filter and average isotropically
    filter_ind(:,iElem) = .true.
    average_ind(:,iElem) = .true.
    !MATTEO: DEBUG
    filtdir_out(iElem) = 8.0
  ELSE ! in the boundary layer
    ! Filter and average only the plane normal to the boundary                    
    !i, average on the four edges                                       
    vec = Elem_xgp(:,PP_N,0,0,iElem) - Elem_xgp(:,0,0,0,iElem)          
    vec = vec + (Elem_xgp(:,PP_N,PP_N,0,iElem) - Elem_xgp(:,0,PP_N,0,iElem))
    vec = vec + (Elem_xgp(:,PP_N,0,PP_N,iElem) - Elem_xgp(:,0,0,PP_N,iElem))
    vec = vec + (Elem_xgp(:,PP_N,PP_N,PP_N,iElem) - Elem_xgp(:,0,PP_N,PP_N,iElem))
    vec = vec/4.0                                                       
    vec = vec/sqrt(sum(vec**2)) !Unit vector of i direction             
    !If i is normal to the distance filter and average in i             
    isn = ISNORMAL(vec,EwallDist(:,iElem))                              
    filter_ind(1,iElem) = isn                                           
    average_ind(1,iElem) = isn                                          
    !MATTEO: DEBUG                                                      
    IF(isn) filtdir_out(iElem) = filtdir_out(iElem)+1                   
    !j                                                                  
    vec = Elem_xgp(:,0,PP_N,0,iElem) - Elem_xgp(:,0,0,0,iElem)          
    vec = vec + (Elem_xgp(:,PP_N,PP_N,0,iElem) - Elem_xgp(:,PP_N,0,0,iElem))
    vec = vec + (Elem_xgp(:,0,PP_N,PP_N,iElem) - Elem_xgp(:,0,0,PP_N,iElem))
    vec = vec + (Elem_xgp(:,PP_N,PP_N,PP_N,iElem) - Elem_xgp(:,PP_N,0,PP_N,iElem))
    vec = vec/4.0                                                       
    vec = vec/sqrt(sum(vec**2))                                         
    !If j is normal to the distance filter and average in j             
    isn = ISNORMAL(vec,EwallDist(:,iElem))                              
    filter_ind(2,iElem) =  isn                                          
    average_ind(2,iElem) = isn                                          
    !MATTEO: DEBUG                                                      
    IF(isn) filtdir_out(iElem) = filtdir_out(iElem)+2                   
    !k                                                                  
    vec = Elem_xgp(:,0,0,PP_N,iElem) - Elem_xgp(:,0,0,0,iElem)          
    vec = vec + (Elem_xgp(:,PP_N,0,PP_N,iElem) - Elem_xgp(:,PP_N,0,0,iElem))
    vec = vec + (Elem_xgp(:,0,PP_N,PP_N,iElem) - Elem_xgp(:,0,PP_N,0,iElem))
    vec = vec + (Elem_xgp(:,PP_N,PP_N,PP_N,iElem) - Elem_xgp(:,PP_N,PP_N,0,iElem))
    vec = vec/4.0                                                       
    vec = vec/sqrt(sum(vec**2))                                         
    !If k is normal to the distance filter and average in k             
    isn = ISNORMAL(vec,EwallDist(:,iElem))                                                                           
    filter_ind(3,iElem) = isn                                           
    average_ind(3,iElem) = isn                                          
    !MATTEO: DEBUG                                                      
    IF(isn) filtdir_out(iElem) = filtdir_out(iElem)+4
  END IF !in or out boundary layer
END DO !iElem

average_type = 0
!build integration weights for each of possibly 7 average types
!sort average_ind to average_type for select cases
CALL GetVandermonde(    PP_N   , NodeTypeCL  , PP_N    , NodeType,   Vdm_CLN_N         , modal=.FALSE.)
DO iElem=1,nElems
  IF     (average_ind(1,iElem).AND.average_ind(2,iElem).AND.average_ind(3,iElem)) THEN !volume
    IntElem(:,:,:,iElem) = 1./sJ(:,:,:,iElem,0)
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(i)*wGP(j)*wGP(k)
    END DO; END DO; END DO
    average_type(iElem) = 1
  ELSEIF ((.NOT.average_ind(1,iElem)).AND.average_ind(2,iElem).AND.average_ind(3,iElem)) THEN ! I-Plane 
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      IntElem(i,j,k,iElem) = SQRT(SUM(Metrics_ftilde(:,i,j,k,iElem,0)**2))
      IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(j)*wGP(k)
    END DO; END DO; END DO
    average_type(iElem) = 2
  ELSEIF (average_ind(1,iElem).AND.(.NOT.average_ind(2,iElem)).AND.average_ind(3,iElem)) THEN ! J-Plane 
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      IntElem(i,j,k,iElem) = SQRT(SUM(Metrics_gtilde(:,i,j,k,iElem,0)**2))
      IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(i)*wGP(k)
    END DO; END DO; END DO
    average_type(iElem) = 3
  ELSEIF (average_ind(1,iElem).AND.average_ind(2,iElem).AND.(.NOT.average_ind(3,iElem))) THEN ! K-Plane 
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      IntElem(i,j,k,iElem) = SQRT(SUM(Metrics_htilde(:,i,j,k,iElem,0)**2))
      IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(i)*wGP(j)
    END DO; END DO; END DO
    average_type(iElem) = 4
  ELSEIF (average_ind(1,iElem).AND.(.NOT.average_ind(2,iElem)).AND.(.NOT.average_ind(3,iElem))) THEN ! I-Line 
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,dXCL_N(1,1:3,:,:,:,iElem),dX_N(:,:,:,:))
    IntElem(:,:,:,iElem) = (dX_N(1,:,:,:)**2+dX_N(2,:,:,:)**2+dX_N(3,:,:,:)**2)**0.5
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(i)
    END DO; END DO; END DO
    average_type(iElem) = 5
  ELSEIF ((.NOT.average_ind(1,iElem)).AND.average_ind(2,iElem).AND.(.NOT.average_ind(3,iElem))) THEN ! J-Line 
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,dXCL_N(2,1:3,:,:,:,iElem),dX_N(:,:,:,:))
    IntElem(:,:,:,iElem) = (dX_N(1,:,:,:)**2+dX_N(2,:,:,:)**2+dX_N(3,:,:,:)**2)**0.5
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(j)
    END DO; END DO; END DO
    average_type(iElem) = 6
  ELSEIF ((.NOT.average_ind(1,iElem)).AND.(.NOT.average_ind(2,iElem)).AND.average_ind(3,iElem)) THEN ! K-Line
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,dXCL_N(3,1:3,:,:,:,iElem),dX_N(:,:,:,:))
    IntElem(:,:,:,iElem) = (dX_N(1,:,:,:)**2+dX_N(2,:,:,:)**2+dX_N(3,:,:,:)**2)**0.5
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(k)
    END DO; END DO; END DO
    average_type(iElem) = 7
  END IF
END DO

N_Testfilter = GETINT('N_Testfilter')
DO i=0,N_Testfilter
  FilterMat_Testfilter(i,i) = 1.
END DO
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
USE MOD_EddyVisc_Vars,     ONLY:SGS_Ind,muSGSmax,S_en_out
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem,i,j,k
REAL,INTENT(IN)                           :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho
REAL,INTENT(INOUT)                        :: muSGS 
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
S_eN = SQRT(S_eN)
! dynsmag model
muSGS = S_eN*rho*SGS_Ind(1,i,j,k,iElem) 
S_en_out(1,i,j,k,iElem) = S_eN
SGS_Ind(2,i,j,k,iElem) = muSGS 
muSGSmax(iElem) = MAX(muSGS,muSGSmax(iElem))
END SUBROUTINE dynsmag

SUBROUTINE compute_cd(U_in)
!===============================================================================================================================
!compute TKE indicator needed for the modification of Smagorinskys viscosity.
!===============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EddyVisc_Vars,          ONLY:SGS_Ind,FilterMat_testfilter, IntElem
USE MOD_EddyVisc_Vars,          ONLY:filter_ind, average_type 
USE MOD_Filter,                 ONLY:Filter_Selective
USE MOD_Lifting_Vars,           ONLY:gradUx,gradUy,gradUz
USE MOD_Mesh_Vars,              ONLY:nElems

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
REAL,DIMENSION(3,3,0:PP_N,0:PP_N,0:PP_N)     :: gradv_all, S_ik, M_ik, L_ik
REAL,DIMENSION(3,0:PP_N,0:PP_N,0:PP_N)       :: V_filtered
REAL, DIMENSION(0:PP_N,0:PP_N,0:PP_N)        :: MM,ML,S_eN_filtered,S_eN,divV
REAL                                         :: D_Ratio
REAL                                         :: MMa,MLa
!===============================================================================================================================
DO iElem=1,nElems
  !store gradients in matrix for readability and filtering
  gradv_all(1,1,:,:,:) = gradUx(2,:,:,:,iElem)
  gradv_all(1,2,:,:,:) = gradUy(2,:,:,:,iElem)
  gradv_all(1,3,:,:,:) = gradUz(2,:,:,:,iElem)
  gradv_all(2,1,:,:,:) = gradUx(3,:,:,:,iElem)
  gradv_all(2,2,:,:,:) = gradUy(3,:,:,:,iElem)
  gradv_all(2,3,:,:,:) = gradUz(3,:,:,:,iElem)
  gradv_all(3,1,:,:,:) = gradUx(4,:,:,:,iElem)
  gradv_all(3,2,:,:,:) = gradUy(4,:,:,:,iElem)
  gradv_all(3,3,:,:,:) = gradUz(4,:,:,:,iElem)
  divV(:,:,:) = 2./3.*(gradv_all(1,1,:,:,:)+gradv_all(2,2,:,:,:)+gradv_all(3,3,:,:,:))

  ! dynamic Smagorinsky model
  DO i=1,3
    V_Filtered(i,:,:,:) = U_in(i+1,:,:,:,iElem)/U_in(1,:,:,:,iElem)
  END DO !i

  !Filter velocities
  CALL Filter_Selective(3,FilterMat_testfilter,V_filtered,filter_ind(:,iElem)) 
  !           _ _   __
  !Compute L=-u u + uu  !!TENSOR ik
  DO i=1,3
    DO k=1,3
      L_ik(i,k,:,:,:) = U_in(i+1,:,:,:,iElem)*U_in(k+1,:,:,:,iElem)/U_in(1,:,:,:,iElem)**2
    END DO ! k
    CALL Filter_Selective(3,FilterMat_testfilter,L_ik(i,1:3,:,:,:),filter_ind(:,iElem)) 
  END DO ! i
  DO i=1,3
    DO k=1,3
      L_ik(i,k,:,:,:) = L_ik(i,k,:,:,:) - V_filtered(i,:,:,:)*V_filtered(k,:,:,:)
    END DO ! k
  END DO ! i
  ! Compute shear rate tensor
  DO i=1,3
    DO j=1,3
      S_ik(i,j,:,:,:) = (gradv_all(i,j,:,:,:)+gradv_all(j,i,:,:,:))*0.5
    END DO
  END DO
  !correct for compressibility
  S_ik(1,1,:,:,:) = S_ik(1,1,:,:,:) -divV(:,:,:)
  S_ik(2,2,:,:,:) = S_ik(2,2,:,:,:) -divV(:,:,:)
  S_ik(3,3,:,:,:) = S_ik(3,3,:,:,:) -divV(:,:,:)
  !Comnpute norm of S_ik
  S_eN(:,:,:)= (S_ik(1,1,:,:,:)**2 + S_ik(2,2,:,:,:)**2. + S_ik(3,3,:,:,:)**2.)
  S_eN(:,:,:)= S_eN(:,:,:) +  2*(S_ik(2,1,:,:,:) )**2.
  S_eN(:,:,:)= S_eN(:,:,:) +  2*(S_ik(3,1,:,:,:) )**2.
  S_eN(:,:,:)= S_eN(:,:,:) +  2*(S_ik(3,2,:,:,:) )**2.
  S_eN(:,:,:)= sqrt(2* S_eN(:,:,:) )

  !            _____            _ _   ____
  !Compute M=-(Delta/Delta)**2*|S|S + |S|S  !!TENSOR ik
  DO i=1,3
    DO k=1,3
      M_ik(i,k,:,:,:) = S_eN(:,:,:)*S_ik(i,k,:,:,:)
    END DO ! k
    CALL Filter_Selective(3,FilterMat_testfilter,M_ik(i,1:3,:,:,:),filter_ind(:,iElem)) 
  END DO ! i

  !filter gradients
  CALL Filter_Selective(3,FilterMat_Testfilter,gradv_all(:,1,:,:,:),filter_ind(:,iElem))
  CALL Filter_Selective(3,FilterMat_Testfilter,gradv_all(:,2,:,:,:),filter_ind(:,iElem))
  CALL Filter_Selective(3,FilterMat_Testfilter,gradv_all(:,3,:,:,:),filter_ind(:,iElem))
  divV(:,:,:) = 2./3.*(gradv_all(1,1,:,:,:)+gradv_all(2,2,:,:,:)+gradv_all(3,3,:,:,:))
  ! Compute filtered shear rate tensor, unfiltered not needed anymore
  DO i=1,3
    DO j=1,3
      S_ik(i,j,:,:,:) = (gradv_all(i,j,:,:,:)+gradv_all(j,i,:,:,:))*0.5
    END DO
  END DO
  S_ik(1,1,:,:,:) = S_ik(1,1,:,:,:) -divV(:,:,:)
  S_ik(2,2,:,:,:) = S_ik(2,2,:,:,:) -divV(:,:,:)
  S_ik(3,3,:,:,:) = S_ik(3,3,:,:,:) -divV(:,:,:)
  !Comnpute norm of S_ik filtered
  S_eN_filtered(:,:,:)= (S_ik(1,1,:,:,:)**2 + S_ik(2,2,:,:,:)**2. + S_ik(3,3,:,:,:)**2.)
  S_eN_filtered(:,:,:)= S_eN_filtered(:,:,:) +  2*(S_ik(2,1,:,:,:) )**2.
  S_eN_filtered(:,:,:)= S_eN_filtered(:,:,:) +  2*(S_ik(3,1,:,:,:) )**2.
  S_eN_filtered(:,:,:)= S_eN_filtered(:,:,:) +  2*(S_ik(3,2,:,:,:) )**2.
  S_eN_filtered(:,:,:)= sqrt(2* S_eN_filtered(:,:,:) )
  !                                 _____
  !D_ratio: square ratio of filter widhts (Delta/Delta)**2, reciprocal of polynomial degrees
  D_Ratio=(REAL(PP_N+1)/REAL(N_testfilter+1))**2
  DO i=1,3
    DO k=1,3
      M_ik(i,k,:,:,:) = M_ik(i,k,:,:,:) - D_Ratio * S_eN_filtered(:,:,:)*S_ik(i,k,:,:,:)
    END DO ! k
  END DO ! i
  !
  !Compute C_d**2 =1/2* M_ik*L_ik / M_ik*M_ik  !!SUM OVER ik !!!for muSGS=C_d**2|S|
  !contract with M_ik according to least square approach of Lilly
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

  ! Selective cell average (on selected directions)
  SELECT CASE (average_type(iElem))
  CASE (1) ! Volume
    MMa = 0.; MLa = 0.
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
      MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
    END DO; END DO; END DO !ijk
    SGS_Ind(1,:,:,:,iElem) = MLa/MMa 
  CASE (2) ! I-Plane
    DO i=0,PP_N
      MMa = 0.; MLa = 0.
      DO j=0,PP_N; DO k=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO; END DO !jk
      SGS_Ind(1,i,:,:,iElem) = MLa/MMa 
    END DO !i
  CASE (3) ! J-Plane
    DO j=0,PP_N
      MMa = 0.; MLa = 0.
      DO i=0,PP_N; DO k=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO; END DO !jk
      SGS_Ind(1,:,j,:,iElem) = MLa/MMa 
    END DO !j
  CASE (4) ! K-Plane
    DO k=0,PP_N
      MMa = 0.; MLa = 0.
      DO i=0,PP_N; DO j=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO; END DO !jk
      SGS_Ind(1,:,:,k,iElem) = MLa/MMa 
    END DO !k
  CASE (5) ! I-Line
    DO j=0,PP_N; DO k=0,PP_N
      MMa = 0.; MLa = 0.
      DO i=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO ! i
      SGS_Ind(1,:,j,k,iElem) = MLa/MMa 
    END DO; END DO !jk
  CASE (6) ! J-Line
    DO i=0,PP_N; DO k=0,PP_N
      MMa = 0.; MLa = 0.
      DO j=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO ! i
      SGS_Ind(1,i,:,k,iElem) = MLa/MMa 
    END DO; END DO !jk
  CASE (7) ! K-Line
    DO i=0,PP_N; DO j=0,PP_N
      MMa = 0.; MLa = 0.
      DO k=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO ! i
      SGS_Ind(1,i,j,:,iElem) = MLa/MMa 
    END DO; END DO!jk
  END SELECT
END DO
END SUBROUTINE compute_cd 

!===============================================================================================================================
! Is vector vec1 "more normal" to vec2 "than parallel"?
!===============================================================================================================================
FUNCTION ISNORMAL(vec1,vec2) result(norm)
REAL, INTENT(IN) :: vec1(:), vec2(:)
LOGICAL          :: norm 

REAL :: vec2Norm(size(vec2,1)), parcomp, normvec(size(vec1,1)), normcomp

! norm of vec2
vec2Norm = vec2/SQRT(SUM(vec2**2))
! component of vec1 parallel to vec2
parcomp = SUM(vec1*vec2Norm)
! part of vec1 normal to vec2
normvec = vec1 - parcomp*vec2Norm
! component of vec1 normal to vec2
normcomp = SQRT(SUM(normvec**2)) 
IF (normcomp .GE. ABS(parcomp)) THEN
  norm = .TRUE.
ELSE
  norm = .FALSE.
ENDIF
END FUNCTION ISNORMAL


!===============================================================================================================================
!> Deallocate arrays and finalize variables used by Smagorinsky SGS model
!===============================================================================================================================
SUBROUTINE Finalizedynsmag()
! MODULES
USE MOD_EddyVisc_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!===============================================================================================================================
DEALLOCATE(DeltaS)
DEALLOCATE(DeltaS_master)
DEALLOCATE(DeltaS_slave)
DEALLOCATE(SGS_Ind)
DEALLOCATE(SGS_Ind_master)
DEALLOCATE(SGS_Ind_slave)
DEALLOCATE(muSGS)
DEALLOCATE(muSGSmax)
DEALLOCATE(FilterMat_Testfilter)
DEALLOCATE(filter_ind)
DEALLOCATE(average_ind)
DEALLOCATE(average_type)
DEALLOCATE(IntELem)
dynsmagInitIsDone = .FALSE.
END SUBROUTINE Finalizedynsmag





#endif
END MODULE MOD_dynsmag
