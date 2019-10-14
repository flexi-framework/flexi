MODULE MOD_Jac_Ex_Vars
!===================================================================================================================================
! Contains global variables used by the Jac_Ex module.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                               :: Jac_Ex_InitIsDone=.FALSE.
REAL,ALLOCATABLE                      :: LL_plus(:,:)              !< LL_plus(i,j) = L^hat_plus(i)*L_plus(j)
REAL,ALLOCATABLE                      :: LL_minus(:,:)             !< LL_minus(i,j) = L^hat_minus(i)*L_minus(j)
#if PARABOLIC
REAL,ALLOCATABLE                      :: L_mp(:,:)                 !< L_mp(i,iLocSide)=either L_minus(i) or L_plus(i) 
REAL,ALLOCATABLE                      :: R_minus(:,:,:,:)   ! BR2 lifting surface term
REAL,ALLOCATABLE                      :: R_plus(:,:,:,:)    ! BR2 lifting surface term
REAL,ALLOCATABLE                      :: JacLiftingFlux(:,:,:,:,:)
!REAL,ALLOCATABLE                      :: PrimConsJac(:,:,:,:,:)  ! Derivative of U_Prim with respect to U_cons
#endif /*PARABOLIC*/
#if FV_ENABLED && FV_RECONSTRUCT
REAL,ALLOCATABLE                      :: UPrim_extended(:,:,:,:,:) !< extended primitive solution array containing additional
                                                                   !< first layers of neighbouring elements
REAL,ALLOCATABLE                      :: FV_sdx_XI_extended(:,:,:,:)   !< extended inverse of distance between neighboring dofs in x
REAL,ALLOCATABLE                      :: FV_sdx_ETA_extended(:,:,:,:)  !< extended inverse of distance between neighboring dofs in y
REAL,ALLOCATABLE                      :: FV_sdx_ZETA_extended(:,:,:,:) !< extended inverse of distance between neighboring dofs in z
#endif
!===================================================================================================================================
END MODULE MOD_Jac_Ex_Vars
