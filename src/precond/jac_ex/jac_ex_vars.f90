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
REAL,ALLOCATABLE                      :: LL_plus(:,:)
REAL,ALLOCATABLE                      :: LL_minus(:,:)
REAL,ALLOCATABLE                      :: LL_mp(:,:)         !LL_mp(i,iLocSide)=LL_minus/plus(i,i)
REAL,ALLOCATABLE                      :: l_mp(:,:)          !l_mp(i,iLocSide)=l_minus/plus(i) 
REAL,ALLOCATABLE                      :: nVec(:,:,:,:,:)    !local normal vectors sorted in the ijk System
REAL,ALLOCATABLE                      :: Surf(:,:,:,:)      !surf element sorted in the ijk system
#if PARABOLIC
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
