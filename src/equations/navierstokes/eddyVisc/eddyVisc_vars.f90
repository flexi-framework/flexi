!===================================================================================================================================
!> Contains the parameters needed for eddy viscosity models
!===================================================================================================================================
MODULE MOD_EddyVisc_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

ABSTRACT INTERFACE
  SUBROUTINE EddyViscInt(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,iElem,i,j,k,muSGS)
  INTEGER,INTENT(IN)  :: iElem  !< index of current element
  !> indices of the c
  INTEGER,INTENT(IN)  :: i,j,k
  !> gradients of the directions
  REAL,INTENT(IN)     :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32
  REAL,INTENT(IN)     :: rho    !< Density
  REAL,INTENT(INOUT)  :: muSGS  !< local SGS viscosity
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE EddyVisc_surfInt(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,DeltaSS,SGS_Ind,muSGS,Face_xGP)
  !> gradients of the velocities w.r.t. all directions
  REAL,INTENT(IN)   :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32
  REAL,INTENT(IN)   :: rho      !< Density
  REAL,INTENT(IN)   :: DeltaSS  !< Filter width
  REAL,INTENT(IN)   :: SGS_Ind  !< Indicator for SGS model
  REAL,INTENT(IN)   :: Face_xGP !< Coordinate for van-Driest damping
  REAL,INTENT(OUT)  :: muSGS    !< local SGS viscosity
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE testfilterInt(U_in)
    USE MOD_PreProc
    USE MOD_Mesh_Vars,              ONLY:nElems
    REAL,INTENT(IN)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
  END SUBROUTINE
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                             :: eddyViscType           !< type of eddy viscosity
CHARACTER(LEN=255)                  :: WallDistFile
PROCEDURE(EddyViscInt)     ,POINTER :: eddyViscosity          !< pointer to routine for computing volume eddy viscosity
PROCEDURE(EddyVisc_surfInt),POINTER :: eddyViscosity_surf     !< pointer to routine for computing surface eddy viscosity
PROCEDURE(testfilterInt),POINTER       :: testfilter     !< pointer to routine for computing test filter

!Smagosinsky Standard
REAL,ALLOCATABLE  :: DeltaS(:)         !< filter width, used by Smagorinsky modell
REAL,ALLOCATABLE  :: DeltaS_master(:)
REAL,ALLOCATABLE  :: DeltaS_slave(:)
REAL,ALLOCATABLE  :: muSGS(:,:,:,:,:)  !< Viscosity for the sub-grid
REAL,ALLOCATABLE  :: muSGSmax(:)       !< Viscosity for the sub-grid
REAL              :: CS                !< Smagorinsky constant, LES
REAL              :: PrSGS             !< Prandtl number for the sub-grid scales
!TKECUbed
REAL,ALLOCATABLE  :: FilterMat_Testfilter(:,:) 
REAL,ALLOCATABLE  :: SGS_Ind(:,:,:,:,:) 
REAL,ALLOCATABLE  :: MM_Avg(:,:,:,:) 
REAL,ALLOCATABLE  :: ML_Avg(:,:,:,:) 
REAL,ALLOCATABLE  :: SGS_Ind_master(:,:,:,:) 
REAL,ALLOCATABLE  :: SGS_Ind_slave(:,:,:,:) 
! Dynamic Smagorinsky
REAL,ALLOCATABLE  :: lineElem(:,:,:,:,:) 
LOGICAL, ALLOCATABLE :: filter_ind(:,:) !< Do filter along i,j,k index?
LOGICAL, ALLOCATABLE :: average_ind(:,:) !< Do average along i,j,k index?
!MATTEO:debug output
REAL,ALLOCATABLE  :: S_en_out(:,:,:,:,:)  !< Debug output of |S|
REAL,ALLOCATABLE  :: filtdir_out(:)  !< Debug output of filtering directions
REAL,ALLOCATABLE  :: walldist_out(:)  !< Debug output of wall distance
REAL,ALLOCATABLE  :: walldist_x(:)  !< Debug output of wall distance
REAL,ALLOCATABLE  :: walldist_y(:)  !< Debug output of wall distance
REAL,ALLOCATABLE  :: walldist_z(:)  !< Debug output of wall distance


LOGICAL           :: VanDriest=.FALSE.
LOGICAL           :: SmagorinskyInitIsDone=.FALSE.
LOGICAL           :: DynSmagInitIsDone=.FALSE.

END MODULE MOD_EddyVisc_Vars
