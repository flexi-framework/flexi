!==================================================================================================================================
!> Contains the global variables needed by the ideal gas equation of state.
!==================================================================================================================================
MODULE MOD_EOS_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
#if PARABOLIC
REAL              :: mu0               !< dynamic viscosity \f$\mu\f$
REAL              :: Pr                !< Prandtl number
REAL              :: KappasPr          !< \f$\kappa\f$/Pr
REAL              :: lambda            !< thermal conductivity
#if PP_VISC==1
REAL              :: Ts                !< Sutherland temperature
REAL              :: cSuth             !< Parameters used in muSuth
#endif
#if (PP_VISC==1) || (PP_VISC==2)
REAL              :: Tref,ExpoSuth     !< Parameters used in muSuth and power law
#endif
#endif /*PARABOLIC*/
REAL              :: cp                !< specific heat at constant pressure
REAL              :: cv                !< specific heat at constant volume
REAL              :: Kappa             !< heat capacity ratio / isentropic exponent
REAL              :: KappaM1           !< = \f$\kappa - 1\f$
REAL              :: sKappaM1          !< = \f$1/(\kappa -1)\f$
REAL              :: KappaP1           !< = \f$\kappa + 1\f$
REAL              :: sKappaP1          !< = \f$1/(\kappa +1)\f$
REAL              :: R                 !< specific gas constant

LOGICAL           :: EosInitIsDone=.FALSE.
!==================================================================================================================================
END MODULE MOD_EOS_Vars
