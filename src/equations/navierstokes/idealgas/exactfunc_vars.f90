MODULE MOD_Exactfunc_Vars
!==================================================================================================================================
! Contains the parameters needed for the Navier Stokes calculation
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
REAL              :: MachShock         !< Shoch Mach speed for ExactFunction = 6 (shock)
REAL              :: PreShockDens      !< Pre-shock density for ExactFunction = 6 (shock)
REAL              :: AdvVel(3)         !< Advection Velocity for the test cases
REAL              :: IniCenter(3)      !< parameter used for Shu vortex
REAL              :: IniAxis(3)        !< parameter used for Shu vortex
REAL              :: IniFrequency      !< parameter used for Shu vortex
REAL              :: IniAmplitude      !< parameter used for Shu vortex
REAL              :: IniHalfwidth      !< parameter used for Shu vortex
REAL              :: P_Parameter       !< parameter for Couette-Poiseuille flow
REAL              :: U_Parameter       !< parameter for Couette-Poiseuille flow


!==================================================================================================================================

END MODULE MOD_Exactfunc_Vars
