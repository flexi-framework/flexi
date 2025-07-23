!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
#include "eos.h"

!==================================================================================================================================
!> Routines to provide boundary conditions for the domain. Fills the boundary part of the fluxes list.
!>
!> Available boundary conditions are:
!>  * 1   : Periodic boundary conditions (no work to be done here, are already filled due to mesh connection)
!>  DIRICHLET BCs:
!>  * 2   : Use the initial exact function (if BC state = 0) or a refstate as dirichlet boundary conditions
!>  * 12  : Read in dirichlet boundary conditions from a HDF5 file
!>  * 121 : Similar to 2, but pre-compute and store the evaluation of an exact func prescribed by BCState
!>  * 22  : Similar to 2, but BCState specifies exact function to be used
!>  WALL BCs:
!>  * 3   : Adiabatic wall
!>  * 4   : Isothermal wall (Temperature specified by refstate)
!>  * 9   : Slip wall
!>  * 91  : Slip wall with correct gradient calculation (expensive)
!>  OUTFLOW BCs:
!>  * 23  : Outflow BC where the second entry of the refstate specifies the desired Mach number at the outflow
!>  * 24  : Pressure outflow BC (pressure specified by refstate)
!>  * 25  : Subsonic outflow BC
!>  INFLOW BCs:
!>  * 27  : Subsonic inflow BC, WARNING: REFSTATE is different: Tt,alpha,beta,<empty>,pT (4th entry ignored), angles in DEG
!==================================================================================================================================
MODULE MOD_GetBoundaryFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: InitBC
PUBLIC:: GetBoundaryFlux
PUBLIC:: GetBoundaryState
PUBLIC:: FinalizeBC
#if FV_ENABLED && FV_RECONSTRUCT
PUBLIC:: GetBoundaryFVgradient
#endif /*FV_ENABLED && FV_RECONSTRUCT*/
#if PARABOLIC
PUBLIC:: Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Initialize boundary conditions. Read parameters and sort boundary conditions by types.
!> Call boundary condition specific init routines.
!==================================================================================================================================
SUBROUTINE InitBC()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Viscosity
USE MOD_ReadInTools
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: nRefState,BCData,BCDataPrim,nBCByType,BCSideID
USE MOD_Equation_Vars     ,ONLY: BCStateFile,RefStatePrim
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,BoundaryType,nBCs,Face_xGP
#if PARABOLIC
USE MOD_Exactfunc_Vars    ,ONLY: delta99_in,x_in,BlasiusInitDone
#endif
USE MOD_EOS               ,ONLY: ConsToPrim
USE MOD_ExactFunc         ,ONLY: ExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: locType,locState
INTEGER :: MaxBCState,MaxBCStateGlobal
LOGICAL :: readBCdone
REAL    :: talpha,tbeta
INTEGER :: p,q
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
   CALL CollectiveStop(__STAMP__,&
     "InitBC not ready to be called or already called.")
END IF
! determine globally max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)

  ! Check for max. Refstate used if current BC requires Refstate
  IF((locType.EQ. 2).OR.(locType.EQ. 4).OR. &
     (locType.EQ.23).OR.(locType.EQ.24).OR. &
     (locType.EQ.25).OR.(locType.EQ.27)     ) MaxBCState = MAX(MaxBCState,locState)

  ! If required, check if Refstate available
  IF (locState.LT.1) THEN
    SELECT CASE (locType)
    CASE(4)
      CALL Abort(__STAMP__,'No refstate (rho,x,x,x,p) defined to compute temperature from density and pressure for BC_TYPE',locType)
    CASE(23)
      CALL Abort(__STAMP__,'No outflow Mach number in refstate (x,Ma,x,x,x) defined for BC_TYPE',locType)
    CASE(24,25)
      CALL Abort(__STAMP__,'No outflow pressure in refstate (x,x,x,x,p) defined for BC_TYPE',locType)
    CASE(27)
      CALL Abort(__STAMP__,'No inflow refstate (Tt,alpha,beta,<empty>,pT) in refstate defined for BC_TYPE',locType)
    CASE(121,22)
      CALL Abort(__STAMP__,'No exactfunc defined for BC_TYPE',locType) ! Technically not a missing refstate but exactfunc
    END SELECT
  END IF
#if FV_RECONSTRUCT
  IF((locType.EQ.3).OR.(locType.EQ.4))THEN
    ASSOCIATE(prim => RefStatePrim(:,locState))
#if PARABOLIC
    IF(VISCOSITY_PRIM(prim).LE.0.) &
#endif
    CALL Abort(__STAMP__,'No-slip BCs cannot be used without viscosity in case of FV-reconstruction!')
    END ASSOCIATE
  END IF
#endif
END DO
MaxBCStateGLobal=MaxBCState
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

! Sanity check for BCs
IF(MaxBCState.GT.nRefState)THEN
  CALL Abort(__STAMP__,&
    'ERROR: Boundary RefState not defined! (MaxBCState,nRefState):',MaxBCState,REAL(nRefState))
END IF

#if PARABOLIC
! Check for Blasius BCs and read parameters if this has not happened in the equation init
IF (.NOT.BlasiusInitDone) THEN
   DO i=1,nBCs
     locType =BoundaryType(i,BC_TYPE)
     locState=BoundaryType(i,BC_STATE)
     IF ((locType.EQ.121).AND.(locState.EQ.1338)) THEN
       delta99_in      = GETREAL('delta99_in')
       x_in            = GETREALARRAY('x_in',2,'(/0.,0./)')
       BlasiusInitDone = .TRUE.
       EXIT
     END IF
   END DO
END IF
#endif

! Allocate buffer array to store temp data for all BC sides
ALLOCATE(BCData(PP_nVar,        0:PP_N,0:PP_NZ,nBCSides))
ALLOCATE(BCDataPrim(PP_nVarPrim,0:PP_N,0:PP_NZ,nBCSides))
BCData=0.
BCDataPrim=0.

! Initialize boundary conditions
readBCdone=.FALSE.
DO i=1,nBCs
  locType =BoundaryType(i,BC_TYPE)
  locState=BoundaryType(i,BC_STATE)
  SELECT CASE (locType)
  CASE(12) ! State File Boundary condition
    IF(.NOT.readBCdone) CALL ReadBCFlow(BCStateFile)
    readBCdone=.TRUE.
  CASE(27) ! Subsonic inflow
    ! Compute normalized direction vector a(1:3) from paper:
    !   "Verification Assessment of Flow Boundary Conditions for CFD", John W. Slater, AIAA 3882, 2021.
    ! to later compute the projection of the velocity direction normal (prescribed with alpha and beta) to the local face normal.
    talpha=TAN(PP_PI/180.*RefStatePrim(2,locState)) ! Convert alpha from degree to radian and compute tan
    tbeta =TAN(PP_PI/180.*RefStatePrim(3,locState)) ! Convert beta  from degree to radian and compute tan
    RefStatePrim(VEL1,locState)=1.    /SQRT((1.+talpha**2+tbeta**2)) ! (8a)
    RefStatePrim(VEL2,locState)=talpha/SQRT((1.+talpha**2+tbeta**2)) ! (8b)
    RefStatePrim(VEL3,locState)=tbeta /SQRT((1.+talpha**2+tbeta**2)) ! (8c)
  END SELECT
END DO

! Initialize Dirichlet BCs that use a pre-computed and then stored evaluation of an exact func
DO iSide=1,nBCSides
  IF (Boundarytype(BC(iSide),BC_TYPE).EQ.121) THEN
    DO q=0,PP_NZ; DO p=0,PP_N
      CALL ExactFunc(Boundarytype(BC(iSide),BC_STATE),0.,Face_xGP(:,p,q,0,iSide),BCData(:,p,q,iSide))
      CALL ConsToPrim(BCDataPrim(:,p,q,iSide),BCData(:,p,q,iSide))
    END DO; END DO ! p,q=0,PP_N
  END IF
END DO

! Count number of sides of each boundary
ALLOCATE(nBCByType(nBCs))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i) nBCByType(i)=nBCByType(i)+1
  END DO
END DO

! Sort BCs by type, store SideIDs
ALLOCATE(BCSideID(nBCs,MAXVAL(nBCByType)))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i)THEN
      nBCByType(i)=nBCByType(i)+1
      BCSideID(i,nBCByType(i))=iSide
    END IF
  END DO
END DO

END SUBROUTINE InitBC


!==================================================================================================================================
!> Computes the boundary state for the different boundary conditions.
!==================================================================================================================================
SUBROUTINE GetBoundaryState(SideID,t,Nloc,UPrim_boundary,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: BoundaryType,BC
USE MOD_EOS          ,ONLY: ConsToPrim,PrimtoCons
USE MOD_EOS          ,ONLY: PRESSURE_RIEMANN
USE MOD_EOS_Vars     ,ONLY: sKappaM1,Kappa,KappaM1,R
USE MOD_ExactFunc    ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc,BCDataPrim,RefStatePrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: SideID                                   !< ID of current side
REAL,INTENT(IN)         :: t                                        !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)      :: Nloc                                     !< polynomial degree
REAL,INTENT(IN)         :: UPrim_master(  PRIM,0:Nloc,0:ZDIM(Nloc)) !< inner surface solution
REAL,INTENT(IN)         :: NormVec(          3,0:Nloc,0:ZDIM(Nloc)) !< normal surface vectors
REAL,INTENT(IN)         :: TangVec1(         3,0:Nloc,0:ZDIM(Nloc)) !< tangent surface vectors 1
REAL,INTENT(IN)         :: TangVec2(         3,0:Nloc,0:ZDIM(Nloc)) !< tangent surface vectors 2
REAL,INTENT(IN)         :: Face_xGP(         3,0:Nloc,0:ZDIM(Nloc)) !< positions of surface flux points
REAL,INTENT(OUT)        :: UPrim_boundary(PRIM,0:Nloc,0:ZDIM(Nloc)) !< resulting boundary state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: p,q
INTEGER                 :: BCType,BCState
REAL,DIMENSION(PP_nVar) :: Cons
REAL                    :: Ma,MaOut,c,cb,pt,pb ! for BCType==23,24,25.27
REAL                    :: U,Tb,Tt,tmp1,tmp2,tmp3,A,Rplus,nv(3) ! for BCType==27
!===================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(2) ! Exact function or refstate
  IF(BCState.EQ.0)THEN
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q),Cons)
      CALL ConsToPrim(UPrim_boundary(:,p,q),Cons)
    END DO; END DO
  ELSE
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      UPrim_boundary(:,p,q) = RefStatePrim(:,BCState)
    END DO; END DO
  END IF

CASE(12)  ! Dirichlet-type: BCState from readin state
  UPrim_boundary(:,:,:) = BCDataPrim(:,:,:,SideID)

CASE(121) ! Dirichlet-type:BCState from exact function computed once at the beginning of the simulation
  UPrim_boundary(:,:,:) = BCDataPrim(:,:,:,SideID)

CASE(22)  ! Dirichlet-type: BCState specifies exactfunc to be used
  DO q=0,ZDIM(Nloc); DO p=0,Nloc
    CALL ExactFunc(BCState,t,Face_xGP(:,p,q),Cons)
    CALL ConsToPrim(UPrim_boundary(:,p,q),Cons)
  END DO; END DO

CASE(3,4,9,91,23,24,25,27)
  ! Initialize boundary state with rotated inner state
  DO q=0,ZDIM(Nloc); DO p=0,Nloc
    ! transform state into normal system
    UPrim_boundary(DENS,p,q)     = UPrim_master(DENS,p,q)
    UPrim_boundary(VEL1,p,q)     = DOT_PRODUCT(UPrim_master(VELV,p,q),NormVec( :,p,q))
    UPrim_boundary(VEL2,p,q)     = DOT_PRODUCT(UPrim_master(VELV,p,q),TangVec1(:,p,q))
    UPrim_boundary(VEL3,p,q)     = DOT_PRODUCT(UPrim_master(VELV,p,q),TangVec2(:,p,q))
    UPrim_boundary(PRES,p,q)     = UPrim_master(PRES,p,q)
    UPrim_boundary(TEMP,p,q)     = UPrim_master(TEMP,p,q)
  END DO; END DO !p,q

  SELECT CASE(BCType)
  CASE(3) ! Adiabatic wall
    ! For adiabatic wall all gradients are 0
    ! We reconstruct the BC State, rho=rho_L, velocity=0, rhoE_wall = p_Riemann/(Kappa-1)
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      UPrim_boundary(PRES,p,q) = PRESSURE_RIEMANN(UPrim_boundary(:,p,q)) ! pressure from solving local Riemann problem
      UPrim_boundary(VELV,p,q) = 0.                                      ! no slip
      UPrim_boundary(TEMP,p,q) = UPrim_master(TEMP,p,q)                  ! adiabatic => temperature from the inside
      ! set density via ideal gas equation, consistent to pressure and temperature
      UPrim_boundary(DENS,p,q) = UPrim_boundary(PRES,p,q)/(UPrim_boundary(TEMP,p,q)*R)
    END DO; END DO ! q,p

  CASE(4) ! Isothermal wall
    ! For isothermal wall, all gradients are from interior
    ! We reconstruct the BC State, rho=rho_L, velocity=0, rhoE_wall =  rho_L*C_v*Twall
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! Set pressure by solving local Riemann problem
      UPrim_boundary(PRES,p,q) = PRESSURE_RIEMANN(UPrim_boundary(:,p,q)) ! pressure from solving local Riemann problem
      UPrim_boundary(VELV,p,q) = 0.                                      ! no slip
      UPrim_boundary(TEMP,p,q) = RefStatePrim(TEMP,BCState)              ! temperature from RefState
      ! set density via ideal gas equation, consistent to pressure and temperature
      UPrim_boundary(DENS,p,q) = UPrim_boundary(PRES,p,q)/(UPrim_boundary(TEMP,p,q)*R)
    END DO; END DO ! q,p

  CASE(9,91) ! Euler (slip) wall
    ! vel=(0,v_in,w_in)
    ! NOTE: from this state ONLY the velocities should actually be used for the diffusive flux
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! Set pressure by solving local Riemann problem
      UPrim_boundary(PRES,p,q) = PRESSURE_RIEMANN(UPrim_boundary(:,p,q)) ! pressure from solving local Riemann problem
      UPrim_boundary(VEL1,p,q) = 0.                                      ! slip in tangential directions
      ! Density is chosen from the inside, following
      ! "Riemann Solvers and Numerical Methods for Fluid Dynamics", Toro (Chapter 6.3.3 Boundary Conditions)
      UPrim_boundary(DENS,p,q) = UPrim_master(DENS,p,q) ! density from inside
      ! set temperature via ideal gas equation, consistent to density and pressure
      UPrim_boundary(TEMP,p,q) = UPrim_boundary(PRES,p,q)/(UPrim_boundary(DENS,p,q)*R)
    END DO; END DO ! q,p

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Cases 21-29 are taken from NASA report:
  ! "Inflow/Outflow Boundary Conditions with Application to FUN3D", Jan-Reneé Carlson, NASA/TM–2011-217181, 2011.
  ! and correspond to the BCs 2.1 to 2.9.
  ! NOTE: Quantities in paper are non-dimensional such that T=c^2.
  !---------------------------------------------------------------------------------------------------------------------------------
  CASE(23) ! Outflow Mach number BC
    ! NOTE: Should not be used with adjacent walls (destroys boundary layer profile, like exact function)
    ! Refstate for this case is special, VelocityX specifies outlet mach number
    ! State: (/dummy,MaOut,dummy,dummy,dummy/)
    MaOut=RefStatePrim(2,BCState) ! Mach number prescribed by user. Corresponds to M_set in paper
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      c  = SQRT(kappa*UPrim_boundary(PRES,p,q)/UPrim_boundary(DENS,p,q)) ! (19) local speed of sound from inside
      Ma = UPrim_Boundary(VEL1,p,q)/c                                    ! (20) Mach number based on (inner) side-normal component
      ! (23) set pressure depending on subsonic or supersonic case
      IF(Ma<1) THEN ! subsonic
        ! Compute local total pressure pt based on local (inner) Mach number and local (inner) pressure with isentropic relation
        pt = UPrim_boundary(PRES,p,q)*((1+0.5*(kappa-1)*Ma   *Ma)   **( kappa*sKappaM1)) ! (21)
        ! Compute local boundary pressure based on local total pressure pt and prescribed boundary Mach number MaOut
        pb =                       pt*((1+0.5*(kappa-1)*MaOut*MaOut)**(-kappa*sKappaM1)) ! (22)
      ELSE
        ! Supersonic: Use local (inner) total pressure instead
        pb = UPrim_boundary(PRES,p,q)+0.5*UPrim_boundary(DENS,p,q)*DOT_PRODUCT(UPrim_Boundary(VELV,p,q),UPrim_Boundary(VELV,p,q))
      END IF
      ! (24) Set boundary state
      UPrim_boundary(DENS,p,q) = kappa*pb/(c*c)           ! Density based on inner speed of sound and boundary pressure
      UPrim_boundary(VELV,p,q) = UPrim_boundary(VELV,p,q) ! Velocity from inner state
      UPrim_boundary(PRES,p,q) = pb                       ! Computed boundary pressure
      ! set temperature via ideal gas equation, consistent to density and pressure
      UPrim_boundary(TEMP,p,q) = UPrim_boundary(PRES,p,q)/(R*UPrim_boundary(DENS,p,q))
    END DO; END DO !p,q

  CASE(24) ! Pressure outflow BC
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! check if sub- or supersonic
      c  = SQRT(kappa*UPrim_boundary(PRES,p,q)/UPrim_boundary(DENS,p,q)) ! (19) local speed of sound from inside
      Ma = UPrim_Boundary(VEL1,p,q)/c                                    ! (20) Mach number based on (inner) side-normal component
      ! (25) set pressure depending on subsonic or supersonic case
      IF(Ma<1) THEN ! subsonic
        ! (26) Set boundary state
        pb = RefStatePrim(5,BCState)                        ! Pressure prescribed at boundary by user
        UPrim_boundary(DENS,p,q) = kappa*pb/(c*c)           ! Density based on inner speed of sound and boundary pressure
        UPrim_boundary(VELV,p,q) = UPrim_boundary(VELV,p,q) ! Velocity from inner state
        UPrim_boundary(PRES,p,q) = pb                       ! Pressure
        ! set temperature via ideal gas equation, consistent to density and pressure
        UPrim_boundary(TEMP,p,q) = UPrim_boundary(PRES,p,q)/(R*UPrim_boundary(DENS,p,q))
      ELSE
        ! Supersonic: State corresponds to pure inner state, which has already been written to UPrim_Boundary.
        !             Hence, nothing to do here!
      ENDIF
    END DO; END DO !p,q

  CASE(25) ! Subsonic outflow BC
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! check if sub- or supersonic
      c  = SQRT(kappa*UPrim_boundary(PRES,p,q)/UPrim_boundary(DENS,p,q)) ! (19) local speed of sound from inside
      Ma = UPrim_Boundary(VEL1,p,q)/c                                    ! (20) Mach number based on (inner) side-normal component
      ! (27) set pressure depending on subsonic or supersonic case
      IF(Ma<1) THEN
        ! Subsonic: pressure prescribed at boundary
        pb = RefStatePrim(5,BCState)
      ELSE
        ! Supersonic: set local (inner) total pressure
        pb = UPrim_boundary(PRES,p,q)+0.5*UPrim_boundary(DENS,p,q)*DOT_PRODUCT(UPrim_Boundary(VELV,p,q),UPrim_Boundary(VELV,p,q))
      ENDIF
      ! (28) set velocity depending on local flow direction (inflow/outflow), i.e. force outflow by setting normal velocity
      !      always to point outwards.
      IF (UPrim_boundary(VEL1,p,q)<0.) THEN
        UPrim_boundary(VEL1,p,q) = ABS(UPrim_boundary(VEL1,p,q)) ! Multiplication with normal vector of side happens
        UPrim_boundary(VEL2,p,q) = 0.                            ! below by rotating back into global coordinate system
        UPrim_boundary(VEL3,p,q) = 0.
      END IF
      ! (29) Set boundary state
      UPrim_boundary(DENS,p,q) = kappa*pb/(c*c)
      UPrim_boundary(PRES,p,q) = RefStatePrim(5,BCState) ! always outflow pressure
      ! set temperature via ideal gas equation, consistent to density and pressure
      UPrim_boundary(TEMP,p,q) = UPrim_boundary(PRES,p,q)/(R*UPrim_boundary(DENS,p,q))
    END DO; END DO !p,q

  CASE(27) ! Subsonic inflow BC
    ! via stagnation temperature Tt, stag. pressure pt, angle of attack alpha and yaw angle beta
    ! Refstate is different: (Tt,alpha,beta,<empty>,pt) (4th entry ignored)  (angles in DEG not RAD)
    ! WARNING: Computation of the speed of sound at boundary is wrong in Carlsen paper, since it does not account for the angles
    ! alpha and beta in (40), (41) and (43). Correction is proposed below based on the paper
    !   "Verification Assessment of Flow Boundary Conditions for CFD", John W. Slater, AIAA 3882, 2001.

    Tt=RefStatePrim(1  ,BCState) ! Prescribed stagnation temperature
    nv=RefStatePrim(2:4,BCState) ! Vector a(1:3) from AIAA paper. Was precomputed in ini routine. NOT the alpha and beta angle
    pt=RefStatePrim(5  ,BCState) ! Prescribed stagnation pressure

    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! Term A from AIAA paper. Describes projection of BC velocity in global coordinates to side-normal component.
      ! Multiply with -1, since normal vector in AIAA paper is defined INTO the domain
      A=-1.*DOT_PRODUCT(nv(:),NormVec(:,p,q)) ! (7) in AIAA paper

      c=SQRT(kappa*UPrim_boundary(PRES,p,q)/UPrim_boundary(DENS,p,q)) ! (19) local speed of sound from inside
      ! 1D Riemann invariant: Rplus = U_i+2*c_i/kappaM1, Rminus = U_b-2c_b/kappaM1, normal component only!
      Rplus=-UPrim_boundary(VEL1,p,q)-2.*c/KappaM1 ! (37) compute outward propagating invariant based on inner state

      ! ATTENTION: Extrapolated Riemann invariant (39) must be computed based on NORMAL component of velocity U_b only.
      !            However, total enthalpy (38) must be computed based on overall MAGNITUDE of velocity U_b.
      !            Hence, projection A has to be applied to transform between normal component and magnitude of velocity when
      !            inserting (39) into (38), as also shown in Eq. 9 in AIAA paper detailed above.
      !            This yields the correct form of equation (40):
      !                  H_t = c_b^2/(gamma-1) + 1/2*[1/A*(R^+ +2*c_b/(gamma-1))]^2
      !            This is identical to rewriting (5) from AIAA paper in terms of sound speed and solving resulting quad. equation.
      !            Re-arranging above equation as a quadratic equation in c_b yields the correct coefficients as
      tmp1 = A**2 + 2./KappaM1                     ! (43) a
      tmp2 = 2*Rplus                               ! (43) b
      tmp3 = KappaM1/2.*Rplus**2 - Kappa*R*Tt*A**2 ! (43) c NOTE: for ideal gas: H = c_p*T = Kappa/(Kappa-1)*R*T

      ! (44) The max. of the two solutions is the physical one
      cb=MAX( (-tmp2+SQRT(tmp2**2-4*tmp1*tmp3))/(2*tmp1),&  ! (42) solution 1
              (-tmp2-SQRT(tmp2**2-4*tmp1*tmp3))/(2*tmp1) )  ! (42) solution 2

      ! Compute remaining variables at boundary (different to paper)
      Tb = cb**2/(Kappa*R)                             ! via c^2=kappa*R*T
      Ma = SQRT(2./KappaM1*(Tt/Tb-1.))                 ! (46) yields Tt/Tb = 1+(kappa-1)/2*M^2
      pb = pt*(1.+0.5*KappaM1*Ma**2)**(-kappa/kappam1) ! (46) ATTENTION: in paper wrong sign in exponent for pressure
                                                       ! isentropic pressure relation
      U  = Ma*SQRT(Kappa*R*Tb)                         ! Velocity magnitude via Mach number

      ! (47) Set boundary state
      UPrim_boundary(DENS,p,q) = pb/(R*Tb) ! density based on boundary state
      UPrim_boundary(VEL1,p,q) = U*DOT_PRODUCT(nv(:),Normvec( :,p,q)) ! Contribution of magnitude in side-normal coords
      UPrim_boundary(VEL2,p,q) = U*DOT_PRODUCT(nv(:),Tangvec1(:,p,q)) ! in all three directions. Will then be transformed
      UPrim_boundary(VEL3,p,q) = U*DOT_PRODUCT(nv(:),Tangvec2(:,p,q)) ! correctly into global coordinates below
      UPrim_boundary(PRES,p,q) = pb
      UPrim_boundary(TEMP,p,q) = Tb
    END DO; END DO !p,q
  END SELECT

  ! rotate state back to physical system
  DO q=0,ZDIM(Nloc); DO p=0,Nloc
    UPrim_boundary(VELV,p,q) = UPrim_boundary(VEL1,p,q)*NormVec( :,p,q) &
                             + UPrim_boundary(VEL2,p,q)*TangVec1(:,p,q) &
                             + UPrim_boundary(VEL3,p,q)*TangVec2(:,p,q)
  END DO; END DO

CASE(1) !Periodic already filled!
  CALL Abort(__STAMP__, &
      "GetBoundaryState called for periodic side!")
CASE DEFAULT ! unknown BCType
  CALL Abort(__STAMP__,&
       'no BC defined in navierstokes/getboundaryflux.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryState


!==================================================================================================================================
!> Computes the boundary fluxes for a given face (defined by SideID).
!> Calls GetBoundaryState and directly uses the returned values for all Riemann-type BCs.
!> For other types of BCs, we directly compute the flux on the interface.
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(SideID,t,Nloc,Flux,UPrim_master,                   &
#if PARABOLIC
                           gradUx_master,gradUy_master,gradUz_master,&
#endif
                           NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: BoundaryType,BC
USE MOD_EOS          ,ONLY: PrimToCons,ConsToPrim
USE MOD_ExactFunc    ,ONLY: ExactFunc
#if PARABOLIC
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
USE MOD_Riemann      ,ONLY: ViscousFlux
#endif
USE MOD_Riemann      ,ONLY: Riemann
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS_master
#endif
#if TESTCASE_BC
USE MOD_TestCase     ,ONLY: GetBoundaryFluxTestcase
#endif /*TESTCASE_BC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: SideID                                         !< ID of current side
REAL,INTENT(IN)      :: t                                              !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)   :: Nloc                                           !< polynomial degree
REAL,INTENT(IN)      :: UPrim_master( PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)) !< inner surface solution
#if PARABOLIC
REAL,INTENT(IN)      :: gradUx_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !< inner surface solution gradients in x-direction
REAL,INTENT(IN)      :: gradUy_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !< inner surface solution gradients in y-direction
REAL,INTENT(IN)      :: gradUz_master(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)) !< inner surface solution gradients in z-direction
#endif /*PARABOLIC*/
REAL,INTENT(IN)      :: NormVec (3,0:Nloc,0:ZDIM(Nloc))                !< normal vector on surfaces
REAL,INTENT(IN)      :: TangVec1(3,0:Nloc,0:ZDIM(Nloc))                !< tangential1 vector on surfaces
REAL,INTENT(IN)      :: TangVec2(3,0:Nloc,0:ZDIM(Nloc))                !< tangential2 vector on surfaces
REAL,INTENT(IN)      :: Face_xGP(3,0:Nloc,0:ZDIM(Nloc))                !< positions of surface flux points
REAL,INTENT(OUT)     :: Flux(PP_nVar,0:Nloc,0:ZDIM(Nloc))              !< resulting boundary fluxes
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: p,q
INTEGER                              :: BCType,BCState
REAL                                 :: UCons_boundary(PP_nVar    ,0:Nloc,0:ZDIM(Nloc))
REAL                                 :: UCons_master  (PP_nVar    ,0:Nloc,0:ZDIM(Nloc))
#if PARABOLIC
INTEGER                              :: iVar
REAL                                 :: nv(3),tv1(3),tv2(3)
REAL                                 :: BCGradMat(1:PP_dim,1:PP_dim)
REAL                                 :: Fd_Face_loc(PP_nVar,    0:Nloc,0:ZDIM(Nloc))
REAL                                 :: Gd_Face_loc(PP_nVar,    0:Nloc,0:ZDIM(Nloc))
REAL                                 :: Hd_Face_loc(PP_nVar,    0:Nloc,0:ZDIM(Nloc))
REAL                                 :: gradUx_Face_loc(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc))
REAL                                 :: gradUy_Face_loc(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc))
REAL                                 :: gradUz_Face_loc(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc))
REAL                                 :: gradUx_vNormal,gradUx_vTang1,gradUy_vNormal,gradUy_vTang1
REAL                                 :: gradUn_vNormal,gradUn_vTang1,gradUt1_vNormal,gradUt1_vTang1
#if PP_dim == 3
REAL                                 :: gradUx_vTang2,gradUy_vTang2,gradUz_vNormal,gradUz_vTang1,gradUz_vTang2
REAL                                 :: gradUn_vTang2,gradUt1_vTang2,gradUt2_vNormal,gradUt2_vTang1,gradUt2_vTang2
#endif
#endif /*PARABOLIC*/
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

#if TESTCASE_BC
IF (BCType.LT.0) THEN ! testcase boundary condition
  CALL GetBoundaryFluxTestcase(SideID,t,Nloc,Flux,UPrim_master,              &
#if PARABOLIC
                               gradUx_master,gradUy_master,gradUz_master,&
#endif
                               NormVec,TangVec1,TangVec2,Face_xGP)
ELSE
#endif /*TESTCASE_BC*/
  CALL GetBoundaryState(SideID,t,Nloc,UPrim_boundary,UPrim_master,&
      NormVec,TangVec1,TangVec2,Face_xGP)

  SELECT CASE(BCType)
  CASE(2,12,121,22,23,24,25,27) ! Riemann-Type BCs
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      CALL PrimToCons(UPrim_master(:,p,q),  UCons_master(:,p,q))
      CALL PrimToCons(UPrim_boundary(:,p,q),UCons_boundary(:,p,q))
    END DO; END DO ! p,q=0,PP_N
    CALL Riemann(Nloc,Flux,UCons_master,UCons_boundary,UPrim_master,UPrim_boundary, &
        NormVec,TangVec1,TangVec2,doBC=.TRUE.)
#if PARABOLIC
    CALL ViscousFlux(Nloc,Fd_Face_loc,UPrim_master,UPrim_boundary,&
         gradUx_master,gradUy_master,gradUz_master,&
         gradUx_master,gradUy_master,gradUz_master,&
         NormVec&
#if EDDYVISCOSITY
        ,muSGS_master(:,:,:,SideID),muSGS_master(:,:,:,SideID)&
#endif
    )
    Flux = Flux + Fd_Face_loc
#endif /*PARABOLIC*/

  CASE(3,4,9,91) ! Walls
#if EDDYVISCOSITY
    muSGS_master(:,:,:,SideID)=0.
#endif
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! Now we compute the 1D Euler flux, but use the info that the normal component u=0
      ! we directly tranform the flux back into the Cartesian coords: F=(0,n1*p,n2*p,n3*p,0)^T
      Flux(DENS,p,q) = 0.
      Flux(MOMV,p,q) = UPrim_boundary(PRES,p,q)*NormVec(:,p,q)
      Flux(ENER,p,q) = 0.
    END DO; END DO !p,q
    ! Diffusion
#if PARABOLIC
    SELECT CASE(BCType)
    CASE(3,4)
      ! Evaluate 3D Diffusion Flux with interior state and symmetry gradients
      CALL EvalDiffFlux3D(Nloc,UPrim_boundary,&
                          gradUx_master, gradUy_master, gradUz_master, &
                          Fd_Face_loc,   Gd_Face_loc,   Hd_Face_loc    &
#if EDDYVISCOSITY
                         ,muSGS_master(:,:,:,SideID) &
#endif
                         )
      IF (BCType.EQ.3) THEN
        ! Enforce energy flux is exactly zero at adiabatic wall
        Fd_Face_loc(ENER,:,:)=0.
        Gd_Face_loc(ENER,:,:)=0.
        Hd_Face_loc(ENER,:,:)=0.
      END IF
    CASE(9)
      ! Euler/(full-)slip wall
      ! Version 1: set the normal derivatives to zero
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        nv = NormVec(:,p,q)
        ! BCGradMat = I - n * n^T = (gradient - normal component of gradient)
#if (PP_dim==3)
        BCGradMat(1,1) = 1. - nv(1)*nv(1)
        BCGradMat(2,2) = 1. - nv(2)*nv(2)
        BCGradMat(3,3) = 1. - nv(3)*nv(3)
        BCGradMat(1,2) = -nv(1)*nv(2)
        BCGradMat(1,3) = -nv(1)*nv(3)
        BCGradMat(3,2) = -nv(3)*nv(2)
        BCGradMat(2,1) = BCGradMat(1,2)
        BCGradMat(3,1) = BCGradMat(1,3)
        BCGradMat(2,3) = BCGradMat(3,2)
        gradUx_Face_loc(:,p,q) = BCGradMat(1,1) * gradUx_master(:,p,q) &
                               + BCGradMat(1,2) * gradUy_master(:,p,q) &
                               + BCGradMat(1,3) * gradUz_master(:,p,q)
        gradUy_Face_loc(:,p,q) = BCGradMat(2,1) * gradUx_master(:,p,q) &
                               + BCGradMat(2,2) * gradUy_master(:,p,q) &
                               + BCGradMat(2,3) * gradUz_master(:,p,q)
        gradUz_Face_loc(:,p,q) = BCGradMat(3,1) * gradUx_master(:,p,q) &
                               + BCGradMat(3,2) * gradUy_master(:,p,q) &
                               + BCGradMat(3,3) * gradUz_master(:,p,q)
#else
        BCGradMat(1,1) = 1. - nv(1)*nv(1)
        BCGradMat(2,2) = 1. - nv(2)*nv(2)
        BCGradMat(1,2) = -nv(1)*nv(2)
        BCGradMat(2,1) = BCGradMat(1,2)
        gradUx_Face_loc(:,p,q) = BCGradMat(1,1) * gradUx_master(:,p,q) &
                               + BCGradMat(1,2) * gradUy_master(:,p,q)
        gradUy_Face_loc(:,p,q) = BCGradMat(2,1) * gradUx_master(:,p,q) &
                               + BCGradMat(2,2) * gradUy_master(:,p,q)
        gradUz_Face_loc(:,p,q) = 0.
#endif
      END DO; END DO !p,q

      ! Evaluate 3D Diffusion Flux with interior state (with normalvel=0) and symmetry gradients
      ! Only velocities will be used from state (=inner velocities, except normal vel=0)
      CALL EvalDiffFlux3D(Nloc, UPrim_boundary,                              &
                          gradUx_Face_loc, gradUy_Face_loc, gradUz_Face_loc, &
                          Fd_Face_loc, Gd_Face_loc, Hd_Face_loc              &
#if EDDYVISCOSITY
                         ,muSGS_master(:,:,:,SideID)                         &
#endif
      )
    CASE(91)
      ! Euler/(full-)slip wall
      ! Version 2: For scalars and tangential velocity, set gradients in normal direction to zero.
      ! For velocity in wall-normal direction, set gradients in wall-tangential direction to zero.
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        nv  = NormVec( :,p,q)
        tv1 = TangVec1(:,p,q)
        tv2 = TangVec2(:,p,q)
        ! BCGradMat = I - n * n^T = (gradient - normal component of gradient)
#if (PP_dim==3)
        BCGradMat(1,1) = 1. - nv(1)*nv(1)
        BCGradMat(2,2) = 1. - nv(2)*nv(2)
        BCGradMat(3,3) = 1. - nv(3)*nv(3)
        BCGradMat(1,2) = -nv(1)*nv(2)
        BCGradMat(1,3) = -nv(1)*nv(3)
        BCGradMat(3,2) = -nv(3)*nv(2)
        BCGradMat(2,1) = BCGradMat(1,2)
        BCGradMat(3,1) = BCGradMat(1,3)
        BCGradMat(2,3) = BCGradMat(3,2)
        gradUx_Face_loc(LIFT_TEMP,p,q) = BCGradMat(1,1) * gradUx_master(LIFT_TEMP,p,q) &
                                       + BCGradMat(1,2) * gradUy_master(LIFT_TEMP,p,q) &
                                       + BCGradMat(1,3) * gradUz_master(LIFT_TEMP,p,q)
        gradUy_Face_loc(LIFT_TEMP,p,q) = BCGradMat(2,1) * gradUx_master(LIFT_TEMP,p,q) &
                                       + BCGradMat(2,2) * gradUy_master(LIFT_TEMP,p,q) &
                                       + BCGradMat(2,3) * gradUz_master(LIFT_TEMP,p,q)
        gradUz_Face_loc(LIFT_TEMP,p,q) = BCGradMat(3,1) * gradUx_master(LIFT_TEMP,p,q) &
                                       + BCGradMat(3,2) * gradUy_master(LIFT_TEMP,p,q) &
                                       + BCGradMat(3,3) * gradUz_master(LIFT_TEMP,p,q)
        ! First: Transform to gradients of wall-aligned velocities
        gradUx_vNormal = nv(1 )*gradUx_master(LIFT_VEL1,p,q)+nv(2 )*gradUx_master(LIFT_VEL2,p,q)+nv(3 )*gradUx_master(LIFT_VEL3,p,q)
        gradUx_vTang1  = tv1(1)*gradUx_master(LIFT_VEL1,p,q)+tv1(2)*gradUx_master(LIFT_VEL2,p,q)+tv1(3)*gradUx_master(LIFT_VEL3,p,q)
        gradUx_vTang2  = tv2(1)*gradUx_master(LIFT_VEL1,p,q)+tv2(2)*gradUx_master(LIFT_VEL2,p,q)+tv2(3)*gradUx_master(LIFT_VEL3,p,q)
        gradUy_vNormal = nv(1 )*gradUy_master(LIFT_VEL1,p,q)+nv(2 )*gradUy_master(LIFT_VEL2,p,q)+nv(3 )*gradUy_master(LIFT_VEL3,p,q)
        gradUy_vTang1  = tv1(1)*gradUy_master(LIFT_VEL1,p,q)+tv1(2)*gradUy_master(LIFT_VEL2,p,q)+tv1(3)*gradUy_master(LIFT_VEL3,p,q)
        gradUy_vTang2  = tv2(1)*gradUy_master(LIFT_VEL1,p,q)+tv2(2)*gradUy_master(LIFT_VEL2,p,q)+tv2(3)*gradUy_master(LIFT_VEL3,p,q)
        gradUz_vNormal = nv(1 )*gradUz_master(LIFT_VEL1,p,q)+nv(2 )*gradUz_master(LIFT_VEL2,p,q)+nv(3 )*gradUz_master(LIFT_VEL3,p,q)
        gradUz_vTang1  = tv1(1)*gradUz_master(LIFT_VEL1,p,q)+tv1(2)*gradUz_master(LIFT_VEL2,p,q)+tv1(3)*gradUz_master(LIFT_VEL3,p,q)
        gradUz_vTang2  = tv2(1)*gradUz_master(LIFT_VEL1,p,q)+tv2(2)*gradUz_master(LIFT_VEL2,p,q)+tv2(3)*gradUz_master(LIFT_VEL3,p,q)
        ! Second: Transform to gradients w.r.t. wall-aligned directions, set boundary conditions
        gradUn_vNormal  = nv( 1)*gradUx_vNormal+nv( 2)*gradUy_vNormal+nv( 3)*gradUz_vNormal
        gradUn_vTang1   = 0.!nv( 1)*gradUx_vTang1 +nv( 2)*gradUy_vTang1 +nv( 3)*gradUz_vTang1
        gradUn_vTang2   = 0.!nv( 1)*gradUx_vTang2 +nv( 2)*gradUy_vTang2 +nv( 3)*gradUz_vTang2
        gradUt1_vNormal = 0.!tv1( 1)*gradUx_vNormal+tv1( 2)*gradUy_vNormal+tv1( 3)*gradUz_vNormal
        gradUt1_vTang1  = tv1( 1)*gradUx_vTang1 +tv1( 2)*gradUy_vTang1 +tv1( 3)*gradUz_vTang1
        gradUt1_vTang2  = tv1( 1)*gradUx_vTang2 +tv1( 2)*gradUy_vTang2 +tv1( 3)*gradUz_vTang2
        gradUt2_vNormal = 0.!tv2( 1)*gradUx_vNormal+tv2( 2)*gradUy_vNormal+tv2( 3)*gradUz_vNormal
        gradUt2_vTang1  = tv2( 1)*gradUx_vTang1 +tv2( 2)*gradUy_vTang1 +tv2( 3)*gradUz_vTang1
        gradUt2_vTang2  = tv2( 1)*gradUx_vTang2 +tv2( 2)*gradUy_vTang2 +tv2( 3)*gradUz_vTang2
        ! Third: Transform back to gradients w.r.t. physical x/y/z-coordinates
        gradUx_vNormal  = nv(1)*gradUn_vNormal+tv1(1)*gradUt1_vNormal+tv2(1)*gradUt2_vNormal
        gradUx_vTang1   = nv(1)*gradUn_vTang1+ tv1(1)*gradUt1_vTang1+ tv2(1)*gradUt2_vTang1
        gradUx_vTang2   = nv(1)*gradUn_vTang2+ tv1(1)*gradUt1_vTang2+ tv2(1)*gradUt2_vTang2
        gradUy_vNormal  = nv(2)*gradUn_vNormal+tv1(2)*gradUt1_vNormal+tv2(2)*gradUt2_vNormal
        gradUy_vTang1   = nv(2)*gradUn_vTang1+ tv1(2)*gradUt1_vTang1+ tv2(2)*gradUt2_vTang1
        gradUy_vTang2   = nv(2)*gradUn_vTang2+ tv1(2)*gradUt1_vTang2+ tv2(2)*gradUt2_vTang2
        gradUz_vNormal  = nv(3)*gradUn_vNormal+tv1(3)*gradUt1_vNormal+tv2(3)*gradUt2_vNormal
        gradUz_vTang1   = nv(3)*gradUn_vTang1+ tv1(3)*gradUt1_vTang1+ tv2(3)*gradUt2_vTang1
        gradUz_vTang2   = nv(3)*gradUn_vTang2+ tv1(3)*gradUt1_vTang2+ tv2(3)*gradUt2_vTang2
        ! Forth: Transform back to gradients of velocities in physical x/y/z-coordinates
        gradUx_Face_loc(LIFT_VEL1,p,q) = nv(1)*gradUx_vNormal+tv1(1)*gradUx_vTang1+tv2(1)*gradUx_vTang2
        gradUx_Face_loc(LIFT_VEL2,p,q) = nv(2)*gradUx_vNormal+tv1(2)*gradUx_vTang1+tv2(2)*gradUx_vTang2
        gradUx_Face_loc(LIFT_VEL3,p,q) = nv(3)*gradUx_vNormal+tv1(3)*gradUx_vTang1+tv2(3)*gradUx_vTang2
        gradUy_Face_loc(LIFT_VEL1,p,q) = nv(1)*gradUy_vNormal+tv1(1)*gradUy_vTang1+tv2(1)*gradUy_vTang2
        gradUy_Face_loc(LIFT_VEL2,p,q) = nv(2)*gradUy_vNormal+tv1(2)*gradUy_vTang1+tv2(2)*gradUy_vTang2
        gradUy_Face_loc(LIFT_VEL3,p,q) = nv(3)*gradUy_vNormal+tv1(3)*gradUy_vTang1+tv2(3)*gradUy_vTang2
        gradUz_Face_loc(LIFT_VEL1,p,q) = nv(1)*gradUz_vNormal+tv1(1)*gradUz_vTang1+tv2(1)*gradUz_vTang2
        gradUz_Face_loc(LIFT_VEL2,p,q) = nv(2)*gradUz_vNormal+tv1(2)*gradUz_vTang1+tv2(2)*gradUz_vTang2
        gradUz_Face_loc(LIFT_VEL3,p,q) = nv(3)*gradUz_vNormal+tv1(3)*gradUz_vTang1+tv2(3)*gradUz_vTang2
#else
        BCGradMat(1,1) = 1. - nv(1)*nv(1)
        BCGradMat(2,2) = 1. - nv(2)*nv(2)
        BCGradMat(1,2) = -nv(1)*nv(2)
        BCGradMat(2,1) = BCGradMat(1,2)
        gradUx_Face_loc(LIFT_TEMP,p,q) = BCGradMat(1,1) * gradUx_master(LIFT_TEMP,p,q) &
                                       + BCGradMat(1,2) * gradUy_master(LIFT_TEMP,p,q)
        gradUy_Face_loc(LIFT_TEMP,p,q) = BCGradMat(2,1) * gradUx_master(LIFT_TEMP,p,q) &
                                       + BCGradMat(2,2) * gradUy_master(LIFT_TEMP,p,q)
        gradUz_Face_loc(LIFT_TEMP,p,q) = 0.
        ! First: Transform to gradients of wall-aligned velocities
        gradUx_vNormal = nv(1 )*gradUx_master(LIFT_VEL1,p,q)+nv(2 )*gradUx_master(LIFT_VEL2,p,q)
        gradUx_vTang1  = tv1(1)*gradUx_master(LIFT_VEL1,p,q)+tv1(2)*gradUx_master(LIFT_VEL2,p,q)
        gradUy_vNormal = nv(1 )*gradUy_master(LIFT_VEL1,p,q)+nv(2 )*gradUy_master(LIFT_VEL2,p,q)
        gradUy_vTang1  = tv1(1)*gradUy_master(LIFT_VEL1,p,q)+tv1(2)*gradUy_master(LIFT_VEL2,p,q)
        ! Second: Transform to gradients w.r.t. wall-aligned directions, set boundary conditions
        gradUn_vNormal  = nv( 1)*gradUx_vNormal+nv( 2)*gradUy_vNormal
        gradUn_vTang1   = 0.!nv( 1)*gradUx_vTang1 +nv( 2)*gradUy_vTang1
        gradUt1_vNormal = 0.!tv1( 1)*gradUx_vNormal+tv1( 2)*gradUy_vNormal
        gradUt1_vTang1  = tv1( 1)*gradUx_vTang1 +tv1( 2)*gradUy_vTang1
        ! Third: Transform back to gradients w.r.t. physical x/y-coordinates
        gradUx_vNormal  = nv(1)*gradUn_vNormal+tv1(1)*gradUt1_vNormal
        gradUx_vTang1   = nv(1)*gradUn_vTang1+ tv1(1)*gradUt1_vTang1
        gradUy_vNormal  = nv(2)*gradUn_vNormal+tv1(2)*gradUt1_vNormal
        gradUy_vTang1   = nv(2)*gradUn_vTang1+ tv1(2)*gradUt1_vTang1
        ! Forth: Transform back to gradients of velocities in physical x/y-coordinates
        gradUx_Face_loc(LIFT_VEL1,p,q) = nv(1)*gradUx_vNormal+tv1(1)*gradUx_vTang1
        gradUx_Face_loc(LIFT_VEL2,p,q) = nv(2)*gradUx_vNormal+tv1(2)*gradUx_vTang1
        gradUy_Face_loc(LIFT_VEL1,p,q) = nv(1)*gradUy_vNormal+tv1(1)*gradUy_vTang1
        gradUy_Face_loc(LIFT_VEL2,p,q) = nv(2)*gradUy_vNormal+tv1(2)*gradUy_vTang1
        gradUz_Face_loc(LIFT_VELV,p,q) = 0.
        gradUx_Face_loc(LIFT_VEL3,p,q) = 0.
        gradUy_Face_loc(LIFT_VEL3,p,q) = 0.
#endif
      END DO; END DO !p,q

      ! Evaluate 3D Diffusion Flux with interior state (with normalvel=0) and symmetry gradients
      ! Only velocities will be used from state (=inner velocities, except normal vel=0)
      CALL EvalDiffFlux3D(Nloc,UPrim_boundary,                            &
                          gradUx_Face_loc,gradUy_Face_loc,gradUz_Face_loc, &
                          Fd_Face_loc,Gd_Face_loc,Hd_Face_loc            &
#if EDDYVISCOSITY
                         ,muSGS_master(:,:,:,SideID)&
#endif
      )
    END SELECT

    ! Sum up Euler and Diffusion Flux
    DO iVar=2,PP_nVar
      Flux(iVar,:,:) = Flux(iVar,:,:) &
                     + NormVec(1,:,:)*Fd_Face_loc(iVar,:,:) &
                     + NormVec(2,:,:)*Gd_Face_loc(iVar,:,:) &
                     + NormVec(3,:,:)*Hd_Face_loc(iVar,:,:)
    END DO ! iVar
#endif /*PARABOLIC*/

  CASE(1)
    ! Periodic already filled!
  CASE DEFAULT ! unknown BCType
    CALL Abort(__STAMP__, 'no BC defined in navierstokes/getboundaryflux.f90!')
  END SELECT
#if TESTCASE_BC
END IF ! BCType < 0
#endif /*TESTCASE_BC*/

END SUBROUTINE GetBoundaryFlux


#if FV_ENABLED && FV_RECONSTRUCT
!==================================================================================================================================
!> Computes the gradient at a boundary for FV subcells.
!==================================================================================================================================
SUBROUTINE GetBoundaryFVgradient(SideID,t,gradU,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP,sdx_Face)
! MODULES
USE MOD_PreProc
USE MOD_Globals       ,ONLY: Abort
USE MOD_DG_Vars       ,ONLY: UPrim_Boundary
USE MOD_Mesh_Vars     ,ONLY: BoundaryType,BC
#if TESTCASE_BC
USE MOD_TestCase      ,ONLY: GetBoundaryFVgradientTestcase
#endif /*TESTCASE_BC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID                                      !< ID of current side
REAL,INTENT(IN)   :: t                                           !< current time (provided by time integration scheme)
REAL,INTENT(IN)   :: UPrim_master(       PRIM,0:PP_N,0:PP_NZ)    !< primitive solution from the inside
REAL,INTENT(OUT)  :: gradU       (       PRIM,0:PP_N,0:PP_NZ)    !< gradient at boundary for FV subcells
REAL,INTENT(IN)   :: NormVec (              3,0:PP_N,0:PP_NZ)    !< normal vector on surfaces
REAL,INTENT(IN)   :: TangVec1(              3,0:PP_N,0:PP_NZ)    !< tangential1 vector on surfaces
REAL,INTENT(IN)   :: TangVec2(              3,0:PP_N,0:PP_NZ)    !< tangential2 vector on surfaces
REAL,INTENT(IN)   :: Face_xGP(              3,0:PP_N,0:PP_NZ)    !< positions of surface flux points
REAL,INTENT(IN)   :: sdx_Face(                0:PP_N,0:PP_NZ,3)  !< distance between center of FV-cell and boundary
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: p,q
INTEGER           :: BCType,BCState
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

#if TESTCASE_BC
IF (BCType.LT.0) THEN ! testcase boundary condition
  CALL GetBoundaryFVgradientTestcase(SideID,t,gradU,UPrim_master)
ELSE
#endif /*TESTCASE_BC*/
  CALL GetBoundaryState(SideID,t,PP_N,UPrim_boundary,UPrim_master,&
      NormVec,TangVec1,TangVec2,Face_xGP)
  SELECT CASE(BCType)
  CASE(2,3,4,9,91,12,121,22,23,24,25,27)
    DO q=0,PP_NZ; DO p=0,PP_N
      gradU(:,p,q) = (UPrim_master(:,p,q) - UPrim_boundary(:,p,q)) * sdx_Face(p,q,3)
    END DO; END DO ! p,q=0,PP_N
  CASE(1) !Periodic already filled!
  CASE DEFAULT ! unknown BCType
    CALL Abort(__STAMP__,&
         'no BC defined in navierstokes/getboundaryflux.f90!')
  END SELECT
#if TESTCASE_BC
END IF ! BCType < 0
#endif /*TESTCASE_BC*/

END SUBROUTINE GetBoundaryFVgradient
#endif


#if PARABOLIC
!==================================================================================================================================
!> Computes the boundary fluxes for the lifting procedure for a given Cartesian mesh face (defined by SideID).
!==================================================================================================================================
SUBROUTINE Lifting_GetBoundaryFlux(SideID,t,UPrim_master,Flux,NormVec,TangVec1,TangVec2,Face_xGP,SurfElem)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_DG_Vars      ,ONLY: UPrim_Boundary
USE MOD_Mesh_Vars    ,ONLY: BoundaryType,BC
USE MOD_Lifting_Vars ,ONLY: doWeakLifting
#if TESTCASE_BC
USE MOD_TestCase     ,ONLY: Lifting_GetBoundaryFluxTestcase
#endif /*TESTCASE_BC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID                                   !< ID of current side
REAL,INTENT(IN)   :: t                                        !< current time (provided by time integration scheme)
REAL,INTENT(IN)   :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ) !< primitive solution from the inside
REAL,INTENT(OUT)  :: Flux(     PP_nVarLifting,0:PP_N,0:PP_NZ) !< lifting boundary flux
REAL,INTENT(IN)   :: NormVec (              3,0:PP_N,0:PP_NZ) !< normal vector on surfaces
REAL,INTENT(IN)   :: TangVec1(              3,0:PP_N,0:PP_NZ) !< tangential1 vector on surfaces
REAL,INTENT(IN)   :: TangVec2(              3,0:PP_N,0:PP_NZ) !< tangential2 vector on surfaces
REAL,INTENT(IN)   :: Face_xGP(              3,0:PP_N,0:PP_NZ) !< positions of surface flux points
REAL,INTENT(IN)   :: SurfElem(                0:PP_N,0:PP_NZ) !< surface element to multiply with flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: p,q
INTEGER           :: BCType,BCState
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

#if TESTCASE_BC
IF (BCType.LT.0) THEN ! testcase boundary conditions
  CALL Lifting_GetBoundaryFluxTestcase(SideID,t,UPrim_master,Flux)
ELSE
#endif /*TESTCASE_BC*/
  CALL GetBoundaryState(SideID,t,PP_N,UPrim_boundary,UPrim_master,&
                        NormVec,TangVec1,TangVec2,Face_xGP)
  SELECT CASE(BCType)
  CASE(2,12,121,22,23,24,25,27) ! Riemann solver based BCs
      Flux = 0.5*(UPrim_master(PRIM_LIFT,:,:)+UPrim_boundary(PRIM_LIFT,:,:))
  CASE(3,4) ! No-slip wall BCs
    DO q=0,PP_NZ; DO p=0,PP_N
#if PP_OPTLIFT == 0
      Flux(LIFT_DENS,p,q) = UPrim_Boundary(DENS,p,q)
      Flux(LIFT_VELV,p,q) = 0.
      Flux(LIFT_TEMP,p,q) = UPrim_Boundary(TEMP,p,q)
#else
      Flux(LIFT_VELV,p,q) = 0.
      Flux(LIFT_TEMP,p,q) = UPrim_Boundary(TEMP,p,q)
#endif
    END DO; END DO !p,q
  CASE(9,91)
    ! Euler/(full-)slip wall, symmetry BC
    ! Solution from the inside with velocity normal component set to 0 (done in GetBoundaryState)
    DO q=0,PP_NZ; DO p=0,PP_N
      ! Compute Flux
#if PP_OPTLIFT == 0
      Flux(LIFT_DENS,p,q) = UPrim_master(  DENS,p,q)
      Flux(LIFT_VELV,p,q) = UPrim_boundary(VELV,p,q)
      Flux(LIFT_TEMP,p,q) = UPrim_master(  TEMP,p,q)
#else
      Flux(LIFT_VELV,p,q) = UPrim_boundary(VELV,p,q)
      Flux(LIFT_TEMP,p,q) = UPrim_master(  TEMP,p,q)
#endif
    END DO; END DO !p,q
  CASE(1) !Periodic already filled!
  CASE DEFAULT ! unknown BCType
    CALL Abort(__STAMP__,&
         'no BC defined in navierstokes/getboundaryflux.f90!')
  END SELECT

  ! in case lifting is done in strong form
  IF(.NOT.doWeakLifting) Flux=Flux-UPrim_master(PRIM_LIFT,:,:)

  DO q=0,PP_NZ; DO p=0,PP_N
    Flux(:,p,q)=Flux(:,p,q)*SurfElem(p,q)
  END DO; END DO
#if TESTCASE_BC
END IF
#endif /*TESTCASE_BC*/

END SUBROUTINE Lifting_GetBoundaryFlux
#endif /*PARABOLIC*/


!==================================================================================================================================
!> Read in a HDF5 file containing the state for a boundary. Used in BC Type 12.
!==================================================================================================================================
SUBROUTINE ReadBCFlow(FileName)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY:BCData,BCDataPrim
USE MOD_Mesh_Vars         ,ONLY:offsetElem,nElems,nBCSides,S2V2,SideToElem,nGlobalElems
USE MOD_HDF5_Input        ,ONLY:OpenDataFile,GetDataProps,CloseDataFile,ReadAttribute,ReadArray
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_ProlongToFace     ,ONLY:EvalElemFace
USE MOD_Interpolation_Vars,ONLY:NodeType
#if (PP_NodeType==1)
USE MOD_Interpolation_Vars,ONLY:L_minus,L_plus
#endif
USE MOD_ChangeBasisByDim  ,ONLY:ChangeBasisVolume
USE MOD_EOS               ,ONLY:ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileName       !< name of file BC data is read from
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,POINTER                  :: U_N(:,:,:,:,:)=>NULL()
REAL,ALLOCATABLE,TARGET       :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE              :: Vdm_NHDF5_N(:,:)
REAL                          :: Uface(PP_nVar,0:PP_N,0:PP_NZ)
INTEGER                       :: nVar_HDF5,N_HDF5,nElems_HDF5,N_HDF5Z
INTEGER                       :: p,q,SideID,ElemID,locSide
CHARACTER(LEN=255)            :: NodeType_HDF5
LOGICAL                       :: InterpolateSolution
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A,A)')'  Read BC state from file "',TRIM(FileName)
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5)

IF(nElems_HDF5.NE.nGlobalElems)THEN
  CALL Abort(__STAMP__,&
             'BaseFlow file does not match solution. Elements',nElems_HDF5)
END IF

#if (PP_dim==2)
N_HDF5Z=0
#else
N_HDF5Z=N_HDF5
#endif
ALLOCATE(U_local(PP_nVar,0:N_HDF5,0:N_HDF5,0:N_HDF5Z,nElems))
CALL ReadArray('DG_Solution',5,(/PP_nVar,N_HDF5+1,N_HDF5+1,N_HDF5Z+1,nElems/),OffsetElem,5,RealArray=U_local)
CALL CloseDataFile()

! Read in state
InterpolateSolution=((N_HDF5.NE.PP_N) .OR. (TRIM(NodeType_HDF5).NE.TRIM(NodeType)))
IF(.NOT. InterpolateSolution)THEN
  ! No interpolation needed, read solution directly from file
  U_N=>U_local
ELSE
  ! We need to interpolate the solution to the new computational grid
  ALLOCATE(U_N(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  ALLOCATE(Vdm_NHDF5_N(0:PP_N,0:N_HDF5))
  CALL GetVandermonde(N_HDF5,NodeType_HDF5,PP_N,NodeType,Vdm_NHDF5_N,modal=.TRUE.)

  SWRITE(UNIT_stdOut,*)'Interpolate base flow from restart grid with N=',N_HDF5,' to computational grid with N=',PP_N
  DO ElemID=1,nElems
    CALL ChangeBasisVolume(PP_nVar,N_HDF5,PP_N,Vdm_NHDF5_N,U_local(:,:,:,:,ElemID),U_N(:,:,:,:,ElemID))
  END DO ! ElemID
  DEALLOCATE(Vdm_NHDF5_N)
END IF

! Prolong boundary state
DO SideID=1,nBCSides
  ElemID  = SideToElem(S2E_ELEM_ID    ,SideID)
  locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)

#if (PP_NodeType==1)
  CALL EvalElemFace(PP_nVar,PP_N,U_N(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)
#else
  CALL EvalElemFace(PP_nVar,PP_N,U_N(:,:,:,:,ElemID),Uface,locSide)
#endif
  DO q=0,PP_NZ; DO p=0,PP_N
    BCData(:,p,q,SideID)=Uface(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
    CALL ConsToPrim(BCDataPrim(:,p,q,SideID),BCData(:,p,q,SideID))
  END DO; END DO
END DO

IF(InterpolateSolution) DEALLOCATE(U_N)
DEALLOCATE(U_local)

SWRITE(UNIT_stdOut,'(A)')'  done initializing BC state!'
END SUBROUTINE ReadBCFlow


!==================================================================================================================================
!> Finalize arrays used for boundary conditions.
!==================================================================================================================================
SUBROUTINE FinalizeBC()
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,BCDataPrim,nBCByType,BCSideID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(BCData)
SDEALLOCATE(BCDataPrim)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC

END MODULE MOD_GetBoundaryFlux
