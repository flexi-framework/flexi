!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
!==================================================================================================================================
!> Contains global variables used by the FV modules.
!==================================================================================================================================
MODULE MOD_FV_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
#if FV_ENABLED
LOGICAL                :: FVInitIsDone=.FALSE.
LOGICAL                :: FVInitBasisIsDone=.FALSE.
! parameters
REAL                   :: FV_IndUpperThreshold   !< Upper threshold: Element is switched from DG to FV if indicator
                                                 !< rises above this value
REAL                   :: FV_IndLowerThreshold   !< Lower threshold: Element is switched from FV to DG if indicator
                                                 !< falls below this value

#if FV_ENABLED == 1
LOGICAL                :: FV_toDG_indicator      !< additional Persson indicator applied to DG solution after switch from FV to DG
                                                 !< to check if DG solution is valid
REAL                   :: FV_toDG_limit          !< limit for ^ this indicator: If FV_toDG_indicator is above limit, keep FV
LOGICAL                :: FV_toDGinRK            !< Flag that allows switching of FV elements to DG during Runge Kutta stages.
                                                 !< This may violated the DG timestep restriction of the element.
LOGICAL                :: FV_IniSharp            !< Maintain a sharp interface in the initial solution in the FV region
LOGICAL                :: FV_IniSupersample      !< Supersample initial solution inside each sub-cell
#endif /*FV_ENABLED*/
LOGICAL                :: switchConservative     !< Perform DG/FV switch in reference element

! Limiting
INTEGER                :: LimiterType            !< Readin variable for type of used fv limiter
REAL                   :: FV_sweby_beta          !< parameter for Sweby limiter

! FV/DG Switching
! TODO: The following variables are only need for Switch-Type FV. Could be hidden behind correct preprocessor flags in the future
INTEGER,ALLOCATABLE    :: FV_Elems(:)            !< indicates if DG element (0) or FV subcells (1) for each element
INTEGER,ALLOCATABLE    :: FV_Elems_master(:)     !< prolongate FV_Elems to faces
INTEGER,ALLOCATABLE    :: FV_Elems_slave(:)
INTEGER,ALLOCATABLE    :: FV_Elems_Sum(:)        !< = FV_Elems_master + 2*FV_Elems_slave
                                                 !< array with values from 0..3: 0=both DG, 1=master FV, 2=slave FV, 3=both FV
INTEGER,ALLOCATABLE    :: FV_Elems_counter(:)    !< counts for every element the RK stages it was DG or FV, nullified after
                                                 !< hdf5-output
INTEGER                :: FV_Switch_counter      !< counts how often FV_Switch is called, nullified after hdf5-output
                                                 !< should be identical to nTimesteps * nRkStages

REAL,ALLOCATABLE       :: FV_Elems_Amount(:)     !< counts for every element the RK stages it was DG or FV, nullified after

! FV/DG Blending
#if FV_ENABLED == 2
REAL,ALLOCATABLE       :: FV_alpha(:)            !< Blending coefficient
REAL,ALLOCATABLE       :: FV_alpha_master(:)     !< Prolongated blending coefficient on master sides
REAL,ALLOCATABLE       :: FV_alpha_slave( :)     !< Prolongated blending coefficient on slave  sides
REAL                   :: FV_alpha_min           !< Minimal blending coefficient (all elems below are treated as pure DG)
REAL                   :: FV_alpha_max           !< Maximal blending coefficient
REAL                   :: FV_alpha_fix           !< fixed blending paramter
REAL                   :: FV_alpha_extScale      !< Amount alpha is scaled when extended into neighbouring elements.
LOGICAL                :: FV_doExtendAlpha       !< Flag whether alpha should be extended into neighbouring elements
INTEGER                :: FV_nExtendAlpha        !< Number of times alpha should be passed towards neighboring elements per timestep
#endif /*FV_ENABLED==2*/

! FV variables on reference element
REAL,ALLOCATABLE       :: FV_X(:)                !< positions of 'midpoints' of FV subcells in [-1,1]
REAL,ALLOCATABLE       :: FV_BdryX(:)            !< positions of boundaries of FV subcells in [-1,1]
REAL,ALLOCATABLE       :: FV_w(:)                !< weights of FV subcells (lenght of subcell)
REAL,ALLOCATABLE       :: FV_w_inv(:)            !< 1/FV_w
REAL,ALLOCATABLE       :: FV_Vdm(:,:)            !< Vandermonde to switch from DG to FV
REAL,ALLOCATABLE       :: FV_sVdm(:,:)           !< Vandermonde to switch from FV to DG
INTEGER                :: FV_CellType            !< Type of FV Cell: -1 = SAME              ,0 = EQUIDISTANT
                                                 !<                   1 = LEGENDRE_GAUSS    ,2 = LEGENDRE_LOBATTO  
                                                 !<                   3 = CHEBYSHEV_LOBATTO

#if FV_RECONSTRUCT
REAL,ALLOCATABLE,TARGET:: FV_surf_gradU(:,:,:,:) !< FD over DG interface
REAL,ALLOCATABLE,TARGET:: FV_multi_master(:,:,:,:)      !< multipurpose array: contains:
REAL,ALLOCATABLE,TARGET:: FV_multi_slave(:,:,:,:)       !< - DG: first inner value of U next to the face
                                                        !< - FV: first inner gradient from points next and second next to face

REAL,ALLOCATABLE       :: FV_sdx_Face (:,:,:,:)  !< 1 / FV_dx_Face
REAL,ALLOCATABLE       :: FV_sdx_XI   (:,:,:,:)  !< 1 / FV_dx_XI
REAL,ALLOCATABLE       :: FV_sdx_ETA  (:,:,:,:)  !< 1 / FV_dx_ETA
REAL,ALLOCATABLE       :: FV_sdx_ZETA (:,:,:,:)  !< 1 / FV_dx_ZETA
REAL,ALLOCATABLE,TARGET:: FV_dx_XI_L  (:,:,:,:)  !< distance between support point of FV sub-cell and the left boundary of that
                                                 !< sub-cell in XI direction
REAL,ALLOCATABLE,TARGET:: FV_dx_ETA_L (:,:,:,:)  !< distance between support point of FV sub-cell and the left boundary of that
                                                 !< sub-cell in ETA direction
REAL,ALLOCATABLE,TARGET:: FV_dx_ZETA_L(:,:,:,:)  !< distance between support point of FV sub-cell and the left boundary of that
                                                 !< sub-cell in ZETA direction
REAL,ALLOCATABLE,TARGET:: FV_dx_XI_R  (:,:,:,:)  !< distance between support point of FV sub-cell and the right boundary of that
                                                 !< sub-cell in XI direction
REAL,ALLOCATABLE,TARGET:: FV_dx_ETA_R (:,:,:,:)  !< distance between support point of FV sub-cell and the right boundary of that
                                                 !< sub-cell in ETA direction
REAL,ALLOCATABLE,TARGET:: FV_dx_ZETA_R(:,:,:,:)  !< distance between support point of FV sub-cell and the right boundary of that
                                                 !< sub-cell in ZETA direction
REAL,ALLOCATABLE       :: FV_dx_slave (:,:,:,:)  !< contains FV_dx_XI_L/FV_dx_XI_R or ETA or ZETA of the slave side element
                                                 !< but in face coordinate system
REAL,ALLOCATABLE       :: FV_dx_master(:,:,:,:)  !< contains FV_dx_XI_L/FV_dx_XI_R or ETA or ZETA of the master side element
                                                 !< but in face coordinate system

! arrays for gradients and metrics
REAL,ALLOCATABLE,TARGET:: gradUxi  (:,:,:,:,:)        !< FD in XI direction (PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL,ALLOCATABLE,TARGET:: gradUeta (:,:,:,:,:)        !< FD in ETA direction
REAL,ALLOCATABLE,TARGET:: gradUzeta(:,:,:,:,:)        !< FD in ZETA direction
#if PARABOLIC
REAL,ALLOCATABLE       :: gradUxi_central  (:,:,:,:,:)!< FD in XI direction (central limited for viscous fluxes)
REAL,ALLOCATABLE       :: gradUeta_central (:,:,:,:,:)!< FD in ETA direction (central limited for viscous fluxes)
REAL,ALLOCATABLE       :: gradUzeta_central(:,:,:,:,:)!< FD in ZETA direction (central limited for viscous fluxes)
#endif
#endif

! Metric terms for FV Volint
REAL,ALLOCATABLE       :: FV_SurfElemXi_sw(:,:,:,:)  !< SurfElem for inner FV faces in XI direction  (1:PP_N,0:PP_N,0:PP_N,PP_nElems)
REAL,ALLOCATABLE       :: FV_SurfElemEta_sw(:,:,:,:) !< SurfElem for inner FV faces in ETA direction (0:PP_N,1:PP_N,0:PP_N,PP_nElems)
REAL,ALLOCATABLE       :: FV_SurfElemZeta_sw(:,:,:,:)!< SurfElem for inner FV faces in ZETA direction(0:PP_N,0:PP_N,1:PP_N,PP_nElems)
REAL,ALLOCATABLE       :: FV_NormVecXi (:,:,:,:,:)   !< Normal vector for inner FV faces in XI direction
REAL,ALLOCATABLE       :: FV_TangVec1Xi(:,:,:,:,:)   !< Tangent1 vector for inner FV faces in XI direction
REAL,ALLOCATABLE       :: FV_TangVec2Xi(:,:,:,:,:)   !< Tangent2 vector for inner FV faces in XI direction
REAL,ALLOCATABLE       :: FV_NormVecEta (:,:,:,:,:)  !< Normal vector for inner FV faces in ETA direction
REAL,ALLOCATABLE       :: FV_TangVec1Eta(:,:,:,:,:)  !< Tangent1 vector for inner FV faces in ETA direction
REAL,ALLOCATABLE       :: FV_TangVec2Eta(:,:,:,:,:)  !< Tangent2 vector for inner FV faces in ETA direction
REAL,ALLOCATABLE       :: FV_NormVecZeta (:,:,:,:,:) !< Normal vector for inner FV faces in ZETA direction
REAL,ALLOCATABLE       :: FV_TangVec1Zeta(:,:,:,:,:) !< Tangent1 vector for inner FV faces in ZETA direction
REAL,ALLOCATABLE       :: FV_TangVec2Zeta(:,:,:,:,:) !< Tangent2 vector for inner FV faces in ZETA direction

#if PARABOLIC
REAL,ALLOCATABLE       :: FV_Metrics_fTilde_sJ(:,:,:,:,:) !< Metrics for FV subcells multiplied with sJ
REAL,ALLOCATABLE       :: FV_Metrics_gTilde_sJ(:,:,:,:,:) !< Metrics for FV subcells multiplied with sJ
REAL,ALLOCATABLE       :: FV_Metrics_hTilde_sJ(:,:,:,:,:) !< Metrics for FV subcells multiplied with sJ
#endif

#endif /* FV_ENABLED */
!==================================================================================================================================
END MODULE MOD_FV_Vars
