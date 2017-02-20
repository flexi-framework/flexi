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
! Define variables for normal and extended state vector
! Normal   U(1:5)  with conservative variables
! Extended U(1:11) with conservative and primitive variables

#define PP_2Var 11

#define CONS 1:5  /* conservative variables */
#define PRIM 6:11 /* primitive variables */

! conservative variables
#define DENS  1   /* density */
#define MOM1  2   /* momentum x */
#define MOM2  3   /* momentum y */
#define MOM3  4   /* momentum z */
#define MOMV  2:4 /* momentum vector */
#define ENER  5   /* energy */

! primitive (exteded) variables
#define SRHO  6   /* specific volume (1./density) */
#define VEL1  7   /* velocity x */
#define VEL2  8   /* velocity y */
#define VEL3  9   /* velocity z */
#define VELV  7:9 /* velocity range */
#define PRES  10  /* pressure */
#define TEMP  11  /* temperature */


! routines to compute physical quantities 
#define KAPPASPR_MAX_TIMESTEP_H()      (MAX(4./3.,KappasPr))
#define THERMAL_CONDUCTIVITY_H(mu)     (mu*cp/Pr)
#define TOTAL_TEMPERATURE_H(T,Mach)    (T*(1+0.5*(kappa-1)*Mach**2))
#define TOTAL_PRESSURE_H(p,Mach)       (p/((1+0.5*(kappa-1)*Mach**2)**(-kappa/(kappa-1.))))
#define BETA_RIEMANN_H()               (SQRT(0.5*kappaM1/kappa))
#define ROEC_RIEMANN_H(RoeH,RoeVel)    (SQRT(kappaM1*(RoeH-0.5*DOT_PRODUCT(RoeVel,RoeVel))))
#define ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)     (kappaM1/(Roec*Roec) * (Delta_U(1)*(RoeH-RoeVel(1)*RoeVel(1)) - Delta_U(6) + RoeVel(1)*Delta_U(2)))

! routines to compute physical quantities from conservative variables or extended variables
! conservative
#define VELOCITY_H(U,sRho)             (U(MOMV)*sRho)
#define SPEEDOFSOUND_H(p,sRho)         (SQRT(Kappa*p*sRho))
#define TOTALENERGY_H(U,sRho,Vel)      (U(ENER)/U(DENS))
#define TOTALENTHALPY_H(U,p,sRho)      ((U(ENER)+p)*sRho)
#define ENTROPY_H(U,T)                 (R*(sKappaM1*LOG(T)-LOG(U(DENS))))
#define TEMPERATURE_H(U)               ((U(ENER)-0.5*DOT_PRODUCT(U(MOMV),U(MOMV))/U(DENS))/(U(DENS)*cv))

! extended (NOTE: compute from cons. When computing derived (neither prim or cons) variables
! assume that both prim and cons vars are filled
#define VELOCITY_HE(UE)                (UE(MOMV)*UE(SRHO))
#define PRESSURE_HE(UE)                (KappaM1*(UE(ENER)-0.5*DOT_PRODUCT(UE(VELV),UE(MOMV))))
#define SPEEDOFSOUND_HE(UE)            (SQRT(Kappa*UE(PRES)*UE(SRHO)))
#define TOTALENERGY_HE(UE)             (UE(ENER)*UE(SRHO))
#define TOTALENTHALPY_HE(UE)           ((UE(ENER)+UE(PRES))*UE(SRHO))
#define TEMPERATURE_HE(UE)             (UE(PRES)*UE(SRHO)/R)
#define ENERGY_HE(UE)                  (sKappaM1*UE(PRES)+0.5*DOT_PRODUCT(UE(MOMV),UE(VELV)))

#if PP_VISC == 0
#define VISCOSITY_PRIM(U)              mu0
#elif PP_VISC == 1
#define VISCOSITY_PRIM(U)              muSuth(U(6))
#elif PP_VISC == 2
#define VISCOSITY_PRIM(U)              mu0*U(6)**ExpoSuth
#endif

#if PP_VISC == 0
#define VISCOSITY_TEMPERATURE(T)       mu0
#elif PP_VISC == 1
#define VISCOSITY_TEMPERATURE(T)       muSuth(T)
#elif PP_VISC == 2
#define VISCOSITY_TEMPERATURE(T)       mu0*T**ExpoSuth
#endif
