/*
!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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
*/

#ifndef VISU3D_PLUGIN_H
#define VISU3D_PLUGIN_H

#include "posti_pluginTypes.h"

extern "C" {
  extern void __mod_visu3d_MOD_visu3d_requestinformation(int* mpi_comm_IN, int* str_len, const char* state_file, struct CharARRAY* varnames);
}

extern "C" {
  extern void __mod_visu3d_MOD_visu3d_cwrapper(int* mpi_comm_IN, 
        int* strlen_prm,   const char* prmfile_IN, 
        int* strlen_posti, const char* postifile_IN, 
        int* strlen_state, const char* statefile_IN,
        struct DoubleARRAY* coords_out,  struct DoubleARRAY* values_out, 
        struct IntARRAY* nodeids_out,    struct DoubleARRAY* coordsFV_out,
        struct DoubleARRAY* valuesFV_out,struct IntARRAY* nodeidsFV_out,
        struct CharARRAY* varnames_out,  struct IntARRAY* components_out);
}

extern "C" {
  extern void __mod_visu3d_MOD_visu3d_dealloc_nodeids();
}

extern "C" {
  extern void __mod_visu3d_MOD_visu3d_finalize();
}

#endif /* VISU3D_PLUGIN_H */
