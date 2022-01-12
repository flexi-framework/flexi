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

#ifndef VISU_PLUGIN_H
#define VISU_PLUGIN_H

#include "pluginTypes_visu.h"

extern "C" {
  extern void __mod_visu_cwrapper_MOD_visu_requestinformation(int* mpi_comm_IN,
        int* str_len, const char* state_file, int* str_len_mesh, const char* mesh_file,  struct CharARRAY* varnames, struct CharARRAY* bcnames);
}

extern "C" {
  extern void __mod_visu_cwrapper_MOD_visu_cwrapper(int* mpi_comm_IN,
        int* HighOrder,
        int* strlen_prm,   const char* prmfile_IN,
        int* strlen_posti, const char* postifile_IN,
        int* strlen_state, const char* statefile_IN,
        struct DoubleARRAY* coords_out,  struct DoubleARRAY* values_out,
        struct IntARRAY* nodeids_out,    struct DoubleARRAY* coordsFV_out,
        struct DoubleARRAY* valuesFV_out,struct IntARRAY* nodeidsFV_out,
        struct CharARRAY* varnames_out,
        struct DoubleARRAY* coordsSurf_out,  struct DoubleARRAY* valuesSurf_out,
        struct IntARRAY* nodeidsSurf_out,    struct DoubleARRAY* coordsSurfFV_out,
        struct DoubleARRAY* valuesSurfFV_out,struct IntARRAY* nodeidsSurfFV_out,
        struct CharARRAY* varnamesSurf_out);
}

extern "C" {
  extern void __mod_visu_cwrapper_MOD_visu_dealloc_nodeids();
}

extern "C" {
  extern void __mod_visu_MOD_finalizevisu();
}

#endif /* VISU_PLUGIN_H */
