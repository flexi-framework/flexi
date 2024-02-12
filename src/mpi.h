!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
!===================================================================================================================================
! Here, MPI variables should be defined
!===================================================================================================================================

#if USE_MPI
#ifdef MPI_F08
! MPI implementation
#define __MPI__ mpi_f08

! MPI data types
#define MPI_TYPE_COMM    TYPE(mpi_comm)
#define MPI_TYPE_INFO    TYPE(mpi_info)
#define MPI_TYPE_REQUEST TYPE(mpi_request)
#define MPI_TYPE_STATUS  TYPE(mpi_status)

! HDF5 still uses integers to handle communicators
#define H5PSET_FAPL_MPIO_F_(list,comm,info,error) H5PSET_FAPL_MPIO_F(list,comm%mpi_val,info%mpi_val,error)
#else
! MPI implementation
#define __MPI__ mpi

! MPI data types
#define MPI_TYPE_COMM    INTEGER
#define MPI_TYPE_INFO    INTEGER
#define MPI_TYPE_REQUEST INTEGER
#define MPI_TYPE_STATUS  INTEGER

! HDF5 still uses integers to handle communicators
#define H5PSET_FAPL_MPIO_F_(list,comm,info,error) H5PSET_FAPL_MPIO_F(list,comm        ,info        ,error)
#endif /*MPI_F08*/
#endif /*USE_MPI*/
