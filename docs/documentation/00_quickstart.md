# Quick Start Guide

This quick start guide allows for a fast installation and setup of **FLEXI** without diving into the general details of the framework, compile options and features. Further information and detailed descriptions are available by following the indicated references.

## Installation and Setup

**FLEXI** is free and open source (GPLv3). The current version of the code-framework is available online and can be acquired from the [GitHub repository](https://github.com/flexi-framework/flexi.git) by either cloning it or downloading the compressed folder, see {ref}`sec:obtain_flexi`.

**FLEXI** requires the following [packages](sec:installing_dependencies) to be installed on the system:

- git
- CMake
- Fortran and C/C++ compilers (GNU compilers recommended)
- MPI libraries (OpenMPI recommended)
- LAPACK/OpenBLAS
- HDF5
- FFTW

LAPACK/OpenBLAS, HDF5, and FFTW can be automatically installed by enabling the corresponding compiler options in the CMake configuration.

**FLEXI** is compiled using CMake. For a compilation in default configuration use:
```bash
mkdir build
cmake -B build
cmake --build build
```
Custom configurations can be generated with
```bash
ccmake -B build
```
including the installation of the third-party packages mentioned above (LAPACK/OpenBLAS, HDF5, FFTW), see section [compiler options](sec:code_options) for a detailed list. The executables will be generated in `./build/bin/`.

## Mesh Generation

For the generation of high-order meshes the standalone mesh generator **HOPR** is required, creating **FLEXI** compatible mesh files in HDF5 format. Simple, structured meshes can be directly generated in **HOPR** using the integrated mesh generator, while the processing of complex geometries can be a based on external meshes in CGNS or GMSH format. In any case, a parameter file for the mesh generation and modification is required. **HOPR** is available on [GitHub](https://github.com/hopr-framework/hopr/releases) and can be compiled using CMake or by simply downloading the provided AppImage.
For an in-depth description we refer to the [**HOPR** documentation](https://hopr.readthedocs.io/en/latest/).

## Running FLEXI

**FLEXI** can be run by executing the generated binary in the build folder `$FLEXIROOT/build/bin/flexi` and providing a parameter file with the simulation-specific definitions. The [feature list](sec:feature_list) provides an overview of the various features implemented in **FLEXI**, while section [parameter file](sec:parameter_file) contains all options and a short description.
For the definition of the initial and boundary conditions we refer to sections {ref}`subsec:ic` and {ref}`subsec:bc`, respectively.
Optionally, simulations can be restarted from an existing state file (i.e. volume solution) by appending it to the argument vector:

```bash
flexi parameter_flexi.ini [Restart_State.h5]
```
Further details concerning the capabilities of **FLEXI** and the application to small testcases, including e.g., the flow around a [NACA0012](NACA0012) airfoil, are included in the [tutorials](Tutorials).

## Tools
**FLEXI** comes with a comprehensive [postprocessing](ToolsOverview) toolchain, such as, e.g., the [interpolation](sec:swap_mesh) between different meshes, the [time averaging](sec:time_averaging) of solution files, and the [animation](sec:animate_tool). Most importantly, it includes the `posti_visu` tool to convert the solution files from the custom *h5* format to *vtu* files readable by **ParaView**, as covered in the [workflow section](subsec:post_processing).
