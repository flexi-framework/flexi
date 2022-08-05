\hypertarget{codeoptions}{}

# Code options \label{chap:codeoptions}
## Compiler options \label{sec:compileroptions}
This section describes the main configuration options which can be set when building **FLEXI** using CMake. 
Some options are dependent on others being enabled (or disabled), so the available ones may change. 

The first set of options describe general CMake behaviour:

* ``CMAKE_BUILD_TYPE``:

    This statically specifies what build type (configuration) will be built in this build tree. Possible values are
    * **Release**
    
        "Normal" execution.
    
    * **Profile**
    
        Performance profiling using gprof.
    
    * **Debug**
    
        Debug compiler for detailed error messages during code development.
    
* ``CMAKE_HOSTNAME``:

    This will display the host name of the machine you are compiling on.

* ``CMAKE_INSTALL_PREFIX``:

    If “make install” is invoked or INSTALL is built, this directory is prepended onto all install directories. This variable defaults to /usr/local on UNIX.

For some external libraries and programs that **FLEXI** uses, the following options apply:

* ``CTAGS_PATH``:

    This variable specifies the Ctags install directory, an optional program used to jump between tags in the source file.

* ``LIBS_BUILD_HDF5``: ON/OFF

    This will be set to ON if no prebuilt HDF5 installation was found on your machine. In this case a HDF5 version will be build and used instead.

* ``HDF5_DIR``:

    If you want to use a prebuilt HDF5 library that has been build using the CMake system, this directory should contain the CMake configuration file for HDF5 (optional).

The following options enable or disable specific features of **FLEXI**. If you want to use a certain feature, make sure to enable it during the build process!

* ``FLEXI_2D``: ON/OFF

    If set to ON the code will run in two-dimensional mode. You have to provide a mesh that consists of only one layer of elements in the third dimension.

* ``FLEXI_EQNSYSNAME``:
    
    This variable defines the equation system, which will be compiled and used for the simulation. Possible values are
    * *Navierstokes*
    * *Linearscalaradvection*
    
* ``FLEXI_FV``:

    Set this to enable the corresponding finite volume subcell shock capturing mechanism. Implemented are
    * OFF    : No finite volume shock capturing
    * SWITCH : Switching DG elements into a finite volume subcell representation [@sonntag2017efficient]
    * BLEND  : Blending the finite volume discretization operator to the DG operator. [@hennemann2021provably]
    
* ``FLEXI_FV_RECONSTRUCTION``:  ON/OFF
    
    Only available if FLEXI_FV is set either to SWITCH or BLEND. Enables the reconstruction of interface values in the finite volume subcells. Needed for calculation of gradients and to use a second order finite volume scheme.
    
* ``FLEXI_LIFTING``:

    Two different lifting methods for the parabolic part of the equation system are implemented. Possible values are
    * **BR1**: First method of Bassi and Rebay [@BR1]
    * **BR2**: Second method of Bassi and Rebay [@BR2]
    
* ``LIBS_USE_MKL``:  ON/OFF
    
    This flag defines, whether Intel's MKL (Math Kernel Library) should be used. This is only meaningful when **FLEXI** is compiled with Intel compiler.    
    
* ``LIBS_USE_MPI``: ON/OFF

    This flag defines, whether **FLEXI** is compiled with MPI (necessary for parallel execution).

* ``FLEXI_NODETYPE``:

    **FLEXI** space discretization is based on a DG method. Here two different basis functions could be used, see [@KoprivaGassner2010] for details.
    
    * GAUSS
    * GAUSS-LOBATTO

* ``LIBS_USE_PAPI``:  ON/OFF
    
    Enable to use the PAPI library to perform performance measurements (e.g. flop counts). 

* ``FLEXI_PARABOLIC``:  ON/OFF
    
    This variable defines, whether the parabolic part of the chosen system should be included or not. Practically, this corresponds to setting the viscosity or diffusivity coefficient to zero, but is much more efficient, since the lifting routines are not necessary.
    
* ``FLEXI_POLYNOMIAL_DEGREE``:

    Since **FLEXI** is a high order CFD solver based on polynomial basis function, the polynomial degree can already be chosen in the compile process.
    If the default value N is chosen, different polynomial degrees can be defined later in the parameter file.
    
    * N
    * \#
    
* ``FLEXI_SPLIT_DG``:  ON/OFF

    Enable to use the split form of the discontinuous Galerkin operator. Allows to use kinetic energy or entropy stable flux functions. Only available for Gauss-Lobatto nodes.
    
* ``FLEXI_TESTCASE``:

    Some specific (and often used) simulation setups are encapsulated in test cases. These include e.g. case-specific initialization, analyze routines, boundary conditions etc. The default test case does not include any additions.
    * default
    * taylorgreenvortex
    * phill
    * channel
    * riemann2d
    
    See section \ref{sec:testcases} for details.

* ``FLEXI_VISCOSITY``:

    There are different modeling approaches for the viscosity in **FLEXI**. You can choose
    * constant
    * sutherland
    * powerlaw

The remaining part of the options deal with the post-processing framework **POSTI**. 
    
* ``POSTI``:  ON/OFF

    Enable to also build the post-processing tools next to the actual simulation software. When this general option is enabled, you will need to enable the specific options for the tools that should be build.
    
    
* ``POSTI_*``:  ON/OFF

    Each of the **POSTI** tools has it's own build option. Enable to build this specific tool.
    
* ``POSTI_VISU_PARAVIEW*``:  ON/OFF

    Enable to build the ParaView plugin for visualization of **FLEXI** simulation data. The ParaView libraries etc. must be available on the system and the environment variable $ParaView_DIR set accordingly.

