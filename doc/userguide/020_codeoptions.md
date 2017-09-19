\hypertarget{codeoptions}{}

# Code options \label{chap:codeoptions}
## Compiler options \label{sec:compileroptions}
This section describes the main configuration options which can be set when building FLEXI. The options should 
work on all supported systems, see section installation. The first set of options specify the components of 
the FLEXI toolkit to compile. Some options are dependent on others being enabled, so the available options may change. 
Components of the FLEXI package can be selected using the following options: \label{missing:compileroptions_not_finalized}

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

* ``CTAGS_PATH``:

    This variable specifies the CTAGS install directory.

* ``FLEXI_2D``: ON/OFF

    If set to ON the code will run in two-dimensional mode. You have to provide a mesh that consists of only one layer of elements in the third dimension.

* ``FLEXI_BUILD_HDF5``: ON/OFF

    This will be set to ON if no rebuilt HDF5 installation was found on your machine. In this case a HDF5 version will be build and used instead.

* ``FLEXI_EQNSYSNAME``:
    
    This variable defines the equation system, which will be compiled and used for the simulation. Possible values are
    * *Navierstokes*
    * *Linearscalaradvection*
    
* ``FLEXI_FV``:  ON/OFF
    
    Set this to ON to enable the usage of the finite volume subcell shock capturing mechanism.
    
* ``FLEXI_FV_RECONSTRUCTION``:  ON/OFF
    
    Only available if FLEXI_FV is set to ON. Enables the reconstruction of interface values in the finite volume subcells. Needed for calculation of gradients and to use a second order finite volume scheme.
    
* ``FLEXI_LIFTING``:

    Two different lifting methods for the parabolic part of the equation system are implemented. Possible values are
    * **BR1**: First method of Bassi and Rebay [@BR1]
    * **BR2**: Second method of Bassi and Rebay [@BR2]
    
* ``FLEXI_MKL``:  ON/OFF
    
    This flag defines, whether Intel's MKL (Math Kernel Library) should be used. This is only meaningful when FLEXI is compiled with Intel compiler.    
    
* ``FLEXI_MPI``: ON/OFF

    This flag defines, whether FLEXI is compiled with MPI (necessary for parallel execution).

* ``FLEXI_NODETYPE``:

    FLEXI space discretization is based on a DG method. Here two different basis functions could be used, see [@KoprivaGassner2010] for details.
    
    * GAUSS
    * GAUSS-LOBATTO

* ``FLEXI_PAPI``:  ON/OFF
    
    Enable to use the PAPI library to perform performance measurements (e.g. flop counts). 

* ``FLEXI_PARABOLIC``:  ON/OFF
    
    This variable defines, whether the parabolic part of the chosen system should be included or not. Practically, this corresponds to setting the viscosity or diffusivity coefficient to zero, but is much more efficient, since the lifting routines are not necessary.
    
* ``FLEXI_POLYNOMIAL_DEGREE``:

    Since FLEXI is a high order CFD solver based on polynomial basis function, the polynomial degree can already be chosen in the compile process.
    If the default value N is chosen, different polynomial degrees can be defined later in the parameter file.
    
    * N
    * \#
    
* ``FLEXI_TESTCASE``:

    Each different CFD simulation needs specific settings, e.g. initialization, source terms, analyze routines, etc. ... . If a new test case is defined in the test case directory, you can find this test case by toggling through this variable. There are some predefined test cases and a default test case, which could be used for code development
    * default
    * taylorgreenvortex
    * phill
    * channel
    
    See section \ref{sec:testcases} for details.
   
* ``FLEXI_TUTORIALS``

    If this variable is set ON, a tutorial folder with different test case will be generated.

* ``FLEXI_USE_PRIMITIVE_GRADIENTS``:  ON/OFF

    FLEXI has two strategies to calculate gradients for the parabolic operator. If FLEXI_USE_PRIMITIVE_GRADIENTS is turned ON, the gradient are calculated based on the primitive variable set, e.g. $(\rho,u,v,w,p)^T$, otherwise they are calculated by the conservative variables $(\rho,\rho u,\rho v,\rho w,\rho E)^T$. Only applicable to the *Navier-Stokes* equations set.

* ``FLEXI_VISCOSITY``:

    There are different modeling approaches for the viscosity in FLEXI. You can choose
    * constant
    * sutherland
    * powerlaw
* ``HDF5_DIR``:

    If you want to use a rebuilt HDF5 library that has been build using the CMake system, this directory should contain the CMake configuration file for HDF5 (optional).

