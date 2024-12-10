# Code Overview

(sec:feature_list)=
## Feature List
**FLEXI** currently has the following features implemented, with most of them being covered in {numref}`Tutorials` {ref}`Tutorials`.

```{topic} Equation Systems
* Compressible Euler equations
* Compressible Navier-Stokes equations
* Linear scalar advection and diffusion equations
* Reynolds-averaged Navier--Stokes equations using Spalart--Allmaras turbulence model
```
```{topic} Space Discretization
* Discontinuous Galerkin Spectral Element Method (DGSEM) {cite}`KoprivaGassner2010,hindenlang2012explicit`
    * Legendre Gauss nodes
    * Legendre Gauss Lobatto nodes
* Finite Volume (FV) shock-capturing by either
    * Switching to finite volume subcells {cite}`sonntag2017efficient` or
    * Blending the finite volume operator {cite}`hennemann2021provably`
    * Several shock indicators available
```
```{topic} Time Discretization
* Explicit Runge-Kutta (RK) schemes
    * Standard RK schemes
    * Low storage RK schemes {cite}`Carpenter1994`
    * Strong stability preserving RK methods {cite}`niegemann2012efficient`
```
```{topic} Computational Domain
* Two- or three-dimensional domains
* Curved Meshes
* Nonconforming Meshes via mortar interfaces {cite}`koprivamortar2002`
* Sponge zone {cite}`flad2014discontinuous`
* Boundary conditions
    * Various subsonic inflow and outflow conditions {cite}`carlson2011inflow`
    * Exact boundaries (Dirichlet)
    * Periodic boundaries
    * Slip wall (Euler wall)
    * Non-slip walls (Navier-Stokes wall): adiabatic / isothermal
```
```{topic} Numerical Scheme
* Classical and split-form DG schemes {cite}`Gassner2016`
* Dealiasing {cite}`gassner2013accuracy`
    * Filtering
    * Overintegration
* Riemann solvers
    * Local Lax-Friedrichs
    * HLL
    * HLLC
    * Roe-Pike
* Lifting methods
    * Bassi Rebay 1 {cite}`BR1`
    * Bassi Rebay 2 {cite}`BR1`
* Time averaging
```

(sec:code_options)=
## Compiler Options
The following table describes the most important configuration options which can be set when building **FLEXI** using CMake.
Some options are dependent on others being enabled (or disabled), such that the available ones may change upon reconfiguring.

```{list-table} Compiler Options
:header-rows: 1
:align: center
:class: longtable
:width: 100%
:widths: 35 20 45

* - Build Option
  - Possible Values
  - Description
* - `CMAKE_BUILD_TYPE`
  - Release / Profile / Debug
  - normal execution / performance profiling using *gprof* / debug compiler for detailed error messages during code development
* - `CTAGS_PATH`
  -
  - install directory of *Ctags*, an optional program used to jump between tags in the source files (see e.g. the implementations [Exuberant Ctags](https://ctags.sourceforge.net/) or [Universal Ctags](https://github.com/universal-ctags/ctags))
* - `LIBS_BUILD_HDF5`
  - on / off
  - will be set to *on* if no pre-built HDF5 installation was found on your machine to build a HDF5 version during compilation
* - `HDF5_DIR`
  -
  - specify the directory of a pre-built HDF5 library that was built using the CMake system,<br/> this directory should contain the CMake configuration files (e.g. *hdf5-config.cmake*)
* - `FLEXI_2D`
  - on / off
  - set to *on* to run two-dimensional simulations, in this case you have to provide a mesh that consists of only one layer of elements in the third dimension
* - `FLEXI_EQNSYSNAME`
  - linearscalaradvection / navierstokes / rans_sa
  - linear scalar advection-diffusion equation / Navier--Stokes equations / Reynolds-averaged Navier--Stokes equations using Spalart--Allmaras turbulence model
* - `FLEXI_FV`
  - off / switch / blend
  - Finite-Volume (FV) shock-capturing: disabled / by *switching* DG elements in FV subcell representation {cite}`sonntag2017efficient` / by *blending* the FV and the DG operator {cite}`hennemann2021provably`
* - `FLEXI_FV_RECONSTRUCTION`
  - on / off
  - only available if `FLEXI_FV` is set either to *switch* or *blend*, enables the linear reconstruction of the solution at the FV subcell faces (second-order FV scheme)<br/> and is needed for the calculation of parabolic gradients
* - `FLEXI_LIFTING`
  - br1 / br2
  - lifting method to compute the DG gradients in the parabolic terms: first {cite}`BR1` / second {cite}`BR2` method of Bassi and Rebay
* - `LIBS_USE_MPI`
  - on / off
  - define whether to compile with MPI (necessary for parallel execution)
* - `FLEXI_NODETYPE`
  - gauss / gauss-lobatto
  - node-set used to define the basis functions of the DG method, see {cite}`KoprivaGassner2010` for details
* - `LIBS_USE_PAPI`
  - on / off
  - enable to use the [PAPI](https://icl.utk.edu/papi/) library to perform performance measurements (e.g. flop counts)
* - `FLEXI_PARABOLIC`
  - on / off
  - define whether the parabolic term of the chosen equation system should be included or not,<br/> more efficient than simply setting the diffusion coefficient to zero since the gradients do not need to be discretized
* - `FLEXI_POLYNOMIAL_DEGREE`
  - N / \{1,2,3,...\}
  - polynomial degree of basis functions: to be set in parameter file / compile with fixed degree for performance (1,2,3,...)
* - `FLEXI_SPLIT_DG`
  - on / off
  - enable the split form of the DG operator, allows to use kinetic energy or entropy stable flux functions
* - `FLEXI_TESTCASE`
  - default / hit / phill / riemann2d / taylorgreenvortex / channel
  - some benchmark simulation setups are encapsulated in test cases (separate sub-folders) with case-specific initialization, analyze routines, boundary conditions, etc.,<br/> while the default test case does not include any additions: see section  {ref}`Tutorials` for more details
* - `FLEXI_VISCOSITY`
  - constant / sutherland / powerlaw
  - modeling approach for the dynamic viscosity: constant / Sutherland's law / power law
* - `POSTI`
  - on / off
  - enable to also build the post-processing tool-set next to the actual simulation code, the specific tools can be selected once this flag is enabled
* - `POSTI_VISU_PARAVIEW`
  - on / off
  - enable to build the **ParaView** plugin for the visualization of **FLEXI** simulation data,<br/> the ParaView libraries must be available on the systems and environment variable `$ParaView_DIR` set accordingly
* - `FLEXI_PERFORMANCE`
  - on / off
  - enables a set of advanced features to improve the performance of FLEXI
* - `FLEXI_PERFORMANCE_OPTILIFT`
  - on / off
  - enable to lift only the gradients of the variables in the flux function of the selected equation system<br/> improves the performance for `FLEXI_PARABOLIC=ON`, but cannot be used if posti tool-set is built (`POSTI=ON`)
* - `FLEXI_PERFORMANCE_PGO`
  - on / off
  - enables profile-guided optimization (PGO) for compilation, currently only supported with GNU compiler<br/> the required two-step compilation process is detailed in section {ref}`sec:tut_ptcf_performance`
```
