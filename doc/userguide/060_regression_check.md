\hypertarget{regressioncheck}{}

# Regression Check \label{chap:regressioncheck}

The regression check (or *Reggie*) is intended to provide continuous code revisal in order to locate bugs.
To achieve this, *Reggie* can automatically execute provided examples and check if they are running correctly in a multitude of ways:

| Type                            | Available Modes       | Description                                         |
| ------------------------------- |:---------------------:| ---------------------------------------------------:|
| No-result                       | -                     | only checks for successful execution                |
| L2 and Linf norms               | 1, 2                  | Error norms calculated via *ExactFunc*              |
| h5diff                          | 3                     | diff multiple HDF5 files                            |

Table: Regression check types.

[comment]: <> (| Record points                   | 1, 2, 3, 4            | time signal of properties recorded over time        |)
[comment]: <> (| Performance                     | 1                     | time measurement of wall time                       |)
[comment]: <> (| analyze-tools                   | 1, 2, 3, 4            | high-level analysis (entropy, FFT, time averaging)  |)
[comment]: <> (| line plot                       | 2, 3, 4               | plot properties over lines                          |)


Some checks can be performed with different modes, e.g., the L2 norm may be compared to pre-calculated values for each equation 
variable or a constant. The following table liste possible comparison types.

| Property                        | Value         | Description                                |
| ------------------------------- |:-------------:| ------------------------------------------:|
| Comparison mode                 | 1             | use constant (pre-defined value)           |
|                                 | 2             | use supplied data file                     |
|                                 | 3             | HDF5 supplied data file for H5-diff        |

Table: Regression check modes.

[comment]: <> (|                                 | 4             | use pre-defined function, e.g., $f(x)=x^2$ |)

Besides running and checking examples using my current build of **FLEXI**, *Reggie* is also able to automatically build **FLEXI** using different combinations of CMake options 
and run examples on each of them. 
This helps to check if changes made to the code don't have negative effects on different equation systems, lifting methods, viscosity types, ... than the 
one currently working on.

The regression check can and should be run by the user when making significant changes to the code to check if no unintendet side effects are introduced. 
See section \ref{sec:runReggi} to learn how to run the regression check manually.
*Reggie* can also be integrated into systems like [GitLab CI](https://about.gitlab.com/gitlab-ci/) to perform checks when somebody pushes to the
repository or more extensive tests on a e.g. weekly basis. This helps to maintain code quality. A example configuration for the GitLab CI ``.gitlab-ci.yml`` is provided
in the **FLEXI** main directory.


## Enable *Reggie* 

By default, *Reggie* is not build automatically with **FLEXI**. If you want to manually use the regressioncheck, run

~~~~~~~~~
make regressioncheck
~~~~~~~~~

in your build directory. If you want to build **FLEXI** and *Reggie* together, you can use


~~~~~~~~~
make all regressioncheck
~~~~~~~~~

## Files

The following files are to be placed within every regression check folder 

| File                            | Value         | Description                                |
| ------------------------------- |:-------------:| ------------------------------------------:|
| configuration.cmake             | required      | compilation flags: EQNSYS, MPI, etc.       |
| parameter_reggie.ini            | required      | Comparison type, general settings          |
| *\*\_mesh.h5*                   | required      | mesh file in HDF5 format                   |
| parameter_flexi.ini             | required      | IC, BC and numerical settings              |
| parameter_hopr.ini              | optional      | for ANSA, ICEM formats etc.                |

Table: Required files for regression check.

### configuration.cmake

The information supplied by this file is used to create cmake builds in every combination possible using the 
compile flags and excluding invalid combinations. The following example (/regressioncheck/examples/run_freestream)

~~~~~~~Bash
! fixed compiler flags
CMAKE_BUILD_TYPE=Release
FLEXI_BUILD_HDF5=OFF
FLEXI_PAPI=OFF
FLEXI_POLYNOMIAL_DEGREE=N
FLEXI_MKL=OFF

! include combinations
FLEXI_EQNSYSNAME=navierstokes,linearscalaradvection
FLEXI_LIFTING=br1,br2
FLEXI_MPI=ON  !,OFF ! MPI=OFF is not supported on pre-compiled runner HDF5
FLEXI_NODETYPE=GAUSS,GAUSS-LOBATTO
FLEXI_PARABOLIC=ON,OFF
FLEXI_VISCOSITY=constant,sutherland,powerlaw

! exclude combinations
EXCLUDE:FLEXI_VISCOSITY=sutherland,FLEXI_PARABOLIC=OFF
EXCLUDE:FLEXI_VISCOSITY=powerlaw,FLEXI_PARABOLIC=OFF
EXCLUDE:FLEXI_LIFTING=br2,FLEXI_EQNSYSNAME=linearscalaradvection
~~~~~~~

results in 48 possible and 24 valid compile combinations. The regression check would proceed and compile to build 24 combinations 
and run them within the considered example.

### parameter_reggie.ini

The parameter input for the regression check in located in this file. The L2/LInf reference data is stored in *referencenorm.txt* and 
the HDF5-diff files for comparison are *ReferenceStateFile* and *CheckedStateFile*. The array name for HDF5-diff is given by 
*ReferenceDataSetName*. If the example is to be restarted from a supplied restart file, the name of the file must be supplied under 
*RestartFileName*.

~~~~~~~Bash
nVar= 5
MPI= T
ReferenceFile= referencenorm.txt
ReferenceStateFile= cavity_reference_State_0000000.200000000.h5
CheckedStateFile= cavity_State_0000000.200000000.h5
ReferenceDataSetName= DG_Solution
RestartFileName=
~~~~~~~

Note: $nVar$ and $MPI$ in their current state are deprecated. This will be changed in a future release.

### parameter_flexi_/*.ini

A parameter file for each *FLEXI_EQNSYSNAME* must be supplied in order to run a regressioncheck with the equation system of choice. 
The regression check automatically detects if a cmake compiled flexi binary is able to be run within the example depending on the 
supplied parameter_flexi_\*.ini file.


## Running the regression check \label{sec:runReggi}

The binary *regressioncheck* (found in the */bin* directory) supplied the user with quick start help information via 

~~~~~~~Bash
./regressioncheck --help
~~~~~~~

The regression check operates three main loops and will terminate accordingly to the errors occuring within a certain level. The 
loops are

* **Building** (level 1): loops different flexi binaries
 
    multiple compiler flag combinations for cmake are supplied and each of them are built. A failiure occuring when 
    a combination is being built will terminate the regression check immediately.
* **Examples** (level 2): loops different folders with setups (mesh+\*.ini files)

    multiple examples may be supplied within the directory */regressioncheck/example/*. The regression check may 
    loop these examples with all builds specified in level 1. A failiure occuring on this level will continue the code and produce 
    an error code for the example.
* **Sub-examples** (level 3): loops different parameter_flexi_\*.ini files

    This level handles all run-time flags, e.g., TimeDisc or finite volume IndicatorType, which are supplied 
    in the parameter_flexi_\*.ini file. A failiure occuring on this level will continue the code and produce an error code for the 
    example.

As the regression check can be executed under multiple modes, multiple input parameter can be supplied. Generally, the regression 
check is performed by running

~~~~~~~Bash
./regressioncheck [RuntimeOption] [RuntimeOptionType]
~~~~~~~

First input argument [RuntimeOption]

| [RuntimeOption]           | mode                                                                      |
|---------------------------|---------------------------------------------------------------------------|
| help                      | prints the information output                                             |
|---------------------------|---------------------------------------------------------------------------|
| run (default)             | runs all examples with prefix *run*, e.g., *run_1* or                     |
|                           | *run_freestream*. These tests all inlcude very                            |
|                           | short simulations with <1sec execution time                               |
|                           | e.g. used for on-check-in tests                                           |
|---------------------------|---------------------------------------------------------------------------|
| build                     | runs the example *run_freestream* on default and                          |
|                           | (requires locally built HDF5 or loaded HDF5 paths)                        |
|                           | compiles all possible compiler flag combinations                          |
|                           | specified in *comfiguration.flexi* and considers                          |
|                           | the specified exclude list for invalid combinations                       |
|                           | e.g. used for nightly tests                                               |
|                           | The cmake output may be shown on-screen by adding *debug* after possible  | 
|                           | [RuntimeOptionType] commands                                              |
|                           | Multi-processor compilation is supported by adding *XX*  after possible   |
|                           | [RuntimeOptionType] commands, where *XX* is the number of processors      |
|---------------------------|---------------------------------------------------------------------------|
| conv_test (ToDo!)         | specific feature test: runs the *conv_test* example                       |
|                           | runs two modes: p-convergence and h-convergence                           |
|                           | e.g. used for weakly tests                                                |
|---------------------------|---------------------------------------------------------------------------|
| performance (ToDo!)       | specific feature test: runs the *performance* example                     |
|                           | automatically checks out specified flexi version tag                      |
|                           | and run the example to acquire the reference                              |
|                           | performance wall time                                                     |
|                           | e.g. used for weakly tests                                                |

Table: First input argument for the regression check.

Second input argument [RuntimeOptionType] depends on the first input argument

|[RuntimeOptionType]        | [RuntimeOption] | mode                                                          |
|---------------------------|-----------------|---------------------------------------------------------------|
|run (default)              |                 | perform all tests for examples beginning with *run_\**        |
|---------------------------|-----------------|---------------------------------------------------------------|
|cavity                     | run             | specific feature test: runs all the the examples that         |
|                           |                 | begin with *cavity*, e.g., *cavity_1* or *cavity_xyz*         |
|                           |                 | the example *cavity*, e.g., tests long time stability,        |
|                           |                 | some BC and ExactFunc ( e.g. used for weakly tests)           |
|                           |                 |                                                               |
|                           | build           | uses the *configurations.flexi* within the *cavity*           |
|                           |                 | directory and builds all compiler flag combinations           |
|                           |                 | that are specified there and runs them all on *cavity*        |


## Practical Examples

1. Print the quick guide and answers the following questions

   * How do I run the regression check?
   * Which input files are needed, e.g., comfiguration.flexi and parameter_reggie.ini
   * What do the Error Codes tell me?

~~~~~~~Bash
./regressioncheck --help
~~~~~~~


2. The default mode is *run* which executes all examples beginning with *run_..* under /regressioncheck/examples/run_..

~~~~~~~Bash
./regressioncheck
~~~~~~~

3. Execute (only run them) all examples beginning with *build_..* under /regressioncheck/examples/build_..

~~~~~~~Bash
./regressioncheck run build
~~~~~~~

4. Execute (only run them) all examples beginning with *feature_..* under /regressioncheck/examples/feature_..

~~~~~~~Bash
./regressioncheck run feature
~~~~~~~

5. Building with camke: 
  uses the configuration.flexi file under /regressioncheck/examples/run_freestream and all compiler flag 
  comibations that are specified there in order to build flexi with cmake using 2 threads/processors and 
  runs this example with each build separately

~~~~~~~Bash
./regressioncheck build 2
~~~~~~~

6. Building with camke with compilation output: 
  uses the configuration.flexi file under /regressioncheck/examples/run_freestream and all compiler flag 
  comibations that are specified there in order to build flexi with cmake and runs this example with each build separately 
  the argument *debug* pipes the complete cmake and flexi build output on screen instead of writing it to the file build_flexi.out

~~~~~~~Bash
./regressioncheck build debug
~~~~~~~

7. Building with camke and special directoreis: 
  uses the configuration.flexi file under /regressioncheck/examples/feature_cavity and all compiler flag 
  comibations that are specified there in order to build flexi with cmake and runs this example with each build separately

~~~~~~~Bash
./regresioncheck build feature_cavity
~~~~~~~



## Error Codes

The regression check prints a summary at the end of its execution with a list of examples that were conducted as well as error 
codes for debugging purposes. These codes are liste below and give an indication of the error source.

| Error Code                 | Description                                            |
| --------------------------:|--------------------------------------------------------|
|                          0 | no error                                               |
|                          1 | failed during build                                    |
|                          2 | computation of example failed                          |
|                          3 | mismatch in norms                                      |
|                          4 | mismatch in dataset                                    |
|                         77 | no flexi executable found for option run               |
|                         99 | fail of execute_system_command                         |

Table: Error codes produced by *./regressioncheck*.

## Creating New Reggies

* if a new *FLEXI_EQNSYSNAME* is to be introduced, the code has to be altered in *regressioncheck_run.f90* with the appropriate *Nvar*
* create a folder under */regressioncheck/examples/*
* supply 
    * the parameter files
    * a mesh file
    * referece states and/or norm files
  


