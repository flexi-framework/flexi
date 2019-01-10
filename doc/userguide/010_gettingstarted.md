\hypertarget{gettingstarted}{}

# Getting started \label{chap:gettingstarted}

## Installation

### Prerequisites
**FLEXI** has been tested for various Linux distributions. This includes Ubuntu 14.04 LTS, 16.04 LTS and 18.04 LTS, OpenSUSE 42.1 and CentOS 7.
The suggested packages in this section can of course be replaced by self compiled versions.

The required packages for the Ubuntu Linux distributions are listed in table \ref{tab:installation_prereqs_ubuntu}. Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git

| Package          | Ubuntu 14.04    | Ubuntu 16.04    | Ubuntu 18.04    |
|:----------------:|:---------------:|:---------------:|:---------------:|
| git              | x               |      x          |      x          |
| cmake            | x               |      x          |      x          |
| cmake-curses-gui | o               |      o          |      o          |
| liblapack3       | x               |      x          |      x          |
| liblapack-dev    | x               |      x          |      x          |
| gfortran         | x               |      x          |      x          |
| g++              | x               |      x          |      x          |
| mpi-default-dev  | x               |      x          |      x          |
| zlib1g-dev       | -               |      x          |      x          |
| exuberant-ctags  | o               |      o          |      o          |

Table: Debian/Ubuntu packages.\label{tab:installation_prereqs_ubuntu}
x: required, o: optional, -: not available

The required packages for OpenSUSE and CentOS are listed in table \ref{tab:installation_prereqs_redhat}.

Under OpenSUSE, packages are installed by the following command.

    sudo zypper install git

The `PATH` variable must be extended by the openmpi path

    export PATH=$PATH:/usr/lib64/mpi/gcc/openmpi/bin

Under CentOS, packages are installed by the following command.

    sudo yum install git

Additionally, the `PATH` variable must be extended by the openmpi path

    export PATH=$PATH:/usr/lib64/openmpi/bin



| Package          | OpenSUSE 42.1 | CentOS 7 |
|:----------------:|:-------------:|:--------:|
| git              |      x        |    x     |
| cmake            |      x        |    x     |
| lapack-devel     |      x        |    x     |
| openmpi          |      x        |    x     |
| openmpi-devel    |      x        |    x     |
| zlib-devel       |      x        |    x     |
| gcc-fortran      |      x        |    x     |
| gcc              |      x        |    -     |
| gcc-c++          |      x        |    x     |
| ctags-etags      |      -        |    o     |

Table: OpenSUSE/CentOS packages.\label{tab:installation_prereqs_redhat}
x: required, o: optional, -: not available


On some systems it may be necessary to increase the size of the stack (part of the memory used to store information about active subroutines) in order to execute **FLEXI** correctly. This is done using the command

~~~~~~~
ulimit -s unlimited
~~~~~~~

from the command line. For convenience, you can add this line to your `.bashrc`.

### Obtaining the source

The **FLEXI** repository is available at GitHub. To obtain the most recent version you have two possibilities:

* Clone the **FLEXI** repository from Github

        git clone https://github.com/flexi-framework/flexi.git

* Download **FLEXI** from Github:

        wget https://github.com/flexi-framework/flexi/archive/master.tar.gz
        tar xzf master.tar.gz

Note that cloning **FLEXI** from GitHub may not be possible on some machines, as e.g. the HLRS at the University of Stuttgart restricts internet access. Please refer to section \ref{sec:cloninghlrs} of this user guide.

### Compiling the code \label{sec:compilingthecode}

* Open a terminal
* Change into the **FLEXI** directory
* Create a new subdirectory and use CMake to configure and compile the code

        mkdir build; cd build
        cmake ../
        make

The executables **flexi** and **posti_visu** are contained in your **FLEXI** directory in `build/bin/`.

Custom configuration of compiler options may be done using

    ccmake ../

For a list of all compiler options see Section \ref{sec:compileroptions}.

#### Directory paths

In the following, we write `$FLEXIROOT` as a substitute for the path to the **FLEXI** repository. Please replace `$FLEXIROOT` in all following commands with the path to your **FLEXI** repository *or* add an environment variable `$FLEXIROOT`. 

Furthermore, the path to executables is omitted in the following, so for example, we write `flexi` instead of `$FLEXIROOT/build/bin/flexi`. 

Here is some explanation for Linux beginners:

In order to execute a file, you have to enter the full path to it in the terminal. There are two different ways to enable typing `flexi` instead of the whole path (do not use both at the same time!)

1. You can add an alias for the path to your executable. Add a command of the form

~~~~~~~
alias flexi='$FLEXIROOT/build/bin/flexi'
~~~~~~~

to the bottom of the file `~/.bashrc`. Source your `~/.bashrc` afterwards with

~~~~~~~
. ~/.bashrc
~~~~~~~

2. You can add the **FLEXI** binary directory to your `$PATH` environment variable by adding

~~~~~~~
export PATH=$PATH:$FLEXIROOT/build/bin
~~~~~~~

to the bottom of the file `~/.bashrc` and sourcing your `~/.bashrc` afterwards.


### Running the code

For a first minimal **FLEXI** run, do the following:

* Open a terminal
* Navigate to a directory, in this case *temp* 

        cd temp

* Copy the *cavity* tutorial folder

        cp -r $FLEXIROOT/tutorials/cavity/Basic_Re100 .
        cd Basic_Re100

* Run flexi

        flexi parameter_flexi.ini

* Convert the output files to the vtu format

        posti_visu cavity_State_0000000.200000000.h5

* Visualize using e.g. ParaView.

## Basic Usage

For a basic overview of the framework and the single components of the flow solver a flowchart is given in Figure \ref{fig:modules_flowchart}.

![Flowchart: Basic modules and files used by **FLEXI**\label{fig:modules_flowchart}](figures/flowchart.pdf)

### HOPR {-}

A standalone high-order preprocessor **HOPR** has been developed to generate high-order meshes from input data from external linear mesh generators. Different file formats are supported. HOPR has been recently made open source under the GPLv3 license. It generates a **FLEXI** conform mesh format in HDF5 for efficient parallel initialization. For a complete overview of HOPR, see [https://www.hopr-project.org](https://www.hopr-project.org).

**HOPR** can be compiled in the same way as **FLEXI** (see section \ref{sec:compilingthecode}). The basic command to run HOPR is

~~~~~~~
hopr parameter.ini
~~~~~~~

Note that the path to the **HOPR** executable is omitted in the command (see \ref{sec:compilingthecode}).

### FLEXI {-}

**FLEXI**, a high order DGSEM based CFD solver, is the core module in the tool chain. Generally **FLEXI** requires two main files as input, a mesh file in HDF5 format generated by **HOPR** and a parameter file where the main settings for the CFD simulation are set. The results files generated by **FLEXI** are also HDF5 files.

The basic command to run **FLEXI** is

~~~~~~~
mpirun -np [no. processors] flexi parameter.ini
~~~~~~~

Note: Adding ```mpirun -np [no. processors]``` before the **FLEXI** executable starts **FLEXI** in parallel with the specified number of threads. If it is omitted, **FLEXI** is run on one processor without MPI. This also applies to all other tools mentioned below.


### Parameter file {-}

The `parameter.ini` file contains the main settings for the CFD simulation. The parameter file defines e.g.

* CFL (Courant-Friedrichs-Lewy) number
* polynomial degree,
* simulation end time and dump/analyze intervals
* boundary conditions
* initial and boundary states

A complete list of all runtime options that can be set in the parameter file can be obtained with the command

~~~~~~~
flexi --help
~~~~~~~

It is also supplied in Section \ref{sec:parameterfile}. The ```--help``` option also works for most other **FLEXI** tools.

### posti_visu tool {-}

To visualize the results e.g. with ParaView, a converter tool is provided. The posti_visu tool takes the HDF5 files generated  by **FLEXI**.

The basic command to run the posti_visu tool is

~~~~~~~
mpirun -np [no. processors] posti_visu parameter.ini [flexi_outputfile.h5]
~~~~~~~

In this case a parameter file is specified in which options like the type and amount of the visualization nodes and mesh options are defined - see Section \ref{sec:postiVisu} for all available options. You can also omit the parameter file argument:

~~~~~~~
mpirun -np [no. processors] posti_visu [flexi_outputfile.h5]
~~~~~~~

This runs posti_visu using only standard options, i.e.

* equidistant visualization nodes
* amount of visualizaion nodes equals number of collocation points per element
* allowing for curved meshes
* visualizing the conservative variables

### HDF5 {-}

HDF5 is a data model, library, and file format for storing and managing data. It supports an unlimited variety of datatypes, and is designed for flexible and efficient I/O and for high volume and complex data. For further information and to download the software, visit the HDF5 website at [https://www.hdfgroup.org](https://www.hdfgroup.org).


## Feature list

The currently implemented features of **FLEXI** include

* Equation systems:
    * compressible Euler equations
    * compressible Navier-Stokes equations
    * linear scalar advection and diffusion
* Space discretization: DGSEM method [@KoprivaGassner2010;@hindenlang2012explicit]
    * Legendre Gauss
    * Legendre Gauss Lobatto
* Time discretization: explicit Runge-Kutta methods
    * standard RK methods
    * low storage RK methods [@Carpenter1994]
    * strong stability preserving RK methods [@niegemann2012efficient]
* Two- or three-dimensional domains
* Riemann solvers:
    * local Lax-Friedrichs
    * HLL
    * HLLC
    * Roe-Pike
* Curved Meshes
* Nonconforming Meshes via mortar interfaces [@koprivamortar2002]
* Shock capturing
    * Employing finite volume subcells
    * Several shock indicators available
* Boundary conditions
    * Various subsonic inflow and outflow conditions [@carlson2011inflow]
    * exact boundaries (Dirichlet)
    * periodic boundaries
    * slip wall (Euler wall)
    * non-slip walls (Navier-Stokes wall)
        * adiabatic
        * isothermal
* Splitform discontinuous Galerkin schemes [@Gassner2016]
* Dealiasing [@gassner2013accuracy]
    * filtering
    * overintegration
* Lifting methods
    * Bassi Rebay 1 [@BR1]
    * Bassi Rebay 2 [@BR1]
* Sponge zone [@flad2014discontinuous]
* Time averaging
