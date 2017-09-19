\hypertarget{gettingstarted}{}

# Getting started \label{chap:gettingstarted}

## Installation

### Prerequisites
**FLEXI** has been tested for various Linux distributions. This includes Ubuntu 14.04 LTS and 16.04 LTS, OpenSUSE 42.1 and CentOS 7. \label{missing:prerequisites_not_finalized}


The required packages for the Ubuntu Linux distributions are listed in table \ref{tab:installation_prereqs_ubuntu}. Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git
    

| Package          | Ubuntu 14.04    | Ubuntu 16.04    |
|:----------------:|:---------------:|:---------------:|
| git              | x               |      x          |
| cmake            | x               |      x          |
| cmake-curses-gui | x               |      x          |
| liblapack3       | x               |      x          |
| liblapack-dev    | x               |      x          |
| gfortran         | x               |      x          |
| g++              | x               |      x          |
|  mpi-default-dev | x               |      x          |
| zlib1g-dev       | -               |     x           |

Table: Required debian packages under Ubuntu.\label{tab:installation_prereqs_ubuntu}

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
| gcc-fortran      |       x       |    x     |
| gcc              |      x        |    -     |
| gcc-c++          |      x        |    x     |

Table: Required Red Hat packages under OpenSUSE and CentOS.\label{tab:installation_prereqs_redhat}

On some systems it may be necessary to increase the size of the stack (part of the memory used to store information about active subroutines) in order to execute FLEXI correctly. This is done using the command

        ulimit -s unlimited

from the command line. For convenience, you can add this line to your `.bashrc`.

### Obtaining the source

The **FLEXI** repository is available at GitHub. To obtain the most recent version you have two possibilities:

* Clone the FLEXI repository from Github

        git clone https://github.com/flexi-framework/flexi.git

* Download **FLEXI** from Github:

        wget https://github.com/flexi-framework/flexi/archive/master.tar.gz
        tar xzf master.tar.gz

Note that cloning FLEXI from GitHub may not be possible on some machines, as e.g. the machines at the HLRS at the University of Stuttgart due to restricted internet access. Please refer to section \ref{sec:cloninghlrs} of this user guide.

### Compiling the code \label{sec:compilingthecode}

* Open a terminal
* Change into the **FLEXI** directory
* Create a new subdirectory and use CMake to configure and compile the code

        mkdir build; cd build
        cmake ../
        make

The executables **flexi** and **flexi2vtk** are contained in your **FLEXI** directory in `build/bin/`. If desired, environment variables and aliases for the executables can be made available by adding the following line to your `.bashrc` \label{missing:doesthisworkwithalllinux}

    . ~/.flexi
and sourcing your `.bashrc` afterwards in the terminal

    . ~/.bashrc
    
Note that this enables the usage of `$FLEXI_TUTORIALS_DIR`.

Custom configuration of compiler options may be done using

    ccmake ../
    
For a list of all compiler options see Section \ref{sec:compileroptions}.


### Running the code

The following examples assume that the **FLEXI** environment variables have been made available (see Section \ref{sec:compilingthecode}).

* Open a terminal
* Navigate to a directory and copy a case folder 

        cd temp

* Copy the *cavity* tutorial folder \label{missing:aliases_run_code}

        cp -r $FLEXI_TUTORIALS_DIR/cavity .
        cd cavity

* Run flexi

        $FLEXI_DIR/flexi parameter_flexi.ini

* Converting the output files to the vtu format 

        $FLEXI_DIR/flexi2vtk parameter_flexi.ini cavity_State_0000000.200000000.h5

* Visualize using e.g. ParaView.

## Basic Usage

For a basic overview of the framework and the single components of the flow solver a flowchart is given in Figure \ref{mylabel}.

![Flowchart: Basic modules and files used by **FLEXI**\label{mylabel}](figures/flowchart.pdf)

### HOPR {-}

A standalone high-order preprocessor HOPR has been developed to generate high-order meshes from input data from external linear mesh generators. Different file formats are supported. HOPR has been recently made open source under the GPLv3 license. It generates a **FLEXI** conform mesh format in HDF5 for efficient parallel initialization. For a complete overview of HOPR, see [https://www.hopr-project.org](https://www.hopr-project.org).

The basic command to run HOPR is

    hopr parameter.ini


### FLEXI {-}

**FLEXI**, a high order DGSEM based CFD solver, is the core module in the tool chain. Generally **FLEXI** requires two main files as input, a mesh file in HDF5 format generated by HOPR and a parameter file where the main settings for the CFD simulation are set. The results files generated by **FLEXI** are also HDF5 files.

The basic command to run **FLEXI** is \label{missing:flexi_variable_notset2}

~~~~~~~
mpirun -np [no. processors] flexi parameter.ini
~~~~~~~



### Parameter file {-}

The `parameter.ini` file contains the main settings for the CFD simulation. The parameter file defines e.g.

* CFL number (Courant-Friedrichs-Lewy)
* polynomial degree,
* simulation end time and dump/analyze intervals
* boundary conditions
* Initial and boundary states

A complete list of all runtime options that can be set in the parameter file is supplied in Section \ref{sec:parameterfile}.

### FLEXI2VTK tool {-}

To visualize the results e.g. with ParaView, a converter tool is provided. The FLEXI2VTK tool takes the HDF5 files generated  by **FLEXI**.

The basic command to run the FLEXI2VTK tool is \label{missing:convert_variable_notset2}

~~~~~~~
mpirun -np [no. processors] $FLEXI_DIR/flexi2vtk parameter.ini [flexi_outputfile.h5]
~~~~~~~

In this case a parameter file is specified in which options like the type and amount of the visualization nodes and mesh options are defined - see Section \ref{sec:convert_tool}
for all available options. If you want to visualize some files using only standard options (equidistant visualization nodes and allowing for curved meshes), you can directly specify
the degree of the visualization basis on the command line and don't need to supply a parameter file. The syntax for this command line mode is 

~~~~~~~
mpirun -np [no. processors] $FLEXI_DIR/flexi2vtk --NVisu=INTEGER [flexi_outputfile.h5]
~~~~~~~

where INTEGER is substituted by the desired degree of the visualization basis. Usually the degree of your visualization basis is chosen higher than in the computation to avoid 
interpolation errors by supersampling the solution.

### HDF5 {-}

HDF5 is a data model, library, and file format for storing and managing data. It supports an unlimited variety of datatypes, and is designed for flexible and efficient I/O and for high volume and complex data. For further information and to download the software, visit the HDF5 website at [https://www.hdfgroup.org](https://www.hdfgroup.org).


## Feature list

The currently implemented features of **FLEXI** include \label{missing:featurelist_notfinal}

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
* Dealiasing [@gassner2013accuracy]
    * filtering
    * overintegration
* Lifting methods
    * Bassi Rebay 1 [@BR1]
    * Bassi Rebay 2 [@BR1]
* Sponge zone [@flad2014discontinuous]
* Time averaging
