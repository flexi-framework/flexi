# Installation

## Prerequisites

Generally, **FLEXI** requires the following packages:
- git
- CMake
- Fortran and C/C++ compilers (GNU compilers recommended)
- MPI libraries (OpenMPI recommended)
- ParaView
- LAPACK/OpenBLAS[^1]
- HDF5[^1]
- FFTW[^1]

(sec:installing_dependencies)=
### Installing the Dependencies from the Package Repositories
**FLEXI** has been tested on various Linux distributions, including Ubuntu and Debian, as well as OpenSUSE, CentOS and Fedora.
The required packages for DEB-based and RPM-based Linux distributions can be obtained from the `apt` and the `dnf` environment, respectively. Refer to {numref}`tab-deb` for the package names.

```{list-table} Package names for Linux distributions
:name: tab-deb
:header-rows: 2
:align: center
:class: longtable
:width: 100%
:widths: 30 40 40

* - Package
  - Debian / Ubuntu
  - RHEL / Fedora
* - Installation Command
  - `sudo apt-get install`
  - `sudo dnf install`
* - git
  - `git`
  - `git`
* - CMake
  - `cmake-extras cmake-curses-gui`
  - `cmake`
* - C/C++/Fortran
  - `g++ gfortran`
  - `gcc-c++ gcc-gfortran`
* - MPI
  - `mpi-default-dev`
  - `mpich-devel`
* - ZLIB
  - `zlib1g-dev`
  - `zlib-ng-devel`
* - ParaView
  - `paraview-dev`
  - `paraview-devel`
* - LAPACK/OpenBLAS[^1]
  - `libopenblas-dev`
  - `openblas-devel`
* - HDF5[^1]
  - `libhdf5-mpi-dev`
  - `hdf5-mpich-devel`
* - FFTW[^1]
  - `libfftw3-dev`
  - `fftw-devel`
```

```{tip}
On RPM-based distributions, you might need to load the MPI module using the command `module load mpi`.
```


### Additional Configuration

On some systems it may be necessary to increase the size of the stack (part of the memory used to store information about active subroutines) in order to execute **FLEXI** correctly. This is done by entering the following command.
```bash
ulimit -s unlimited
```


(sec:obtain_flexi)=
## Obtaining the Source Code

The **FLEXI** repository is available at GitHub. To obtain the most recent version you have three possibilities:

* Clone the **FLEXI** repository from GitHub
  ```bash
  git clone https://github.com/flexi-framework/flexi.git
  ```
* Download **FLEXI** from GitHub:
  ```bash
  wget https://github.com/flexi-framework/flexi/archive/master.tar.gz
  tar xzf master.tar.gz
  ```
* Download a release **FLEXI** repository from GitHub
  
  <https://github.com/flexi-framework/flexi/tags>

```{attention}
Cloning **FLEXI** from GitHub may not be possible on some HPC clusters due to restricted internet access. Please refer to the cluster's user instructions for possible remedies, such as establishing a SOCKS proxy on a machine with unlimited internet access, as documented [here](https://kb.hlrs.de/platforms/index.php/Secure_Shell_ssh#HTTP(S)) for the HLRS at the University of Stuttgart.
```


(sec:compile_flexi)=
## Compiling the Code

In order to compile the code, change into the **FLEXI** root directory, create a new sub-folder, and use CMake to configure and compile the code
```bash
mkdir build
cmake -B build
cmake --build build
```
Custom configuration of compiler options may be done using
```bash
ccmake -B build
```
For a list of all compiler options see section {ref}`sec:code_options`.

The executables `flexi` and `posti_visu` (if enabled) are generated in the sub-directory `build/bin/`.

```{note}
In the remainder of this user guide, we omit the path to the **FLEXI** executable (and related tools), but assume assume it can be executed directly by typing `flexi`. This can be achieved by defining an *alias* or *symbolic link*, for example.
```
For Linux beginners, we provide a short explanation on how to achieve this usage behavior. In general, in order to execute a file, the command either has to be in the `PATH` environment variable or you have to enter the full path to it in the terminal. To enable typing only `flexi`, you can add a symbolic link to the **FLEXI** executable in the current directory, e.g., test case folder, by entering
```bash
ln -s [FLEXI_ROOT]/build/bin/flexi
```
Among the files in the current directory, this symbolic link will be listed as
```bash
flexi -> [FLEXI_ROOT]/build/bin/flexi
```


## Running the Code

For a first minimal **FLEXI** simulation, navigate to the *cavity* tutorial folder and run **FLEXI**:
```bash
cd [FLEXI_ROOT]/tutorials/cavity/Basic_Re100
flexi parameter_flexi.ini
```
Convert the output files to the *vtu* format by entering
```bash
posti_visu cavity_State_0000000.200000000.h5
```
and visualize the generated files using, e.g., **ParaView**. Note that this conversion step requires enabling the `posti_visu` tool by toggling the `POSTI` flag in the CMake configuration (see section {ref}`sec:compile_flexi` above).


[^1]: Package can be automatically installed through **FLEXI** as compiler option.
