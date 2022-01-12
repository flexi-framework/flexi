## Flexi installation procedure

## Prerequisites

Flexi has been tested for various Linux distributions. This includes Ubuntu 14.04 LTS and 16.04 LTS, OpenSUSE 13.2 and CentOS 7. In addition Paraview or Tecplot can be used for visualization.

The required packages for the Ubuntu Linux distributions are listed in the tables below. Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git


| Package          | Ubuntu 14.04    | Ubuntu 16.04    |
|:----------------:|:---------------:|:---------------:|
| git              | x               |      x          |
| cmake            | x               |      x          |
| cmake-curses-gui | x               |      x          |
| liblapack3       | x               |      x          |
| liplapack-dev    | x               |      x          |
| gfortran         | x               |      x          |
| g++              | x               |      x          |
|  mpi-default-dev | x               |      x          |
| zlib1g-dev       | -               |     x           |

Table: Required debian packages under Ubuntu.


Under OpenSUSE, packages are installed by the following comand.

    sudo zypper install git

The `PATH` variable must be extended by the openmpi path

    export PATH=$PATH:/usr/lib64/mpi/gcc/openmpi/bin

Under CentOS, packages are installed by the following comand.

    sudo yum install git

Additionally, the `PATH` variable must be extended by the openmpi path

    export PATH=$PATH:/usr/lib64/openmpi/bin

| Package          | OpenSUSE 13.2 | CentOS 7 |
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

Table: Required Red Hat packages under OpenSUSE and CentOS.

## Compiling the code

* Open a terminal
* Change into the Flexi directory
* Create a new subdirectory and use CMake to configure and compile the code

        mkdir build; cd build
        cmake ../
        make

The executables **flexi** and **convert** are contained in your Flexi directory in `build/bin/`.

Custom configuration of compiler options may be done using

    ccmake ../

## Running the code

* Open a terminal
* Navigate to a directory and copy a case folder

        cd temp

* Copy the *cavity* tutorial folder

        cp -r $FLEXI_TUTORIALS/cavity .
        cd cavity

* Run flexi

        $flexi parameter_flexi.ini
