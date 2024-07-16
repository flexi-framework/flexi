[![logo](https://numericsresearchgroup.org/images/icons/flexi.svg "FLEXI")][flexi]


[![license](https://img.shields.io/github/license/flexi-framework/flexi.svg?maxAge=2592000 "GPL-3.0 License")](LICENSE.md)
[![doi](https://img.shields.io/badge/DOI-10.1016/j.camwa.2020.05.004-blue "DOI")](https://doi.org/10.1016/j.camwa.2020.05.004)
[![youtube](https://img.shields.io/badge/YouTube-red?logo=youtube "YouTube")](https://www.youtube.com/@nrgiag8633)
[![userguide](https://img.shields.io/badge/Userguide-silver "Userguide")][userguide]
[![gallery](https://img.shields.io/badge/Gallery-teal "Gallery")][gallery]

# About

[FLEXI][flexi] is a high-order numerical framework for solving PDEs, with a special focus on Computational Fluid Dynamics.
[FLEXI][flexi] is based on the Discontinuous Galerkin Spectral Element Method (DGSEM), which allows for high-order of accuracy 
and fully unstructured hexahedral meshes. The solver is parallelized very efficiently for large-scale applications and
scales to 500,000+ cores. Moreover, [FLEXI][flexi] comes with a capable pre- and postprocessing suite that enables complex
simulation setups up to the finished visualization.

[FLEXI][flexi] has been developed by the [Numerics Research Group (NRG)][nrg] founded by Prof. Claus-Dieter Munz and currently
lead by Prof. Andrea Beck at the Institute of Aerodynamics and Gasdynamics at the University of Stuttgart, Germany.

You can find detailed installation instructions, the extensive documentation and
several tutorial cases for FLEXI [here][flexi].

FLEXI is Copyright (C) 2016, Prof. Claus-Dieter Munz and is released under the **GNU General Public License v3.0**.
For the full license terms see the included [license file](LICENSE.md).

Numerous people have worked on and with FLEXI over the last years.
We would like to thank all these [contributors](CONTRIBUTORS.md) for their efforts they spent on building FLEXI.
 
In case you have questions regarding FLEXI or want to contribute yourself
by either reporting bugs, requesting features or adding somthing
different to the project, feel free to open an issue or pull request.

# Cite
FLEXI is a scientific project. If you use FLEXI for publications or
presentations in science, please support the project by citing it.
As general reference, please cite
```
Krais, N., Beck, A., Bolemann, T., Frank, H., Flad, D., Gassner, G., Hindenlang, F., Hoffmann, M., Kuhn, T., Sonntag, M., & Munz, C.-D. (2021).
FLEXI: A high order discontinuous Galerkin framework for hyperbolic–parabolic conservation laws,
Computers & Mathematics with Applications, 81, 186-219.
```
or use the following Bibtex entry

    @article{flexi,
      title = {{FLEXI}: {A} high order discontinuous {G}alerkin framework for hyperbolic-parabolic conservation laws},
      journal = {Computers \& Mathematics with Applications},
      volume = {81},
      pages = {186-219},
      year = {2021},
      doi = {https://doi.org/10.1016/j.camwa.2020.05.004},
      author = {Nico Krais and Andrea Beck and Thomas Bolemann and Hannes Frank and David Flad and Gregor Gassner and Florian Hindenlang and Malte Hoffmann and Thomas Kuhn and Matthias Sonntag and Claus-Dieter Munz},
    }

To refer to specific applications and features, you can also cite the appropriate paper from this [list][publications].

# Quick Start Guide
For a more detailed installation instructions, please see the documention [here][userguide].

FLEXI is tested for various Linux distributions including Ubuntu, OpenSUSE, CentOS or Arch and also runs on MacOS.
For installation you require the following dependencies:

| Package          | Required | Installed by FLEXI |
|:-----------------|:--------:|:------------------:|
| Git              |      x   |                    |
| CMake            |      x   |                    |
| C/C++ Compiler   |      x   |                    |
| Fortran Compiler |      x   |                    |
| LAPACK           |      x   |      x             |
| HDF5             |      x   |      x             |
| MPI              |     (x)  |                    |

The MPI library is only required for running parallel simulations on multiple ranks and the HDF5 and LAPACK libraries
can be installed automatically during the FLEXI build process.
The names of the packages and the package manager might differ depending on the specific distribution used.

### Getting the code
Open a terminal, download FLEXI via git and optionally export the FLEXI directory:

    git clone https://github.com/flexi-framework/flexi.git
    export FLEXI_DIR="$(pwd)/flexi"

### Compiling the code
Enter the FLEXI directory, create a build directory and use CMake to configure and compile the code

    cd $FLEXI_DIR
    mkdir build; cd build
    cmake ../
    make

The executable `flexi` is now contained in the FLEXI directory in `build/bin/`.
Custom configurations of the compiler options, dependencies and code features can be set using

    ccmake ../

### Running the code
Navigate to the directory of the tutorial **cavity** and run FLEXI

    cd $FLEXI_DIR/tutorials/cavity
    $FLEXI_DIR/build/bin/flexi parameter_flexi.ini

# Used libraries
FLEXI uses several external libraries as well as auxiliary functions from open source projects, including:
* [HDF5](https://www.hdfgroup.org/)
* [MPI](https://www.mcs.anl.gov/research/projects/mpi/)
* [LAPACK](https://www.netlib.org/lapack/)
* [OpenMP](https://www.openmp.org/)
* [FFTW](https://www.fftw.org/)
* [CMake](https://cmake.org/)
* [Reggie2.0](https://github.com/reggie-framework/reggie2.0/)
* [PAPI](https://icl.cs.utk.edu/papi/)

[nrg]:           https://numericsresearchgroup.org/index.html
[flexi]:         https://numericsresearchgroup.org/flexi_index.html
[publications]:  https://numericsresearchgroup.org/publications.html#services
[userguide]:     https://numericsresearchgroup.org/userguide/userguide.pdf
[gallery]:       https://numericsresearchgroup.org/gallery.html#portfolio
[youtube]:       https://www.youtube.com/@nrgiag8633 
