[![logo](https://numericsresearchgroup.org/images/icons/flexi.svg "FLEXI")][flexi]


[![license](https://img.shields.io/github/license/flexi-framework/flexi.svg?maxAge=2592000 "GPL-3.0 License")](LICENSE.md)
[![doi](https://img.shields.io/badge/DOI-10.1016/j.camwa.2020.05.004-blue "DOI")](https://doi.org/10.1016/j.camwa.2020.05.004)
[![youtube](https://img.shields.io/badge/YouTube-red?logo=youtube "YouTube")](https://www.youtube.com/@nrgiag8633)
[![userguide](https://img.shields.io/badge/Userguide-silver "Userguide")][userguide]
[![readthedocs](https://img.shields.io/badge/ReadTheDocs-2980b9 "ReadTheDocs")][readthedocs]
[![gallery](https://img.shields.io/badge/Gallery-teal "Gallery")][gallery]

# About

[FLEXI][flexi] is a high-order numerical framework for solving PDEs, with a special focus on Computational Fluid Dynamics. [FLEXI][flexi] is based on the Discontinuous Galerkin Spectral Element Method (DGSEM), which allows for high-order of accuracy and fully unstructured hexahedral meshes. The solver is parallelized very efficiently for large-scale applications and scales to 500,000+ cores. Moreover, [FLEXI][flexi] comes with a capable pre- and post-processing suite that enables complex simulation setups up to the finished visualization.

[FLEXI][flexi] has been developed by the [Numerics Research Group (NRG)][nrg] founded by Prof. Claus-Dieter Munz and currently lead by Prof. Andrea Beck at the Institute of Aerodynamics and Gasdynamics at the University of Stuttgart, Germany.

You can find detailed installation instructions, the extensive documentation and several tutorial cases for FLEXI [here][flexi].

[FLEXI][flexi] is Copyright (C) 2010-2022 Prof. Claus-Dieter Munz, Copyright (C) 2022-2024 Prof. Andrea Beck, and is released under the **GNU General Public License v3.0**. For the full license terms see the included [license file](LICENSE.md).

Numerous people have worked on and with [FLEXI][flexi] over the last years. We would like to thank all these [contributors](CONTRIBUTORS.md) for their efforts they spent on building [FLEXI][flexi].
 
In case you have questions regarding [FLEXI][flexi] or want to contribute yourself by either reporting bugs, requesting features or adding somthing different to the project, feel free to open an issue or pull request.

# Cite
[FLEXI][flexi] is a scientific project. If you use FLEXI for publications or presentations in science, please support the project by citing it. As general reference, please cite
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
For a more detailed installation instructions, please see the [online documentation][readthedocs] or the [userguide][userguide].

[FLEXI][flexi] is tested for various Linux distributions including Ubuntu, OpenSUSE, CentOS, or Arch. ƎLexi also runs on macOS. For the installation, you require the following dependencies:

| Package          | Required | Installed by FLEXI |
|:-----------------|:--------:|:------------------:|
| Git              |      x   |                    |
| CMake            |      x   |                    |
| C/C++ Compiler   |      x   |                    |
| Fortran Compiler |      x   |                    |
| LAPACK           |      x   |      x             |
| HDF5             |      x   |      x             |
| MPI              |     (x)  |                    |

The MPI library is only required for running parallel simulations on multiple ranks. The HDF5 and LAPACK libraries can are optionally built and locally installed during the [FLEXI][flexi] build process. The names of the packages and the package manager might differ depending on the specific distribution used.

### Getting the code
Open a terminal, download [FLEXI][flexi] via git

    git clone https://github.com/flexi-framework/flexi.git

### Compiling the code
Enter the [FLEXI][flexi] directory, create a build directory and use CMake to configure and compile the code

    cd flexi
    cmake -B build
    cmake --build build

The executable `flexi` is now contained in the [FLEXI][flexi] directory in `build/bin/`. Custom configurations of the compiler options, dependencies, and code features can be set using

    ccmake -B build

### Running the code
Navigate to the directory of the tutorial **cavity** and run [FLEXI][flexi]

    cd tutorials/cavity
    flexi parameter_flexi.ini

# Used libraries
[FLEXI][flexi] uses several external libraries as well as auxiliary functions from open source projects, including:
* [CMake](https://cmake.org)
* [FFTW](https://www.fftw.org)
* [HDF5](https://www.hdfgroup.org)
* [LAPACK](https://www.netlib.org/lapack)
* [MPI](https://www.mcs.anl.gov/research/projects/mpi)
* [OpenMP](https://www.openmp.org)
* [OpenBLAS](https://www.openblas.net)
* [PAPI](https://icl.cs.utk.edu/papi)
* [Reggie2.0](https://github.com/reggie-framework/reggie2.0)

[nrg]:           https://numericsresearchgroup.org/index.html
[flexi]:         https://numericsresearchgroup.org/flexi_index.html
[publications]:  https://numericsresearchgroup.org/publications.html#services
[userguide]:     https://numericsresearchgroup.org/userguide/pdf/userguide.pdf
[readthedocs]:   https://numericsresearchgroup.org/userguide/html/index.html
[gallery]:       https://numericsresearchgroup.org/gallery.html#portfolio
[youtube]:       https://www.youtube.com/@nrgiag8633 
