# FLEXI

[![license](https://img.shields.io/github/license/flexi-framework/flexi.svg?maxAge=2592000)]()

FLEXI is a high-order numerical framework for solving PDEs,
with a focus on Computational Fluid Dynamics.
FLEXI based on the Discontinuous Galerkin Spectral Element
Method (DGSEM), which allows for high-order of accuracy 
and fully unstructured hexahedral meshes.
The solver is parallelized very efficiently and scales up
to hundreds of thousand cores.

FLEXI has been developed by the [Numerics Research Group (NRG)][nrg]
lead by Prof. Claus-Dieter Munz at the Institute of Aerodynamics
and Gasdynamics at the University of Stuttgart, Germany.

This is a scientific project. If you use FLEXI for publications or
presentations in science, please support the project by citing
our publications given in [references](REFERENCE.md).

## Installation / Documentation

For installation instruction see [install](INSTALL.md).

See the full documentation including usage instructions and
turorial for FLEXI [here][flexi].
 
In case you have question regarding FLEXI, want to report bugs
or contribute to the project you can use the mailing list
<flexi@listserv.uni-stuttgart.de>.
You can also subscribe to the mailing list [here][list].

## License
FLEXI is Copyright (C) 2016, Prof. Claus-Dieter Munz and is 
released under the terms of the
GNU General Public License v3.0. For the full license terms see
the included license file [license](LICENSE.md).

## List of Contributors
Numerous people have worked on and with FLEXI over the last years.
We would like to thank all these [contributors](CONTRIBUTORS.md)
for their efforts they spent on building FLEXI.

## Used libraries

FLEXI uses several external libraries as well as auxilliary functions from open source projects, including:
* [HDF5](https://www.hdfgroup.org/)
* [MPI](http://www.mcs.anl.gov/research/projects/mpi/)
* [LAPACK](http://www.netlib.org/lapack/)
* [PAPI](http://icl.cs.utk.edu/papi/)
* [FOUL](http://foul.sourceforge.net/)
* [OpenMP](http://www.openmp.org/)
* [FFTW](http://www.fftw.org/)

[nrg]:  https://www.iag.uni-stuttgart.de/arbeitsgruppen/numerische-methoden/
[flexi]: https://www.flexi-project.org/
[list]: https://listserv.uni-stuttgart.de/mailman/listinfo/flexi

## Regressioncheck

For informations about the regressioncheck see [reggie](REGGIE.md).
