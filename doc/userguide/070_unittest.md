\hypertarget{unittest}{}

# Unit tests \label{chap:unittest}

Unit tests are used to test individual key units of the source code. Currently these key routines include:

* Calculation of node positions and integration weights.
* Calculation of Vandermonde matrices.
* Calculation of derivative matrices.
* Algorithms to interpolate data from one set of nodes to another.
* Algorithm to prolong volume data to sides.
* Algorithm to perform a surface integral.
* Functionality of Read-In tools.

## Integration of unit test with CTest

These unit tests are integrated into the **FLEXI** build process using the [CTest](https://cmake.org/Wiki/CMake/Testing_With_CTest) tool. Usually CTest will be run every time you
build **FLEXI** and give you an overview on the exit status of each test that looks something like this:

~~~~
Run unit tests
Test project /home/FLEXI/build
    Start 1: NodesAndWeights
1/6 Test #1: NodesAndWeights ..............   Passed    0.01 sec
    Start 2: Vandermonde
2/6 Test #2: Vandermonde ..................   Passed    0.01 sec
    Start 3: DerivativeMatrix
3/6 Test #3: DerivativeMatrix .............   Passed    0.00 sec
    Start 4: ChangeBasis
4/6 Test #4: ChangeBasis ..................   Passed    0.00 sec
    Start 5: SurfInt
5/6 Test #5: SurfInt ......................   Passed    0.00 sec
    Start 6: ProlongToFace
6/6 Test #6: ProlongToFace ................   Passed    0.00 sec

100% tests passed, 0 tests failed out of 6

Total Test time (real) =   0.05 sec
~~~~

To manually run the tests after a build use the CTest command

~~~~
ctest
~~~~

in your build directory. The manual page of CTest can give you an overview of all available options.

If you don't want to run the test after each build there is a CMake option called FLEXI_UNITTESTS that can be used to turn the tests on and off.
This is an advanced option that CCMake will only show if you enter the advanced mode by pressing the t key.

## Implementation of unit tests

All unit tests are implemented in FORTRAN and can be found in the subdirectory unitTests in your **FLEXI** directory alongside a seperate CMakeLists.txt and some binary input and reference files.

### CMakeLists.txt

The CMakeLists.txt defines a custom function called add_unit_test which can be used in the CMakeLists.txt to add a single test to the CTest tool. The syntax is

~~~~
add_unit_test(NAME SOURCEFILE.F90)
~~~~

All tests are defined using this function. At the end of the CMakeLists.txt a custom target all_tests is defined which includes all unit tests and will run the ctest command after it has been build.

The whole CMakeLists.txt content is included in the main CMakeLists.txt if the option FLEXI_UNITTESTS is set to ON (default) by CMake.

### General unit test structure

The general structure of the unit tests is the same in all cases. They are implemented as FORTRAN programs. The unit test will call a function or subroutine from the **FLEXI** framework with input either set in the program itself or read from a binary file.
The output of this call will then be compared to some precomputed reference results (also stored as binary files) with a certain tolerance to account for differences in e.g. compiler versions and system architecture. If the results are within the given tolerance,
the test will be passed, otherwise it will fail by returning a value other than 0.

The programs usually also contain a command line option that can be uses to generate the reference solution from a code version that is known to work correctly.

Have a look at the source code of one of the already implemented unit tests if you want to have a more detailed idea about how to implement your own tests.

### Generation of reference mesh data

Some of the unit tests require parts of the mesh data structure to be able to call the funtions to be tested. For this purpose, a curved single element is created and all the mesh data stored as a binary file called ``UnittestElementData.bin``. This binary file can then be read during runtime
by the unit test programs.

To generate the curved single element mesh, run **HOPR** with the parameter file provided in the ``unitTest`` subdirectory of **FLEXI**. To generate the binary file, run **FLEXI** with the following command line argument and the parameter file
provided in the ``unitTest`` subdirectory:


~~~~
flexi --generateUnittestReferenceData parameter.ini
~~~~
