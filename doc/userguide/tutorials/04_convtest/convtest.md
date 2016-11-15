## Convergence Test

In this tutorial the order of convergence for **FLEXI** is computed. The procedure is fully scripted, such that in the end a number of runs have been performed on a variation of grids or a variation of polynomial degrees and the order of convergence is computed automatically. A plot of the corresponding L2 error norms is produced and copied into the directory where the convergence test is executed from. Note, the script is written in *Python 2.7*. 

The convergence test is separated into two parts, first an inviscid test where the order of convergence for the advective part, i.e. the Euler equation without viscous fluxes, is computed. Then a viscous convergence test calculates the order of convergence with consideration of the viscosity.

In this tutorial, the following topics are dealt with:

* Calling **FLEXI** using a Python script
* Compute the order of convergence
* Switch from a Navier-Stokes simulation to an Euler simulation
* Treatment of source terms in **FLEXI**
* Execution of **FLEXI** with Mortar meshes

### Manufactured Solution

The basic principle of computing the order of convergence of **FLEXI** is to apply a benchmark test where an exact analytical solution is known. We apply here the method of manufactured solutions. A detailed description is found in Roache [@roache_code_2001]. In general, a smooth function is proposed as solution to the equation system. Since the proposed function is generally not a solution to the equation system, a source term is calculated to force the corresponding solution. Therefore, the function is derived analytically and inserted into the equation system, cf. for the continuity equation

\begin{equation}
\rho_t + (\rho u)_x = Q(x,t) \quad \text{with} \quad \rho = A + \sin(B(x,t)) \quad \text{and} \quad u = \text{const}
.
\end{equation}

The source term has to be added accordingly within the time integration loop of the flow solver. In **FLEXI**, sine waves in the density are advected on a constant velocity field. For the different equation systems, i.e. Euler or Navier-Stokes, different source terms have to be considered.

### Inviscid Convergence Test

Copy the *convtest* tutorial folder \label{missing:aliases_tutorial_convtest}

        cp -r $FLEXI_TUTORIALS/convtest .

#### Compiler Options

Since this test case is valid for the Euler equations only, the code has to be compiled without parabolic fluxes. Hence, make sure
to compile **FLEXI** with the cmake options in the following table


| Option                          | Value         | Comment      |
| ------------------------------- |:-------------:| ------------:|
| CMAKE_BUILD_TYPE                | Release       |              |
| FLEXI_EQYNSYSNAME               | navierstokes  |              |
| FLEXI_PARABOLIC                 | OFF           |              |
| FLEXI_MPI                       | ON            |  optional    |

Table: Cmake options for the convergence test simulation. \label{tab:convtest_cmakeoptions}

#### Mesh Generation with **HOPR**

The mesh files for this tutorial can be found in the tutorial directory. Parameter files are provided to produce the meshes using
**HOPR**. With *parameter_hopr.ini* conform meshes are created and the file *parameter_mortar_hopr.ini* generates non-conform
meshes.

#### Manufactured Solution For The Inviscid Case

The manufactured solution is

\begin{equation}
\rho = A*(1 + B*\sin(\Omega*|\underline{x} - \underline{v}t|)) \quad \text{with} \quad A,B,\Omega,\underline{v} = const
.
\end{equation}

The advantage is that for the Euler equation, the resulting source term is zero, $Q(x,t) \equiv 0$.

As mesh, a Cartesian box is used with periodic boundaries. The mesh and the corresponding solution is shown in \ref{fig:convtest_mesh_and_result}.


![](tutorials/04_convtest/convtest_mesh.png)   ![](tutorials/04_convtest/convtest_result.png) 
Figure: Mesh and flow field solution of the density. View in $x$-$y$-plane.\label{fig:convtest_mesh_and_result}

To investigate the order of convergence of a given polynomial degree $N$, the mesh resolution has to increase. We provide meshes with 1, 2, 4 and 8 elements. They are provided in the tutorial directory with an according parameter file for the preprocessing tool HOPR.

#### Flow Simulation with FLEXI

The inviscid convergence test is run from the parameter file *parameter_convtest_flexi.ini*. Essentially, any valid parameter file can be used, as a manufactured solution is simulated. This allows to test the various methods and features of the code and investigate their order of convergence. However, for this tutorial we restrict the parameter file to a very simple baseline test case.  The following entries have to be made

**IniExactFunc = 2** : The manufactured solution is applied by means of the exact function. This function is used to initialize **FLEXI** and can in general be used for Dirichlet boundary conditions. Here, we apply only periodic boundary conditions. 

**N\_Analyze = 10** : The number of interpolation nodes for the analyze routines, needed for calculating the error norms. We suggest at least $N\_Analyze = 2N$.

**AdvVel = (/0.3,0.,0./)** : The constant velocity vector $\underline{v}$

**CalcErrorNorms = T** : Flag to turn on the calculation of L2 and L_inf error norms

**ProjectName = ConvTest** : Project name used to name the output files of **FLEXI** and the convergence test scripts.


#### Numerical settings

The default settings for the time integration are displayed in table \ref{tab:convtest_num_set}.

Table: Numerical settings for time integration \label{tab:freestream_num_set}

| Variable        | Description                            | Value         |
| --------------- |:--------------------------------------:|:-------------:|
| tend            |                                        | 0.5           |
| Analyze_dt      |                                        | 0.5           |
| nWriteData      |                                        | 1             |
| CFLscale        |                                        | 0.9           |
| DFLscale        |                                        | 0.9           |


The remaining numerical settings necessary, e.g. the polynomial degree and the mesh filename are set via the script file. The script can be found in the directory 

       $FLEXIROOT/tools/convergence_test

Two scripts are provided, the file *convergence_grid* calculates the order of grid convergence for a given polynomial degree $N$ with increasing mesh resolution. The file *convergence* calculates spectral convergence on a given mesh with increasing polynomial degree. Both scripts are written in Python 2.7. In the first case of grid convergence, the polynomial degree and the set of meshes can be adjusted. Here we choose a polynomial degree of $3$, i.e. the theoretical order of convergence is $N+1 = 4$. The spectral convergence is calculated for polynomials of degree $N \in [1,10]$ on a mesh with 4 elements in each spatial direction.

The command

~~~~~~~
$FLEXIROOT/tools/convergence_test/convergence_grid $FLEXIDIR/bin/flexi parameter_convtest_flexi.ini
~~~~~~~

runs the code. The standard output of **FLEXI** is written into the logfile *ConvTest.log*. An ASCII file *ConvTest_convfile_grid.csv* is written that includes all L2 and L_inf error norms for the state vector $U$ for all meshes and the corresponding orders of convergence. Furthermore, a PDF file *ConvTest_convtest_grid.pdf* is generated that plots the L2 error of the momentum in $x$-direction against the number of elements of the meshes. Another curve represents the theoretical order of convergence for the chosen polynomial degree.

Spectral convergence can be investigated using the command

~~~~~~~
$FLEXIROOT/tools/convergence_test/convergence $FLEXIDIR/bin/flexi parameter_convtest_flexi.ini
~~~~~~~

Corresponding files are produced, where *\_grid* is replaced by *\_N*. Figure \ref{fig:convtest_convergenceplots} shows the result for grid convergence.

![Plot of grid convergence \label{fig:convtest_convergence_grid}](tutorials/04_convtest/ConvTest_convtest_grid.pdf)

Figure \ref{fig:convtest_convergence_N} shows the result for spectral convergence.

![Plot of spectral convergence \label{fig:convtest_convergence_N}](tutorials/04_convtest/ConvTest_convtest_N.pdf) 


### Viscous Convergence Test

The second convergence test includes the viscous terms. Make sure, **FLEXI** is compiled with the setting listed in the table \ref{tab:convtest_visc_cmakeoptions} to include parabolic
fluxes.

| Option                          | Value         | Comment      |
| ------------------------------- |:-------------:| ------------:|
| CMAKE_BUILD_TYPE                | Release       |              |
| FLEXI_EQYNSYSNAME               | navierstokes  |              |
| FLEXI_PARABOLIC                 | ON            |              |
| FLEXI_MPI                       | ON            |  optional    |

Table: Cmake options for the convergence test simulation. \label{tab:convtest_visc_cmakeoptions}


For this case, another manufactured solution is chosen

\begin{equation}
\rho = 2 + A*\sin(\Omega*|\underline{x}| - v\pi t) \quad \text{with} \quad A,\Omega,v = const
.
\end{equation}

The same function is applied to the momentum in all spatial directions. The mass specific total energy in this case is $\rho e = \rho \rho$.
This manufactured solution has a non-zero source term. In **FLEXI**, this source term is added in the routine *CalcSource* in the file 

         $FLEXIROOT/src/equations/navierstokes/equations.f90


Note that this manufacture can also be solved without considering the viscous terms. In this case the source term does not vanish.

The parameter file for this case is *parameter_convtestvisc_flexi.ini*. Some modifications have to be made, compared to the inviscid convergence test. The correct exact function for this test case is *IniExactFunc = 4*. Note, that a viscosity has to be prescribed. Our default value for this case is $Mu0 = 0.03$.

Execution of the convergence tests is analogously to the inviscid case, e.g.

~~~~~~~
$FLEXIROOT/tools/convergence_test/convergence_grid $FLEXIDIR/bin/flexi parameter_convtestvisc_flexi.ini
~~~~~~~


### Execution of the Convergence Test with non-conform meshes


**FLEXI** is capable of using non-conform meshes with mortar interfaces. To illustrate the execution of **FLEXI** with these meshes,
four mortar meshes are provided for the convergene test. Their file names include the term *MORTAR*. To use the convergence test
script, just open the script file at *$FLEXIROOT/tools/convergence_test/convergence_grid* and exchange the file names of the meshes.
The script can then be run similar as described above.

### Execution of Convergence Test Outside of Tutorials

The convergence test scripts are provided in the directory

       $FLEXIROOT/tools/convergence_test

including the Python script to execute **FLEXI**

       $FLEXIROOT/tools/convergence_test/execute_flexi.py

The mesh files and the parameter files are also available in the directory

       $FLEXIROOT/ini/convtest

