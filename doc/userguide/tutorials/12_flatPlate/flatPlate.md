## Flow over a flat plate with RANS equations
\label{sec:tut_flatPlate}

**FLEXI** is mainly designed for scale-resolving simulations of turbulent flow, i.e. direct numerical simulation (DNS)
 or large eddy simulation (LES). Unfortunately, the required resolution of the grids in both space and time for those
simulation approaches exceeds the available computational power by several orders of magnitude for many applications that
are relevant to e.g. engineering applications.

An alternative for those situations is to solve the so-called Reynolds Averaged Navier-Stokes or RANS equations. The
RANS equations are derived by applying a temporal filter to the Navier-Stokes equations, which must be wide enough to filter
out any fluctuating components of the solution variables. One is thus left with a set of equations that describe a
stationary solution. Formally, those filtered equations will look like the original Navier-Stokes equations, except for
closure terms that appear due the filtering of the non-linear fluxes. Those closure terms depend on the unfiltered
variables, and thus need to be modelled. For RANS, the closure terms must contain the effects of all scales except for the
stationary flow field. The modelling of those is a challenge, and many different approaches have developed over the years.
In **FLEXI**, a rather simple but widely used model after Spalart and Allmars (SA model, see [@allmaras2012modifications])
is implemented. It belongs to the class of one-equation models, and solves a single additional equation for a turbulent viscosity
$\tilde{\mu}$. All effects of the unresolved scales are thus modeled by assuming an increased viscosity of
$\mu+\mu_{turb}(\tilde{\mu})$.

In this tutorial, we will simulate the flow over a flat plate using the RANS equations. For the presented setup, the Reynolds number at the
streamwise position $x=1m$ along the plate will be equal to $Re_x = \frac{\rho_{\infty} u_{\infty}x}{\mu}=5 \cdot 10^6$. The mach
number is set to $Ma=0.2$.

First, we explain how the parameters for the RANS equations have to be set. Next, it is shown how implicit time discretization
can help to save computational time.

Copy the ``naca0012`` tutorial folder to your desired working directory.

        cp -r $FLEXI_TUTORIALS/flatPlate .

### Compiler options

Make sure that **FLEXI** is compiled with the cmake options listed in the following table. Note that since only a stationary
flow over the plate is considered, the example is purely two-dimensional, which additionally saves computational resources.


| Option                          | Value         | Comment      |
| ------------------------------- |:-------------:| ------------:|
| CMAKE_BUILD_TYPE                | Release       |              |
| FLEXI_2D                        | ON            |              |
| FLEXI_EQYNSYSNAME               | rans_sa       |              |
| FLEXI_FV                        | OFF           |              |
| FLEXI_LIFTING                   | br1           |              |
| FLEXI_PARABOLIC                 | ON            |              |
| FLEXI_MPI                       | ON            |  optional    |
| FLEXI_NODETYPE                  | GAUSS         |              |

Table: CMake options for the RANS flat plate simulation. \label{tab:flatplate_cmakeoptions}

~~~~~~~~~~~
ccmake [flexi root directory]
~~~~~~~~~~~

Set the above options and then compile the code by issuing

~~~~~~~~~~~
make
~~~~~~~~~~~


### Flow Simulation with FLEXI

The simulation setup is defined in ``parameter.ini``.

``IniRefState = 1`` : the initial condition uses ``RefState 1`` for the initial flow field solution.

``IniExactFunc = 1`` : exact function routine for initialization, case 1 initializes a freestream state based on ``IniRefState``.

Since we solve an additional equation for the SA model, the ``RefState`` will also have an additional entry compared to the Navier-Stokes
simulations. In the sixth entry, a value for $\tilde{\nu}=\frac{\tilde{\mu}}{\rho}$ must be set. For fully turbulent simulations, a value
of $\tilde{\mu}=\frac{\mu}{10}$ is appropriate.

Material properties are given in table \ref{tab:flatplate_materialproperties}.

| Property                        | Variable      | Value       |
| ------------------------------- |:-------------:| -----------:|
| dynamic viscosity $\mu$         | mu0           | 0.0000171   |
| isentropic coefficient $\kappa$ | kappa         |  1.4        |

Table: Material properties set in the parameter file \label{tab:flatplate_materialproperties}

### Numerical settings

The DG solution on the mesh is represented by piecewise polynomials and the polynomial degree in this tutorial is chosen as $N=4$.

The main code settings are displayed in table \ref{tab:flatplate_num_set}.


| Variable        | Description                            | Value         |
| --------------- |:---------------------------------------|--------------:|
| N               | Polynomial degree                      | 4             |
| MeshFile        | Mesh file to be used                   |meshes/CART_HEX2D_FlatPlateC_mesh.h5|
| tend            | end time of the simulation             | 0.25          |
| Analyze_dt      | time interval for analysis             | 0.001         |
| nWriteData      | dump solution every n'th Analyze_dt    | 50            |
| CFLscale        |                                        | 0.9           |
| DFLscale        |                                        | 0.9           |

Table: Numerical settings \label{tab:flatplate_num_set}

### Mesh and boundary conditions

The mesh itself is of a simple rectangular shape. On the left site and on the top, the free-stream value given by the
initialization is prescribed directly. On the right, a pressure-based outflow condition is used. Since the flow
is stationary, no additional treatment of the outflow is necessary. On the bottom of the domain, for $0.2m$ a symmetry condition is applied,
since we don't want to directly start with the flat plate at the inflow. After $x=0.2m$, an adiabatic wall begins. The mesh is stretched, such
that the nodes in wall-normal direction are clustered towards the wall. The gradients in (even the stationary) turbulent boundary layer are very
large, such that we need a lot of resolution in that direction. The grid is also stretched in the wall-parallel direction, and nodes are clustered
at the leading edge of the plate. Here, larger gradients are expected since the boundary layer is developing very rapidly at the start. At later
$x$-positions, we can use very large cells in the wall-parallel direction, since not a lot of change is going on there.

Several versions of the mesh are included in the subfolder ``mesh``. They only differ in the amount of cells used in the wall-normal direction.
The main influence is the size of the cell closest to the wall, which needs to resolve the largest gradient. Mesh *A* is the finest,
with a cell size of $y^+ \approx 4$, expressed in [wall-units](https://en.wikipedia.org/wiki/Law_of_the_wall) in the middle of the flat plate.
 This does not yet take into account that there are multiple nodes inside of a cell for DG methods! Mesh *C*, which is used in the default setup,
in combination  with the polynomial degree of $N=4$ leads to $y^+ \approx 4$ for the first **node**. Though this is much larger than the
usual requirement of $y^+ \approx 1$, due to the high order nature of the simulation this is enough to gain results that are in reasonable agreement
with experimental and other numerical reference solutions.

### Wall distance

The additional equation for the SA model uses the distance from the closest solid wall as an input. Since this is a static quantity, we can
compute this in a pre-processing step and then simply read that information if we start our actual simulation. For the default setup,
the file containing the information about the wall distance is already included, in the file ``meshes/CART_HEX2D_FlatPlateC_walldistance.h5``.
If you want to use other meshes or a different polynomial degree, you will need to generate the file yourself. For this purpose, the pre-processing
tool ``walldistance`` is included in **FLEXI**. You need to change some compile options so it will be available,  which are listed
in the following table.

| Option                          | Value         | Comment      |
| ------------------------------- |:-------------:| ------------:|
| FLEXI_BUILD_POSTI               | ON            |              |
| POSTI_BUILD_WALLDISTANCE        | ON            |              |

Table: CMake options for the walldistance tool. \label{tab:flatplate_cmakeoptions}

Once the tool is build, you can run it by issuing a command like

~~~~~~~~~~~
posti_walldistance walldistance.ini
~~~~~~~~~~~

The tool only takes a parameter file as input. A sample ``.ini`` can be found in the ``meshes`` subfolder, you need to specify the required
polynomial degree and the mesh to use. The program will then use a gradient descent technique to find the distance to the closest wall
for each solution point. Since for Gauss-Lobatto nodes the boundaries and thus walls are included in the set of nodes, we would have
points where the distance is zero. Unfortunately, the SA equation divide by that distance, which does not work. Thus, only use the RANS
equation system with Gauss nodes!

### Running the code
We proceed by running the code in parallel. For example using 4 processors, use the following command

~~~~~~~
mpirun -np 4 flexi parameter.ini
~~~~~~~

On a 2017 laptop with core i5 processor, this simulation takes about 20 hours until the default end-time is reached. Note
that the time here has no physical meaning, since we are interested in the stationary solution. The analyze routines will report
the remaining residual for each quantity ($L_2$ norm of the spatial DG operator). These are indicators for how much the
solution is still changing. We call our solution converged if a certain threshold for the residuals is reached, or if they
do not get smaller anymore.

### Implicit time discretization

This simulation setup is now used to illustrate how to use implicit time discretization instead of explicit methods. Since the mesh contains very small cells
to resolve the gradient at the wall, the explicit time step is really small. Also, we don't require a time-accurate solution of the problem, instead we want to converge to
a stationary solution as fast as possible. This makes our simulation perfect for implicit time discretization.
 In FLEXI several explicit and implicit Runge-Kutta methods are implemented. To get a list of the available options run the command

~~~~~~~
flexi --help
~~~~~~~

Here is a list how the result could look like

~~~~~~~
TimeDiscMethod        =      CarpenterRK4-5 ! Specifies the type of time-discretization to be used, e.g. the name of a specific Runge-Kutta scheme. Possible values:
                                     !   * standardrk3-3
                                     !   * carpenterrk4-5
                                     !   * niegemannrk4-14
                                     !   * toulorgerk4-8c
                                     !   * toulorgerk3-7c
                                     !   * toulorgerk4-8f
                                     !   * ketchesonrk4-20
                                     !   * ketchesonrk4-18
                                     !   * eulerimplicit
                                     !   * cranknicolson2-2
                                     !   * esdirk2-3
                                     !   * esdirk3-4
                                     !   * esdirk4-6
~~~~~~~

The list starts with explicit methods which are followed by the implicit ones starting from ``eulerimplicit``. Implicit time discretization methods allow for a higher CFL number than the explicit ones as they are typically unconditionally stable. The drawback of implicit methods is that they require the solution of large non-linear equation systems in each timestep.
In FLEXI the arising non-linear equation system is solved with Newton's method. This necessitates a linear solver for each Newton step. For this the iterative GMRES method is used.
Iterative methods such as Newton's method and GMRES require the definition of convergence criteria and iteration parameters. Below, they are shortly summarized:

~~~~~~~
!===========================================
! Implicit
!===========================================
adaptepsNewton        =                   F ! Adaptive Newton convergence criterion by Runge-Kutta error estimation
EpsNewton             =             0.1E-02 ! Newton tolerance, only used if adaptepsNewton=F
nNewtonIter           =                  50 ! Maximum amount of Newton iterations

EisenstatWalker       =                   F ! Adaptive abort criterion for GMRES
gammaEW               =                 0.9 ! Parameter for Eisenstat Walker adaptation
EpsGMRES              =             0.1E-02 ! GMRES Tolerance, only used if EisenstatWalker=F
nRestarts             =                  10 ! Maximum number of GMRES restarts
nKDim                 =                  30 ! Maximum number of Krylov subspaces for GMRES, after that a restart is performed
~~~~~~~

Note that the adaptive Newton tolerance is only available for ``esdirk2-3``,``esdirk3-4`` and ``esdirk4-6``. When this option is true, Newton's convergence criterion is smaller the smaller the timestep is, hence this option has to be used with care. An adaptive GRMES convergence criterion means that a coarser tolerance is chosen for the first Newton steps, which is decreased for increasing Newton steps. If Newton's method or GMRES do not converge during calculation the abort criteria and the maximum iterations have to be adjusted in the parameter file.

To accelerate both, Newton's method and GMRES, different options are included in FLEXI. We start with Newton's method:
There are different options to choose as the convergence property of Newton's method depends highly on the start value. Note that the method 3 (dense output extrapolation) is only available for ``esdirk3-4`` and ``esdirk4-6``.

~~~~~~~
!===========================================
! Implicit
!===========================================
PredictorType         =                   0 ! Type of predictor to be used, 0: use current U, 1: use right hand side, 2: polynomial extrapolation, 3: dense output formula of RK scheme
PredictorOrder        =                   0 ! Order of predictor to be used (PredictorType=2)
~~~~~~~
Typically a predictor is beneficial for 'small' CFL numbers as here the extrapolation gives good results. For 'larger' timesteps, a predictor which does not use the current U can increase the required amount of iterations.

The solution of the linear system which arises from Newton's method requires the inversion of the large system matrix. As FLEXI is designed for large scale parallel simulations, it is not desirable to build up this large matrix and multiply it with the state vector but to approximate it via a finite difference. The parameters for the finite difference are summarized below.

~~~~~~~
!===========================================
! Implicit
!===========================================
FD_Order              =                   2 ! Order of FD approximation (1/2)
Eps_Method            =                   2 ! Method to determine the step size of the FD approximation of A*v in GMRES, 1: sqrt(machineAccuracy)*scaleps, 2: take norm of solution into account
scaleps               =                 1.0 ! Scaling factor for step size in FD, mainly used in Eps_Method=1
~~~~~~~

In general, the order of the finite difference approximation can be kept at 1. But in flows with very low Mach numbers or when very little change in a stationary solution is happening, the degradation of the accuracy of the finite difference by machine precision can become an issue. In such cases the order of the finite difference can be increased, as it is done in this tutorial. The method of how to calculate the step size of the finite difference sometimes has to be adjusted when using nondimensional equations where e.g. the Mach or Reynolds number explicitly occur.

To accelerate the convergence of GMRES a preconditioner can be applied. The preconditioner approximates the inverse of the system matrix (Jacobian). Again, several options are available:

~~~~~~~
!===========================================
! Preconditioner
!===========================================
PrecondType           =                   1 ! Preconditioner Type (0: no Preconditioner, 1: analytic, 2: finite difference, 3: compute both and compare)
PrecondIter           =                   1 ! Defines how often preconditioner is built
SolveSystem           =                   1 ! Solver of the preconditioned system (0: exact LU inversion, 1: inexact ILU(0) inversion, always with NoFillIn=T)
DebugMatrix           =                   0 ! Write Jacobians to file for debug purposes (0: no output, 1: non-inverted matrix,2: additionally inverted matrix, 3: additionally check inversion accuracy)
HyperbolicPrecond     =                   F ! Preconditioner only for the hyperbolic flux
NoFillIn              =                   F ! Precond for parabolic system forced to have the same sparsity as the Euler Precond
DoDisplayPrecond      =                   F ! Display building time of preconditioner
~~~~~~~

The parameters ``PrecondType``,``DebugMatrix`` and  ``DoDisplayPrecond``only have to be changed if analysis concerning the preconditioner shall be conducted. In an application case those parameters should remain unchanged. An important parameter is ``SolveSystem``: Using the exact LU inversion is sometimes necessary to obtain convergence of the linear solver. Nevertheless the application of the inexact ILU(0) is typically faster, especially for higher orders and higher dimensions (3d).

### Runing the simulation with implicit time discretization

The implicit method FLEXI uses (Newton-GMRES) is well suited for unsteady calculations as here the initial guess for Newton's method is relatively good. Starting a simulation from an unphysical homogeneous state causes very poor convergence properties of Newton's method. Hence, it is beneficial to start the current simulation with an explicit time discretization method and than do a restart after several timesteps have been calculated.
For the current simulation we suggest to simulate until ``tend=0.02`` with an explicit scheme and than use this state for a restart with ``eulerimplicit``. The CFL number can then be massively increased, e.g. to $100$. The following plot shows how the residuals of the different conservative variables then decrease over time. The residuals are the $L_2$-Norm of the spatial DG operator at a stage of the time-discretization method. They decrease by several orders of magnitude before reaching a plateau. Here, the Newton algorithm will detect that nearly no change is happening in the iterative procedure, and no iterations will be performed any more. A stationary solution has been reached.

![Residuals over time for the implicit part of the flat plate computation](tutorials/12_flatPlate/residuals.png)
