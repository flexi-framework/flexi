## Flow over a flat plate with RANS equations
\label{sec:tut_flatPlate}

In this tutorial, the simulation with the RANS equations of a flow over a flat plat is considered. First, we explain how the parameters for the RANS equations have to be set. Next, it is shown how implicit time discretization can help to save computational time.

Copy the ``naca0012`` tutorial folder to your desired working directory.

        cp -r $FLEXI_TUTORIALS/flatPlate .

### Compiler options
        
Make sure that **FLEXI** is compiled with the cmake options listed in the following table.


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

Table: Cmake options for the RANS flat plate simulation. \label{tab:flatplate_cmakeoptions}

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

Material properties are given in table \ref{tab:flatplate_materialproperties}. The RefState and viscosity are chosen such that the Reynolds number based on the position along the plate is equal to 5 million at x=1, and the Mach number based on the free stream velocity is 0.2.

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

### Boundary conditions

todo!

### Running the code 
We proceed by running the code in parallel. For example using 4 processors, use the following command

~~~~~~~
mpirun -np 4 flexi parameter.ini
~~~~~~~

On a 2017 laptop with core i5 processor, this simulation takes about 20 hours.

### Implicit time discretization

This simulation setup is now used to illustrate how to use implicit time discretization instead of explicit methods. In FLEXI several explicit and implicit Runge-Kutta methods are implemented. To get a list of the available options run the command

~~~~~~~
flexi --help
~~~~~~~

Here is a list how the result could look like

~~~~~~~
TimeDiscMethod        =      CarpenterRK4-5 ! Specifies the type of time-discretization to be used, e.g. the name of a  
                                            ! specific Runge-Kutta scheme. Possible values:  
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

The list starts with explicit methods which are followed by the implicit ones starting from ``eulerimplicit``. Implicit time discretization methods allow for a higher CFL number than the explicit ones as they are typically unconditionally stable. The drawback of implicit methods is that they require the solution of large non-equation systems in each timestep.
In FLEXI the arising non-linear equation system is solved with Newton's method. This leeds to the need of a linear solver for each Newton's step. For this the iterative GMRES method is used.
The use of iterative methods such as Newton's method and GMRES requires the definition of convergence criteria and iteration parameters. Below, they are summarized shortly:
~~~~~~~
!============================================================================================================================
! Implicit
!============================================================================================================================
adaptepsNewton        =                   F ! Adaptive Newton convergence criterion by Runge-Kutta error estimation  
EpsNewton             =             0.1E-02 ! Newton tolerance, only used if adaptepsNewton=F  
nNewtonIter           =                  50 ! Maximum amount of Newton iterations  

EisenstatWalker       =                   F ! Adaptive abort criterion for GMRES  
gammaEW               =                 0.9 ! Parameter for Eisenstat Walker adaptation  
EpsGMRES              =             0.1E-02 ! GMRES Tolerance, only used of EisenstatWalker=F  
nRestarts             =                  10 ! Maximum number of GMRES Restarts  
nKDim                 =                  30 ! Maxmim number of Krylov subspaces for GMRES, after that a restart is performed  
~~~~~~~

Note that the adaptive Newton tolerance is only available for ``esdirk2-3``,``esdirk3-4`` and ``esdirk4-6``. When this option is true, Newton's convergence criterion is smaller the smaller the timestep is, hence this option has to be used with care. An adaptive GRMES convergence criterion means that a coarser tolerance is chosen for the first Newton steps and is decreased for increasing Newton steps. If Newton's method or GMRES do not converge during calculation the abort criteria and the maximum iterations have to be adjusted in the parameter file.
 
To accelerate the both, Newton's method and GMRES, different options are included in FLEXI. We start with Newton's method:
As the convergece property of Newton's method highly depends on the starting value, there are different options to choose. Note that the method 3 (dense output extrapolation) is only available for ``esdirk3-4`` and ``esdirk4-6``.
~~~~~~~
!============================================================================================================================
! Implicit
!============================================================================================================================
PredictorType         =                   0 ! Type of predictor to be used, 0: use current U, 1: use right hand side, 2:  
                                            ! polynomial extrapolation, 3: dense output formula of RK scheme  
PredictorOrder        =                   0 ! Order of predictor to be used (PredictorType=2)  
~~~~~~~
Typically chosing a predictor is beneficial for 'small' CFL numbers as here the extrapolation gives good results. For 'larger' timesteps chosing a predictor other than the current U can increase the required amount of iterations.

The solution of the linear system arising from Newton's method requires the inversion of the large system matrix. As FLEXI is designed for large scale parallel simulations it is not desirable to build up this large matrix and multiply it with the state vector but to approximate it via a finite difference. Parameters how this finite difference is built are summarized below.

~~~~~~~
!============================================================================================================================
! Implicit
!============================================================================================================================
FD_Order              =                   1 ! Order of FD approximation (1/2)  
Eps_Method            =                   2 ! Method of determining the step size of FD approximation of A*v in GMRES, 1:  
                                            ! sqrt(machineAccuracy)*scaleps, 2: take norm of solution into account  
scaleps               =                 1.0 ! Scaling factor for step size in FD, mainly used in Eps_Method=1  
~~~~~~~

Typically those parameters do not have to be changed. When simulating flows with very low Mach numbers the degradation of the accuracy of the finite difference by machine precission can become an issue. In such cases the order of the finite difference can be increased. The method of how to calculate the step size of the finite difference sometimes has to be adjusted when using nondimensional equations where e.g. the Mach or Reynolds number explicitly occur.

To accelerate the convergence of GMRES a preconditioner can be applied. The preconditioner approximates the inverse of the system matrix (Jacobian). Again, several options are available:
~~~~~~~
!============================================================================================================================
! Preconditioner
!============================================================================================================================
PrecondType           =                   1 ! Preconditioner Type (0: no Preconditioner, 1: analytic, 2: finite difference, 3:  
                                            ! compute both and compare)  
PrecondIter           =                   1 ! Defines how often preconditioner is built  
SolveSystem           =                   1 ! Solver of the preconditioned system (0: exact LU inversion, 1: inexact ILU(0)  
                                            ! inversion, always with NoFillIn=T)  
DebugMatrix           =                   0 ! Write Jacobians to file for debug purposes (0: no output, 1: non-inverted  
                                            ! matrix,2: additionally inverted matrix, 3: additionally check inversion  
                                            ! accuracy)  
HyperbolicPrecond     =                   F ! Preconditioner only for the hyperbolic flux  
NoFillIn              =                   F ! Precond for parabolic system forced to have the same sparsity as the Euler  
                                            ! Precond  
DoDisplayPrecond      =                   F ! Display building time of preconditioner  
~~~~~~~

The parameters ``PrecondType``,``DebugMatrix`` and  ``DoDisplayPrecond``only have to be changed if analysis concerning the preconditioner shall be conducted. In an application case those parameters should remain unchanged. An important parameter is ``SolveSystem``: Using the exact LU inversion is sometimes necesarry to obtain convergence of the linear solver. Nevertheless the application of the inexact ILU(0) is typically faster, especially for higher orders and higer dimensions (3d).

### Runing the simulation with implicit time discretization

The implicit method FLEXI uses (Newton-GMRES) is well suited for unsteady calculations as here the initial guess for Newton's method is relatively good. Starting a simulation from an unphysical homogenous state causes very poor convergence properties of Newton's method. Hence, it is beneficial to start the current simulation with an explicit time discretization method and than do a restart after several timesteps have been calculated.
For the current simulation we suggest to simulate until ``tend=0.02`` with an explicit scheme and than use this state for a restart with ``eulerimplicit``.
