## Taylor Green Vortex 
\label{sec:tut_tgv}

This tutorial describes how to set up and run the basic test case for turbulent flows, the Taylor-Green-Vortex (TGV) - see e.g. [@gassner2013accuracy]. We will learn how to avoid catastrophic failure of the code due to non-linear instabilities. This is done by using
polynomial de-aliasing. In a second step we add the sub grid scale model of Smagorinsky. The tutorial assumes that you are familiar with the general **FLEXI** and **HOPR** work flow (please finish the previous tutorials first if this sounds strange to you).

### Flow description

The initial condition to the (TGV) is a sinus distribution in the u and v velocity components. This leads to rapid production of turbulent structures, after a short initial laminar phase. While the test case is incompressible in principle,
we solve it here in a compressible setting. The chosen Mach number with respect to the highest velocity in the field is $0.1$. The Reynolds number of the flow is defined as $1/\nu$. The domain is set up as a triple periodic box with edge length $2\pi$.

![](tutorials/05_taylorGreenVortex/dns_reference.png)

### Compiler options
        
Make sure that **FLEXI** is compiled with the cmake options listed in the following table.


| Option                          | Value              | Comment      |
| ------------------------------- |:-------------:     | ------------:|
| CMAKE_BUILD_TYPE                | Release            |              |
| FLEXI_EQYNSYSNAME               | navierstokes       |              |
| FLEXI_PARABOLIC                 | ON                 |              |
| LIBS_USE_MPI                    | ON                 |  optional    |
| FLEXI_EDDYVISCOSITY             | ON                 |  optional    |
| FLEXI_TESTCASE                  | taylorgreenvortex  |  optional    |

Table: CMake options for the cavity simulation. \label{tab:tgv_cmakeoptions}

For others you may keep the default values. Compile the code.

#### Mesh Generation with HOPR

We use a mesh with 4 cells per direction for the tutorial. In case you want to generate other meshes the parameter file for **HOPR** is included in the tutorial directory (parameter_hopr.ini),
the default mesh is included. Using 4 cells with a polynomial degree of $N=7$, means we use the typical large eddy setup of $32$ DOF per direction.

### Tutorial - Flow at $Re=1600$

Copy the ``tgv`` tutorial folder to your working directory

        cp -r $FLEXI_TUTORIALS/taylorgreenvortex .
        
        
Step into the folder. In case you do not want to generate the mesh files yourselves, the meshes have already been provided. 

#### Preparing the Flow Simulation with FLEXI

The simulation setup is defined in *parameter_flexi.ini*. To get help on any of the parameters listed therein, you can run **FLEXI** from the command line by typing

     flexi --help
     

The parameters in the file are grouped thematically, however, this is not mandatory. All lines starting with a "!" are comments. 

##### Output 

The test case has its own analyze output (PROJECTNAME_TGVAnalysis.csv), that we will use. We don't look at flow visualization in this tutorial. Besides other interesting quantities, the file contains the incompressible dissipation rate. This is the resolved dissipation of the gradient field, computed as the integral over the domain of the strain rate tensor norm $S_{ij}S_{ij}$, times viscosity times $2$. It is stored in the second column of the file. We will use this quantity in the tutorial to verify your results.

##### Interpolation / Discretization parameters

~~~~~~~

    ! ================================================ !
    ! INTERPOLATION
    ! ================================================ !
    N             = 7  

~~~~~~~

The parameter *N* sets the degree of the solution polynomial, e.g. in this example, the solution is approximated by a polynomial of degree $7$ in each spatial direction. This results in $(N+1)^3$ degrees of freedom for each (3D) element. In general, *N* can be chosen to be any integer greater or equal to $1$, however, the discretization and the timestep calculation has not extensively been tested beyond $N\approx 23$. Usually, a good compromise of performance and accuracy is found for $N\in[3,..,9]$.

To apply polynomial de-aliasing there are the following options:

~~~~~~~

    ! ================================================ !
    ! OVERINTEGRATION (ADVECTION PART ONLY)
    ! ================================================ !
    OverintegrationType=0  ! 0:off 1:cut-off filter 
                           ! 2: conservative cut-off 
    NUnder        = 7      ! specifies effective polydeg 
                           ! (modes > NUnder are thrown away)
                           ! for types 1 and 2

~~~~~~~

**FLEXI** has three ways of doing polynomial de-aliasing. Mode 0: don't do it. Mode 1: a filter is applied to the time-update $(J*U_t)$. The filter is formulated as a Galerkin projection of degree $N$ to $NUnder$, the effective resolution is thus $NUnder$. Mode 2: in principle identical to Mode 1, but takes into account non-linear metric terms. For the linear mesh of this tutorial the result is identical, while Mode 2 is slightly more computational expensive, so we omit it.

For **FLEXI** we can run under-resolved computations without sub grid scale model. The only artifical dissipation is then provided by the Riemann solver used for the inter-cell fluxes. You can change the Riemann solver to see the effect with the following parameters:

~~~~~~~
    ! ================================================ !
    ! Riemann
    ! ================================================ !
    Riemann=  RoeEntropyFix ! Riemann solver to be used: 
                            ! LF, HLLC, Roe,  
                            ! RoeEntropyFix, HLL, HLLE, HLLEM  
~~~~~~~
To add Smagorinsky's model set the following parameter to $1$, here CS is the Smagorinsky constant usually chosen around $0.1$ for isotropic turbulence (such as TGV). 

~~~~~~~
    ! ================================================ !
    ! LES MODEL
    ! ================================================ !
    eddyViscType = 0       ! Choose LES model, 1:Smagorinsky
    CS = 0.02              ! Smagorinsky constant
    PrSGS = 0.6            ! turbulent Prandtl number
~~~~~~~

#### Running the Simulation and Results

The command

~~~~~~~
flexi parameter_flexi.ini > std.out
~~~~~~~

runs the code and dumps all output into the file *std.out*. If you wish to run the code in parallel using MPI, the standard command is

~~~~~~~
mpirun -np XX flexi parameter_flexi.ini > std.out
~~~~~~~

where $XX$ is an integer denoting the number of processes to be used in parallel. Note that **FLEXI** uses an element-based parallelization strategy, so the minimum load per process/core is *one* grid element, i.e. do not use more cores than cells in the grid! 
**FLEXI** writes a TGV analyze file (PROJECTNAME_TGVAnalysis.csv). You can use yor favorite plotting programm to visualize the ASCII data. Using gnuplot you can create plots with the following syntax:

~~~~~~
gnuplot
set key autotitle columnhead
set datafile separator ","
plot "PROJECTNAME_TGVAnalysis.csv" u 1:2 w l
~~~~~~
Where you replace PROJECTNAME with the projectname you defined in the parameter_flexi.ini file.

#### Part I: crashing simulation

First we run **FLEXI** without overintegration/de-aliasing. We will find that the code crashes, once scale production becomes relevant. You can compare your result to the plot in the tutorial folder 

![](tutorials/05_taylorGreenVortex/crash_no_dealiasing.png)
      


#### Part II: Overintegration
We now use overintegration by changing the respective settings in the parameter_flexi.ini file as described above. For ``Overintegration==1`` set $N=11$ and $NUnder=7$. You can compare your result to the plot below.

![](tutorials/05_taylorGreenVortex/les_dealiasing.png)



#### Part III: Explicit LES model

To see the effect of adding explicit eddy viscosity we activate the LES model (Smagorinsky) as described above. To obtain the reference result of the following plot set $CS=0.1$. Don't forget to switch overintegration of again and set polynomial degree to $N=7$. Feel free to play around with the constant, have fun!

![](tutorials/05_taylorGreenVortex/les_smago_oi.png)
