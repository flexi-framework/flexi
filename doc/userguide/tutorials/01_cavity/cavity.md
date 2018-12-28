## Lid-driven Cavity
\label{sec:tut_cavity}

This tutorial describes how to set up and run the first non-trivial flow problem. The lid-driven cavity flow is a standard test case for numerical schemes, and a number of results have been published in literature, see e.g. [@ghia_cavity], [@gao_cavity]. This tutorial assumes that you have completed the previous tutorial, know how to edit files and postprocess the solution with your favorite visualization tool, e.g. ParaView. Also, the later parts of the tutorial assume that you have access to a computer with an MPI-based parallelization with at least 4 computing cores - otherwise, it will just take a lot longer :). 

The tutorial is split into two parts: The basic part will teach you about setting up the code and running the simulations. The advanced part will build on this and give you a glimpse on how to make modification to code to accommodate more complex simulation and add features you might need. If you are just interested in running the code as is, you may skip the advanced part, or just complete parts of it. 

### Flow description
The flow under consideration is essentially incompressible and two-dimensional, but we will use the three-dimensional code for the compressible Navier-Stokes equations to solve it here. This is not the most efficient way to compute this flow, but it works well to show you how to set up and run a simulation in this tutorial. The computation will be conducted in a three-dimensional, square domain with periodic boundary conditions in the "third" direction. The wall of the cavity are modelled as isothermal walls, and a fixed flow is prescribed at the upper boundary, i.e. the lid of the domain. For the Reynolds numbers investigated here, this generates a steady, vortical flow field in the cavity. 

The following picture (figure \ref{fig:cavity_re400_velmag}) shows the resulting velocity field and streamlines for $Re=400$.

![Contours of velocity magnitude for the $Re=400$ lid-driven cavity case. \label{fig:cavity_re400_velmag}](tutorials/01_cavity/re400_velcontours.png)

### Compiler options
        
Make sure that **FLEXI** is compiled with the cmake options listed in the following table.


| Option                          | Value         | Comment      |
| ------------------------------- |:-------------:| ------------:|
| CMAKE_BUILD_TYPE                | Release       |              |
| FLEXI_EQYNSYSNAME               | navierstokes  |              |
| FLEXI_PARABOLIC                 | ON            |              |
| FLEXI_MPI                       | ON            |  optional    |

Table: Cmake options for the cavity simulation. \label{tab:cavity_cmakeoptions}

The standard settings are sufficient for this example. To check whether they are set, change to your ``build`` folder and open the cmake GUI 

~~~~~~~~~~~
ccmake [flexi root directory]
~~~~~~~~~~~

If necessary, set the above options and then compile the code by issuing

~~~~~~~~~~~
make
~~~~~~~~~~~


### Part I: Basic Tutorial - Flow at $Re=100$

Copy the ``cavity`` tutorial folder \label{missing:aliases_tutorial_cavity1}

        cp -r $FLEXI_TUTORIALS/cavity .
        
        
Step into the ``Basic_Re100`` subfolder. The contents of the folder should look like in picture (figure \ref{fig:cavity_dirls}). In case you do not want to generate the mesh files yourselves, the meshes have already been provided. 

![Files for the basic lid-driven cavity tutorial. \label{fig:cavity_dirls}](tutorials/01_cavity/dirls.png)        
        


#### Mesh Generation with HOPR
The domain of interest consists of a square 2D geometry. Although the flow field is two-dimensional, we will create a three-dimensional domain here and apply periodic boundary conditions in the $z-$direction. Also, we will only use *one* element in that direction to save computational costs. 

Like before, the mesh files used by **FLEXI** are created by supplying an input file *parameter_hopr.ini* with the appropriate information to **HOPR**.

    ./hopr parameter_hopr.ini

A full tutorial on how to run **HOPR** is available at the project site [HOPR-project](http://hopr-project.org), so we will just review the basics here. In *parameter_hopr.ini*, the following lines create a cubical domain from $[0,0,0]\rightarrow[1,1,1]$, discretized by $4\times 4$ elements in the $x-y$ plane and $1$ element in the $z-$direction. 

~~~~~~~
    Corner =(/0.,0.,0. ,,1.,0.,0. ,,1.,1.,0. ,,  0.,1.,0.,, 0.,0.,1. ,,1.,0.,1. ,,1.,1.,1. ,,  0.,1.,1. /)
    ! Corner node positions: (/ x_1,y_1,z_1, x_2,y_2,z_2,..... , x_8,y_8,z_8/)
    nElems =(/4,4,1/)              ! number of elements in each direction
~~~~~~~

The next line maps the 6 faces of the cube to the boundary conditions following the CGNS notation, i.e. it indicates which of the boundary conditions defined below belong to which face. 

~~~~~~~
      BCIndex      =(/1,3,6,4,5,2/)        ! Indices of Boundary Conditions for  six Boundary Faces (z-,y-,x+,y+,x-,z+)
~~~~~~~

For example, the first boundary condition defined will be applied to the $z-$ face of the cube, while the third one will be applied to the $y-$ face and so on. The order $(z-,y-,x+,y+,x-,z+)$ is given by the CGNS standard and thus fixed, so if you keep this line as it is, the boundary conditions below can be defined in the (less confusing) order of  $(z-,z+,y-,y+,x-,x+)$.

![BC Indices for the cube \label{fig:cavity_bcindex_front}](tutorials/01_cavity/bcindex_front.png)![BC Indices for the cube \label{fig:cavity_bcindex_back}](tutorials/01_cavity/bcindex_back.png)

Figure: BC Indices for the cube, front and back view. \label{fig:cavity_bcindex}

The boundary conditions are now assigned to faces in the following lines: 

~~~~~~~
    BoundaryName=BC_zminus             ! BC index 1 
    BoundaryType=(/1,0,0,1/)           ! (/ Type, curveIndex, State, alpha /)
    BoundaryName=BC_zplus              ! BC index 2
    BoundaryType=(/1,0,0,-1/)          ! 
    vv=(/0.,0.,1./)                    ! 

    BoundaryName=BC_wall_lower         ! BC index 3
    BoundaryType=(/4,0,1,0/)
    BoundaryName=BC_free               ! BC index 4
    BoundaryType=(/2,0,0,0/)

    BoundaryName=BC_wall_left          ! BC index 5
    BoundaryType=(/4,0,1,0/)
    BoundaryName=BC_wall_right         ! BC index 6
    BoundaryType=(/4,0,1,0/)

~~~~~~~

The identifier *BoundaryName* can be chosen freely, it is good practice to start it with *BC_*. The identifier *BoundaryType* specifies what kind of boundary is to be applied to the face. For the $z-$oriented faces, periodic boundary conditions are chosen, which are denoted by type $1$, the first index in the array. The second index indicates if the face is *curved* (not for this test case), and the third indicates which *reference state* will be used at the boundary. This will be discussed in more detail later, it is important to note that the state can be set here when generating the grid, but can also be overwritten in *parameter_flexi.ini* later. The fourth entry is used to match periodic boundary condition sides, i.e. these need to form a pair of the form $[a,-a]$, where $a$ is a unique integer for each periodic pair and the sign specifies the orientation. In this example, the $z-$ and $z+$ sides are connected and the associated connection vector *vv* is in the $z-$direction. In summary, to specify a periodic side pair, the first index in *BoundaryType* must be $1$, the fourth index must be a matching pair of integers and the connection vector must be specified. 

In the $y-$direction, we choose an isothermal wall boundary condition at the lower side (type $4$) and a Dirichlet / fixed state at the top of the cavity (type $2$). Again, the corresponding states can be set here, but overwritten later. 

In the $x-$direction, the left and right walls are again chosen to be isothermal wall boundary conditions. 

![BC Types for the lid driven cavity \label{fig:cavity_type_front}](tutorials/01_cavity/bctype_front.png)![BC Types for the lid driven cavity \label{fig:cavity_type_back}](tutorials/01_cavity/bctype_back.png)

Figure: BC Types for the lid driven cavity, front and back view. \label{fig:cavity_bctype}

For this basic tutorial, the simple meshes shown in picture (figure \ref{fig:cavity_meshes}) will be used. 

![Meshes for the basic lid-driven cavity tutorial. \label{fig:cavity_meshes}](tutorials/01_cavity/meshes_re100.png)   


#### Preparing the Flow Simulation with FLEXI

The simulation setup is defined in *parameter_flexi.ini*. To get help on any of the parameters listed therein, you can run **FLEXI** from the command line by typing

     ./flexi --help
     or
     ./flexi --help SECTION
     or
     ./flexi --help PARAMETER
     
The first command lists all parameters grouped thematically, however, this is not mandatory.
The second command lists all parameters of a certain section and the last command only gives help about a specific option. 
All lines starting with a "!" are comments. 

##### Output parameters

~~~~~~~

    ! ================================================ !
    ! OUTPUT 
    ! ================================================ !
    ProjectName   = Tutorial_Cavity_Re100
    outputFormat  = 0  

~~~~~~~

The parameter *ProjectName* is used to name all the files generated by the simulation, e.g. the state at time $0.3$ (the solution vector of conserved variables at each node)  is written to the file *Tutorial_Cavity_Re100_State_0000000.300000000.h5* and so on. The parameter *outputFormat* is set to $0$ here, indicating that the output is *off*.  Turning it on will result in an on-the-fly **visualization** output of the solution states, but this is *not* recommended as good practice - it will slow down the code considerably, in particular in parallel mode. State-files in the HDF5-format are written always independent of this parameter. The recommended way is to use the *posti_visu* tool after the simulation to generate Paraview-readable files. Note that currently, the only output format available is *vtk* (by setting outputFormat  = 3), since Tecplot output is not available due to GPL licensing issues. 

##### Interpolation / Discretization parameters

~~~~~~~

    ! ================================================ !
    ! INTERPOLATION
    ! ================================================ !
    N             = 3  
    NAnalyze      = 10 

~~~~~~~

The parameter *N* sets the degree of the solution polynomial, e.g. in this example, the solution is approximated by a polynomial of degree $3$ in each spatial direction. This results in $(N+1)^3$ degrees of freedom for each (3D) element. In general, *N* can be chosen to be any integer greater or equal to $1$, however, the discretization and the timestep calculation has not extensively been tested beyond $N\approx 23$. Usually, for a good compromise of performance and accuracy is found for $N\in[3,..,9]$.

*NAnalyze* determines the polynomial degree for the analysis routines, e.g. the accuracy of the calculation of error norms or testcase specific integrals during the computation. A good rule of thumb is to set $NAnalyze=2\times N$.

##### Mesh parameters

~~~~~~~

    ! ================================================ !
    ! MESH  
    ! ================================================ !
    MeshFile      = cavity2x2_mesh.h5 
    useCurveds    = F

    !BoundaryName=BC_wall_lower             
    !BoundaryType=(/4,1/)           
    BoundaryName=BC_wall_left         
    BoundaryType=(/4,1/)                
    BoundaryName=BC_wall_right       
    BoundaryType=(/4,1/) 
    BoundaryName=BC_free
    BoundaryType=(/2,1/)

~~~~~~~

The parameter *MeshFile* contains the name of the **HOPR** mesh file in HDF5 format (and/or the full path to it). *UseCurveds* indicates whether the mesh is considered to be curved, i.e. if high-order mesh information should be used. Setting this to *F* can be used to discard high order information in the mesh file and treat it as a linear mesh. For the current tutorial, the meshes are linear by design. 

The boundary conditions are set via the *BoundaryName* identifier, which specifies a name that must be present in the mesh file (see section on mesh generation above). Each line containing the boundary name must be followed by a line containing the *BoundaryType*. A list of types available for the Navier-Stokes equations can be found in table \ref{tab:boundaryconditions}. For types that require additional information (like Dirichlet boundaries), the second index in *BoundaryType* refers to the *RefState* (short for reference state) which is used to determine the unknowns / quantities from the outside for this boundary condition. For example, for a Dirichlet boundary (Type 2), the full reference state is set at the boundary, so the lines 

~~~~~~~

    BoundaryName=BC_free
    BoundaryType=(/2,1/)

~~~~~~~

indicate that the *first* reference state vector listed below is set at this boundary (the lid part of the cavity). Note that the reference vectors are always in primitive variables, i.e. $(\rho, u, v, w, p)^T$ *unless* specified otherwise. Also, for the boundaries *BC_wall_left* and *BC_wall_right*, the same reference state is chosen. These boundaries are isothermal walls, so a wall temperature needs to be specified - this is computed from the primitive variables in the associated reference state. 

Note that the lines for the lower wall boundary are commented out. As has been mentioned above, the boundary conditions can be set in **HOPR** directly, and then be overwritten here. This is just an example to show how both approaches work. Later, when running **FLEXI** with these settings, it is good practice to inspect the boundary condition information as understood by **FLEXI**. In this case, the output of **FLEXI** to the console will look like this:

~~~~~~~
    ----------------------------------------------------------------
    |          BoundaryName |     BC_wall_left | *CUSTOM | 
    |          BoundaryType |       (/ 4, 1 /) | *CUSTOM | 
    |          BoundaryName |    BC_wall_right | *CUSTOM | 
    |          BoundaryType |       (/ 4, 1 /) | *CUSTOM | 
    |          BoundaryName |          BC_free | *CUSTOM | 
    |          BoundaryType |       (/ 2, 1 /) | *CUSTOM | 
    |     Boundary in HDF file found |  BC_free
    |                            was |  2 0
    |                      is set to |  2 1
    |     Boundary in HDF file found |  BC_wall_left
    |                            was |  4 1
    |                      is set to |  4 1
    |     Boundary in HDF file found |  BC_wall_right
    |                            was |  4 1
    |                      is set to |  4 1
    .............................................................
    BOUNDARY CONDITIONS|           Name   Type  State     Alpha
    |                  |      BC_zminus      1      0         1
    |                  |       BC_zplus      1      0        -1
    |                  |  BC_wall_lower      4      1         0
    |                  |        BC_free      2      1         0
    |                  |   BC_wall_left      4      1         0
    |                  |  BC_wall_right      4      1         0
    .............................................................

~~~~~~~



##### Equation parameters


In this section of the parameter file, the settings associated with the equation to be solved are set, including initial conditions and reference data for the boundary conditions. The parameter *IniExactFunc* specifies which solution or function should be used to fill the initial solution vector, i.e. it specifies what the starting flow field looks like. Setting this to $1$ selects a uniform initial state in the whole domain. Note that this solution is also used to compute the errors norms and can also be time-dependent. The state itself is defined by the parameter *IniRefstate*, in this case, the second one is used. 



~~~~~~~
    ! ============================================== !
    ! EQUATION
    ! ============================================== !
    IniExactFunc  = 1
    IniRefState   = 2
    RefState      = (/1.0,1.,0.,0.,71.4285714286/)
    RefState      = (/1.0,0.,0.,0.,71.4285714286/)
    mu0           = 0.01
    R             = 1
    Pr            = 0.72
    kappa         = 1.4
~~~~~~~

As described above, the RefStates are given in primitive variables. The second one describes a fluid at rest and is used to initialize a resting fluid in the cavity. The first state is used to determine the driving flow at the top of the cavity (and to compute the wall temperatures for the boundaries). Constant flow properties like the gas constant are given in table \ref{tab:cavity_flow_prop} 
and define the gas behavior in combination with the ideal gas law $p=\rho R T$. Note that the code itself does not distinguish between dimensional and non-dimensional quantities and it is the user's responsibility to set all data consistently. For anything other than an ideal gas with constant viscosity and heat conductivity, physically meaningful quantities should be set. 

Table: Numerical settings \label{tab:cavity_flow_prop}

| Property                        | Variable      | Value       |
| ------------------------------- |:-------------:| -----------:|
| dynamic viscosity $\mu$         | mu0           |  0.1        |
| ideal gas constant $R$          | R             |  1          |
| isentropic coefficient $\kappa$ | kappa         |  1.4        |
| Prandtl number            $\Pr$ | Pr            |  0.72       |

From these settings, the Mach- and Reynolds number can be computed as follows, taking into account a reference cavity length of $1$ and the magnitude of the driving velocity:  
\begin{equation}
Mach=u/c=1.0/\sqrt{\kappa \frac{p}{\rho}}=1.0/10.0=0.1
\end{equation}
\begin{equation}
Re=\frac{u L \rho}{\mu_0}=\frac{1.0}{0.01}=100
\end{equation}

Since we are comparing against an incompressible reference solution, setting the Mach number to $0.1$ is a good compromise between accuracy and efficiency. 


#### Temporal discretization parameters
~~~~~~~
	! ================================================ !
	! TIMEDISC
	! ================================================ !
	tend          = 5.0
	Analyze_dt    = 0.1
	nWriteData    = 1
	CFLscale      = 0.9
	DFLscale      = 0.4

~~~~~~~

The parameter *tend* determines the end time of the solution, *Analyze_dt* the interval at which the analysis routines (like error computation, checking of wall velocities etc. see below) are called. The multiplier *nWriteData* determines the interval at which the solution state files are written to the file system, e.g. in this case *nWriteData* $\times$ *Analyze_dt* = 0.1 is the interval for writing to disc. The CFL and DFL numbers determine the explicit time step restriction for the advective and viscous parts. Note that these values should always be chosen to be $<1$, but since the determination of the timestep includes some heuristics, both values should be chosen conservatively. 


#### Analysis parameters

~~~~~~~
	! ============================================================== !
	! ANALYZE
	! ============================================================== !
	CalcErrorNorms=   T   ! Calculate error norms
	CalcBodyForces=   T   ! Calculate body forces
	CalcWallVelocity= T   ! Calculate wall velocities
	CalcMeanFlux=     T   ! Calculate mean flux through boundaries

~~~~~~~

These switches trigger the output of analysis files. 


Now that all the necessary preparations have been made, the simulation can be started. 


#### Running the Simulation and Results
\label{sec:tut_cavity_running_and_results}

The command

~~~~~~~
./flexi parameter_flexi.ini > std.out
~~~~~~~

runs the code and dumps all output into the file *std.out*. If you wish to run the code in parallel using MPI, the standard command is

~~~~~~~
mpirun -np XX ./flexi parameter_flexi.ini > std.out
~~~~~~~

where $XX$ is an integer denoting the number of processes to be used in parallel. Note that **FLEXI** uses an element-based parallelization strategy, so the minimum load per process/core is *one* grid element, i.e. do not use more cores than cells in the grid! 

If the simulation runs without error, the screen output should look similar to figure \ref{fig:cavity_re100_stdout}, showing information on the timestep and producing the analysis results selected via the parameter file. Once the run has completed, your working folder will contain a number of *State* files in HDF5 format and * *.dat* ASCII files for the analysis data, which can be opened in Tecplot, a spreadsheet tool etc. 

![Sample screen output of the cavity simulation. \label{fig:cavity_re100_stdout}](tutorials/01_cavity/re100_stdout.png)

Since we start the simulation from a fluid at rest, it will take some iterations / timesteps to achieve a steady state solution. One way to check if the solution has converged to a steady state is to check some characteristic quantities, e.g. the velocities at the walls, body forces etc. In figure \ref{fig:cavity_re100_wallvel}, the temporal evolution of the velocities at the lower wall are plotted over time for 4 different simulations, and for all cases, at $t_{end}=5$, a sufficiently stationary solution has been achieved. Note that since the boundary conditions are applied *weakly* in a DG setting, a velocity slip at walls can occur, with its magnitude depending on the local wall resolution. 


![Time evolution of wall velocity at lower wall for the $Re=100$ lid-driven cavity case. \label{fig:cavity_re100_wallvel}](tutorials/01_cavity/re100_wallvel.png)

For comparison with published data, 4 simulations where run on the provided meshes. Details are listed in table \ref{tab:cavity_sims}.

Table: $Re=100$ lid driven cavity: Simulations \label{tab:cavity_sims}

| Simulation    | Mesh          |     N          |
| ------------- |:-------------:| --------------:|
| 1             | cavity_2x2    |     3          |
| 2             | cavity_4x4    |     3          |
| 3             | cavity_8x8    |     3          |
| 4             | cavity_2x2    |     7          |

![Contours of velocity magnitude for the $Re=100$ lid-driven cavity case. \label{fig:cavity_re100_velcontours}](tutorials/01_cavity/re100_velcontours.png)


\pagebreak 

A contour plot of the velocity magnitude at the end time is given in figre \ref{fig:cavity_re100_velcontours}. To generate this plot, convert the *State* files to a Paraview format by invoking


~~~~~~~
posti_visu parameter_postiVisu.ini parameter_flexi.ini your_Statefiles.h5
~~~~~~~


For a more quantitative comparison with published data, you can generate a plot of the $u-$velocity on the centerline ($x=0$) of the cavity. Figure \ref{fig:cavity_re100_u_over_y} shows the results for the 4 simulations run here, along with published data available in [@ghia_cavity], [@gao_cavity]. 

![Comparison of centerline velocities for the $Re=100$ lid-driven cavity case with published results. \label{fig:cavity_re100_u_over_y}](tutorials/01_cavity/re100_u_over_y.png)

For simulation 1, the agreement with literature results is fair. This is due to the coarse resolution with $2\times (3+1) = 8$ degrees of freedom per x- and y-direction. Doubling the grid elements results in an improved match, doubling it again (simulation 3), the agreement with the published data is excellent. It should be noted that the same accuracy can be achieved by increasing $N$ and keeping the coarse grid. Simulation 3 and 4 have nearly identical results, although the number of degrees of freedom differs by a factor of 2. This is a feature of high order schemes for smooth problems. 

### Part II: Advanced Tutorial - Flow at $Re=400$

In this part of the tutorial, we will expand upon what has been learned in the basic section above. The general setup of the computation remains the same, but the Reynolds number is increased, which will be accounted for by a new mesh and a higher resolution. Also, we will give a first glimpse at how to make your own changes to the code and show how to add a new function to be used as a custom initial or boundary condition. It is recommend that you have completed the basic tutorial, have access to at least 4 cores for the computation (or just have the patience to wait longer) and are familiar with **FORTRAN** syntax. 

If not done so already,copy the ``cavity`` tutorial folder \label{missing:aliases_tutorial_cavity1}

        cp -r $FLEXI_TUTORIALS/cavity .
        
        
and step into the ``Advanced_Re400`` subfolder.

      


#### Mesh Generation with HOPR
To account for the increased Reynolds number, the number of elements in the $x-y$-plane is increased to $12\times 12$. Also, a stretching in the $y-$direction is introduced, as can be seen in picture \ref{fig:cavity_re400_stretchmesh}. In *parameter_hopr.ini*, the following line facilitates an exponential stretching:


~~~~~~~

	factor       =(/1.0,-1.2,1./)                ! element stretching a constant growth factor (+/- changes direction)

~~~~~~~

You can generate your own mesh or re-use the provided one, labeled *cavity12x12_stretch_mesh.h5*.

![Stretched mesh for the $Re=400$ lid-driven cavity case. \label{fig:cavity_re400_stretchmesh}](tutorials/01_cavity/re400_stretchmesh.png)

#### Adding a custom initial or boundary function
Before running the simulation, we will include a custom function to be used as boundary condition for the top boundary, i.e. the velocity driving the flow in the cavity. In this, we follow the suggestions from [@gao_cavity], where the $u(x)$ velocity at the lid is given as

\begin{equation}\label{eq:cavity:code}
    u(x)= 
\begin{cases}
    c_1 x^4 + c_2 x^3 + c_3 x^2 +c_4 x,& \text{if } 0\le x <0.2\\
    d_1 x^4 + d_2 x^3 + d_3 x^2 +d_4 x + d_5,& \text{if } 0.8 < x \le 1.0\\
    1,              & \text{otherwise}
\end{cases}
\end{equation}
with 
\begin{equation}
\begin{split}
[c_1,c_2,c_3,c_4]=& 1000\times[4.9333,-1.4267,0.1297,-0.0033]\\
[d_1,d_2,d_3,d_4,d_5]=& 10000\times[0.4933,-1.8307,2.5450,1-.5709,0.3633]\\
\end{split}
\end{equation}

This assumes the top bondary to be from $0 \ldots 1$, as is the case in our domain. To add this new function to **FLEXI**, locate the file "$FLEXI/src/equations/navierstokes/idealgas/exactfunc.f90" and open it in the text editor of your choice. Locate the *SUBROUTINE* *ExactFunc*, which provides functions to the boundary, initial condition and analysis routines of the code. The header of the routine you are looking for is shown in figure \ref{fig:cavity_code1}.

![Header for the subroutine *ExactFunc* \label{fig:cavity_code1}](tutorials/01_cavity/code1.png)


To add equation \ref{eq:cavity:code} to the code, add a new *CASE* to the routine, in which you define the state vector for the primitive variables, and then convert them to conservative ones (see e.g. case 8 for how this is done). You might also need to introduce some new local variables for this routine. To check if your changes are syntactically correct, compile the code with your changes by invoking *make* from your build directory. If the compilation process was not successful, check the output on the screen that will give you some hints on what might be wrong. 

To benefit from this tutorial, it is recommend that you do try to complete this programming task. For reference, figure \ref{fig:cavity_code2} shows a sample piece of code with a correct implementation as *CASE 9*. From now on, we will refer to *CASE 9* as the case number in question, yours might be different of course. 

![Sample code for the custom equation from [@gao_cavity]. \label{fig:cavity_code2}](tutorials/01_cavity/code2.png)

Note that in the sample code, we do not specify the full primitive state vector within the routine, but re-use the *IniRefState* and just overwrite the $u$-velocity. This is a matter of choice, but it allows to set the Mach number by setting the *IniRefState* accordingly.  


#### Preparing the Flow Simulation with FLEXI
To setup the simulation, you can either modify the *parameter_flexi.ini* files from the basic tutorial or use the provided ones. We will conduct 2 simulations: one with the constant driving flow boundary conditions as before, and one with the new custom equation specified as equation *CASE 9*. Two parameter files for these cases are provided, the one with the custom boundary condition named *parameter_flexi_custombc.ini*. For both cases, the changes in the parameterfile are:

~~~~~~~
		MeshFile      = cavity12x12_stretch_mesh.h5
		...

		mu0           = 0.0025
                ...

		tend          = 100.0

~~~~~~~


The mesh file needs to be adjusted as well as the viscosity to set the Reynolds number. Also, the end time is set to $100$ to account for the longer accomondation period the simulation will need to achieve a steady state. Setting it to $100$ is very conservative, as we will see later, feel free to lower it as you see fit. 

For the custom boundary condition case, the parameter file needs to be adjusted to address the new boundary settings. According to table \ref{tab:boundaryconditions}, Dirichlet boundary conditions with a specified reference equation (instead of a reference _state_) are of type $22$. The second index in the entry refers thus to the equation (*CASE 9*) we have programmed above. 

~~~~~~~
	BoundaryName=BC_free
	BoundaryType=(/22,9/)
~~~~~~~

Now we are ready to run both simulations and compare the results. 


#### Running the Simulation and Results

The command

~~~~~~~
./flexi parameter_flexi.ini > std.out
~~~~~~~

runs the code and dumps all output into the file *std.out*. If you wish to run the code in parallel using MPI, the standard command is

~~~~~~~
mpirun -np XX ./flexi parameter_flexi.ini > std.out
~~~~~~~

![Evolution of wall velocities at the lower wall for $Re=400$ lid driven cavity simulations. \label{fig:cavity_re400_wallvel}](tutorials/01_cavity/re400_wallvel.png)

Change the parameter file to *parameter_flexi_custombc.ini* for the second run. Note that you can adust the end time if you wish. From figure \ref{fig:cavity_re400_wallvel}, which shows the evolution of the mean velocity at the lower wall, the flow reaches steady state after about $t=35$. The following figures show the flow field and comparison of the centerline velocities with published results.






![](tutorials/01_cavity/re400_velmag_const.png)     ![](tutorials/01_cavity/re400_velmag_custom.png)

Figure: Steady state solution of velocity magnitude of $Re=400$ lid driven cavity. Left: Constant boundary condition, Right: Custom boundary condition.


![Evolution of wall velocities at the lower wall for $Re=400$ lid driven cavity simulations. \label{fig:cavity_re400_u_over_y}](tutorials/01_cavity/re400_u_over_y.png)
