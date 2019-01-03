## Plane Turbulent Channel Flow
\label{sec:tut_ptcf}

This tutorial describes how to set up and run the Plane-Turbulent-Channel-Flow test case. We will learn how to use the split form DG method to guarantee non-linear stability of the turbulent channel flow. In a second step, we add the sub grid scale model of Smagorinsky combined with vanDriest damping to run stable wall-bounded turbulent flows with explicit small scale dissipation. The tutorial assumes that you are familiar with the general FLEXI and HOPR work flow (please finish the previous tutorials first if this sounds strange to you).

### Flow description

The initial condition of the channel flow test case is an analytical mean u velocity profile superimposed by a sinus distribution in the u, v and w velocity components. This leads to rapid production of turbulent channel flow structures. To maintain a constant mass flux, a constant pressure source term is superimposed. While the test case is incompressible in principle, we solve it here in a compressible setting. The chosen Mach number with respect to the highest velocity in the field is $Ma=0.1$ according to the Moser channel test case. The Reynolds number of the flow is defined as $Re_{\tau}=1/\nu$. The domain is set up as a Cartesian box with length $2\pi$, height $2$ and span width $\pi$ with periodic boundaries in the x- and z-directions as well as no-slip walls at the top and the bottom of the domain.   

### Compiler options
        
Make sure that **FLEXI** is compiled with the CMake options listed in the following table.


| Option                          | Value              | Comment      |
| ------------------------------- |:-------------:     | ------------:|
| CMAKE_BUILD_TYPE                | Release            |              |
| FLEXI_EQYNSYSNAME               | navierstokes       |              |
| FLEXI_PARABOLIC                 | ON                 |              |
| FLEXI_MPI                       | ON                 |  optional    |
| FLEXI_EDDYVISCOSITY             | ON                 |  optional    |
| FLEXI_NODETYPE                  | GAUSS-LOBATTO      |              |
| FLEXI_SPLIT_DG                  | ON                 |              |
| FLEXI_TESTCASE                  | channel            |              |
| FLEXI_BUILDPOSTI                | ON                 |              |
| POSTI_BUILD_CHANNEL_FFT         | ON                 |              |

Table: CMake options for the plane turbulent channel flow test case simulation. \label{tab:ptcf_cmakeoptions}

For all other CMake options you may keep the default values. Compile the code.

#### Mesh Generation with HOPR

We use a mesh with 4 cells per direction for the tutorial. In case you want to generate other meshes the parameter file for HOPR is included in the tutorial directory (*parameter_hopr.ini*),
the default mesh is included. Using 4 cells with a polynomial degree of $N=5$, means we use a large eddy simulation setup of $24$ DOFs per direction.

### Tutorial - Flow at $Re_{\tau}=180$

Copy the ``plane_turbulent_channel_flow`` tutorial folder to your working directory.

        cp -r $FLEXI_TUTORIALS/plane_turbulent_channel_flow .
        
Step into the folder. In case you do not want to generate the mesh files yourself, a default mesh has already been provided. 

#### Preparing the Flow Simulation with FLEXI

The simulation setup is defined in *parameter_flexi.ini*. To get help on any of the parameters listed therein, you can run **FLEXI** from the command line by typing

     ./flexi --help
     
The parameters in the file are grouped thematically, however, this is not mandatory. All lines starting with a "!" are comments. 

##### Output 

In this tutorial we don't look at the flow visualization of the instantaneous state files. Here, we will rather post process consecutive, instantaneous state files with the ``posti_channel_fft`` tool. As an output, we receive mean velocity and Reynolds stress profiles as well as turbulent energy spectra at different locations normal to the channel wall.

##### Interpolation / Discretization parameters

In this tutorial we use the split form DG method to guarantee non-linear stability of the turbulent channel flow simulation. As already specified in the CMake options table \ref{tab:ptcf_cmakeoptions}, the ``FLEXI_SPLIT_DG`` option has to be switched ON in combination with the ``FLEXI_NODETYPE`` ``GAUSS-LOBATTO``. **FLEXI** has five split flux formulations implemented. Therefore, a specific split flux formulation has to be set in the *parameter_flexi.ini* file. In this tutorial the pre-defined split flux formulation by Pirozzoli is used, which results in a kinetic energy preserving DG scheme. 

~~~~~~~

    ! ================================================ !
    ! SplitDG
    ! ================================================ !
    SplitDG       = PI     ! SplitDG formulation to be used: SD, MO, DU, KG, PI     

~~~~~~~

To switch on Smagorinsky's model set the eddyViscType to $1$ in the *paramerter_flexi.ini* file. In addition, the following parameters have to be set. CS is the Smagorinsky constant usually chosen around $0.11$ for wall bounded turbulent flows and the turbulent Prandtl number is commonly set to $0.6$. To guarantee a stable wall-bounded turbulent flow simulation in combination with Smagorinsky's model VanDriest damping has to be switched on to ensure zero eddy viscosity close to the channel walls.  

~~~~~~~
    ! ================================================ !
    ! LES MODEL
    ! ================================================ !
    eddyViscType = 0       ! Choose LES model, 1:Smagorinsky
    VanDriest = T          ! Van Driest damping for LES viscosity (channel flow only)
    CS = 0.11              ! Smagorinsky constant
    PrSGS = 0.6            ! turbulent Prandtl number
~~~~~~~

#### Running the Simulation and Results

The command

~~~~~~~
./flexi parameter_flexi.ini > std.out
~~~~~~~

runs the code and dumps all output into the file *std.out*. If you wish to run the code in parallel using MPI, the standard command is

~~~~~~~
mpirun -np XX ./flexi parameter_flexi.ini > std.out
~~~~~~~

where $XX$ is an integer denoting the number of processes to be used in parallel. Note that **FLEXI** uses an element-based parallelization strategy, so the minimum load per process/core is *one* grid element, i.e. do not use more cores than cells in the grid! 
Once the simulation finished state files can be post processed by the ``posti_channel_fft`` tool which was build by the ``POSTI_BUILD_CHANNEL_FFT`` CMake option. To run the postprocessing, the standard command is 

~~~~~~~
./posti_channel_fft parameter_channel_fft.ini [State1 State2 ...]
~~~~~~~

where the *parameter_channel_fft.ini* file is given in the tutorial folder and the amount of statefiles is specified by the user. In this tutorial we use all state files with a timestamp between $t=10.0$ and $t=15.0$. As an output you receive three different files. One containing the mean velocity profiles as well as the Reynolds stress profiles and the other two files contain turbulent erergy spectra. 
To visualize those files you can run the python script ``plotChannelFFT.py`` in the ``tools/testcases`` folder with the following command in your simulation directory

~~~~~~~
python $TOOLDIR/plotChannelFFT.py -p $PROJECTNAME -t $POSTITIME
~~~~~~~

where ``$TOOLDIR`` specifies the path to the ``tools/testcases`` folder, ``$PROJECTNAME`` the projectname specified in the *parameter_flexi.ini* file and ``$POSTITIME`` the timestamp of your output files from the ``posti_channel_fft`` tool.

#### Part I: SplitDG iLES
First, we run **FLEXI** without Smagorinsky's model which we call an implicit LES (iLES), as no explicit sub-grid scale dissipation model is added. The resulting mean velocity and Reynolds stress profiles as well as turbulent energy spectra close to the centre of the channel are given in Figure \ref{fig:Re180_turbulentChannel}.

![Mean velocity and Reynolds stress profiles (left) as well as turbulent energy spectra close to the centre of the channel (right} of an implicit LES at $Re_{\tau}=180$. \label{fig:Re180_turbulentChannel}](tutorials/10_planeTurbulentChannelFlow/Re180_turbulentChannel.png)

#### Part II: SplitDG with explicit LES model 
In a second step, we run **FLEXI** with Smagorinsky's model and VanDriest damping which needs to be switched on in the parameter file as described above. The resulting mean velocity and Reynolds stress profiles as well as turbulent energy spectra close to the centre of the channel are given in Figure \ref{fig:Re180_turbulentChannel_Smag}. In comparison to the previous simulation you might recognize the effect of the explicit damping on the Reynolds stress profile $\overline{u'u'}$ close to the maximum, most. To further study the influence of Smagorinsky's model play around with the spatial resolution both in terms of grid resolution as well as the polynomial degree N. You can also increase the Reynoldsnumber to $Re_{\tau}=395$ or $Re_{\tau}=590$ and compare the results to DNS results from Moser et al. [@moser1999direct]. 

![Mean velocity and Reynolds stress profiles (left) as well as turbulent energy spectra close to the centre of the channel (right} of a LES with Smagorinsky's model and vanDriest damping at $Re_{\tau}=180$. \label{fig:Re180_turbulentChannel_Smag}](tutorials/10_planeTurbulentChannelFlow/Re180_turbulentChannel_Smag.png)

