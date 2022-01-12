## Flow around a cylinder

In this tutorial, the simulation around a two-dimensional circular cylinder at $Re_D=200$ and
$Ma=0.2$ is considered. The goal of this tutorial is to introduce the usage of a new
functionality of the **posti_visualizerecordpoints** tool and the new tool **posti_dmd**.

### Compiler options

Make sure that **FLEXI** is compiled with the CMake options listed in the following table.


| Option                             | Value         | Comment      |
| -----------------------------------|:-------------:| ------------:|
| CMAKE_BUILD_TYPE                   | Release       |              |
| FLEXI_2D                           | ON            |              |
| FLEXI_EQYNSYSNAME                  | navierstokes  |              |
| FLEXI_PARABOLIC                    | ON            |              |
| LIBS_USE_MPI                       | ON            |  optional    |
| POSTI_RP_VISUALIZE                 | ON            |              |
| POSTI_DMD                          | ON            |              |

Table: Cmake options for the cylinder simulation. \label{tab:cylinder_cmakeoptions}

To check whether they are set, change to your ``build`` folder and open the CMake GUI

~~~~~~~~~~~
ccmake [flexi root directory]
~~~~~~~~~~~

If necessary, set the above options and then compile the code by issuing

~~~~~~~~~~~
make
~~~~~~~~~~~

\pagebreak


### Mesh Generation with HOPR

The mesh file used by **FLEXI** is created by **HOPR**

    ./hopr parameter_hopr.ini

This creates the mesh file *Cylinder_Re200_mesh.h5* in HDF5 format. If **HOPR** is not available, the mesh file is supplied in this tutorial.


### Flow Simulation with FLEXI

The simulation setup is defined in ``parameter_flexi.ini``. The initial condition is selected via the variable vector ``RefState=(/1.,1.0,0.,0.,17.857/)`` which represents the vector of primitive solution variables $(\rho, u, v, w, p)^T$.


Material properties are given in table \ref{tab:cylinder_materialproperties}. Based on the ideal gas law, we get

$$Ma=1/\sqrt{\kappa p/\rho}=0.2$$

Note that in this non-dimensional setup the mesh is scaled such that the reference length is unity, i.e. $D=1$. Then to arrive at $Re=\rho u D / \mu = 200$, the viscosity is set to

$$ \mu = \rho u D / Re = 1/Re = 0.005 $$


| Property                        | Variable      | Value       |
| ------------------------------- |:-------------:| -----------:|
| dynamic viscosity $\mu$         | mu0           |  0.005      |
| ideal gas constant $R$          | R             |  17.857     |
| Prandtl number                  | Pr            |  0.72       |
| isentropic coefficient $\kappa$ | kappa         |  1.4        |

Table: Material properties set in the parameter file \label{tab:cylinder_materialproperties}

### Numerical settings

The DG solution on the mesh is represented by piecewise polynomials and the polynomial degree in this tutorial is chosen as $N=4$.

The main code settings are shown in table \ref{tab:cylinder_num_set}.


| Variable        | Description                            | Value         |
| --------------- |:---------------------------------------|--------------:|
| N               | Polynomial degree                      | 4             |
| MeshFile        | Mesh file to be used                   |Cylinder_Re200_mesh.h5|
| tend            | end time of the simulation             | 300           |
| Analyze_dt      | time interval for analysis             | 0.01          |
| nWriteData      | dump solution every n'th Analyze_dt    | 500            |
| CFLscale        |                                        | 0.9           |
| DFLscale        |                                        | 0.9           |

Table: Numerical settings \label{tab:cylinder_num_set}

### Boundary conditions

The boundary conditions were already set in the mesh file by **HOPR**. Thus, the simulation runs without specifying the boundary conditions in the **FLEXI** parameter file. The freestream boundaries of the mesh are Dirichlet boundaries using the same state as the initialization, the wall is modeled as an adiabatic wall. The boundary conditions in $z$ direction are not relevant for this 2D example, but would be realized as periodic boundaries for a 3D simulation. All boundary conditions used are listed below.

~~~~~~~
         Name      Type     State     Alpha
   BC_cylinder         3         0         0
   BC_farfield         2         0         0
~~~~~~~


### Running the code
We proceed by running the code in parallel. For example, using 4 processors, use the following command

~~~~~~~
mpirun -np 4 flexi parameter_flexi.ini
~~~~~~~

The simulation runs for 300 convective time units to achieve periodic vortex
shedding, thus the simulation can take up to one to two hours.

### Evaluation of Strouhal number

The Strouhal number (which is a non-dimensional frequency, $Sr=\frac{f \cdot D}{u}$, describing the oscillatory motion of the flow) is estimated using the forces acting on the cylinder induced by the vortex shedding. The forces are calculated on the fly during runtime. The associated flags in the parameter file are

~~~~~~~~~~~~
CalcBodyForces=T
WriteBodyForces=T
~~~~~~~~~~~~

The first line activates the calculation of the forces at each ``Analyze_dt``, the second line enforces output of the forces to a file.  In figure \ref{fig:cylinder_bodyforce} the force in y-direction is plotted. By measuring the time from peak to peak over several periods the Strouhal number can be estimated to $0.1959$ which is close to the expected value from literature.


![Resulting forces on the airfoil up to $t=10$.\label{fig:cylinder_bodyforce}](tutorials/11_cylinder/cylinder_fy.png)


### Evaluation of the separation angle

The mean separation angle is evaluated using the **record points**-tool as
introduced in \ref{sec:postiRecordpoints}.
The simulation setup already contains the **record points** set and the record points
are written during the simulation. The **record points** set contains a plane
within the boundary layer of the upper cylinder side.
This time we want to use the **Plane_doBLProps** functionality within the **posti_visualizerecordpoints** tool. With this tool we want to analyze the boundary layer properties such as the wall friction to estimate the separation point.
The parameter needed are already set in the parameter_visualizeRecordpoints.ini
file.

You can run the tool using

~~~~~~~
posti_visualizerecordpoints parameter_visualizeRecordpoints.ini Cylinder_Re200_RP_*
~~~~~~~

After executing the tool, you will get a file named
**Cylinder_RP_BLProps_upperSide_BLPla000001.vts**
which can be visualized with
ParaView. The Data we want to visualize is one dimensional, so you won't be able
to see the data in the render view. To visualize it you need the apply the
plot over time filter. ParaView should automatically apply the correct range to
plot on. By plotting  **tau_w** over the circumference you can estimate the separation
angle to $113~deg$ (the intersection with the $y=0$ line).



### Dynamic Mode Decomposition

In this part of the tutorial we want to introduce the posti tool **posti_dmd**.
The dynamic mode decomposition is an algorithm divide a temporal series into a
set of modes which are associated with a frequency and grow/decay rate. With
this tool we are also capable to determine the Strouhal frequency. The dynamic
mode decomposition (DMD) is
implemented according to Schmid et al. [@schmid2009dynamic].

To use this tool, we need a higher temporal resolution of the written state files.
Thus, we change the time **tend** to $310$ and **nWriteData** to $1$. We restart
the simulation from the latest state file:

~~~~~~~
flexi parameter_flexi.ini Cylinder_Re200_State_0000300.000000000.h5
~~~~~~~


To execute the DMD on the density run the following command:

~~~~~~~
posti_dmd parameter_dmd.ini Cylinder_Re200_State_00003*
~~~~~~~

Depending on the available memory you might have to decrease the number of input
state files.
After execution you will see two additional files **Cylinder_Re200_DMD_0000300.000000000.h5** and **Cylinder_Re200_DMD_Spec_0000300.000000000.dat**.
The first file contains the field representation of the different modes and the
second file contains the ritz spectrum of the modes.

To visualize the field run the following command:

~~~~~~~
posti_visu parameter_postivisu.ini Cylinder_Re200_DMD_0000300.000000000.h5
~~~~~~~

The new file **Cylinder_Re200_Solution_0000300.000000000.vtu** now contains five
modes to visualize. Figure \ref{fig:cylinder_modes} shows the steady, the global, the first and the
second harmonic mode. The global mode is the mode of the considered Strouhal
number.

![DMD modes of the density field. Top left steady mode, top right global mode, bottom left first harmonic, bottom right second harmonic. \label{fig:cylinder_modes}](tutorials/11_cylinder/dmd_modes.png)

With the python script *plot_RitzSpectrum.py* the Ritz spectrum of the DMD can
be plotted. The script is placed in the tools folder of **FLEXI**. To plot the
spectrum execute:

~~~~~~~
python plot_RitzSpectrum.py -d Cylinder_Re200_DMD_Spec_0000300.000000000.dat
~~~~~~~

The result is a Ritz spectrum as shown in fig \ref{fig:cylinder_spec}. On the
*x-axis* the frequency of the modes and on the *y-axis* the growth/decay factor
is plotted, whereat modes with $\omega_r<0$ are damped. The modes placed
directly on the x-axis are the already discussed modes, from left to right the global, the first, the second harmonic mode and so on. The color and size of the plotted modes represent the Euclidian norm of the mode which can be interpreted as an energy norm of the mode.

![Ritz spectrum. \label{fig:cylinder_spec}](tutorials/11_cylinder/RitzSpec.png)



