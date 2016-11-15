## Freestream

The setup considers a freestream scenario with constant pressure $p=101325.0$ Pa, density $\rho=1.225$ kg/$m^3$ and flow vector $\textbf{U}=(1,1,1)^T$ m/s.

![](tutorials/00_freestream/freestream_mesh.png)     ![](tutorials/00_freestream/freestream_result.png)

Figure: Mesh and flow field solution with velocity vector plot of the velocity field. View in $x$-$y$-plane.\label{fig:freestream_mesh_and_result}


Copy the ``freestream`` tutorial folder \label{missing:aliases_tutorial_freestreem}

        cp -r $FLEXI_TUTORIALS/freestream .

### Mesh Generation with HOPR

The mesh files used by **FLEXI** are created by supplying an input file *parameter_hopr.ini* with the appropriate information.

    ./hopr parameter_hopr.ini

This creates the mesh file *cartbox_mesh.h5* in HDF5 format. Alternatively, if you do not want to run **hopr**, you can also use the provided mesh.

### Flow Simulation with FLEXI

The simulation setup is defined in *parameter_flexi.ini*. The initial condition is selected via the variable vector
**RefState=(/1.225,1.,1.,1.,101325./)** which represents the solution vector $(\rho, u, v, w, p)^T$. **FLEXI** allows for multiple 
**RefState** vectors and numerates them respectively for them to be used by different functions. In this example a single **RefState** is 
supplied and therefore is given the number **1**.


**IniRefState = 1** : the initial condition uses **RefState 1** for the initial flow field solution.

**IniExactFunc = 1** : the used exact function routine uses **RefState 1**, e.g., for the calculation of the $L_2$ error norm.

Constant flow properties like the gas constant are given in table \ref{tab:freestream_flow_prop} 
and define the gas behavior in combination with the ideal gas law $p=\rho R T$.

Table: Numerical settings \label{tab:freestream_flow_prop}

| Property                        | Variable      | Value       |
| ------------------------------- |:-------------:| -----------:|
| dynamic viscosity $\mu$         | mu0           | 0.000018547 |
| ideal gas constant $R$          | R             |  276        |
| isentropic coefficient $\kappa$ | kappa         |  1.4        |

### Numerical settings

The DG solution on the mesh is represented by piecewise polynomials and the polynomial degree in this tutorial is chosen as $N=3$.

The default settings for these properties are displayed in table \ref{tab:freestream_num_set}. 

Table: Numerical settings \label{tab:freestream_num_set}

| Variable        | Description                            | Value         |
| --------------- |:--------------------------------------:|:-------------:|
| N               | polynomial degree                      | 3             |
| MeshFile        |                                        |cartbox_mesh.h5|
| tend            |                                        | 1e-6          |
| Analyze_dt      |                                        | 1e-6          |
| nWriteData      |                                        | 1             |
| CFLscale        |                                        | 0.99          |
| DFLscale        |                                        | 0.4           |


The command

~~~~~~~
flexi parameter_flexi.ini > std.out
~~~~~~~

runs the code and dumps all output into the file *std.out*. 
If the run has completed successfully, which should take only a brief moment, the contents of the working folder should look like in figure \ref{fig:freestream_folder}

![The folder contents after a successful run\label{fig:freestream_folder}](tutorials/00_freestream/freestream_folder.png)

Two additional files have been created, which are are named  **Projectname_State_Timestamp.h5**. They contain the solution vector of the conserved variables at each interpolation node at the given time, which corresponds to multiplies of **Analyze_dt**. If these files are not present, something went wrong during the execution of **FLEXI**. In that case, check the _std.out_ file for an error message. 

After a successful completion, the last lines in this files should look like in figure \ref{fig:freestream_stdout}

![The _std.out_ file after a successful run\label{fig:freestream_stdout}](tutorials/00_freestream/freestream_stdout.png)

To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**. Issue the command 

~~~~~~~
flexi2vtk parameter_flexi2vtk.ini cartbox_State_000000*.h5
~~~~~~~
to generate the corresponding *vtu*-files, which can then be loaded into **ParaView**. 
