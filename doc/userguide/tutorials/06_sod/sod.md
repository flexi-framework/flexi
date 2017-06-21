## SOD Shock tube

The Sod shock tube example is one of the most basic test cases to investigate the shock capturing capabilities of a CFD code. A initial discontinuity is located in the middle of the one dimensional domain $[0,1]$. The left and right states are given by $\rho=1, v=0, p=1$ and $\rho=0.125, v=0, p=0.1$.
These states are already set as ``RefState`` in the the *parameter_flexi.ini* file.

Copy the ``sod`` tutorial folder 

        cp -r $FLEXI_TUTORIALS/sod .

### Mesh Generation with HOPR

The mesh files used by **FLEXI** are created by supplying an input file *parameter_hopr.ini* with the appropriate information.

    ./hopr parameter_hopr.ini

This creates the mesh file *SOD_100_mesh.h5* in HDF5 format.

### Flow Simulation with FLEXI

This example requires the Finite Volume shock capturing. Therefore turn the option ``FLEXI_FV`` in the cmake configuration on. Additionally you should switch ``FLEXI_PARABOLIC`` off and recompile the **FLEXI** code.
The simulation setup is defined in *parameter_flexi.ini* and now includes parameters for the Finite Volume shock capturing.  


| Parameter                     | Value       | Description                                                  |
| ----------------------------- | ----------- | -------------------------------------------------------------|
| IndicatorType                 | Persson     |                                                              |
| IndVar                        | 1           | first conservative (density) used for indicator evaluation   |
| IndStartTime                  | 0.001       | until this time FV is used in the whole domain               |
| FV_LimiterType                | MinMod      |                                                              |
| FV_IndUpperThreshold          | -3.         | upper threshold (if IndValue above this value, switch to FV) |
| FV_IndLowerThreshold          | -4.         | lower threshold (if IndValue below this value, switch to DG) |



The command

~~~~~~~
flexi parameter_flexi.ini 
~~~~~~~

runs the code and generates 5 state files **sod_State_TIMESTAMP.h5** for $t=0.0, 0.05, 0.10, 0.15, 0.20$.
To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**. Execute the command 

~~~~~~~
flexi2vtk parameter_flexi.ini sod_State_*.h5
~~~~~~~
to generate the corresponding *vtu*- and *vtm*-files, which can then be loaded into **ParaView**. 
There are two types of *vtu*-files, which contain either the DG or the FV part of the solution. 
The *vtm*-files combine the DG and FV *vtu*-file of every timestamp. Load the *vtm*-files into **ParaView**.

Since this example is a 1D test case use the ``Plot Over Line``-Filter in **ParaView** and select the ``X Axis`` as line.
The result shoul look like in figure \ref{fig:sod_result}.

![Solution of the Sod shock tube example at $t=0.2$.\label{fig:sod_result}](tutorials/06_sod/sod_paraview_visualization.png)

