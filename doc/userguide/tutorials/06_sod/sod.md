## SOD Shock tube

The Sod shock tube example [@Sod1978] is one of the most basic test cases to investigate the shock capturing capabilities of a CFD code. A initial discontinuity is located in the middle of the one dimensional domain $[0,1]$. The left and right states are given by $\rho=1, v=0, p=1$ and $\rho=0.125, v=0, p=0.1$.
These states are already set as ``RefState`` in the the *parameter_flexi.ini* file.

Copy the ``sod`` tutorial folder 

        cp -r $FLEXI_TUTORIALS/sod .

### Mesh Generation with HOPR

The mesh files used by **FLEXI** are created by supplying an input file *parameter_hopr.ini* with the appropriate information.

    ./hopr parameter_hopr.ini

This creates the mesh file *SOD_100_mesh.h5* in HDF5 format.

### Flow Simulation with FLEXI

This example requires the Finite Volume shock capturing and the Euler equations. Therefore set the following options in cmake:

~~~~~~
FLEXI_FV          ON 
FLEXI_PARABOLIC   OFF
~~~~~~

and recompile the **FLEXI** code.
The simulation setup is defined in *parameter_flexi.ini* and includes options for the Finite Volume shock capturing.  


| Option                        | Value       | Description                                                  |
| ----------------------------- | ----------- | -------------------------------------------------------------|
| IndicatorType                 | Persson     |                                                              |
| IndVar                        | 1           | first conservative (density) used for indicator evaluation   |
| IndStartTime                  | 0.001       | until this time FV is used in the whole domain               |
| FV_LimiterType                | MinMod      |                                                              |
| FV_IndUpperThreshold          | -3.         | upper threshold (if IndValue above this value, switch to FV) |
| FV_IndLowerThreshold          | -4.         | lower threshold (if IndValue below this value, switch to DG) |

Explanation of the Finite Volume specific options:
\label{sec:fv_options}

* ``IndicatorType``: specifies the indicator function that is used to detect DG elements that contain a discontinuity. The ``Persson`` indicator [@Persson06shock] is an element local indicator, which compares the different modes of the DG polynomial. If the amount of solution in the highest mode compared to the amount in the lower modes is high the DG polynomial may oscillate. All indicator functions return a high value for trouble elements and a low value for smooth elements.
* ``IndVar``: Variable that is used to evaluate the indicator function. Here the density is used. In general the pressure (6) is a good choice.
* ``IndStartTime``: Until this time the actual indicator function is overwritten by a very high value to force the use of FV elements in the hole domain. This is necessary if initial discontinuities (like in this example) are placed perfectly at the element boundaries. In this case an element local indicator (like Persson) can not detect the discontinuity.
* ``FV_LimiterType``: Limiter of the second order reconstruction. 
* ``FV_IndUpperThreshold`` and ``FV_IndLowerThreshold``: These two threshold values are used to decide in which elements the DG method and where the FV sub-cell scheme should be used. In general a single threshold could be sufficient. If the indicator value raises above this threshold the element used the FV scheme, and if the value is below the DG scheme is used. Practice shows that this could lead to an ongoing switching between the two schemes. Therefore a DG element is switched only to FV if the indicator value gets greater than the upper threshold, but it does not switch back to DG immediately if the value falls below this threshold. This only happens if the indicator value falls below the lower threshold. 

The command

~~~~~~~
flexi parameter_flexi.ini 
~~~~~~~

runs the code and generates 5 state files **sod_State_TIMESTAMP.h5** for $t=0.0, 0.05, 0.10, 0.15, 0.20$.
To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**. Execute the command 

~~~~~~~
posti_visu parameter_postiVisu.ini parameter_flexi.ini sod_State_0000000.*
~~~~~~~

to generate the corresponding *vtu*- and *vtm*-files, which can then be loaded into **ParaView**. 
There are two types of *vtu*-files, which contain either the DG or the FV part of the solution. 
The *vtm*-files combine the DG and FV *vtu*-file of every timestamp. Load the *vtm*-files into **ParaView**.

Since this example is a 1D test case use the ``Plot Over Line``-Filter in **ParaView** and select the ``X Axis`` as line.
The result should look like in figure \ref{fig:sod_result}.

![Solution of the Sod shock tube example at $t=0.2$.\label{fig:sod_result}](tutorials/06_sod/sod_paraview_visualization.png)

