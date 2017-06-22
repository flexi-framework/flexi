## Double Mach Reflection

The Double Mach Reflection is a classical test case to investigate the abilities of a numerical scheme to represent shock and contact discontinuities. 
It was invented by Woodward and Colella [@Woodward1984].
A Mach 10 oblique shock wave hits a reflecting wall and the initial conditions are given by the Rankine-Hugoniot conditions
\begin{equation}
(\rho, v_1, v_2, p) =\\ \begin{cases} 
   \left(8.0, 8.25\cdot cos(30^\circ), -8.25 \cdot sin(30^\circ), 116.5 \right) & x < x_0 + \sqrt{\frac{1}{3}} y \\
   \left(1.4, 0.0, 0.0, 1.0\right) & x \ge x_0 + \sqrt{\frac{1}{3}} y 
\end{cases},
\end{equation}
where $x_0= \frac{1}{6}$ is the start of the wall and the computational domain is $\Omega = [0,4] \times [0,1]$, which is discretized by an equidistant cartesian mesh.


Copy the ``dmr`` tutorial folder 

        cp -r $FLEXI_TUTORIALS/dmr .

### Mesh Generation with HOPR

The mesh files used by **FLEXI** are created by supplying an input file *parameter_hopr.ini* with the appropriate information.

    ./hopr parameter_hopr.ini

This creates the mesh file *DMR_mesh.h5* in HDF5 format.

### Flow Simulation with FLEXI

This example requires the Finite Volume shock capturing. Therefore turn the option ``FLEXI_FV`` in the cmake configuration on. Additionally you should switch ``FLEXI_PARABOLIC`` off and recompile the **FLEXI** code.
The simulation setup is defined in *parameter_flexi.ini* and includes options for the Finite Volume shock capturing.  


| Option                        | Value       | Description                                                  |
| ----------------------------- | ----------- | -------------------------------------------------------------|
| IndicatorType                 | Jameson     |                                                              |
| IndVar                        | 6           | Pressure is used for indicator evaluation                    |
| FV_LimiterType                | MinMod      |                                                              |
| FV_IndUpperThreshold          | 0.010       | upper threshold (if IndValue above this value, switch to FV) |
| FV_IndLowerThreshold          | 0.005       | lower threshold (if IndValue below this value, switch to DG) |
| FV_toDG_indicator             | T           | enable an additional Persson indicator to check switch from FV to DG |
| FV_toDG_limit                 | -5.5        | threshold for additional Persson indicator                   |
| FV_IniSupersample             | T           | supersample initial solution for every FV sub-cell, since DG polynomial oscillates for initial solution |

Explanation of the Finite Volume specific options (read also the explanations for the Sod shock tube in chapter \ref{sec:fv_options}):

* ``IndicatorType``: The ``Jameson`` indicator is used to detect troubled cells, that should use the FV operator. This indicator is an adaption of the switching function of the Jameson-Schmidt-Turkel scheme [@Jameson1981] to Finite Volume sub-cells. In contrast to the Persson indicator this indicator is not element local and therefore more robust for travelling discontinuities.
* ``FV_toDG_indicator``: if set to ``T`` (true) an additional Persson indicator for the switch from FV to DG is used. A FV element designated for the switch to DG is converted to DG and the Persson indicator is evaluated for this DG polynomial to test if the polynomial is oscillating. Only if this test is passed the element becomes a DG element. Otherwise it remains a FV element. This helps especially for all indicator functions acting on the FV sub-cells. Here the Jameson indicator only acts on the direct adjacent DOFs of a specific DOF, which does not allow to capture all high frequencies of a polynomial.
* ``FV_toDG_limit``: Threshold for the additional Persson indicator. An element can only switch to FV if the indicator value is below this threshold.
* ``FV_IniSupersample``: If this option is ``F`` (false) the solution is initialized as DG polynomials for all elements. The indicator function is evaluated to find all troubled elements, which are than converted to FV elements. This causes major trouble if discontinuities lay inside DG elements which leads to heavy oscillating polynomials. Converting these oscillations to FV may produce invalid solutions (i.e. negative density). Switching this option on (``T``) enables a super sampling of the initial solution for every FV sub-cell, which removes the mentioned problems with oscillating polynomials. The mean value of every FV sub-cell is computed by evaluating the initial solution in $(N+1)^3$ equidistant points inside the sub-cell and then taking the arithmetic mean value.



The command

~~~~~~~
flexi parameter_flexi.ini 
~~~~~~~

runs the code and generates 11 state files **dmr_State_TIMESTAMP.h5** for $t=0.0, 0.02, \ldots, 0.20$.
To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**. Execute the command 

~~~~~~~
flexi2vtk parameter_flexi.ini dmr_State_*.h5
~~~~~~~
to generate the corresponding *vtu*- and *vtm*-files, which can then be loaded into **ParaView**. 
There are two types of *vtu*-files, which contain either the DG or the FV part of the solution. 
The *vtm*-files combine the DG and FV *vtu*-file of every timestamp. Load the *vtm*-files into **ParaView**.

The result at $t=0.2$ should look like in figure \ref{fig:dmr_result}.

![Distribution of DG and FV elements (top) and density (bottom) of Double Mach Reflection at $t=0.2$.\label{fig:dmr_result}](tutorials/07_dmr/dmr_paraview_visualization.png)

