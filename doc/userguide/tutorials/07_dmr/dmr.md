## Double Mach Reflection
\label{sec:tut_dmr}

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

    hopr parameter_hopr.ini

This creates the mesh file *DMR_mesh.h5* in HDF5 format.

### Flow Simulation with FLEXI
\label{sec:tut_dmr_simulation}

This example requires the finite volume (FV) shock capturing.
Two variants of the shock capturing are implemented in **FLEXI**, which are both based on the FV sub-cell technique.
For this, the DG element is subdivided into finite volume sub-cells, where each sub-cell corresponds to one degree of freedom of the initial DG element.
This allows to keep the total number of degrees of freedom in the simulation constant and to switch between the DG and the FV operator for each individual element.
The first approach switches DG elements into the FV sub-cell representation such that each element is either represented by its DG solution and DG discretization operator or the element is represented by FV sub-cells, i.e. a piecewise constant representation, as detailed in [@sonntag2017efficient].
The alternative is to blend the FV operator to the DG operator, while the solution itself always remains a DG polynomial.
Please see [@hennemann2021provably] for more details.
In the following, both approaches will be applied to simulate the Double Mach Reflection.

#### Using the Finite Volume Switching
First, set the option ``FLEXI_FV=SWITCH`` in the cmake configuration to enable the switching-based shock capturing. Additionally you should switch ``FLEXI_PARABOLIC`` off and recompile the **FLEXI** code.
The simulation setup is defined in *parameter_flexi_switch.ini* and includes options for the shock capturing via switching.


| Option                        | Value       | Description                                                  |
| ----------------------------- | ----------- | -------------------------------------------------------------|
| FV_LimiterType                | MinMod      |                                                              |
| IndicatorType                 | Jameson     |                                                              |
| IndVar                        | 6           | Pressure is used for indicator evaluation                    |
| FV_IndUpperThreshold          | 0.010       | upper threshold (if IndValue above this value, switch to FV) |
| FV_IndLowerThreshold          | 0.005       | lower threshold (if IndValue below this value, switch to DG) |
| FV_toDG_indicator             | T           | enable an additional Persson indicator to check switch from FV to DG |
| FV_toDG_limit                 | -5.5        | threshold for additional Persson indicator                   |
| FV_IniSupersample             | T           | supersample initial solution for every FV sub-cell, since DG polynomial oscillates for initial solution |

Explanation of the Finite Volume specific options (read also the explanations for the Sod shock tube in chapter \ref{sec:fv_options}):

* ``IndicatorType``: The ``Jameson`` indicator is used to detect troubled cells, that should use the FV operator. This indicator is an adaption of the switching function of the Jameson-Schmidt-Turkel scheme [@Jameson1981] to Finite Volume sub-cells. In contrast to the Persson indicator this indicator is not element local and therefore more robust for travelling discontinuities.
* ``FV_toDG_indicator``: if set to ``T`` (true) an additional Persson indicator for the switch from FV to DG is used. A FV element designated for the switch to DG is converted to DG and the Persson indicator is evaluated for this DG polynomial to test if the polynomial is oscillating. Only if this test is passed the element becomes a DG element. Otherwise it remains a FV element. This helps especially for all indicator functions acting on the FV sub-cells. Here the Jameson indicator only acts on the direct adjacent DOFs of a specific DOF, which does not allow to capture all high frequencies of a polynomial.
* ``FV_toDG_limit``: Threshold for the additional Persson indicator. An element can only switch back to DG from FV if the indicator value on the DG representation is below this threshold.
* ``FV_IniSupersample``: If this option is ``F`` (false) the solution is initialized as DG polynomials for all elements. The indicator function is evaluated to find all troubled elements, which are than converted to FV elements. This causes major trouble if discontinuities lay inside DG elements which leads to heavy oscillating polynomials. Converting these oscillations to FV may produce invalid solutions (i.e. negative density). Switching this option on (``T``) enables a super sampling of the initial solution for every FV sub-cell, which removes the mentioned problems with oscillating polynomials. The mean value of every FV sub-cell is computed by evaluating the initial solution in $(N+1)^3$ equidistant points inside the sub-cell and then taking the arithmetic mean value.



The command

~~~~~~~
flexi parameter_flexi_switch.ini
~~~~~~~

runs the code and generates 11 state files **dmr_SWITCH_State_TIMESTAMP.h5** for $t=0.0, 0.02, \ldots, 0.20$.
To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**. Execute the command 

~~~~~~~
posti_visu parameter_postiVisu.ini parameter_flexi_switch.ini dmr_SWITCH_State_0000000.0*
~~~~~~~
to generate the corresponding *vtu*- and *vtm*-files, which can then be loaded into **ParaView**. 
There are two types of *vtu*-files, which contain either the DG or the FV part of the solution. 
The *vtm*-files combine the DG and FV *vtu*-file of every timestamp. Load the *vtm*-files into **ParaView**.

The result at $t=0.2$ should look like in figure \ref{fig:dmr_result_switch}.

![Distribution of DG and FV elements (top) and density (bottom) of Double Mach Reflection at $t=0.2$.\label{fig:dmr_result_switch}](tutorials/07_dmr/dmr_paraview_visualization.png)

#### Using the Finite Volume Blending
Next, we will investigate the blending approach.
For this, set the option ``FLEXI_FV=BLEND`` in the cmake configuration to enable the blending-based shock capturing.
Moreover, the FV blending requires to set ``FLEXI_NODETYPE=GAUSS-LOBATTO``.
Also enable the split-form DG with ``FLEXI_SPLIT_DG=ON`` and recompile **FLEXI**.
More details on the split formulation are given later in Section \ref{sec:tut_ptcf}.
In this tutorial we use an entropy-stable split formulation to ensure the stability of the FV blending approach.
For the FV sub-cell blending, the elements are not switched completely to the FV operator, but instead the DG operator $R_{DG}$ and FV operator $R_{FV}$ are blended as
$$ R = \alpha R_{FV} + (1-\alpha) R_{DG}$$
with the blending coefficient $\alpha$.
Instead of switching between a DG and a FV discretization, the blending allows a continuous transition between the DG and FV operators.
The blending factor is computed based on the indicator proposed by [@hennemann2021provably], which is parameter-free and does not require any parameters to be tuned by the user.

There are some specific options to for the FV blending that can be set in the parameter file.

| Option                        | Value       | Description                                                  |
| ----------------------------- | ----------- | -------------------------------------------------------------|
| FV_LimiterType                | MinMod      |                                                              |
| FV_alpha_min                  | 0.01        | All elements with a blending factor 'alpha' below this threshold are treated as pure DG elements. This improves the computational efficiency, since the FV discretization operator does not have to be computed for these element.                                                             |
| FV_alpha_max                  | 0.5         | Maximum blending factor 'alpha'. All blending coefficients exceeding this value are clipped to FV_alpha_max |
| FV_doExtendAlpha              | T           | If true, the blending factors are propagated to the neighboring elements |
| FV_alpha_ExtScale             | 0.5         | Scaling factor by which the blending factor is scaled when propagated to its neighbor. Has to be between ``0`` and ``1`` |
| FV_nExtendAlpha               | 1           | Number of times this propagation of the blending should be performed. Higher values correspond to a wider sphere of influence |

Similarly to the switching-based approach, the command

~~~~~~~
flexi parameter_flexi_blend.ini
~~~~~~~

runs the code and generates 11 state files **dmr_BLEND_State_TIMESTAMP.h5** for $t=0.0, 0.02, \ldots, 0.20$.
To visualize the solution, execute the command

~~~~~~~
posti_visu parameter_postiVisu.ini parameter_flexi_blend.ini dmr_BLEND_State_0000000.0*
~~~~~~~
to generate the corresponding *vtu*- and *vtm*-files, which can then be loaded into **ParaView**.

The result at $t=0.2$ should look like in figure \ref{fig:dmr_result_blend}.

![Distribution of the blending factor between the DG and FV operators ``FV_alpha`` (top) and density (bottom) of Double Mach Reflection at $t=0.2$.\label{fig:dmr_result_blend}](tutorials/07_dmr/dmr_paraview_visualization_blend.png)
