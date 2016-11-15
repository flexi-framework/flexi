## Linear Scalar Advection-Diffusion Equation

Besides the Navier-Stokes equations, FLEXI provides another equation system, the three-dimensional linear scalar advection-diffusion (LinAdvDiff for short) equation:

$$  \frac{\partial \Phi}{\partial t} + \nabla \cdot (\mathbf{u}\Phi)=d \nabla^2 \Phi, $$

where a scalar solution $\Phi$ is advected with the constant (three-dimensional) velocity $\mathbf{u}$ and is subjected to diffusion 
with a constant scalar diffusion coefficient $d$.

The LinAdvDiff equation is useful to test and develop features in a simple and computationally cheap environment as well as to investigate basic properties of the DG operator, as we will
do in this tutorial.

### Theoretical Background

The dispersion and dissipation properties of the DGSEM operator can be analysed [@gassner2011disp] by looking at the
evolution of solutions to the one-dimensional linear scalar advection equation (without physical dissipation, $d=0$) when a wave with a specific wavenumber is initialised. 
You sould have a look in the paper if you are interested in the details of the analysis, here only a short summary is presented.

Under the mentioned simplifications, the equation to consider simply reads

$$ \frac{\partial \Phi}{\partial t} + u\frac{\partial \Phi}{\partial x}=0 $$

and we consider a wave-like analytical solution on an infinite domain

$$ \Phi(x,t)=e^{i(kx-\omega t)} $$

with the constant scalar transport velocity $u$, the angular frequency $\omega=ku$ and the wavenumber $k$. 

We now assume to have a uniform mesh with a mesh size of $\Delta x$ and seek numerical solutions of the form

$$ \mathbf{\Phi}^l = \hat{\mathbf{\Phi}} e^{i(kl\Delta x-\omega t)}  $$

where $\mathbf{\Phi}^l$ is a vector containing the degrees of freedom in cell $l$ and $\hat{\mathbf{\Phi}}$ is a complex amplitude vector, both of size $N+1$. If we look at the matrix
notation of the semi-discrete system (meaning only spatial discretization) in a single element, we can formulate an algebraic eigenvalue problem

$$ \underline{\underline{A}} \hat{\mathbf{\Phi}} = \Omega \hat{\mathbf{\Phi}},$$

with $\Omega=\frac{\omega \Delta x}{a}$. The matrix $\underline{\underline{A}}$ contains the spatial discretization and depends on the non-dimensional wavenumber $K=k\Delta x$.

By looking at the solutions of this eigenvalue problem, one can obtain relations for the dissipation and dispersion behavior 
inherent of the (spatial) numerical scheme for different wavenumbers. In figure \ref{fig:linadv_dispdiss} we plot these relationships depending on the polynomial degree $N$ using DGSEM with Gauss nodes for the so called physical mode.
This mode is associated with the eigenvalue that follows the exact dispersion relation for the largest range of wavenumbers and also has the biggest influence on the overall numerical solution, at least for rather well-resolved waves.

![Dispersion and dissipation relationship for $N=1-10$ over the modified wavenumber $K^*$. For the dispersion, the dashed line gives the exact relation.\label{fig:linadv_dispdiss}](tutorials/09_linadv/dispdiss.jpg)

The quantities are normalized by the number of grid points, $K^*=\frac{K}{N+1}$ and $\Omega^*=\frac{\Omega}{N+1}$ to give a fair comparison between different polynomial degrees. This means
that if $K^*=\pi$, we have two points per wavelength which is the theoretical minimum to resolve a wave given by the Nyquist-theorem while a normalized wavenumber of $0$ indicates a constant solution.

As we can see, the dissipation properties of higher-order approximations are significantly improved compared to low-order schemes. They are able to preserve waves with higher frequencies without dissipating them. 
We also observe a sharp increase of the dissipation error in the high modes associated with higher-order approximations, which is one of the reasons the timestep is scaled down for larger $N$.

We will now try to show these properties by conducting some numerical experiments using the linear scalar advection equation.

### Compiling FLEXI with LinAdvDiff

If you want to use the LinAdvDiff equations, you need to specify the equation system during the configuration. We are going to create a second ``build``-folder to keep any existing FLEXI binaries that might have been compiled using the Navier-Stokes equation system.

~~~~~~~~~~~~
cd $FLEXI_DIR
mkdir buildLinAdv && cd buildLinAdv 
cmake -DFLEXI_EQNSYSNAME=linearscalaradvection ../
make
~~~~~~~~~~~~

Since we don't want to have diffusion in this tutorial we can either turn of the parabolic terms in general during configuration by using the following ``cmake`` command instead of the previous one:

~~~~~~~~~~~~
cmake -DFLEXI_EQNSYSNAME=linearscalaradvection -DFLEXI_PARABOLIC=OFF ../
~~~~~~~~~~~~

or simply specify a diffusion coefficient of $0$ in the parameter file.

Of course you can also use ``ccmake ../`` and set the parameters in the interface. The ones necessary are listed in the table below.

| Option                          | Value                  | Comment      |
| ------------------------------- |:----------------------:| ------------:|
| FLEXI_EQYNSYSNAME               | linearscalaradvection  |              |
| FLEXI_PARABOLIC                 | OFF                    | optional     |


### Setup

We want to run numerical experiments and solve the one-dimensional linear scalar advection equation. As our initialization condition we choose a single wave with a specific 
angular frequency and observe the behaviour of this wave depending on the normalised non-dimensional wavenumber $K^*$ and the polynomial degree $N$. Since the analysis in the previous section
has been done on a infinite domain, we are going to make our computational domain large enough to neglect influences from the boundaries. An alternative would be to use periodic boundary conditions,
but this would limit our choice of wavelengths since they have to fit in the domain and could also negatively influence the stability of the time discretization. 
When we conduct the experiments, we are going to focus our attention to a single element in the middle of the domain and observe the evolution of the wave in there.

#### Mesh

For this tutorial we are going to create a (quasi-)one-dimensional equidistant grid. We choose $\Delta x=2$ which is equal to the reference element. We are going to make the computational domain rather long, 
so the boundary conditions will have no influence on our solution. Recall that for LinAdv information propagates with velocity $u$, so if our boundaries are a distance of $L$ away from the cell we are looking at,
we can compute solutions up to $t=uL$ without having any influence of the boundary conditions. 

In our case we are creating a cartesian grid with $x \in [-61,61]$ and 61 cells to reach our desired $\Delta x$. We will later focus on the single element in $x \in [-1,1]$. To achieve a (quasi-)one-dimensional simulation,
 we choose periodic boundary conditions in $y-$ and $z-$direction. The boundary conditions in the $x$-direction will not matter, so we simply use Dirichlet-type BCs (BC type 2) with the analytical wave function 
as the boundary state. We specify this by setting the BC state to zero, which means that the initialization function will be used instead of a seperate function.

We provide a parameter file for HOPR to generate this mesh. It is located in the subdirectory ``mesh`` within the directory that belongs to this tutorial. Either run ``hopr`` with this parameter file

~~~~~~~~~~~~
cd $FLEXI_TUTORIALS_DIR/linadv/mesh
hopr parameter_hopr.ini
~~~~~~~~~~~~

or use the mesh file ``CART_1D_mesh.h5`` which is provided in this directory.

#### Flexi parameters

A ``parameter_flexi.ini`` is provided in the directory ``$FLEXI_TUTORIALS_DIR/linadv/``. We are going to discuss the parameters that are specific to the LinAdvDiff equation system and the ones needed for this tutorial. 
In the section ``EQUATION`` in the parameter file you will find the options that directly influence the behaviour of the equation:

~~~~~~~~~~~~
! ============================================================ !
! EQUATION
! ============================================================ !
AdvVel        = (/1.,0.,0./)
DiffC         = 0.
~~~~~~~~~~~~

The parameter ``AdvVel`` sets the advection velocity in all three space directions. Since we are simulating a one-dimensional case, we only set the first entry to something other than zero.
Here we choose one, although the special initialization that we will be using will set the advection velocity to this value no matter what we put in here.
Since we don't want any physical diffusion, we set the diffusion coefficient ``DiffC`` to zero. If you did compile without the ``FLEXI_PARABOLIC`` option set to ``ON``, you don't need to set this.

We also need to set the initial condition and we choose a special one here to match the initial conditions of the analysis we want to compare against. Recall that the exact solution was 

$$ \Phi(x,t)=e^{i(kx-\omega t)}, $$

where the amplitude of the wave is given by the real part of this complex expression. This real part is then given by

$$ A(x,t) = cos(kx-\omega t).$$

This function is implemented as ``ExactFunc 6`` in the LinAdvDiff equation system, so we are specifying this in the parameter file:

~~~~~~~~~~~~
IniExactFunc  = 6
~~~~~~~~~~~~

What is left is to specify the angular frequency $\omega$, since the wavenumber is given by $k=\frac{\omega}{u}$. The angular frequency can be set during runtime using the parameter

~~~~~~~~~~~~
OmegaRef      = 2.
~~~~~~~~~~~~

We will adopt this parameter when we run our simulations.

To reduce the possibility that the error introduced by the time discretization will alter our results (we are focused here on the behaviour of the spacial discretization), we set our CFL number to a rather small value:

~~~~~~~~~~~~
CFLscale      = 0.1
~~~~~~~~~~~~

Since this is a really small and fast computation, we are also going to use visualization routines during runtime, so we don't need to convert the state files using the ``flexi2vtk`` program. To do this, we set the following parameters in
our parameter file:

~~~~~~~~~~~~
outputFormat  = 3
...
NVisu         = 30
~~~~~~~~~~~~

which will enable the ``.vtu`` output. Whenever the analyze routines are called, there will also be a call to the visualization routines. Since we did not set a ``Analyze_dt`` in the parameter file, the analyse routines will only
be called at the beginning and the end of the simulation. The parameter ``NVisu`` sets the polynomial degree of the output basis.

**Remark**: There is also a parameter that allows us to only call the visualize routines every $n$-th call of the analyze functions. This parameter is called ``nWriteData``. So if you set 

~~~~~~~~~~~~
nWriteData  = 10
~~~~~~~~~~~~

in the parameter file, the visualization file will only be written every 10th analyze routine call. The default for this value is simply ``1``.

### Simulation and results

We are going to start with the simulation of a well-resolved wave and a moderate polynomial degree of $N=4$ which is already set in the parameter file. For a modified wavenumer of $K^*=\frac{1}{4}\pi$, the dissipation and dispersion relations tell us that we
can expect very little errors.
To get the desired modified wavenumber, we need to calculate the corresponding angular frequency. Using the relations given above, the angular frequency can be computed as 

$$ \omega = \frac{K^* (N+1)}{\Delta x} u.$$

This means we need to set the following frequency in the parameter file:

~~~~~~~~~~~~
OmegaRef      = 1.96349540849
~~~~~~~~~~~~

All other options are already set correctly. The simulation will be performed up to a time of $t=5$. We can now run the simulation using

~~~~~~~~~~~~
cd $FLEXI_TUTORIALS_DIR/linadv
$FLEXI_DIR/buildLinAdv/bin/flexi parameter_flexi.ini
~~~~~~~~~~~~

Even on a single processor this calculation should only take a few seconds. 

Now we can have a look at the result by opening the file ``LinAdvCosineWave_Solution_0000005.000000000.vtu`` with ParaView. We are only interested in the central element in $x \in [-1,1]$, so we can clip everything else.
Also since we are performing a one-dimensional simulation, we can extract the solution along the $x$-axis and create a simple line-plot to look at the results. 
To see how dispersion and dissipation introduced by the numerics influence our results, we will also plot the analytical solution to the wave transportation problem, given by the equation above. The result can be seen in figure \ref{fig:N4_wellresolved}. 
As was expected, only very little deviation can be found from the analytical solution, since we are considering a very well resolved wave with a low value for $K^*$.

![Result for $N=4$ and a well resolved wave with $K^*=\frac{1}{4}\pi$ at $t=5$. Comparison against exact wave transport.\label{fig:N4_wellresolved}](tutorials/09_linadv/N4wellresolved.png)

Of course it is much more interesting to look at the behaviour of waves that are not well resolved and how the polynomial degree influences the behaviour of these waves. To do this, we are going to perform four simulations with polynomial degrees ranging from
2 to 11 and set the normalised wavenumber to $K^*=1.6$. Table \ref{tab:omegaRef} will give you an overview of the polynomial degrees and the angular frequencies needed to achieve a constant $K^*$.

| N    | OmegaRef   |
| ---- |:-----------|
| 2    | 2.4        |
| 4    | 4.0        |
| 6    | 5.6        |
| 11   | 9.6        |

Table: Angular frequencies needed to get $K^*=1.6$. \label{tab:omegaRef}

Now run these four simulations by adjusting the values for ``N`` and ``omegaRef`` in the parameter file before each run. Remember to give each of them a unique ``ProjectName`` or else the results will be overwritten each time. 
Also make sure that ``NVisu`` is set to at least three times the value of $N$ to get a meaningful visualization.

We now compare the results as they can be found in the next figure. Again we plot the solution at $t=5$ for the considered polynomial degrees against the analytical solution.
Since we compare using a constant normalized wavenumber, the actual frequency of the wave is increased with the polynomial degree. If we take a look at the result for $N=2$ (which should not be considered high-order), we observe both a significant amount of dissipation 
(a drop in amplitude by about 85%) as well as severe dispersion, e.g. a change in angular frequency and phase angle.  

\newpage

 For a polynomial degree of $N=4$ (at the lower end of what can be considered high-order), we already observe some improvements.
The amplitude only drops by around 65%, although the phase shift is still clearly visible. This trend continues with higher polynomial degrees. For the highest value of $N$ tested here, we observe only small changes in phase angle or amplitude.

These results show how higher-order schemes are able to capture waves with fewer points per wavelength than the lower-order approximations. This is one of the central aspects why we use such schemes.

![](tutorials/09_linadv/comparisonN_1.png)  
![](tutorials/09_linadv/comparisonN_2.png)  
Figure: Result for  a under-resolved wave with $K^*=1.6$ at $t=5$ for $N=2$ (top left), $N=4$ (top right), $N=6$ (bottom left) and $N=11$ (bottom right).\label{fig:Results_underresolved}

