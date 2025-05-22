(sec:tut_freestream)=
# Freestream
Unlike the previous tutorial, which dealt with simpler equations, we will now consider the behavior of a compressible fluid flow described by the Navier-Stokes equations. The simplest valid flow solution imaginable is a freestream scenario under conditions of pressure, density, and velocity. For the current setup, we specify a freestream scenario with constant pressure $p=101325.0\ \mathrm{Pa}$, density $\rho=1.225\ \mathrm{kg/m^3}$ and velocity vector $\boldsymbol{U}=(1,1,1)^T\ \mathrm{m/s}$. This configuration provides a baseline for analyzing how the solver handles simple flow conditions and sets the stage for more challenging simulations.

## Mesh Generation
In the tutorial directory, we provide the necessary mesh file, `cartbox_mesh.h5`, along with a parameter file for **HOPR** to generate this mesh. You can recreate the mesh by running the following command.
```bash
hopr parameter_hopr.ini
```

## Build Configuration
**FLEXI** should be compiled with the `freestream` preset using the following commands.
```bash
cmake -B build --preset freestream
cmake --build build
```

## Simulation Parameters
The parameter file to run the simulation is supplied as `parameter_flexi.ini`.  The parameters specific to the Navier-Stokes equation system can be found in the `EQUATION` section of the file.
```ini
! ============================================================ !
! EQUATION
! ============================================================ !
IniExactFunc  = 1
IniRefState   = 1
RefState      = (/1.225,1.0,1.,1.,101325./)
```
The initial condition is set via the variable vector `RefState` which represents the solution vector $(\rho, u, v, w, p)^T$. **FLEXI** permits for multiple `RefState` vectors, allowing each to be referenced by its corresponding cardinal number for its order number within the .ini file. In this example, only a single `RefState` vector is defined, referenced as $1$. Thus, the selection process for the solution vector looks like the following.

- `IniRefState  = 1`: the initial condition uses `RefState 1` for the initial flow field solution.
- `IniExactFunc = 1`: the employed exact function routine uses `RefState 1`, such as for calculating the $L_2$ error norm.

The material properties of the fluid medium, such as the ideal gas constant, are given in {numref}`tbl:freestream_flow_prop` and define the gas behavior in combination with the ideal gas law $p=\rho R T$.
```{list-table} Material properties for the freestream tutorial.
:header-rows: 1
:name: tbl:freestream_flow_prop
:width: 100%
:widths: 30 30 40

* - Variable
  - Value
  - Description
* - mu0
  - 1.8547e-5
  - Dynamic viscosity $\mu$ 
* - R
  - 276
  - Ideal gas constant $R$
* - kappa
  - 1.4
  - Isentropic coefficient $\kappa$
```
The Discontinuous Galerkin (DG) solution is represented by piecewise polynomials on the computational mesh. In this tutorial, the polynomial degree $N$ is chosen as $N=3$. The remaining numerical parameters are outlined in {numref}`tbl:freestream_num_set`.
```{list-table} Numerical settings for the freestream tutorial.
:header-rows: 1
:name: tbl:freestream_num_set
:width: 100%
:widths: 30 30 40

* - Variable
  - Value
  - Description
* - N
  - 3
  - Polynomial degree
* - MeshFile
  - cartbox_mesh.h5
  - Mesh file to be used                     
* - tend
  - 1e-6
  - End time of the simulation             
* - Analyze_dt
  - 1e-6
  - Time interval for analysis             
* - CFLscale
  - 0.99
  - Scaling for the theoretical CFL number
* - DFLscale
  - 0.4
  - Scaling for the theoretical DFL number
```

## Simulation and Results
We proceed by running the code with the following command.
```bash
flexi parameter_flexi.ini
```
Running the code prints all output to `STDOUT`. If the run completes successfully, the last lines should appear similar to the following (condensed) output.
```
==============================================================
 INITIALIZATION DONE! [ 0.01 sec ]
==============================================================
--------------------------------------------------------------
 Sys date  :    06.11.2024 14:08:54
#GridCells :    8.0000000E+00
#DOFs      :    5.1200000E+02
#Procs     :    1.0000000E+00
#DOFs/Proc :    5.1200000E+02
 WRITING INITIAL SOLUTION:
--------------------------------------------------------------
 Initial Timestep  :    1.0000000E-06
--------------------------------------------------------------
 Errors of initial solution:
 Sim time   :    0.000E+00
 L_2        :    2.902E-16   2.902E-16   2.902E-16   2.902E-16   3.815E-11
 L_inf      :    6.661E-16   6.661E-16   6.661E-16   6.661E-16   1.164E-10
--------------------------------------------------------------
   Time = 0.0E+00  dt = 0.1E-05  ETA [d:h:m]   0:00:00:00 |>  | [  0.00%] 
==============================================================
 FLEXI RUNNING cartbox... [ 0.01 sec ] [ 0:00:00:00 ]
==============================================================
 CALCULATION RUNNING...

--------------------------------------------------------------
 Sys date  :    06.11.2024 14:08:54
 CALCULATION TIME PER STAGE/DOF:            [ 7.05792E-07 sec ]
 EFFICIENCY: CALCULATION TIME [s]/[Core-h]: [ 5.90480E-01 sec/h ]
 Timestep   :    1.063E-04
#Timesteps  :    1.000E+00
 Sim time   :    1.000E-06
 L_2        :    2.902E-16   2.521E-15   2.477E-15   2.513E-15   3.815E-11
 L_inf      :    6.661E-16   1.643E-14   1.421E-14   1.443E-14   1.164E-10
--------------------------------------------------------------
   Time = 0.1000E-05  dt = 0.1000E-05  ETA [d:h:m]:<1 min |==>| [100.00%] 
==============================================================
 FLEXI RUNNING cartbox... [ 0.02 sec ] [ 0:00:00:00 ]
==============================================================
==============================================================
 FLEXI FINISHED! [ 0.02 sec ] [ 0:00:00:00 ]
==============================================================
```
```{error}
If the output does not look like the one above, check for any error messages to diagnose the issue.
```

After completing the simulation, examine the contents of the working directory. For a successful run, the directory should contain additional generated files named `<PROJECTNAME>_State_<TIMESTAMP>.h5`. Each file stores the solution vector of the conserved variables at each interpolation node at a specific time, corresponding to multiples of `Analyze_dt`.
```{code-block} bash
freestream
├── cartbox_mesh.h5
├── cartbox_State_0000000.000000000.h5
├── cartbox_State_0000000.000001000.h5
├── parameter_convert.ini
├── parameter_flexi.ini
├── parameter_hopr.ini
├── parameter_postiVisu.ini
```

## Visualization
**FLEXI** relies on [ParaView](https://www.paraview.org) for visualization. In order to visualize the **FLEXI** solution, its format has to be converted from the HDF5 format into another format suitable for **Paraview**. **FLEXI** provides a post-processing tool [posti_visu](tools-visualization) which generates files in VTK format with the following command.
```bash
posti_visu parameter_postiVisu.ini parameter_flexi.ini cartbox_State_0*
```
