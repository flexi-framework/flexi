\hypertarget{toolsoverview}{}

# Tools Overview \label{chap:toolsoverview}

This section gives an overview over the tools contained in the `FLEXI` repository. It also provides references to the tutorials where their usage is explained. There are two different kinds of tools: 

* `POSTI`-tools can be compiled together with `FELXI` given the according `cmake` options.
* In the `tools` folder, a collection of mostly shell or python scripts can be found. They are mostly used to manage `FLEXI` runs and `FLEXI` output files. 

## POSTI 

`POSTI` tools are mostly documented in section \ref{sec:postiRecordpoints} and \ref{sec:postiVisu}. Here, an overview is given together with references to the respective tutorials. 

### Visualization

---------------------------------------------------------------------------------------------
**posti_visu**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------
                                            
---------------------------------------------------------------------------------------------
**Paraview plugin**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------

### Swap meshes

---------------------------------------------------------------------------------------------
**posti_swapmesh**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------

### Record points

----------------------------------------------------------------------------------------------
**posti_evaluaterecordpoints**
--------------------------------- -------------------------------------------------------------
Brief description                  No description so far 

Basic usage                        `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in               No tutorials so far
-----------------------------------------------------------------------------------------------
                                            
----------------------------------------------------------------------------------------------
**posti_preparerecordpoints**
-------------------------------- -------------------------------------------------------------
Brief description                 No description so far 

Basic usage                       `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in              No tutorials so far
----------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------
**posti_visualizerecordpoints**
--------------------------------- ------------------------------------------------------------
Brief description                  No description so far 

Basic usage                        `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in               No tutorials so far
----------------------------------------------------------------------------------------------
                                            

### Time averaging

---------------------------------------------------------------------------------------------
**posti_mergetimeaverages**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------
                                            
---------------------------------------------------------------------------------------------
**posti_calcfluctuations**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------
                                            

### Test cases

---------------------------------------------------------------------------------------------
**posti_channel_fft**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------
**posti_init_hit**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------




<!---
==================================================================================================================================
-->

## `tools` folder 

The scripts provided in the `tools` folder are generally not part of the tutorials.
They are briefly described below.

For most python tools, possible arguments and syntax can be shown with the `-h` argument: 

~~~~~~~
python2 [toolname.py] -h
~~~~~~~

### `animate`

The python script **animate_paraview.py** creates movies from a series of state files using `PvBatch`, a GUI-less interface to ParaView. 
You need ParaView installed on your system (details can be found in the ParaView <!---TODO--> section) and the directory containing the PvBatch executable needs to be a part of your `$PATH`. Before running this script, you have to visualize your `FLEXI` state file with ParaView and save the current view via `Save State...`, e.g. under the name `pvstate.pvsm`. You also need the `mencoder` tool installed. The basic command to run the script is 

~~~~~~~
python2 animate_paraview.py -l [pvstate.pvsm] -r [path_to_posti_paraview_plugin] [statefile1.h5 statefile2.h5 ...]
~~~~~~~

Apart from the movie file, the script also outputs a `.png`-file for each HDF5 file given as input. 
In order to visualize a set of `.vtu`-files, e.g. from `POSTI_VISU` output, omit the `-r` argument and pass `.vtu`-files instead of `.h5`-files. Further options can be shown with the `-h` argument. 

There are further tools for image handling in this folder: 

The tool **concatenatepics.py** stitches several pairs of images (e.g. creates a time series of stitched images from two time series of images). A possible command could look like this (*Further options can be shown with the `-h` argument*):

~~~~~~~
python concatenatepics.py -d e -p left*.png  -a right*.png
~~~~~~~


The tool **crop.py** crops several images to the same size. Simply pass all images as arguments: 

~~~~~~~
python crop.py [image*.png]
~~~~~~~

The script **pics2movie.py** creates a movie from several images using the `mencoder` tool (which is also done as part of the `animate_paraview.py` script. Basic usage is again 

~~~~~~~
python pics2movie.py [image*.png]
~~~~~~~

and further options can again be shown with the `-h` argument. 


<!---
..................................................................................................................................
-->

### `convergence_test`

The python scripts `convergence` and `convergence_grid` provide automated convergence tests for p- and h-convergence, respectively. 
The basic command is

~~~~~~~
python2 convergence $FLEXI_DIR/flexi [parameter.ini]
~~~~~~~

where `convergence` can be replaced by `convergence_grid` for h-convergence. Further options can again be shown with the `-h` option. 

Note that for h convergence, the mesh names are hard-coded to the form `CART_HEX_PERIODIC_MORTAR_XXX_2D_mesh.h5`, where `XXX` is the number of elements in each direction, and `MORTAR` and `2D` are optional. The polynomial degree in the parameter file is *always* overwritten by that passed to the script as an optional argument, with a default value of 3, if no such argument is passed. <!---TODO: change this in the script--> 


<!---
..................................................................................................................................
-->

### `userblock`

The `userblock` contains complete information about a `FLEXI` run (git branch of the repository, differences to that branch, `cmake` configuration and parameter file) and is prepended to every `.h5` state file. The parameter file is prepended in ASCII format, the rest is binary and is generated automatically during the build process with the `generateuserblock.sh` script. It can be extracted and printed using the `extract_userblock.py` script. Its basic usage is 

~~~~~~~
python2 extract_userblock.py -XXX [statefile.h5]
~~~~~~~

where `-XXX` can be replaced by 

* `-s` to show all available parts of the userblock (such as `CMAKE` or `GIT BRANCH`)
* `-a` to print the complete userblock
* `-p [part]` to print one of the parts listed with the `-s` command.

The second python tool in this folder is `rebuild.py`. It extracts the userblock from a state file and builds a `FLEXI` repository and binary identical to the one that the state file was created with. In order to do so, it clones a `FLEXI` git repository, checks out the given branch, applies the stored changes to the git `HEAD` and builds `FLEXI` with the stored `cmake` options. If run with the parameter file given in the `INIFILE` part of the userblock, this binary should reproduce the same results/behaviour (possible remaining sources of different output are for example differences in restart files, compilers, linked libraries or machines). The basic usage is 

~~~~~~~
python2 rebuild.py [dir] [statefile.h5]
~~~~~~~

where `dir` is an empty directory that the repository is cloned into and where the `FLEXI` executable is built. `statefile.h5` is the state file whose userblock is used to rebuild the `FLEXI` executable. Help can be shown via `-h` for both userblock scripts.


<!---
..................................................................................................................................
-->

### `sortfiles.sh`

The `sortfiles.sh` script sorts all `.h5`-files in subfolders `State`, `BaseFlow`, `TimeAvg` and `RP`, while keeping the last time instance at the upper level. It also copies `Log.*.sdb`, `.log` and `.out` files into a `logs` subdirectory. The project name is hard-coded in the script and has to be adapted there, the directory that is to be sorted is passed as an argument.

### `getload.py`

This script is specific to runs on HPC systems. It calculates a suitable number of nodes and cores to achieve 

* a number of degrees of freedom per core which is close to a target
* an average number of elements per core which is just below a close integer, such that parallel efficiency is not impaired by a few cores with higher load that the others have to wait for.

No arguments are passed to this script, all input values are hard-coded and have to be adjusted in the script.
