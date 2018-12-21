\hypertarget{toolsoverview}{}

# Tools Overview \label{chap:toolsoverview}

This section gives an overview over the tools contained in the flexi repository. It also provides references to the tutorials where their usage is explained. There are two different kinds of tools: 

* `POSTI`-tools can be compiled together with `FELXI` given the according cmake options.
* In the `tools` folder, a collection of mostly shell or python scripts can be found. They are mostly used to manage `FLEXI` runs and `FLEXI` output files. 

## POSTI 

`POSTI` tools are mostly documented in the tutorials in section \ref{sec:tutorials}. Here, an overview is given together with references to the respective tutorials. 

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

---------------------------------------------------------------------------------------------
**posti_evaluaterecordpoints**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------
                                            
---------------------------------------------------------------------------------------------
**posti_preparerecordpoints**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------------
**posti_visualizerecordpoints**
--------------------------- -----------------------------------------------------------------
Brief description            No description so far 

Basic usage                  `posti_visu [parameter.ini] [statefilde.h5]`

Further explained in         No tutorials so far
---------------------------------------------------------------------------------------------
                                            

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




## `tools` folder 

The scripts provided in the `tools` folder are generally not part of the tutorials.
They are briefly described below.

### `animate`

The python script **animate_paraview.py** creates movies from a series of state files using `PvBatch`, a GUI-less interface to Paraview. 
You need ParaView installed on your system (details can be found in the ParaView <!---TODO--> section) and the directory containing the PvBatch executable needs to be a part of your `$PATH$`. Before running this script, you have to visualize your FLEXI state file with paraview and save the current view via `Save State...`, e.g. under the name `pvstate.pvsm`. You also need the mencoder tool installed. The basic command to run the script is 

~~~~~~~
python2 animate_paraview.py -l [pvstate.pvsm] -r [path_to_posti_paraview_plugin] [statefile1.h5 statefile2.h5 ...]
~~~~~~~

Apart from the movie file, the script also outputs a `.png`-file for each HDF5 file given as input. 
Further options can be shown via 

~~~~~~~
python2 animate_paraview.py -h
~~~~~~~

There are further tools for image handling in this folder: 

The tool **concatenatepics.py** stitches several pairs of images (e.g. creates a time series of stitched images from two time series of images). A possible command could look like this:

~~~~~~~
python concatenatepics.py -d e -p left*.png  -a right*.png
~~~~~~~

Further options can be shown via 

~~~~~~~
python concatenatepics.py -h
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

### `convergence_test`



### `userblock`

### `sortfiles.sh`

### `getload.py`
