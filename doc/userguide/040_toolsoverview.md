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
They are breifly described below.

### animate

The python tool `animate_paraview.py` creates movies from a series of state files using a paraview 

