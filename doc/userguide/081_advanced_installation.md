\hypertarget{advancedinstallation}{}

# Advanced installation using Environment Modules \label{chap:modules}

If multiple compilers and thus versions of libraries and tools
are in use the management of all these dependencies can quickly
become complicated and time consuming. Environment modules are
frequently used to manage such setups, but are themselves
hard to set up. Flexi therefore comes with an own set of
environment modules and script to automatically build, update
and generate modules for dependencies such as MPI, HDF5 or
Paraview. This is advantageous, as most tools of the Flexi
framework (HOPR, postprocessing tools) require the same set
of libraries, which only have to be built once.

The `envbuild.sh` script in the `tools/envbuild` folder will
ease the installation of Flexi and build other related tools
including their environment modules. It is tested with Ubuntu
LTS 16.04.

## Environment Modules

First the module environment needs to be installed and populated
with some basic modules, the other modules will later be
generated automatically by the `envbuild` script.
  
Install the packages for environment modules, e.g.

       sudo apt-get install environment-modules (Ubuntu)
       sudo zypper install Modules              (OpenSUSE)

You need to specify the folder, where the modules will
be placed in the file
    
       /etc/environment-modules/modulespath        (Ubuntu)
       /usr/share/Modules/3.2.10/init/.moduleshome (OpenSUSE)

Here add the desired module location to the **first** position 
in the list (e.g. `/opt/Modules`) or use one of the predefined
locations. The environment variable $MODULEPATH should
contain all the selected paths and is set by the init command
below. The `envbuild.sh` script will install its module files
(not the actual libraries) in the **first** path found in this
variable.

Copy the contents of the `tools/envbuild/Modules` folder to the
path you specified (`/opt/Modules`). It already comes with few
modules for a set of common tools. The only required modules
are `env` and `compiler`.
Feel free to remove the other modules if you don't need them.

For non-ssh sessions (non-login shells, i.e. when using a GUI
on a local desktop) the module environment is often not loaded
by default. In this case you need to load it by placing these
commands in your .bashrc (usuming bash is used)

       . /usr/share/modules/init/bash           (Ubuntu)
       . /usr/share/Modules/3.2.10/init/bash    (OpenSUSE)
                                                           
If these paths are not present the environment variable
$MODULESHOME should point to the correct path with the module
files. The variable $MODULEPATH may only be available after
initializing the environment as described above.

If the configuration was successfull the modules should show up
by typing `module avail`.

## Using the envbuild script

Edit the script file `envbuild.sh` and follow the instructions
in the file header, where you can select the libraries and tools
to be built. Then run the script using 

       bash envbuild.sh

It will download and install all required libraries for the
Flexi framework and and some additional post-processing tools.
It will also create module files and place them in the module
environment. Note that depending on the path of the environment
modules and the installation path of the libraries, root
permissions may be required for the script to run.

Note that the script will not install the prerequisites required
for building the tools and libs, which is subject to the user.
It is intended to allow a more simple organisation of multiple
library versions for different compilers and simple updating and
compiling new versions of some libraries.

The main modules for different compilers are the env modules.
Using the command

       module load env/gnu

all relevant modules for the GNU compiler are loaded. Switching
to another compiler environment then can be performed by
   
       module switch env env/intel
