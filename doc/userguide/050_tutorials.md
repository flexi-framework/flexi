\hypertarget{tutorials}{}

# Tutorials \label{chap:tutorials}

This chapter will give a detailed overview of flow simulations with **FLEXI**. It assumes that you are familiar with how to set the compiler options and how to compile the code. It also assumes that you either have set the environment variables for **flexi** and **flexi2vtk** as described in \ref{sec:compilingthecode} _or_ that you copy or link these binaries to the folder you are working from. 

Each tutorial is equipped with .ini files *parameter_hopr.ini* and *parameter_flexi.ini* as well as mesh file *\*\_mesh.h5* in HDF5 format (created with **HOPR**).

~~~~~~
parameter_hopr.ini
parameter_flexi.ini
mesh.h5
~~~~~~
 
We suggest to copy each folder to a new directory, where you can run and modify the "INI-files". Also, link the **flexi** and **flexi2vtk** executables into that folder, if you have not set the appropriate environment variables. 
