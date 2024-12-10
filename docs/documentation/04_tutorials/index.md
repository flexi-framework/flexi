(Tutorials)=
# Tutorials

This chapter provides a detailed overview of flow simulations with **FLEXI**, assuming familiarity with setting compiler options and code compilation. The path to all executables is omitted here. We assume you have either symlinked **flexi**, **hopr**, and all posti tools into the runtime directory or call these executables at their relative location.

Each tutorial directory contains the necessary .ini files - `parameter_hopr.ini`, `parameter_flexi.ini`, `parameter_postiVisu.ini` - as well as the mesh file `*_mesh.h5` in HDF5 format (generated with **HOPR**).
```{code-block} bash
---
caption: Directory tree for a tutorial.
---
tutorial
├── mesh.h5
├── parameter_flexi.ini
├── parameter_hopr.ini
└── parameter_postiVisu.ini
```

```{tip}
While each tutorial can be run directly in its own directory, we recommend copying each folder to a new directory. This way, you can run the simulations and freely modify the .ini files without altering the original setup.
```

```{toctree}
---
maxdepth: 1
caption: Table of Contents
---
01_linadv/index.md
02_freestream/index.md
03_convtest/index.md
04_cavity/index.md
05_tgv/index.md
06_sod/index.md
07_dmr/index.md
08_channel/index.md
09_cylinder/index.md
10_naca0012/index.md
```
