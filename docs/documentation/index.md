% ```{include} ../../README.md
% ---
% relative-docs: docs/
% relative-images:
% ---
% ```

% to prevent the landing page from being labeled as "<no title> - FLEXI Documentation" (displayed in the browser tab, for example),
% we provide a title, but hide it (suggested workaround in the related feature request: https://github.com/sphinx-doc/sphinx/issues/8356)

<div style="visibility: hidden;">

FLEXI
========

</div>

```{image} ./figures/flexi_logo.jpg
:alt: logo
:width: 600px
:align: center
```

[FLEXI][flexi] is a high-order numerical framework for solving PDEs, with a special focus on Computational Fluid Dynamics. [FLEXI][flexi] is based on the Discontinuous Galerkin Spectral Element Method (DGSEM), which allows for high-order of accuracy and fully unstructured hexahedral meshes. The solver is parallelized very efficiently for large-scale applications and scales to 500,000+ cores. Moreover, [FLEXI][flexi] comes with a capable pre- and post-processing suite that enables complex simulation setups up to the finished visualization.

[FLEXI][flexi] has been developed by the [Numerics Research Group (NRG)][nrg] founded by Prof. Claus-Dieter Munz and currently lead by Prof. Andrea Beck at the Institute of Aerodynamics and Gas Dynamics at the University of Stuttgart, Germany.

You can find detailed installation instructions, the extensive documentation and several tutorial cases for FLEXI [here][flexi].

FLEXI is Copyright (C) 2016, Prof. Claus-Dieter Munz and is released under the **GNU General Public License v3.0**. For the full license terms see the included license file.

Numerous people have worked on and with FLEXI over the last years. We would like to thank all these contributors for their efforts they spent on building FLEXI.

In case you have questions regarding FLEXI or want to contribute yourself by either reporting bugs, requesting features or adding something different to the project, feel free to open an issue or pull request.

```{toctree}
---
maxdepth: 2
caption: User Guide
numbered:
---
00_quickstart.md
01_installation.md
02_codeoverview.md
03_workflow.md
04_tutorials/index.md
05_toolsoverview.md
A1_parameter.md
```

```{toctree}
---
maxdepth: 1
caption: References
---
references.rst
```

[flexi]:         https://numericsresearchgroup.org/flexi_index.html
[nrg]:           https://numericsresearchgroup.org/index.html
