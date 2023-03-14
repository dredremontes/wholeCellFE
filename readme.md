
wholeCellFE - cellBP
===============================================================================

Please see ../eltopo3d/readme.txt for information on building the El Topo 
library.

This is an application based on El Topo library to simulate invagination.  This readme describes:

- how to build the executable
- how to run the executable
- an overview of the executable source code

Note: This is research-grade code -- full of hacks, bugs, dead code, etc.

Pre-requisites
=====================
C compiler (gcc 8.3.0 suggested)
OpenBLAS
Netlib-lapack

Building cellBP:
=====================

1. Create Makefile.local_defs in ./topoCell/cellBP/

The included Makefile reads a file called Makefile.local_defs, which contains 
platform-specific definitions.  You must create this file!  
Makefile.example_defs includes suggested settings for Linux and OS/X.

If you want to use the command-line interface, define NO_GUI in your build 
(e.g. in Makefile.local_defs).  If you want to use the GUI, be sure to link 
against OpenGL and GLUT.  

2. Generate dependencies and compile

Building the executable is done by running "make depend" followed by 
"make release" or "make debug".  This should also automatically build the
El Topo library.

Example:
$> make depend
$> make release

This will create the cellBP_release executable which can be used with input scripts to run the whole cell model.

Using cellBP:
=====================

Launching the executable requires two command line parameters:
1. A path to a text "script" file specifying the simulation type and initial 
geometry
2. A path specifying where output files will be written

Example with cell and substrate example input scripts:
$> ./cellBP_release scripts/cell_example.txt scripts/subs_example.txt

If you are running the GUI, this should pop up a GLUT window with a view of the 
triangle mesh surface (see main.cpp to learn what the keyboard does).  If you 
are running with no GUI, it will immediately start running the simulation 
defined by the script, outputting one binary mesh file per frame in /var/tmp/.

Post-processing:
=====================

For the example runs, simulation frames (.vtk files) will be saved in the ./wholeCellFE/cellBP/cell_subs_example/ 

Save locations can be changed in the cell and substrate input scripts by changing relative_output_path 

Simulation frames can be read by paraview

Updating integrin stiffness run:
=====================
The three integrin stiffness settings (baseline, 31 pN/nm, MD-driven) require separate main.cpp and UL_growth_explicit.cpp files. In ./wholeCellFE/backup_main there are .cpp files for each simulation setting which would have to replace the main.cpp and UL_growth_explicit.cpp files in ./wholeCellFE/cellBP/. The .cpp files used must be named main.cpp and UL_growth_explicit.cpp.

Alternatively, for constant integrin stiffness values only, the variable, "kint", can be updated in pN/um within UL_growth_explicit.cpp.

Code base:
=====================

Summaries of the important source files used by El Topo are found in the readme
in that directory.

Drivers:
---------------------

In cellBP, the dynamic surface is assigned a linear velocity per vertex by a 
"mesh driver".

References
=====================

[Bridson et al. 2007]: R. Bridson, J. Hourihan, and M. Nordenstam,  Curl noise 
for procedural fluid flow, Proc. ACM SIGGRAPH 2007

[Brochu and Bridson 2009]: Tyson Brochu and Robert Bridson, Robust Topological 
Operations for Dynamic Explicit Surfaces, SIAM J. Sci. Comput., vol. 31, no. 4 
(2009), pp. 2472-2493 

[Enright et al. 2002]: D. Enright, R. Fedkiw, J. Ferziger, and I. Mitchell, A 
hybrid particle level set method for improved interface capturing, J. Comput. 
Phys., 183 (2002), pp. 83–116.

[Jiao 2007]: X. Jiao, Face offsetting: A unified framework for explicit moving 
interfaces, J. Comput. Phys., 220 (2007), pp. 612–625.


