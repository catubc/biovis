# biovis
Biovis is an OpenGL-based visualization tool built for the Allen Institute (Seattle) biophysically detailed simulator and network builder (forthcoming).

The code runs dynamically in jupyter and can be used to visualize single cell morphologies, synapses and connectivity matrices in networks of 10000 cells or more. Future implementations will visualize activated networks (example of previous simulations of 20,000 cell networks: https://www.dropbox.com/s/jcjj3owtbylw3w7/20kcells_movies.mov?dl=0).

**DEPENDENCIES**

There are a number of basic dependencies (e.g. numpy, PyQt4).  QtOpenGL and OpenGL can be installed via apt-get and pip in linux operating systems.

QtOpenGL: "sudo apt-get install python-qt4-gl"

OpenGL: "sudo pip install PyOpenGL PyOpenGL_accelerate"
