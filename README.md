# Cosmotools
Collection of tools for editing and visualizing HEALPix maps, and other useful things:

* map_editor:
  A Fortran compiled set of tools for manipulating HEALPix maps, with functionalities including, 
  but not limited to, smoothing with a beam, ud_grading, returning map statistics (mean, min, max, rms),
  scaline, converting QU maps to P, threshold masking, etc.

* python tools:
  List is small for now, but there currently exist three importable scripts which allow for:
       * Converting units for frequencies or over bandpasses
       * creating noise realizations from rms maps
       * correlating two maps

More tools will be added in the future.

---------------------------------------------------

Compiling map_editor:

If you already have compiled Commander, all you have to do is change the LOCAL variable in map_editor/Makefile to your own commander/build/install directory, e.g.:
export LOCAL=/mn/stornext/u3/(username)/git/commander/build/install

and make sure to have the right modules loaded by for instance having this in .bashrc:
module load intel/oneapi mpi/latest compiler-rt/latest compiler/latest mkl/latest

