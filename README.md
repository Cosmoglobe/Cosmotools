# Cosmotools
Collection of tools for editing and visualizing HEALPix maps, and other useful things:

- map_editor:
  A Fortran compiled set of tools for manipulating HEALPix maps, with functionalities including, 
  but not limited to, smoothing with a beam, ud_grading, returning map statistics (mean, min, max, rms),
  scaline, converting QU maps to P, threshold masking, etc.

- python tools:
  List is small for now, but there currently exist three importable scripts which allow for:
       - Converting units for frequencies or over bandpasses
       - creating noise realizations from rms maps
       - correlating two maps

More tools will be added in the future.

######################################################

For those of you on the owls in Oslo, map_editor should compile by simply typing 'make' in the map_editor 
directory. If it does not, try running `export LD_LIBRARY_PATH=/mn/stornext/u3/hke/owl/local/lib:$LD_LIBRARY_PATH`
in the terminal and then try to make again.

A list of map_editor dependencies is on the to-do list.

- Daniel

