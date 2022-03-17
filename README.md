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

---------------------------------------------------

In order to make sure everything works well on the owls, we need to point to the compiled libraries.
In your bashrc, add the following line to ensure that you are always pointing to the right directories:
`export LD_LIBRARY_PATH=/mn/stornext/u3/hke/owl/local/lib:$LD_LIBRARY_PATH`

A list of map_editor dependencies is on the to-do list.

- Daniel

