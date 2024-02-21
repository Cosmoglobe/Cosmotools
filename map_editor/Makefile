# The ranlib command
export RANLIB := ranlib

# ================= Compilers ====================
# Serial Fortran 90 compiler.  Default is "f90".
export F90 := ifort

# =========== Compiler Optimizations =============
# Fortran 90 compiler flags
export F90FLAGS := -g -O3 -assume byterecl
#export F90FLAGS := -g -C -assume byterecl -traceback

# Extra flags used for linking
export LDFLAGS := -cxxlib -qopenmp
export OPENMP := -qopenmp

# =============== Main path ================
# Let this point to your own Commander directory
export LOCAL=/mn/stornext/u3/ingunnkw/git/commander/build/install

# =============== CFITSIO Linking ================
# The include and linking commands for cfitsio.
export CFITSIO_INCLUDE := 
export CFITSIO_LINK := -L$(LOCAL)/lib -lcfitsio 

# =============== SHARP Linking ================
export SHARP_INCLUDE :=
export SHARP_LINK := -L$(LOCAL)/lib/ -lsharp

# =============== LAPACK Linking =================
# The include and linking commands for LAPACK.
MKLPATH := $(MKLROOT)
export LAPACK_INCLUDE := 
export LAPACK_LINK := -shared-intel -Wl,-rpath,$(MKLPATH)/lib/intel64 -L$(MKLPATH)/lib/intel64\
  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread

# =============== HEALPix Linking =================
# The include and linking commands for HEALPix.
export HEALPIX := $(LOCAL)/healpix
export HEALPIX_INCLUDE := -I$(HEALPIX)/include
export HEALPIX_LINK := -L$(HEALPIX)/lib -lhealpix

# =============== FFTW3 ===========================
export FFTW_INCLUDE := -I$(LOCAL)/include
export FFTW_LINK := -L$(LOCAL)/lib -lfftw3 -lfftw3f -lfftw3_threads -lfftw3f_threads
export PNG_LINK := -lpng

# =============== HDF ============================                              
export HDF_LINK := -L$(LOCAL)/lib -lhdf5_fortran -lhdf5
export HDF_LINK_CPP := -L$(LOCAL)/lib -lhdf5_cpp -lhdf5
export HDF_INCLUDE := -I$(LOCAL)/include

# ======== collecting it all =================================
export F90COMP := $(F90FLAGS) -I. $(HEALPIX_INCLUDE) $(LAPACK_INCLUDE) $(SHARP_INCLUDE) $(CFITSIO_INCLUDE) $(FFTW_INCLUDE) $(HDF_INCLUDE)
export LINK := $(HEALPIX_LINK) $(LAPACK_LINK) $(SHARP_LINK) $(CFITSIO_LINK) $(FFTW_LINK) $(LDFLAGS) $(HDF_LINK) $(OPENMP)

COBJS  := math_tools.o sort_utils.o locate_mod.o spline_1D_mod.o spline_2D_mod.o \
	map_editor_utils.o map_editor_simple_ops_mod.o map_editor_complex_ops_mod.o 

all : libmap_editor.a map_editor 

map_editor : libmap_editor.a map_editor.o
	$(F90) -o map_editor map_editor.o libmap_editor.a $(LINK) 

map_editor_complex_ops_mod.o    : map_editor_utils.o 

libmap_editor.a : $(COBJS)
	$(AR) $(ARFLAGS) libmap_editor.a $(COBJS)
	$(RANLIB) libmap_editor.a

%.o : %.F90
	$(F90) $(F90COMP) -c $<

%.o : %.f90
	echo "$(F90)  $(F90COMP)";\
	$(F90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

%.o : %.cpp
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ map_editor

