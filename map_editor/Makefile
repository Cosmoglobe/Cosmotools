include config/config.gnu

export F90FLAGS := -g -O3

export F90COMP := $(F90FLAGS) $(HEALPIX_INCLUDE) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(FFTW_INCLUDE) $(HDF_INCLUDE)
export LINK := $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK) $(FFTW_LINK) $(HDF_LINK) $(OPENMP)

COBJS  := math_tools.o sort_utils.o locate_mod.o spline_1D_mod.o spline_2D_mod.o \
	map_editor_utils.o map_editor_simple_ops_mod.o map_editor_complex_ops_mod.o 

all : map_editor libmap_editor.a

map_editor : libmap_editor.a map_editor.o
	$(MPF90) -o map_editor map_editor.o libmap_editor.a $(LINK) 

map_editor_complex_ops_mod.o    : map_editor_utils.o 

libmap_editor.a : $(COBJS)
	$(AR) $(ARFLAGS) libmap_editor.a $(COBJS)
	$(RANLIB) libmap_editor.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

%.o : %.cpp
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ map_editor

