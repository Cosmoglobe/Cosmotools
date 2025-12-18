# ===============================================================
# Makefile for Cosmotools utilities
# Works with Healpix 3.83, FFTW3, HDF5, PNG
# ===============================================================

# Compiler
CXX := g++

# ================================
# External libraries
# ================================

# Healpix
HEALPIX := /home/dwatts/Downloads/Healpix_3.83
HEALPIX_INCLUDE := -I$(HEALPIX)/include -I$(HEALPIX)/include/healpix_cxx
HEALPIX_LIB := -L$(HEALPIX)/lib -lhealpix_cxx -lgomp

# FFTW3
FFTW := /home/dwatts/Downloads/fftw-3.3.10
FFTW_INCLUDE := -I$(FFTW)/include
FFTW_LIB := -L$(FFTW)/lib -lfftw3 -lfftw3_threads -lm -lpthread

# HDF5
HDF5 := /home/dwatts/local/hdf5-2.0.0
HDF5_INCLUDE := -I$(HDF5)/include
HDF5_LIB := -L$(HDF5)/lib -lhdf5 -lz -ldl -lm -pthread  # drop C++ lib for now

# PNG
PNG_INCLUDE := -I/usr/include/libpng16
PNG_LIB := -lpng

# ================================
# Project libraries
# ================================
LIBUTIL := ../libutil
EPHCOM := ../ephcom

LIB_INCLUDE := -I$(LIBUTIL) -I$(EPHCOM) -I.
LIB_LIB := -L$(LIBUTIL) -lutil -L$(EPHCOM) -lephcom

# ================================
# Compiler and linker flags
# ================================
CPPFLAGS := $(LIB_INCLUDE) $(HEALPIX_INCLUDE) $(FFTW_INCLUDE) $(HDF5_INCLUDE) $(PNG_INCLUDE)
CXXFLAGS := -O3 #-Wall -Wextra
LDFLAGS := $(LIB_LIB) $(HEALPIX_LIB) $(FFTW_LIB) $(HDF5_LIB) $(PNG_LIB)

# ================================
# Targets
# ================================
TARGETS := ppm2png powspec fft gp2tab mjd2unix tcat eph idate cut2accept \
           ces_translate date2mjd l1_translate map2png
           #ps2fits ps2ascii patch_ident ces2runlist lxtxt ceslist_merge \
           #ceslistdiff osamize sfor hcat reverse ps2map map2ps ps2fits smap

# ================================
# Build rules
# ================================
all: $(TARGETS)

%: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGETS) *.o

