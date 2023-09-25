# Config file for a generic GCC based compile,
# Using the reference LAPACK/Blas and assuming
# both LAPACK and CFITSIO are in /usr/local.


# =============== Installation ===================

# Final install location.  This defaults to a
# subdirectory of the top-level directory called
# install_$COMMANDER.
#export INSTALL := /usr/local


# =============== System Tools ===================

# The defaults should be fine for any normal
# GNU/Linux operating system...

# Here you can override the make version to use
#export MAKE := gmake

# How to make a tar.gz file
#export CTAR := tar czvf

# The archive command
#export AR := ar

# Archive options to create a library
#export ARFLAGS := crv

# The ranlib command
export RANLIB := ranlib


# ================= Compilers ====================

# Serial Fortran 90 compiler.  Default is "f90".
export F90 := gfortran

# MPI Fortran 90 compiler.  Default is "mpif90".
export MPF90 := mpif90

# MPI Fortran 77 compiler.  Default is "mpif77".
export MPF77 := mpif77

# MPI C compiler.  Default is "mpicc".
export MPCC := mpigcc

# =========== Compiler Optimizations =============

# Fortran 90 compiler flags
export F90FLAGS := -O3 -assume byterecl -qopenmp
#export F90FLAGS := -g -C -assume byterecl -traceback

# Fortran 77 compiler flags
export FFLAGS := -O2

# C compiler flags.
export CFLAGS := -O2

# Extra flags used for linking
export LDFLAGS := -fopenmp
export OPENMP := -fopenmp


COMMON_FLAGS_DEBUG := -O0 -g
COMMON_FLAGS_FORTRAN_DEBUG := $(COMMON_FLAGS_DEBUG) $(BOUNDS) -check bounds -check format -check pointers -check uninit -check output_conversion -assume byterecl -traceback -heap-arrays 16384
#COMMON_FLAGS_FORTRAN_DEBUG := $(COMMON_FLAGS_DEBUG) -check bounds -check format -check pointers -check uninit -check output_conversion -assume byterecl
#export F90FLAGS_DEBUG := -vec_report0 $(COMMON_FLAGS_FORTRAN_DEBUG)
#export FFLAGS_DEBUG := $(COMMON_FLAGS_FORTRAN_DEBUG)
#export CFLAGS_DEBUG := $(COMMON_FLAGS_DEBUG)
#export LDFLAGS_DEBUG := -lm -cxxlib -qopenmp

# =============== Latex and PDF for documentation ================

# The include and linking commands for cfitsio.
export LATEX_COMP   := latex
export DVI2PDF_COMP := dvi2pdf

# =============== Main path ================

export LOCAL := /home/daniher/local

# =============== CFITSIO Linking ================

# The include and linking commands for cfitsio.
export CFITSIO_INCLUDE := 
export CFITSIO_LINK := -L$(LOCAL)/lib -lcfitsio 


# =============== LAPACK Linking =================
# The include and linking commands for LAPACK.
export LAPACK_INCLUDE := 
export LAPACK_LINK := -L/usr/lib -llapack

# =============== HEALPix Linking =================

# The include and linking commands for HEALPix.
export HEALPIX := /home/daniher/libraries/HEALPix/Healpix_3.82
export HEALPIX_INCLUDE := -I$(HEALPIX)/include
export HEALPIX_LINK := -L$(HEALPIX)/lib -lhealpix

# =============== FFTW3 ===========================
export FFTW_INCLUDE := -I$(LOCAL)/include
export FFTW_LINK := -L$(LOCAL)/lib -lfftw3 #-lfftw3f -lfftw3_threads -lfftw3f_threads

export PNG_LINK := -lpng

# =============== NOVAS ============================                            
export PNG_LINK := -lpng

# =============== HDF ============================                              
export HDF_LINK := -L$(LOCAL)/lib -lhdf5_fortran -lhdf5
export HDF_LINK_CPP := -L$(LOCAL)/lib -lhdf5_cpp -lhdf5
export HDF_INCLUDE := -I$(LOCAL)/include
