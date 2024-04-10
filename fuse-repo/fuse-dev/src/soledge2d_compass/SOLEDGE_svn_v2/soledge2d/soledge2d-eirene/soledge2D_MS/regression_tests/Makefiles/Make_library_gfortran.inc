# Fortran compiler and linker
#ORIG FORTCOMP = ifort
FORTCOMP = gfortran

# Compilation options
#ORIG OPT_COMP = -c -O2 -fpp -r8
OPT_COMP =  -c -O2 -ffree-line-length-none -x f95-cpp-input #-fdefault-real-8
#OPT_COMP = -c -O0 -fpp -p -g -traceback -check all -fpe0


# HDF5 include files and library
#ORIG HDF5_DIR=$(HOME)/libs/hdf5
#ORIG SZLIB_DIR=$(HOME)/libs/szip/lib
#ORIG ZLIB_DIR=$(HOME)/libs/zlib/lib
HDF5_DIR=/usr
SZLIB_DIR=/usr/lib
ZLIB_DIR=/usr/lib/x86_64-linux-gnu
INCHDF5 = -I$(HDF5_DIR)/include
#ORIG LIBHDF5 = -L$(HDF5_DIR)/lib -L$(SZLIB_DIR) -L$(ZLIB_DIR) -lhdf5_fortran -lhdf5 #-lz -lsz
LIBHDF5 = -L$(HDF5_DIR)/lib -L$(SZLIB_DIR) -L$(ZLIB_DIR) -lhdf5_fortran -lhdf5 

#MKL FOR LAPACK AND BLAS
#ORIG LIBLAPACK = -L$(MKL_HOME) -Wl,-R$(MKL_HOME) -lmkl_lapack95_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group
LIBLAPACK=-llapack -lblas
