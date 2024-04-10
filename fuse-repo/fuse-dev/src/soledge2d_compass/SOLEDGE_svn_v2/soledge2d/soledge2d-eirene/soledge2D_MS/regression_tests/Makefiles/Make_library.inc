# Fortran compiler and linker
FORTCOMP = ifort

# Compilation options
OPT_COMP = -c -O2 -fpp -r8
#OPT_COMP = -c -O0 -fpp -p -g -traceback -check all


# HDF5 include files and library
INCHDF5 = -I$(HOME)/libs/hdf5/include
LIBHDF5 = -L$(HOME)/libs/hdf5/lib -lhdf5_fortran -lhdf5 $(HOME)/libs/zlib/lib/libz.a $(HOME)/libs/szip/lib/libsz.a 

# Linear solvers include files and libraries
PASTIXDIR = /home/bufferand/libs/pastix
SCOTCHDIR = /home/bufferand/libs/scotch
INCSOLVER = -I$(PASTIXDIR)
LIBSOLVER = -L$(PASTIXDIR) -lpastix -lm -lrt -L$(SCOTCHDIR)/lib  -lscotch -lscotcherrexit -lpthread

#MKL FOR LAPACK AND BLAS
LIBLAPACK = -L$(MKL_HOME) -Wl,-R$(MKL_HOME) -lmkl_lapack95_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group

# OPENMP library
INCOMP = #-openmp
LIBOMP = #-openmp
