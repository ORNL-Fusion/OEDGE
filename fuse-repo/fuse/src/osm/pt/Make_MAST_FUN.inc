#############################
# Platform type
#############################

PLATFORM = MAST_FUN

#############################
# Compilation and linking
#############################

# Fortran compiler and linker
FORTCOMP = ifort

# Compilation options
OPT_COMP = -c -C -O0

# Linking options
OPT_LNK =

# Performance options
OPT_PERF =

#############################
# Libraries
#############################

# MPI include files and library
#INCMPI = -I/usr/local/intellibs/include
#LIBMPI = -L/usr/local/intellibs/lib -lmpich

# BLAS and LAPACK libraries
LIBBLAS = -lblas_LINUX
LIBLAPACK = -L. -llapack_LINUX

# BLACS and SCALAPACK
#LIBBLACS = $(HOME)/f90/libs/lib/blacs_MPI-MAST_FUN-0.a   $(HOME)/f90/libs/lib/blacsF77init_MPI-MAST_FUN-0.a  $(HOME)/f90/libs/lib/blacs_MPI-MAST_FUN-0.a
#LIBSCALAPACK = $(HOME)/f90/libs/lib/libscalapack-MAST_FUN.a

# SPARSKIT library
#LIBSPARSKIT = $(HOME)/f90/libs/lib/libskit-MAST_FUN.a

# Linear solvers include files and libraries
#INCSOLVER = -I$(HOME)/f90/libs/include
#LIBSOLVER = $(HOME)/f90/libs/lib/libpord-MAST_FUN.a $(HOME)/f90/libs/lib/libdmumps-MAST_FUN.a $(HOME)/f90/libs/lib/libzmumps-MAST_FUN.a

# Fast Fourier Transform
INCFFT = -I/usr/local/intellibs/include
LIBFFT = -L/usr/local/intellibs/lib -lfftw3
