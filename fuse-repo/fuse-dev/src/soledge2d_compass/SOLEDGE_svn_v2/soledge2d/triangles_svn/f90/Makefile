#-------------------------
#  Makefile for SOLEDGE-2D
#-------------------------


# Fortran compiler and linker
FORTCOMP = ifort

# Compilation options
OPT_COMP = -c -O2 -cpp

# Linking options
OPT_LNK = 

# Performance options
OPT_PERF =

# Defining output (executable) file name
TARGET = SOLEDGE2D

# HDF5 include files and library
INCHDF5 = -I/home/bufferand/libs/hdf5/include
LIBHDF5 = -L/home/bufferand/libs/hdf5/lib -lhdf5_fortran -lhdf5 /home/bufferand/libs/zlib/lib/libz.a /home/bufferand/libs/szip/lib/libsz.a

# OPENMP library
INCOMP = -openmp
LIBOMP = -openmp


# Defining fortran 90 source files
SRCS_f90 = domain_types.f90 \
	solver_types.f90 \
	RungeKutta.f90 \
	equimag.f90 \
	interp_soledge2d_eirene.f90 \
	jacobs.f90\
	metrix.f90 \
	neutrals.f90 \
	saving.f90 \
	penalisation.f90 \
	multidomain.f90 \
	center.f90 \
	implicit.f90 \
	initialisation.f90 \
	derivatives.f90 \
	sources_perp.f90 \
	sources_para.f90 \
	N_solver.f90 \
	G_solver.f90 \
	Te_solver.f90 \
	Ti_solver.f90 \
	CFL_dt.f90 \
	tridiag.f90 \
	monitors.f90 \
	soledge2m.f90

# Object files corresponding to source files
OBJS = $(SRCS_f90:.f90=.o)

# Defining paths to libraries


# Assembling required librqries and include files
LIB = $(LIBHDF5) $(LIBOMP)
INC = $(INCHDF5) $(INCOMP)

.SUFFIXES: .o .f90

# Linking main program
$(TARGET):	$(OBJS)
		$(FORTCOMP) $(OPT_LNK) -o $(TARGET) $(OBJS) $(LIB)

# Rules for compiling files
.f90.o:
	$(FORTCOMP) $(OPT_COMP) $(INC) $<

# For cleaning compilation created files
clean:	
	rm -f $(OBJS) *~ *.mod *.kmo *.l *.L core $(TARGET)
