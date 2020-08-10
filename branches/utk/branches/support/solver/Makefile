
.SUFFIXES: .o .f90 .f

F90=pgf90
F90GCC=gcc
#
# jdemod - pgf compiler
#
#FFLAGS= -g -c -traceback -Mbounds -Mdepchk -Mchkptr -Mchkfpstk -Mfptrap
FFLAGS= -g -c -Mlist -Mdclchk -Mdepchk -Mfptrap -Mbounds -Msave -traceback -Mchkstk -Mchkfpstk

#FFLAGSGCC= -g -c -fbacktrace -O0 -ffree-line-length-none -march=x86-64 -Wall -fcheck=bounds -fcheck=array-temps -fcheck=do -fcheck=mem -fcheck=pointer 

FFLAGSGCC= -g -c -fbacktrace -O0 -ffree-line-length-none -march=x86-64  -fcheck=bounds -fcheck=array-temps -fcheck=do -fcheck=mem -fcheck=pointer 

#FFLAGS= -c -O2


OPTS=
#OPTS= -Mnolist -Mdepchk -Minform=inform -Minfo -Mchkptr -Mchkfpstk -Mfptrap -Mbounds
#OPTS= -Mnolist -Mdepchk -Minform=inform -Minfo -Mchkfpstk -Mfptrap -Mbounds

DIVMAIN=$(HOME)/divimp-s22
LOCAL_LIB=$(DIVMAIN)/local/lib
LIBS= -L$(LOCAL_LIB) -lghost -lpostsc

IDIR=

#VPATH=../../../comsrc/f90:../:../../code/utils
VPATH=../../comsrc/f90:../../comsrc:../code/utils:./


MODULES = divimp_types.o error_handling.o mod_io_units.o allocate_arrays.o debug_options.o mod_reader.o mod_solparams.o mod_solswitch.o mod_solcommon.o mod_cadas2.o adas.o mod_sol22_debug.o mod_cx.o mod_cadas.o mod_sol22pmom.o mod_sol22pei.o mod_sol22phelpi.o mod_sol22pcx.o mod_solrk.o mod_allocate_sol22_storage.o  mod_sol22_sources.o mod_sol22_utils.o mod_sol22_solver.o mod_sol22_input.o mod_sol22_output.o mod_sol22_plots.o  mod_sol22.o 

SOL22_OBJECTS= utility_com.o unstructured_input.o mod_calcsol_interface.o sol22_test.o 
SOL22_TARG=sol22_test


.f90.o:
	$(F90) $(FFLAGS) $(OPTS) $?

.f.o:
	$(F90) $(FFLAGS) $(OPTS) $?

sol22: $(MODULES) $(SOL22_OBJECTS)
	$(F90) $(MODULES) $(SOL22_OBJECTS) $(OPT) $(LIBS) -o $(SOL22_TARG)

gcc:
	$(MAKE) "F90=$(F90GCC)" "FFLAGS=$(FFLAGSGCC)" "LIBS=-L/home/david/divimp-gf/local/lib -lgfortran -lm -lghost -lpostsc"

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.lst

