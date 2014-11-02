
.SUFFIXES: .o .f90

F90=pgf90

#
# jdemod - pgf compiler
#
FFLAGS= -g -c -traceback -Mbounds -Mdepchk -Mchkptr -Mchkfpstk -Mfptrap

#FFLAGS= -c -O2


OPTS=
#OPTS= -Mnolist -Mdepchk -Minform=inform -Minfo -Mchkptr -Mchkfpstk -Mfptrap -Mbounds
#OPTS= -Mnolist -Mdepchk -Minform=inform -Minfo -Mchkfpstk -Mfptrap -Mbounds

LIBS= 

IDIR=

VPATH=../../../comsrc/f90:../:../../code/utils


MODULES = divimp_types.o error_handling.o allocate_arrays.o utilities.o oedge_plasma_interface.o

SP_OBJECTS= swing_probe.o
SP_TARG=write_sp

TEST_OBJECTS= plasma_test.o
TEST_TARG= test_op

DIMES_OBJECTS= write_dimes_plasma.o
DIMES_TARG= wr_dimes

ERO_OBJECTS= ero_plasma.o write_ero_plasma.o
ERO_TARG = ero_plasma


.f90.o:
	$(F90) $(FFLAGS) $(OPTS) $?

dimes: $(MODULES) $(DIMES_OBJECTS)
	$(F90) $(MODULES) $(DIMES_OBJECTS) $(OPT) $(LIBS) -o $(DIMES_TARG)

ero: $(MODULES) $(ERO_OBJECTS)
	$(F90) $(MODULES) $(ERO_OBJECTS) $(OPT) $(LIBS) -o $(ERO_TARG)

swing: $(MODULES) $(SP_OBJECTS)
	$(F90) $(MODULES) $(SP_OBJECTS) $(OPT) $(LIBS) -o $(SP_TARG)

test: $(MODULES) $(TEST_OBJECTS)
	$(F90) $(MODULES) $(TEST_OBJECTS) $(OPT) $(LIBS) -o $(TEST_TARG)

clean:
	rm *.o
	rm *.mod

