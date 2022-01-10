
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

VPATH=.:../utils

TARG= procts

MODULES = error_handling.o utilities.o array_sort.o filter_elms.o analyse_ts.o 

OBJECTS= processts.o


.f90.o:
	$(F90) $(FFLAGS) $(OPTS) $?

test: $(MODULES) $(OBJECTS)
	$(F90) $(MODULES) $(OBJECTS) $(OPT) $(LIBS) -o $(TARG)

clean:
	rm *.o
	rm *.mod

