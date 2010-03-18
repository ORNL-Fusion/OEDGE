# suffixes
.SUFFIXES: .f6a
# defines
#
# Compiler
#
F77=xlf
#
# Common block directory
#
IDIR=../commons
#
# Optimization and debugging flags
#
OPTG = -g -C -qsource -qextchk -bloadmap:load
OPTO = -O
#
# Fortran compiler options
#
FFLAGS= $(OPT)
#
# Name of target to build
#
TARGOPT= trn6
TARGDBG= trn6dbg
#
# Other source directories
#
DSRC= ../comsrc/
#
# Libraries
#
LIBS= -lxlf
#
# Objects to compile
#
OBJECTS=iotrn.o divtrn.o utility.o sysrs6k.o
#
# rules
#
.f6a.o:
	cp $*.f6a $*.f
	$(F77) $(FFLAGS) -I$(IDIR) -c $*.f
	rm $*.f

$(TARG): $(OBJECTS)
	$(F77) $(OBJECTS) $(FFLAGS) $(LIBS) -o $(TARG)
#
# Various MAKE targets for different platforms and optimization
# levels - need to do an rm *.o to make sure that the whole
# thing is recompiled - it could be set up to do this automatically
# except for the one environment used for development.
#
opt:
	$(MAKE) "TARG=$(TARGOPT)" "OPT=$(OPTO)"
dbg:
	$(MAKE) "TARG=$(TARGDBG)" "OPT=$(OPTG)"
#
# Directions to make files with source in other directories.
#

utility.o: $(DSRC)utility.u6a
	cp $(DSRC)utility.u6a utility.f
	$(F77) $(FFLAGS) -I$(IDIR) -c utility.f
	rm utility.f

sysrs6k.o: $(DSRC)sysrs6k.u6a
	cp $(DSRC)sysrs6k.u6a sysrs6k.f
	$(F77) $(FFLAGS) -I$(IDIR) -c sysrs6k.f
	rm sysrs6k.f

clean:
	rm *.o
	rm *.lst
