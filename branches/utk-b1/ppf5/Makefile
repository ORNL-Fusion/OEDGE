# suffixes
.SUFFIXES: .f5a
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
TARGOPT= ppf5
TARGDBG= ppf5dbg
#
# Other source directories
#
DSRC= ../comsrc/
#
# Libraries
#
LIBS= -lppf -lxlf
#
# Objects to compile
#
OBJECTS=divppf.o wrtdda.o ioppf.o gff.o utlppf.o \
        adas.o harw.o sysrs6k.o
#
# rules
#
.f5a.o:
	cp $*.f5a $*.f
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

adas.o: $(DSRC)adas.u5a
	cp $(DSRC)adas.u5a adas.f
	$(F77) $(FFLAGS) -I$(IDIR) -c adas.f
	rm adas.f

harw.o: $(DSRC)harw.u5a
	cp $(DSRC)harw.u5a harw.f
	$(F77) $(FFLAGS) -I$(IDIR) -c harw.f
	rm harw.f

sysrs6k.o: $(DSRC)sysrs6k.u5a
	cp $(DSRC)sysrs6k.u5a sysrs6k.f
	$(F77) $(FFLAGS) -I$(IDIR) -c sysrs6k.f
	rm sysrs6k.f

clean:
	rm *.o
	rm *.lst





