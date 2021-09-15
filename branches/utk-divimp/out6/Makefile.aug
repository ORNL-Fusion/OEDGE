#
# Definintions
#
# Source filename extensions
#
OUTEXT=.o6a
COMEXT=.u6a
CEXT=.c
#
# VERSION to build
#
# VER=aix
#
# VER=pgi
#
VER=sun
#
# Source file locations
#
DIVMAIN=$(HOME)/divimp
COMMONS=$(DIVMAIN)/commons
#
OUTSRC =$(DIVMAIN)/out6/src
COMSRC =$(DIVMAIN)/comsrc
#
#
VPATH=$(OUTSRC):$(COMSRC)

#
# Suffixes
#
.SUFFIXES: $(OUTEXT) $(COMEXT)

#
# Compiler
#
# FCOMP=xlfnb
#
FCOMP=f90
RM=rm -f
#
# Optimization and debugging flags
#
#
# SUN - F90
#
OPTC= -O
OPTG= -g  -ftrap=%none
OPTO= -O  -ftrap=%none
#
# SET options
#
OPT = $(OPTG)
#
# Fortran compiler options
#
FFLAGS= $(OPT) 
#
# Common block directory
#
IDIR=$(COMMONS)

#
# C-compiler
#
CC= gcc
#
# C-flags
#
CFLAGS= $(OPTC) -I$(IDIR)

#
# Libraries
#
#
# SUN
#
LIBSCOL = -L/usr/local/lib -L$(HOME)/divimp/lib -lc -lghost -lpostcl -ljpeg
LIBSBW = -L/usr/local/lib -L$(HOME)/divimp/lib -lc -lghost -lpostsc -ljpeg
#
LIBS=$(LIBSCOL)
#
# Name of target to build
#
TARG=out6O$(VER)
TARGALT=out6aO$(VER)
TARGCOL=out6OC$(VER)
TARGBW=out6OBW$(VER)

#
# Objects to compile
#

OBJECTS= adas.o amjuel.o datetime$(VER).o write_jpeg$(VER).o harw.o nc.o sys$(VER).o utility.o utility_com.o utility2.o contin.o dummy.o ioout.o outmain.o outbolo.o outcontour.o outlos.o outlos3D.o outplot.o outrcp.o outring.o outxsection.o out000.o out100.o out200.o out300.o out400.o out500.o out600.o out700.o out800.o out900.o outa00.o out966.o out970.o out972.o out974.o out978.o out980.o out981.o out982.o out983.o out990.o plrp.o reiser.o reiser_out.o slmod.o SLoutplot.o SLtrace.o trace.o 

#
#
# Rules
#

$(OUTEXT).o:
	cp $? $*.f
	$(FCOMP) $(FFLAGS) $(ARCH) -I$(IDIR) $(OPTS) -c $*.f
	$(RM) $*.f

$(COMEXT).o:
	cp $? $*.f
	$(FCOMP) $(FFLAGS) $(ARCH) -I$(IDIR) $(OPTS) -c $*.f
	$(RM) $*.f

$(TARG): $(OBJECTS)
	$(FCOMP) $(OBJECTS) $(ARCH) $(OPT) $(LIBS) -o $(TARG)
	-mv $(TARG) ../bin/@sys

#
# Various MAKE targets for different platforms and optimization
# levels - need to do an rm *.o to make sure that the whole
# thing is recompiled - it could be set up to do this automatically
# except for the one environment used for development.  
#

opt:
	$(MAKE)  "OPT=$(OPTO)"

col:
	$(MAKE)  "TARG=$(TARGCOL)" "LIBS=$(LIBSCOL)"

bw:
	$(MAKE)  "TARG=$(TARGBW)" "LIBS=$(LIBSBW)"

alt:
	$(MAKE)  "TARG=$(TARGALT)" "OPT=$(OPTG)"

altopt:
	$(MAKE)  "TARG=$(TARGALT)" "OPT=$(OPTO)"


clean:
	$(RM) *.o
	$(RM) *.lst







