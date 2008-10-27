#
# Definintions
#
# Source filename extensions
#
DIVEXT=.d6a
COMEXT=.u6a
CEXT=.c
#
# VERSION to build
#
# VER=aix
#
VER=sun
#
# Source file locations
#
DIVMAIN=$(HOME)/divimp
COMMONS=$(DIVMAIN)/commons
GHOST=/afs/ipp/common/soft/ghost_v8/@sys/lib
#
DIVSRC =$(DIVMAIN)/div6/src
COMSRC =$(DIVMAIN)/comsrc
#
#
VPATH=$(DIVSRC):$(COMSRC)

#
# Suffixes
#
.SUFFIXES: $(DIVEXT) $(COMEXT)

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
OPTG= -g -ftrap=%none
OPTO= -O -ftrap=%none -xmaxopt -native -xlibmil -xlibmopt 
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
CFLAGS= $(OPTC)

#
# Libraries
#
#
# SUN
#
LIBS= -L/usr/local/lib -L$(HOME)/divimp/lib -lc -lghost -lpostsc
#
# LIBS=-L$(GHOST) -lghost -lgrid_postcl
#
# Name of target to build
#
TARG=div6O$(VER)
TARGALT=div6aO$(VER)

#
# Objects to compile
#

OBJECTS= adas.o adpak.o amjuel.o datetime$(VER).o harw.o nc.o sys$(VER).o utility.o utility_com.o utility2.o bgplasma.o cfd_osm.o cxrec.o div.o divinput.o divstore.o divoutput.o divtrn.o dummy.o eirediv.o eirene.o geier.o grad.o grid.o ion_transport.o ion_crossfield_transport.o ion_parallel_transport.o iztau.o mon.o neut.o neutone.o output.o pindiv.o plasma.o pputil.o redefves.o reiser.o relax.o rundiv.o setup.o sltmp.o sol.o sol23.o solascv0.o solascv1.o solascv2.o solascv3.o soledge.o sputter.o tau.o theta.o walls.o vacuum.o

#
# Rules
#

$(DIVEXT).o:
	cp $? $*.f
	$(FCOMP) $(FFLAGS) $(ARCH) -I$(IDIR) $(OPTS) -c $*.f
	$(RM) $*.f

$(COMEXT).o:
	cp $? $*.f
	$(FCOMP) $(FFLAGS) $(ARCH) -I$(IDIR) $(OPTS) -c $*.f
	$(RM) $*.f

#
# Target
#

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

alt:
	$(MAKE)  "TARG=$(TARGALT)" "OPT=$(OPTG)"

altopt:
	$(MAKE)  "TARG=$(TARGALT)" "OPT=$(OPTO)"

clean:
	$(RM) *.o
	$(RM) *.lst







