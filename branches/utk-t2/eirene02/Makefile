
.SUFFIXES: .o .f

F77=$(PGI)/linux86/bin/pgf90

FFLAGS= -c -fast -tp athlon
#FFLAGS= -c -g -tp athlon
#FFLAGS= -c -g -O -Munroll -Mnoframe -Ktrap=fp -tp athlon
#FFLAGS= -c -C -g -fast -tp athlon 
#FFLAGS= -c -C -g -Ktrap=fp -fast -tp athlon 

OPTS=

LIBS=

IDIR=

VPATH=./Sources:./Sources/Modules 

TARG= eirene

MODULES= PRECISION.o avltree.o PARMMOD.o BRAEIR.o BRASCL.o BRASPOI.o CADGEO.o CAI.o CCONA.o CCOUPL.o CESTIM.o CFPLK.o CGEOM.o CGRID.o CGRPTL.o CINIT.o CLAST.o CLGIN.o CLMSUR.o CLOGAU.o COMNNL.o COMPRT.o COMSIG.o COMSOU.o COMSPL.o COMUSR.o COMXS.o COUTAU.o CPES.o CPL3D.o CPLMSK.o CPLOT.o CPOLYG.o CRAND.o CRECH.o CREF.o CREFMOD.o CSDVI_BGK.o CSDVI_COP.o CSDVI.o CSPEI.o CSPEZ.o CSTEP.o CTETRA.o CTEXT.o CTRCEI.o CTRIG.o CTSURF.o CUPD.o CZT1.o EIRBRA.o

OBJECTS= assistent.o broadcast_dummy.o couple_divimp.o diagno.o file-handling.o geometry.o iterate.o main-routines.o mathematic.o nag-dummy.o output.o particle-tracing.o plot_dummy.o sampling.o scoring.o startup-routines.o surface-processes.o tetraeder.o user.o volume-processes.o system-pgi.o

.f.o:
	$(F77) $(FFLAGS) $(OPTS) $?

eirene: $(MODULES) $(OBJECTS)
	$(F77) $(MODULES) $(OBJECTS) $(OPT) $(LIBS) -o $(TARG)

clean:
	rm *.o
	rm *.lst

