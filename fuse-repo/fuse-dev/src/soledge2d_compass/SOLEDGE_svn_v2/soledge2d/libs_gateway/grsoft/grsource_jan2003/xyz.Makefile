#
#             GR/GR3 Makefile
#
#               Marlene Busch     
#         Institut fuer Mathematik                 
#     Forschungszentrum Juelich GmbH
#	      W-52425 Juelich
#
#             Aug 24, 2000
#
 
      SHELL = /bin/ksh
   GLI_HOME = "/home/boerner/grgli/gligks/gks_intel"
         CC = gcc
     CFLAGS = -g
        F77 = f77
   F77FLAGS = -static
      SEGLD = cc
         LD = ld
    LDFLAGS =
         AR = ar
     RANLIB = ranlib
   BLDFLAG = ar
    GRHHFL = grhhfl.o
    SOFLAGS =
   XLIBPATH =
     XLIBS = $(XLIBPATH) -lXt -lX11
      LIBS = -lm


OBJS =  compress.o \
 gr2isl.o                                   gr3pan.o  gr3ste.o gr3trc.o \
 gr3zbu.o gr90dg.o grarrw.o graxlin.o graxlog.o graxs.o \
 graxsl.o grbld.o grchn.o grchnc.o grchrc.o grclp.o \
 grcrcl.o grdn.o grdrax.o grdrdm.o grdrhs.o grdrlg.o \
 grdrw.o grdrws.o grdsh.o grend.o grfill.o grflls.o \
 grfont.o grftoc.o grgfld.o grgks.o  grhhnl.o \
 grhhtr.o grhpan.o grjmp.o grjmps.o grlgar.o grlgnd.o \
 grln.o grlncn.o grlnlg.o grmrks.o grmskf.o grmskn.o \
 grnwpn.o grnxtf.o grpctr.o grpts.o grrahm.o grscan.o \
 grscax.o grscdl.o grscla.o grsclc.o grsclp.o grsclv.o \
 grshd.o grshow.o grsphr.o grspts.o grstrt.o grtxscn.o \
 grtxt.o grtxtc.o grvar.o grvfld.o grwin.o gzclip1.o \
 gzcol.o  gzinit1.o gzplot1.o gzseed.o gztran1.o \
 gzzz1.o kurvef.o pskurf.o skurf.o  gr3plo1.o gr3plo2.o   gr3dim.o

F90OBJS =  grsoft_interface_block.o grhhfl.o

CINTOBJS =  grsoft_c_interface.o main.o

 
.SUFFIXES: .a .o .f .c  .f90

 
.f.o:
	$(F77) -c  $(F77FLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

.f90.o:
	$(F90) -c $(F90FLAGS) $<

default:
	@make mod
	@make `./config`

mod:
	@chmod 555 config

usage:  @echo "Can't obtain system information."; \
"Usage: make [aix|alpha|cray|sun|sunos|irix|hpux-nagware|pclinux-nageware|pclinux-pgf90|pclinux_intel]"

aix :
	@make targetsaix F77="xlf" \
	CFLAGS="-I$(GLI_HOME) -Daix -g" \
	F77FLAGS="-w -qfloat=hssngl:rsqrt -g"\
	F90="xlf90" \
	F90FLAGS="-w -qfloat=hssngl:rsqrt -g"\
	RANLIB="ar ts" \
	BLDFLAG="ar rv" \
	LIBS="-lgks" \
	SOFLAGS="-bM:SRE -bE:gr.exp -e _nostart" \
	LD="xlf" \
	XLIBS="-lXt -lXext -lX11"

sun4: sun
sun :
	@make targets F77="f77" \
	CC="gcc " \
	CFLAGS="-I$(GLI_HOME)" \
	F77FLAGS="-O3 " \
	F90=" " \
	F90FLAGS=" " \
	F90OBJS=" " \
	RANLIB="ranlib" \
	BLDFLAG="ar cr"  \
	SOFLAGS="" \
	LIBS="-lgks" \
	LDFLAGS="-Bstatic -L/usr/local/fortran/SC1.0" \
	XLIBPATH="-L/usr/openwin/lib" LIBS="-lnsl"

sunos :
	@make targets  F77="f90" \
	CFLAGS="-I$(GLI_HOME)" \
	F77FLAGS="-O3 " \
	F90="f90 " \
	F90FLAGS="-O3 " \
	SHLIBS="libgr.so" SOFLAGS="-G" XLIBPATH="-L/usr/openwin/lib" \
	LIBS="-lgks -lnsl" \
	RANLIB="echo" \
	BLDFLAG="ar cr" 

irix:
	@make targets F77="f90" \
	CFLAGS="-I$(GLI_HOME) -n32" \
	F77FLAGS="-O2 -n32" \
	F90="f90 -n32" \
	F90FLAGS="-O2 " \
	RANLIB="echo" \
	BLDFLAG="ar cr" \
	LD="ld -shared -n32" \
	LIBS="" 

alpha :
	@make targets F77="f90" \
	CFLAGS="-I$(GLI_HOME)" \
	F77FLAGS="-O " \
	F90="f90" \
	F90FLAGS="-O " \
	RANLIB="ar ts" \
	BLDFLAG="ar rv"  \
	LD="ld -shared" \
	XLIBS="-lXt -lXext -lX11" \
	LDFLAGS="-O3" \
	LIBS="-lgks -lUfor -lfor -lFutil -lots  -lm -lc"
 
sn2401: cray
sn7104: cray
sn7119: cray
sn9612: cray
sn1011: cray
#       grsoft_interface.block.f90 jetzt nicht mehr in library
#       liegt als Source unter /usr/local/grsoft und kann
#       dann mit den Benutzer spez. Optionen herangezogen werden
#       Option: -dp alles einfach genau
cray :
	@make   libgr.a F77="f90 " \
	CFLAGS="-Dcray -I$(GLI_HOME)" \
	F77FLAGS="-en -dp  -O3 " \
	F90="f90 " \
	F90FLAGS=" -O3 -dp " \
	F90OBJS=" grhhfl.o" \
	RANLIB="bld tev" \
	BLDFLAG="bld rz"  
# ausprobiert fuer Nag95 Compiler unter AIX
hpux-nagware:
	@make targetsaix  CFLAGS="-I$(GLI_HOME) -Daix" \
	F77="nagf95" F77FLAGS="-f77 -mismatch_all" SEGLD="f95"  \
	LIBS="-lgks" \
	F90="nagf95" F90FLAGS="-f77 -mismatch_all" \
	RANLIB="ranlib" \
	BLDFLAG="ar cr"  \
	SOFLAGS="" LD="nagf95"
i386: pclinux
i486: pclinux
i586: pclinux
i686: pclinux
pclinux-nagware:
	@make targets F77="f95" \
	CFLAGS="-I$(GLI_HOME)" \
	F77FLAGS=" -mismatch_all  " \
	F90="f95" \
	F90FLAGS=" -mismatch_all   " \
	RANLIB="ranlib" \
	BLDFLAG="ar cr" \
	SEGLD="cc" \
	XLIBS="-L/usr/X11R6/lib -lXt -lXext -lX11" \
	LD="pgf90 -shared " \
	LDFLAGS="-L/usr/local/NAGWare_f95/lib -lf95 -u MAIN__" \
	LIBS="-lgks"
pclinux-pgf90:
	@make targetsaix F77="pgf90" \
	CFLAGS="-I$(GLI_HOME)" \
	F77FLAGS="  " \
	F90="pgf90" \
	F90FLAGS="  " \
	RANLIB="ranlib" \
	BLDFLAG="ar cr" \
	SEGLD="pgf90 -M nomain" \
	XLIBS="-L/usr/X11R6/lib -lXt -lXext -lX11" \
	LD="ld  " \
	LDFLAGS=" " \
	LIBS="-lgks"
pclinux-intel:
	@make targetsaix F77="ifc7.1" \
	CFLAGS="-I$(GLI_HOME)" \
	F77FLAGS="  " \
	F90="ifc7.1" \
	F90FLAGS="  " \
	RANLIB="ranlib" \
	BLDFLAG="ar cr" \
	SEGLD="cc" \
	XLIBS="-L/usr/X11R6/lib -lXt -lXext -lX11" \
	LD="ifc7.1" \
	LDFLAGS=" " \
	LIBS="-lgks"

targets: libgr.a libgr.so
targetsaix: libgr.a

libgr.a: $(OBJS) $(F90OBJS) $(CINTOBJS)
	$(BLDFLAG)  $@ $?
	$(RANLIB)  libgr.a
 
libgrf.a: $(OBJS) $(F90OBJS) $(CINTOBJS)
	$(BLDFLAG)  $@ $?
	$(RANLIB)  libgrf.a

libgr.so: $(OBJS) $(F90OBJS) $(CINTOBJS)
	$(LD) -o $@ $(SOFLAGS) $(OBJS) $(F90OBJS) $(CINTOBJS) \
	$(LDFLAGS) $(XLIBS) $(LIBS)

 
clean:
	rm -f libgr.a
	rm -f $(OBJS) $(F90OBJS) $(CINTOBJS) grsoft_interface_block.mod
#rm -f libgrf.a
# bisher bei irix: LIBS="-lgks -lsun -lm" 
