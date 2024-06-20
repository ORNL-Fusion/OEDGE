# 
#  Copyright @ 1984 - 1995   Josef Heinen
# 
#  Permission to use, copy, and distribute this software and its
#  documentation for any purpose with or without fee is hereby granted,
#  provided that the above copyright notice appear in all copies and
#  that both that copyright notice and this permission notice appear
#  in supporting documentation.
# 
#  Permission to modify the software is granted, but not the right to
#  distribute the modified code.  Modifications are to be distributed
#  as patches to released version.
# 
#  This software is provided "as is" without express or implied warranty.
# 
#  Send your comments or suggestions to
#   J.Heinen@kfa-juelich.de.
# 

FOR = f77l3
FFLAGS = /NO /NW /D1 F77L3
AR = 386LIB
LINK = 386LINK
RM = DEL

.SUFFIXES: .obj .for

.for.obj:
	$(FOR) $*.for $(FFLAGS)

all: libgks.lib demo.exe

libgks.lib: gksio.obj gks.obj gksinq.obj gkserror.obj gksroot.obj \
	    gksmisc.obj gksdidd.obj gksdps.obj gksdhpgl.obj gksdpbm.obj \
	    gksdwiss.obj gkswiss.obj gksuns.obj gksafm.obj gksddos.obj
	$(AR) $@ -Create \
	    gksio.obj gks.obj gksinq.obj gkserror.obj gksroot.obj \
	    gksmisc.obj gksdidd.obj gksdps.obj gksdhpgl.obj gksdpbm.obj \
	    gksdwiss.obj gkswiss.obj gksuns.obj gksafm.obj gksddos.obj

clean:
	$(RM) *.exe
	$(RM) *.map
	$(RM) libgks.*
	$(RM) *.obj
	$(RM) *.sld

demo.exe: demo.obj libgks.lib
	$(LINK) demo.obj libgks.lib -lib graph3.lib
