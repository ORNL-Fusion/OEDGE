CFLAGS        = -I$(INCDIR) 

XLIBS         = -lX11 -lXext

OLIBS         = -lm -lbsd

ARFLAGS       = ruv

CC            = cc

DEST	      = .

EXTHDRS	      = /usr/include/stdio.h 

HDRSl	      = com.h \
		defs.h \
		globals.h

INCDIR        = /usr/include

LDFLAGS	      = 

LIBS	      = $(XLIBS) $(OLIBS) -lcurses -ltermcap

LINKER	      = cc

MAKEFILE      = Makefile

OBJS	      = lswap.o

PRINT	      = pr

PROGRAM	      = lswap

RANLIB        = ranlib

SRCS	      = lswap.c

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS)
		@echo -n "Loading $(PROGRAM) ... "
		@$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM) 
		@echo "done"

clean:;		@rm -f $(OBJS)

depend:;	@mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@install -s $(PROGRAM) $(DEST)

print:;		@$(PRINT) $(HDRS) $(SRCS)

program:        $(PROGRAM)

tags:           $(HDRS) $(SRCS); @etags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

###


