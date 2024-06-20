C*
C* Copyright @ 1984 - 1994   Josef Heinen
C*
C* Permission to use, copy, and distribute this software and its
C* documentation for any purpose with or without fee is hereby granted, 
C* provided that the above copyright notice appear in all copies and
C* that both that copyright notice and this permission notice appear
C* in supporting documentation.
C*
C* Permission to modify the software is granted, but not the right to
C* distribute the modified code.  Modifications are to be distributed
C* as patches to released version.
C*
C* This software is provided "as is" without express or implied warranty.
C*
C* Send your comments or suggestions to
C*  J.Heinen@kfa-juelich.de.
C*
C*

$IF DEFINED(F77L3)

	SUBROUTINE GKDGA (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     *      LC, CHARS)
C*  GKS logical device driver for IBM PC w\ Lahey FORTRAN graphics library

	INTEGER LC, LR1, LR2
	INTEGER FCTID, DX, DY, DIMX, IA(3)
	REAL R1(3), R2(3)
	CHARACTER*(*) CHARS

	EXTERNAL GKDPL, GKDSCI, GKDFA

	INTEGER LTYPE, ICOLOR
	INTEGER ERRIND, TNR, STYLE

	INTEGER IARRAY(9)
	REAL RARRAY(7)
	CHARACTER STRING*40

	INTEGER NUMX, NUMY, NUMC
	REAL YRES, MAXX, MAXY

	INTEGER ISTAT, NP
	REAL X, Y
	
	SAVE

C*  include GKS symbol definitions
	INCLUDE 'gksdefs.i'

C*  connection and type
	INTEGER CONID
C*  workstation state
	INTEGER STATE
C*  workstation transformation
	REAL WINDOW(4), VIEWPT(4)

	GOTO (999, 999,   2,   3,   4,   5,   6, 999,   8, 999, 
     *         10, 999,  12,  13,  14,  15,  16, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999,  48, 999, 
     *        999, 999, 999, 999,  54,  55, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999,  81,  82, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999) FCTID+1

	IF (LR1 .LT. 0 .OR. LR2 .LT. 0 .OR. LC .LT. 0)
     *      CALL GKDERR ('GKS: internal inconsistency')
	GOTO 999

C*  open workstation
    2   CONTINUE
	CONID = IA(2)

	I1 = 0
	I2 = 1
	IMODE = 0
	CALL PLOTS (I1, I2, IMODE)
	CALL FACTOR (1.0)

	CALL GRINFO (IARRAY, RARRAY, STRING)
	YRES = 1.0/IARRAY(2)
	MAXX = RARRAY(2)
	MAXY = RARRAY(3)
	NUMX = IARRAY(1)
	NUMY = IARRAY(2)
	NUMC = IARRAY(4)

C*  initialize colormap
	CALL GKDICM (NUMC)

C*  set default workstation window
	WINDOW(1) = 0.
	WINDOW(2) = 1.
	WINDOW(3) = 0.
	WINDOW(4) = 1.

C*  set default workstation viewport
	VIEWPT(1) = 0.032
	VIEWPT(2) = 0.224
	VIEWPT(3) = 0.
	VIEWPT(4) = 0.192

C*  set up device transformation
	CALL GKDSDT (WINDOW, VIEWPT, MAXX, MAXY)
	GOTO 999

C*  close workstation
    3   CONTINUE
	CALL PLOT (0, 0, 999)
	GOTO 999

C*  activate workstation
    4   CONTINUE
	STATE = GACTIV
	GOTO 999

C*  deactivate workstation
    5   CONTINUE
	STATE = GINACT
	GOTO 999

C*  clear workstation
    6   CONTINUE
	CALL PLOT (0, 0, -999)
	GOTO 999

C*  update workstation
    8   CONTINUE
	GOTO 999

C*  message
   10   CONTINUE
	WRITE (CONID, *) CHARS
	GOTO 999

C*  polyline
   12   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
	  CALL GQCNTN (ERRIND, TNR)
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
C*  inquire current linewidth, and color index
	  CALL GQLN (ERRIND, LTYPE)
	  CALL GQPLCI (ERRIND, ICOLOR)
	  CALL GKDSCI (ICOLOR)
	  CALL GKDPL (IA(1), R1, R2, LTYPE, TNR)
	END IF
	GOTO 999

C*  polymarker
   13   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
C*  inquire polymarker color index
	  CALL GQPMCI (ERRIND, ICOLOR)
	  CALL GKDSCI (ICOLOR)
	  CALL GPOLMK (IA(1), R1, R2, GKDPL)
	END IF
	GOTO 999

C*  text
   14   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
	  CALL GQTXCI (ERRIND, ICOLOR)
	  CALL GKDSCI (ICOLOR)
	  CALL GSIMTX (R1, R2, CHARS, GKDPL, GKDFA)
	END IF
	GOTO 999

C*  fill area
   15   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
	  CALL GQCNTN (ERRIND, TNR)
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
	  CALL GQFACI (ERRIND, ICOLOR)
	  CALL GKDSCI (ICOLOR)
C*  inquire current fill area interior style
	  CALL GQFAIS (ERRIND, STYLE)
	  IF (STYLE .NE. GSOLID) THEN
	    CALL GFILLA (IA(1), R1, R2, GKDPL, YRES)
	  ELSE
	    CALL GKDFA (IA(1), R1, R2, TNR)
	  END IF
	END IF
	GOTO 999

C*  cell array
   16   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
	  CALL GCELLA (R1(1), R1(2), R2(1), R2(2), DX, DY, 
     *      DIMX, IA, GKDSCI, GKDFA)
	END IF
	GOTO 999

C*  set color representation
   48   CONTINUE
	I = IA(2)
	IF (I.GE.0 .AND. I.LE.8) THEN
	  R = R1(1)
	  G = R1(2)
	  B = R1(3)
	  CALL GKDSCR (I, R, G, B)
	END IF
	GOTO 999

C*  set workstation window
   54   CONTINUE
	WINDOW(1) = R1(1)
	WINDOW(2) = R1(2)
	WINDOW(3) = R2(1)
	WINDOW(4) = R2(2)
	CALL GKDSDT (WINDOW, VIEWPT, MAXX, MAXY)
	GOTO 999

C*  set workstation viewport
   55   CONTINUE
	VIEWPT(1) = R1(1)
	VIEWPT(2) = R1(2)
	VIEWPT(3) = R2(1)
	VIEWPT(4) = R2(2)
	CALL GKFVP (VIEWPT, 0.256, 0.192)
	CALL GKDSDT (WINDOW, VIEWPT, MAXX, MAXY)
	GOTO 999

C*  request locator
   81   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
	  CALL GKDGIN (NUMX, NUMY, ISTAT, X, Y)
	  IA(1) = ISTAT
	  R1(1) = X
	  R2(1) = Y
	END IF
	GOTO 999

C*  request stroke
   82   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
	  NP = IA(3)
	  DO 8201 I = 1, NP
	    CALL GKDGIN (NUMX, NUMY, ISTAT, X, Y)
	    IF (ISTAT .NE. GOK) GOTO 8202
	    R1(I) = X
	    R2(I) = Y
	    J = I-1
	    IF (I .EQ. 1) THEN
	      J = 1
	      R1(2) = X
	      R2(2) = Y
	    END IF
	    CALL GKDPL (2, R1(J), R2(J), GLSOLI, 0)
 8201     CONTINUE
 8202     IF (I .GT. 2) ISTAT = GOK
	  IA(1) = ISTAT
	  IA(3) = I-1
	END IF

  999   RETURN
	END


	SUBROUTINE GKDPL (N, PX, PY, LTYPE, TNR)
C*  polyline

	REAL PX(1), PY(1)
	INTEGER N, LTYPE, TNR

	EXTERNAL GKDMOV, GKDDRW

C*  switch to the polyline utility
	CALL GPOLIN (N, PX, PY, LTYPE, TNR, GKDMOV, GKDDRW)

	RETURN
	END


	SUBROUTINE GKDMOV (X, Y)
C*  move

	EXTERNAL GKDVEC

	REAL X, Y

C*  switch to the move utility
	CALL GMOVE (X, Y, GKDVEC)

	RETURN
	END


	SUBROUTINE GKDDRW (X, Y)
C*  draw

	EXTERNAL GKDVEC, GKDCOD
	REAL X, Y

C*  switch to the dashed line generator
	CALL GDASH (X, Y, GKDVEC, GKDCOD)

	RETURN
	END


	SUBROUTINE GKDFA (N, PX, PY, TNR)
C*  fill area

	REAL PX(1), PY(1)
	INTEGER N, TNR
	
	EXTERNAL GKDBFA, GKDCFA

C*  switch to the polyline utility
	CALL GPOLIN (N, PX, PY, 0, TNR, GKDBFA, GKDCFA)
	CALL GKDEFA

	RETURN
	END


	SUBROUTINE GKDSDT (WN, VP, MAXX, MAXY)
C*  set up device transformation

	INTEGER MAXBUF
	PARAMETER (MAXBUF = 2000)

	REAL WN(4), VP(4), MAXX, MAXY
	INTEGER NUMX, NUMY, ISTAT
	REAL X, Y

	REAL A, B, C, D, E, F, G, H
	INTEGER INTARY(9), IX, IY
	INTEGER BUTTON, MOUSEX, MOUSEY

	REAL RX, RY, SCALEX, SCALEY
	INTEGER N
	REAL XARRAY(MAXBUF), YARRAY(MAXBUF)

	SAVE

C*  include GKS symbol definitions
	INCLUDE 'gksdefs.i'

	DATA MOUSEX /-1/, MOUSEY /-1/

	SCALEX = MAXX
	SCALEY = MAXY
	
	E = (VP(2)-VP(1))/(WN(2)-WN(1))
	F = MAXX/0.256
	G = (VP(4)-VP(3))/(WN(4)-WN(3))
	H = MAXY/0.192

	A = E*F
	B = F*(VP(1)-WN(1)*E)
	C = G*H
	D = H*(VP(3)-WN(3)*G)

	RETURN


	ENTRY GKDVEC (X, Y)
C*  re-initialize vector drawing sequence

	RX = A*X+B
	RY = C*Y+D
	CALL PLOT (RX, RY, 3)
	
	RETURN


	ENTRY GKDCOD (X, Y)

C*  device transformation
	RX = A*X+B
	RY = C*Y+D
	CALL PLOT (RX, RY, 2)
	
	RETURN


	ENTRY GKDBFA (X, Y)
C*  begin fill area
	
	N = 0

	ENTRY GKDCFA (X, Y)
C*  continue fill area

	IF (N .EQ. MAXBUF)
     *      CALL GKDERR ('GKS: number of points exceeded')

	N = N + 1
	XARRAY(N) = A*X+B
	YARRAY(N) = C*Y+D

	RETURN


	ENTRY GKDEFA
C*  end fill area

	CALL FILL (N, XARRAY, YARRAY)

	RETURN


	ENTRY GKDGIN (NUMX, NUMY, ISTAT, X, Y)
C*  graphic input

C*  check mouse
	INTARY(1) = 0
	CALL INTRUP (INTARY, 51)
	IF (INTARY(1) .NE. 65535)
     *      CALL GKDERR ('GKS: mouse not installed')

C*  show mouse
	INTARY(1) = 1
	CALL INTRUP (INTARY, 51)

	IF (MOUSEX .GE. 0 .AND. MOUSEY .GE. 0) THEN
C*  set mouse position
	    INTARY(1) = 4
	    INTARY(3) = MOUSEX
	    INTARY(4) = MOUSEY
	    CALL INTRUP (INTARY, 51)
	END IF

	BUTTON = 0
C*  get mouse position
   1    INTARY(1) = 3
	INTARY(2) = 1+2+4
	CALL INTRUP (INTARY, 51)
	IF (BUTTON .EQ. 0) THEN
	    BUTTON = INTARY(2)
	    MOUSEX = INTARY(3)
	    MOUSEY = INTARY(4)
	    GOTO 1
	ELSE IF (INTARY(2) .NE. 0) THEN
	    GOTO 1
	END IF

C*  check mouse button
	IF (BUTTON .EQ. 1) THEN
	    ISTAT = GOK
	    IX = MOUSEX
	    IY = MOUSEY
	ELSE
	    ISTAT = GNONE
	END IF

	IF (ISTAT .EQ. GOK) THEN
	    RX = SCALEX * IX / NUMX
	    RY = SCALEY * (NUMY-IY) / NUMY
	    X = (RX-B)/A
	    Y = (RY-D)/C
	END IF

C*  hide mouse
	INTARY(1) = 2
	CALL INTRUP (INTARY, 51)

	RETURN
	END


	SUBROUTINE GKDICM (NUMC)
C*  initialize colormap

	INTEGER NUMC, ICOLOR
	REAL R, G, B

	INTEGER NUMCOL, COLI
	INTEGER ML(0:11), PALETT(0:15), CMAP(0:15)

	INTEGER CI

	SAVE
	
C*  colormap          D  W  R  G  B  C  Y  M  R* G* B* C* Y* M*
	DATA PALETT / 0,15,12,10, 9,11,14,13, 4, 2, 1, 3, 6, 5, 8, 7 /
C*  map location
	DATA ML / 4,10, 7,13, 2, 8, 6,12, 3, 9, 5,11 /
	
	NUMCOL = NUMC
	DO 1, I = 0, 15
	    CMAP(I) = MIN(PALETT(I), NUMCOL-1)
   1    CONTINUE
   
	COLI = 1
	CALL NEWPEN (CMAP(COLI))

	RETURN


	ENTRY GKDSCI (ICOLOR)
C*  set color index

	IF (ICOLOR .GE. 8) THEN
	    CI = ML(MOD(ICOLOR-8,72)/6)
	ELSE
	    CI = ICOLOR
	END IF

	IF (CI .NE. COLI) THEN
	    CALL NEWPEN (CMAP(CI))
	    COLI = CI
	END IF

	RETURN


	ENTRY GKDSCR (ICOLOR, R, G, B)
C*  set color representation

	IF (ICOLOR .LT. 8 .AND. ICOLOR .LT. NUMCOL) THEN
	    CALL GKRMAP (ICOLOR, R, G, B, 8, PALETT, CMAP)
	END IF

	RETURN
	END


	SUBROUTINE GKDERR (ERRMSG)
C*  produce an error message and exit

	CHARACTER*(*) ERRMSG

	CALL PLOT (0, 0, 999)
	WRITE (*,*) ERRMSG
	STOP
	
	END

$ELSE

$IF .NOT. DEFINED(_MSFORTRAN_)
$DEFINE _MSFORTRAN_ = 100
$ENDIF

$IF _MSFORTRAN_ .LT. 300

	INCLUDE 'FGRAPH.FI'

	INTERFACE TO SUBROUTINE INT86 [C] (INTR, INREG, OUTREG)
        STRUCTURE /WREGS/
            INTEGER*2 AX, BX, CX, DX, SI, DI, CFLAG
        END STRUCTURE
        STRUCTURE /BREGS/
            INTEGER*1 AL, AH, BL, BH, CL, CH, DL, DH
        END STRUCTURE
        STRUCTURE /REGS/
            UNION
                MAP
                    RECORD /WREGS/ X
                END MAP
                MAP
                    RECORD /BREGS/ H
                END MAP
            END UNION
        END STRUCTURE
	INTEGER*2 INTR [VALUE]
	RECORD /REGS/ INREG [REFERENCE]
	RECORD /REGS/ OUTREG [REFERENCE]
	END
$ENDIF


	SUBROUTINE GKDGA (FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     *      LC, CHARS)
C*  GKS logical device driver for IBM PC w\ Microsoft FORTRAN graphics library

	INTEGER LC, LR1, LR2
	INTEGER FCTID, DX, DY, DIMX, IA(3)
	REAL R1(3), R2(3)
	CHARACTER*(*) CHARS

	EXTERNAL GKDPL, GKDSCI, GKDFA

	INTEGER I, LTYPE
	INTEGER ICOLOR
	INTEGER ERRIND, TNR
	REAL R, G, B
	INTEGER NUMX, NUMY, NUMC

	INTEGER ISTAT, NP
	REAL X, Y
	
	SAVE

C*  include GKS symbol definitions
	INCLUDE 'gksdefs.i'

	INCLUDE 'FGRAPH.FD'

C*  connection and type
	INTEGER CONID
C*  workstation state
	INTEGER STATE
C*  workstation transformation
	REAL WINDOW(4), VIEWPT(4)

	INTEGER*2 STATUS

C*  video/window configuration
	RECORD /VIDEOCONFIG/ SCREEN

	GOTO (999, 999,   2,   3,   4,   5,   6, 999,   8, 999, 
     *         10, 999,  12,  13,  14,  15,  16, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999,  48, 999, 
     *        999, 999, 999, 999,  54,  55, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999,  81,  82, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999) FCTID+1

	IF (LR1 .LT. 0 .OR. LR2 .LT. 0 .OR. LC .LT. 0)
     *      CALL GKDERR ('GKS: internal inconsistency')
	GOTO 999

C*  open workstation
    2   CONTINUE
	CONID = IA(2)

$IF _MSFORTRAN_ .LT. 300
	STATUS = SETVIDEOMODE ($MAXRESMODE)
	IF (STATUS .EQ. 0) CALL GKDERR ('GKS: cannot set graphics mode')
$ENDIF
	CALL GETVIDEOCONFIG (SCREEN)
	NUMX = SCREEN.NUMXPIXELS
	NUMY = SCREEN.NUMYPIXELS
	NUMC = SCREEN.NUMCOLORS

C*  initialize colormap
	CALL GKDICM (NUMC)

C*  set default workstation window
	WINDOW(1) = 0.
	WINDOW(2) = 1.
	WINDOW(3) = 0.
	WINDOW(4) = 1.

C*  set default workstation viewport
	VIEWPT(1) = 0.032
	VIEWPT(2) = 0.224
	VIEWPT(3) = 0.
	VIEWPT(4) = 0.192

C*  set up device transformation
	CALL GKDSDT (WINDOW, VIEWPT, NUMX, NUMY)
	GOTO 999

C*  close workstation
    3   CONTINUE
	STATUS = SETVIDEOMODE ($DEFAULTMODE)
	GOTO 999

C*  activate workstation
    4   CONTINUE
	STATE = GACTIV
	GOTO 999

C*  deactivate workstation
    5   CONTINUE
	STATE = GINACT
	GOTO 999

C*  clear workstation
    6   CONTINUE
	CALL CLEARSCREEN ($GCLEARSCREEN)
	GOTO 999

C*  update workstation
    8   CONTINUE
	GOTO 999

C*  message
   10   CONTINUE
	WRITE (CONID, *) CHARS
	GOTO 999

C*  polyline
   12   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
	  CALL GQCNTN (ERRIND, TNR)
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
C*  inquire current linewidth, and color index
	  CALL GQLN (ERRIND, LTYPE)
	  CALL GQPLCI (ERRIND, ICOLOR)
	  CALL GKDSCI (ICOLOR)
	  CALL GKDPL (IA(1), R1, R2, LTYPE, TNR)
	END IF
	GOTO 999

C*  polymarker
   13   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
C*  inquire polymarker color index
	  CALL GQPMCI (ERRIND, ICOLOR)
	  CALL GKDSCI (ICOLOR)
	  CALL GPOLMK (IA(1), R1, R2, GKDPL)
	END IF
	GOTO 999

C*  text
   14   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
	  CALL GQTXCI (ERRIND, ICOLOR)
	  CALL GKDSCI (ICOLOR)
	  CALL GSIMTX (R1, R2, CHARS, GKDPL, GKDFA)
	END IF
	GOTO 999

C*  fill area
   15   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
	  CALL GQCNTN (ERRIND, TNR)
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
	  CALL GQFACI (ERRIND, ICOLOR)
	  CALL GKDSCI (ICOLOR)
	  CALL GKDFA (IA(1), R1, R2, TNR)
	END IF
	GOTO 999

C*  cell array
   16   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
	  CALL GCELLA (R1(1), R1(2), R2(1), R2(2), DX, DY, 
     *      DIMX, IA, GKDSCI, GKDFA)
	END IF
	GOTO 999

C*  set color representation
   48   CONTINUE
	I = IA(2)
	IF (I.GE.0 .AND. I.LE.8) THEN
	  R = R1(1)
	  G = R1(2)
	  B = R1(3)
	  CALL GKDSCR (I, R, G, B)
	END IF
	GOTO 999

C*  set workstation window
   54   CONTINUE
	WINDOW(1) = R1(1)
	WINDOW(2) = R1(2)
	WINDOW(3) = R2(1)
	WINDOW(4) = R2(2)
	CALL GKDSDT (WINDOW, VIEWPT, NUMX, NUMY)
	GOTO 999

C*  set workstation viewport
   55   CONTINUE
	VIEWPT(1) = R1(1)
	VIEWPT(2) = R1(2)
	VIEWPT(3) = R2(1)
	VIEWPT(4) = R2(2)
	CALL GKFVP (VIEWPT, 0.256, 0.192)
	CALL GKDSDT (WINDOW, VIEWPT, NUMX, NUMY)
	GOTO 999

$IF _MSFORTRAN_ .LT. 300

C*  request locator
   81   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
	  CALL GKDGIN (ISTAT, X, Y)
	  IA(1) = ISTAT
	  R1(1) = X
	  R2(1) = Y
	END IF
	GOTO 999

C*  request stroke
   82   CONTINUE
	IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
	  CALL GSDT (WINDOW, VIEWPT)
	  NP = IA(3)
	  DO 8201 I = 1, NP
	    CALL GKDGIN (ISTAT, X, Y)
	    IF (ISTAT .NE. GOK) GOTO 8202
	    R1(I) = X
	    R2(I) = Y
	    J = I-1
	    IF (I .EQ. 1) THEN
	      J = 1
	      R1(2) = X
	      R2(2) = Y
	    END IF
	    CALL GKDPL (2, R1(J), R2(J), GLSOLI, 0)
 8201     CONTINUE
 8202     IF (I .GT. 2) ISTAT = GOK
	  IA(1) = ISTAT
	  IA(3) = I-1
	END IF
$ELSE
   81   CONTINUE
   82   CONTINUE
	IA(1) = GNONE
$ENDIF

  999   RETURN
	END


	SUBROUTINE GKDPL (N, PX, PY, LTYPE, TNR)
C*  polyline

	REAL PX(1), PY(1)
	INTEGER N, LTYPE, TNR

	EXTERNAL GKDMOV, GKDDRW

C*  switch to the polyline utility
	CALL GPOLIN (N, PX, PY, LTYPE, TNR, GKDMOV, GKDDRW)

	RETURN
	END


	SUBROUTINE GKDMOV (X, Y)
C*  move

	EXTERNAL GKDVEC

	REAL X, Y

C*  switch to the move utility
	CALL GMOVE (X, Y, GKDVEC)

	RETURN
	END


	SUBROUTINE GKDDRW (X, Y)
C*  draw

	EXTERNAL GKDVEC, GKDCOD
	REAL X, Y

C*  switch to the dashed line generator
	CALL GDASH (X, Y, GKDVEC, GKDCOD)

	RETURN
	END


	SUBROUTINE GKDFA (N, PX, PY, TNR)
C*  fill area

	REAL PX(1), PY(1)
	INTEGER N, TNR
	
	EXTERNAL GKDBFA, GKDCFA

C*  switch to the polyline utility
	CALL GPOLIN (N, PX, PY, 0, TNR, GKDBFA, GKDCFA)
	CALL GKDEFA

	RETURN
	END


	SUBROUTINE GKDSDT (WN, VP, NUMX, NUMY)
C*  set up device transformation

	INTEGER MAXBUF
	PARAMETER (MAXBUF = 2000)

	INTEGER ZERO
	PARAMETER (ZERO = 64)

	REAL WN(4), VP(4)
	INTEGER NUMX, NUMY
	REAL X, Y
	INTEGER ISTAT

	REAL A, B, C, D, E, F, G, H
	INTEGER ERRIND, STYLE, INDEX, PA(0:32)
	INTEGER MAXX, MAXY

	SAVE

C*  include GKS symbol definitions
	INCLUDE 'gksdefs.i'

	INCLUDE 'FGRAPH.FD'

        STRUCTURE /WREGS/
            INTEGER*2 AX, BX, CX, DX, SI, DI, CFLAG
        END STRUCTURE
        STRUCTURE /BREGS/
            INTEGER*1 AL, AH, BL, BH, CL, CH, DL, DH
        END STRUCTURE
        STRUCTURE /REGS/
            UNION
                MAP
                    RECORD /WREGS/ X
                END MAP
                MAP
                    RECORD /BREGS/ H
                END MAP
            END UNION
        END STRUCTURE

	INTEGER*2 INTR, STATUS, IX, IY, NP, CONTRL
	INTEGER*1 OLDSTYLE(8), PATT(8)
	RECORD /XYCOORD/ XY, POINTS(MAXBUF)
	INTEGER*2 BUTTON, MOUSEX, MOUSEY, CROSSX, CROSSY
	INTEGER*2 ASCII, SCANC
	RECORD /REGS/ INREG, OUTREG

	DATA MOUSEX /-1/, MOUSEY /-1/

	MAXX = NUMX
	MAXY = NUMY
	
	E = (VP(2)-VP(1))/(WN(2)-WN(1))
	F = (NUMX-1)/0.256
	G = (VP(4)-VP(3))/(WN(4)-WN(3))
	H = (NUMY-1)/0.192

	A = E*F
	B = F*(VP(1)-WN(1)*E)
	C = G*H
	D = H*(VP(3)-WN(3)*G)

	RETURN


	ENTRY GKDVEC (X, Y)
C*  re-initialize vector drawing sequence

	IX = NINT(A*X+B)
	IY = MAXY-NINT(C*Y+D)
	CALL MOVETO (IX, IY, XY)
	
	RETURN


	ENTRY GKDCOD (X, Y)

C*  device transformation
	IX = NINT(A*X+B)
	IY = MAXY-NINT(C*Y+D)
	STATUS = LINETO (IX, IY)
	
	RETURN


	ENTRY GKDBFA (X, Y)
C*  begin fill area
	
	NP = 0
C*  inquire current fill area interior style
	CALL GQFAIS (ERRIND, STYLE)
C*  inquire current fill area style index
	CALL GQFASI (ERRIND, INDEX)

	CALL GETFILLMASK (OLDSTYLE)
	IF (STYLE .EQ. GHOLLO) THEN
	    CONTRL = $GBORDER
	ELSE
	    CONTRL = $GFILLINTERIOR
	    IF (STYLE .EQ. GSOLID) THEN 
		INDEX = 0
	    ELSE IF (STYLE .EQ. GHATCH) THEN 
		INDEX = INDEX + 108
	    END IF     
	    CALL GKQPA (INDEX, PA)
	    IF (PA(0) .EQ. 4) THEN
	        DO 10, I = 5,8
		    PA(I) = PA(I-4)
   10		CONTINUE
	    END IF
	    DO 20, I = 1,8
		PATT(I) = NOT(PA(I))
   20       CONTINUE
	    CALL SETFILLMASK (PATT)
	END IF


	ENTRY GKDCFA (X, Y)
C*  continue fill area

	IF (NP .EQ. MAXBUF)
     *      CALL GKDERR ('GKS: number of points exceeded')

	NP = NP + 1
	POINTS(NP).XCOORD = NINT(A*X+B)
	POINTS(NP).YCOORD = MAXY-NINT(C*Y+D)

	RETURN


	ENTRY GKDEFA
C*  end fill area

	STATUS = POLYGON (CONTRL, POINTS, NP)
	CALL SETFILLMASK (OLDSTYLE)

	RETURN


$IF _MSFORTRAN_ .LT. 300

	ENTRY GKDGIN (ISTAT, X, Y)
C*  graphic input

C*  check mouse
	INTR = 51
	INREG.X.AX = 0
	CALL INT86(INTR, INREG, OUTREG)
	IF (OUTREG.X.AX .NE. -1)
     *      CALL GKDERR ('GKS: mouse not installed')

	IF (MOUSEX .GE. 0 .AND. MOUSEY .GE. 0) THEN
C*  set mouse position
	    INTR = 51
	    INREG.X.AX = 4
	    INREG.X.CX = MOUSEX
	    INREG.X.DX = MOUSEY
	    CALL INT86(INTR, INREG, OUTREG)
	END IF

	STATUS = SETWRITEMODE ($GXOR)

	CALL MOVETO (MOUSEX, 0, XY)
	STATUS = LINETO (MOUSEX, MAXY)
	CALL MOVETO (0, MOUSEY, XY)
	STATUS = LINETO (MAXX, MOUSEY)
	CROSSX = MOUSEX
	CROSSY = MOUSEY

	BUTTON = 0
   1    IF (MOUSEX .NE. CROSSX .OR. MOUSEY .NE. CROSSY) THEN
	    CALL MOVETO (CROSSX, 0, XY)
	    STATUS = LINETO (CROSSX, MAXY)
	    CALL MOVETO (0, CROSSY, XY)
	    STATUS = LINETO (MAXX, CROSSY)
C*  update crosshair cursor
	    CALL MOVETO (MOUSEX, 0, XY)
	    STATUS = LINETO (MOUSEX, MAXY)
	    CALL MOVETO (0, MOUSEY, XY)
	    STATUS = LINETO (MAXX, MOUSEY)
	    CROSSX = MOUSEX
	    CROSSY = MOUSEY
	END IF

C*  check keyboard buffer
	INTR = 33
        INREG.H.AH = 11
	CALL INT86(INTR, INREG, OUTREG)

	IF (OUTREG.H.AL .NE. 0) THEN
	    INTR = 22
	    INREG.H.AH = 0
	    CALL INT86(INTR, INREG, OUTREG)

	    ASCII = OUTREG.H.AL
	    IF (ASCII .EQ. 0) THEN
C*  non-ASCII character
		SCANC = OUTREG.H.AH
C*  check for cursor key
		IF (SCANC .EQ. 72) THEN
		    MOUSEY = CROSSY - 1
		ELSE IF (SCANC .EQ. 75) THEN
		    MOUSEX = CROSSX - 1
		ELSE IF (SCANC .EQ. 77) THEN
		    MOUSEX = CROSSX + 1
		ELSE IF (SCANC .EQ. 80) THEN
		    MOUSEY = CROSSY + 1
		END IF
C*  set mouse position
	        INTR = 51
		INREG.X.AX = 4
	        INREG.X.CX = MOUSEX
	        INREG.X.DX = MOUSEY
	        CALL INT86(INTR, INREG, OUTREG)
                GOTO 1

	    ELSE IF (ASCII .EQ. 4 .OR. ASCII .EQ. 26) THEN
C*  end-of-file character (CTRL/D, CTRL/Z)
		BUTTON = 2
	    ELSE
		BUTTON = 1
	    END IF
	ELSE
C*  get mouse position
	    INTR = 51
	    INREG.X.AX = 3
	    INREG.X.BX = 1+2+4
	    CALL INT86(INTR, INREG, OUTREG)

            IF (BUTTON .EQ. 0) THEN
                BUTTON = OUTREG.X.BX
                MOUSEX = OUTREG.X.CX
                MOUSEY = OUTREG.X.DX
                GOTO 1
            ELSE IF (OUTREG.X.BX .NE. 0) THEN
                GOTO 1
            END IF
	END IF

C*  check mouse button
        IF (BUTTON .EQ. 1) THEN
	    ISTAT = GOK
	    IX = MOUSEX
	    IY = MOUSEY
	ELSE
	    ISTAT = GNONE
	END IF

	IF (ISTAT .EQ. GOK) THEN
	    X = (IX-B)/A
	    Y = (MAXY-IY-D)/C
	END IF

	CALL MOVETO (CROSSX, 0, XY)
	STATUS = LINETO (CROSSX, MAXY)
	CALL MOVETO (0, CROSSY, XY)
	STATUS = LINETO (MAXX, CROSSY)

	STATUS = SETWRITEMODE ($GPSET)

	RETURN
$ENDIF
	END


	SUBROUTINE GKDICM (NUMC)
C*  initialize colormap

	INTEGER NUMC, ICOLOR
	REAL R, G, B

	INTEGER NUMCOL, COLI
	INTEGER ML(0:11), PALETT(0:15), CMAP(0:15)

	INTEGER CI
	INTEGER*2 STATUS

	SAVE
	
	INCLUDE 'FGRAPH.FD'

C*  colormap          W  D  R  G  B  C  Y  M  R* G* B* C* Y* M*
	DATA PALETT / 0,15,12,10, 9,11,14,13, 4, 2, 1, 3, 6, 5, 8, 7 /
C*  map location
	DATA ML / 4,10, 7,13, 2, 8, 6,12, 3, 9, 5,11 /

	NUMCOL = NUMC        
	IF (NUMCOL .GE. 16) THEN
	    STATUS = REMAPPALETTE ( 0, #3F3F3F)
	    STATUS = REMAPPALETTE (15, #000000)
	END IF

	DO 1, I = 0, 15
	    CMAP(I) = MIN(PALETT(I), NUMCOL-1)
   1    CONTINUE

	COLI = 1
	STATUS = SETCOLOR (CMAP(COLI))

	RETURN


	ENTRY GKDSCI (ICOLOR)
C*  set color index

	IF (ICOLOR .GE. 8) THEN
	    CI = ML(MOD(ICOLOR-8,72)/6)
	ELSE
	    CI = ICOLOR
	END IF

	IF (CI .NE. COLI) THEN
	    STATUS = SETCOLOR (CMAP(CI))
	    COLI = CI
	END IF

	RETURN


	ENTRY GKDSCR (ICOLOR, R, G, B)
C*  set color representation

	IF (ICOLOR .LT. 8 .AND. ICOLOR .LT. NUMCOL) THEN
	    CALL GKRMAP (ICOLOR, R, G, B, 8, PALETT, CMAP)
	END IF

	RETURN
	END


	SUBROUTINE GKDERR (ERRMSG)
C*  produce an error message and exit

	CHARACTER*(*) ERRMSG

	INCLUDE 'FGRAPH.FD'

	INTEGER*2 STATUS

	STATUS = SETVIDEOMODE ($DEFAULTMODE)
	WRITE (*,*) ERRMSG
	STOP
	
	END

$ENDIF
