C*
C* Copyright @ 1984 - 1995   Josef Heinen
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

        SUBROUTINE GKDBM(FCTID, DX, DY, DIMX, IA, LR1, R1, LR2, R2,
     *    LC, CHARS)
C*  GKS logical device driver for bitmap devices

        INTEGER LC, LR1, LR2
        INTEGER FCTID, DX, DY, DIMX, IA(*)
        REAL R1(*), R2(*)
        CHARACTER*(*) CHARS

        EXTERNAL GKBPL, GKBSPL, GKBSFA

        INTEGER BORDER, ICOLOR, LTYPE, NCHARS
        REAL YRES 
        INTEGER ERRIND, TNR, CLEAR

        SAVE

C*  include GKS symbol definitons
        INCLUDE 'gksdefs.i'

C*  connection identifier, workstation type
        INTEGER CONID, WTYPE
C*  workstation state
        INTEGER STATE
C*  workstation transformation
        REAL WINDOW(4), VIEWPT(4)
C*  resolution (dots/inch)
        INTEGER DPI

        GOTO (999, 999,   2,   3,   4,   5,   6, 999, 999, 999, 
     *        999, 999,  12,  13,  14,  15, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999,  54,  55, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 
     *        999, 999, 999, 999, 999, 999, 999) FCTID + 1

        IF (DX .LT. 0 .OR. DY .LT. 0 .OR. DIMX .LT. 0 .OR.
     *      LR1 .LT. 0 .OR. LR2 .LT. 0 .OR. LC .LT. 0)
     *      STOP 'GKS: internal inconsistency'
        GOTO 999

C*  open workstation
    2   CONTINUE
        CONID = IA(2)
        WTYPE = IA(3)

C*  set up connection
        CALL GKBSCO(CONID)

C*  set default workstation window
        WINDOW(1) = 0.
        WINDOW(2) = 1.
        WINDOW(3) = 0.
        WINDOW(4) = 1.

C*  set default workstation viewport
        VIEWPT(1) = 0
        VIEWPT(2) = 0.2032
        VIEWPT(3) = 0.
        VIEWPT(4) = VIEWPT(2)

        IF (WTYPE .EQ. 103) THEN
          DPI = 72
        ELSE
          DPI = 75
        END IF
C*  set up device transformation
        CALL GKBSDT(WINDOW, VIEWPT, DPI)

        CLEAR = 1
        GOTO 999

C*  close workstation
    3   CONTINUE
        IF (CLEAR .NE. 0) THEN
          CALL GKBCBM
          CLEAR = 0
        END IF
        CALL GKBWBM
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
        CLEAR = 1
        GOTO 999

C*  polyline
   12   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (CLEAR .NE. 0) THEN
            CALL GKBCBM
            CLEAR = 0
          END IF
C*  inquire current normalization transformation
          CALL GQCNTN(ERRIND, TNR)
C*  set up device transformation
          CALL GSDT(WINDOW, VIEWPT)
C*  inquire current line type
          CALL GQLN(ERRIND, LTYPE)
C*  inquire polyline color index
          CALL GQPLCI(ERRIND, ICOLOR)
          CALL GKBSCI(ICOLOR)
          CALL GKBPL(IA(1), R1, R2, LTYPE, TNR)
        END IF
        GOTO 999

C*  polymarker
   13   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (CLEAR .NE. 0) THEN
            CALL GKBCBM
            CLEAR = 0
          END IF
C*  set up device transformation
          CALL GSDT(WINDOW, VIEWPT)
C*  inquire polymarker color index
          CALL GQPMCI(ERRIND, ICOLOR)
          CALL GKBSCI(ICOLOR)
          CALL GPOLMK(IA(1), R1, R2, GKBPL)
        END IF
        GOTO 999

C*  text
   14   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (CLEAR .NE. 0) THEN
            CALL GKBCBM
            CLEAR = 0
          END IF
C*  set up device transformation
          CALL GSDT(WINDOW, VIEWPT)
C*  inquire text color index
          CALL GQTXCI(ERRIND, ICOLOR)
          CALL GKBSCI(ICOLOR)
          NCHARS = LEN(CHARS)
          CALL GSIMTX(R1, R2, NCHARS, CHARS, GKBSPL, GKBSFA)
        END IF
        GOTO 999

C*  fill area
   15   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          IF (CLEAR .NE. 0) THEN
            CALL GKBCBM
            CLEAR = 0
          END IF
C*  inquire current normalization transformation
          CALL GQCNTN(ERRIND, TNR)
C*  set up device transformation
          CALL GSDT(WINDOW, VIEWPT)
C*  inquire fill area color index
          CALL GQFACI(ERRIND, ICOLOR)
          CALL GKBSCI(ICOLOR)
          IF (ICOLOR .LE. 1) THEN
            YRES = 1.0 / 600.0
            CALL GFILLA(IA(1), R1, R2, TNR, GKBPL, YRES)
          ELSE
            BORDER = 0
            CALL GKBPL(IA(1), R1, R2, BORDER, TNR)
          END IF
        END IF
        GOTO 999

C*  set workstation window
   54   CONTINUE
        WINDOW(1) = R1(1)
        WINDOW(2) = R1(2)
        WINDOW(3) = R2(1)
        WINDOW(4) = R2(2)
        CALL GKBSDT(WINDOW, VIEWPT, DPI)
        GOTO 999

C*  set workstation viewport
   55   CONTINUE
        VIEWPT(1) = R1(1)
        VIEWPT(2) = R1(2)
        VIEWPT(3) = R2(1)
        VIEWPT(4) = R2(2)
        CALL GKFVP(VIEWPT, 0.254, 0.2032)
        CALL GKBSDT(WINDOW, VIEWPT, DPI)
        GOTO 999

  999   RETURN
        END


        SUBROUTINE GKBPL(N, PX, PY, LTYPE, TNR)
C*  polyline

        REAL PX(N), PY(N)
        INTEGER N, LTYPE, TNR

        EXTERNAL GKBMOV, GKBDRW

C*  switch to the polyline utility
        CALL GPOLIN(N, PX, PY, LTYPE, TNR, GKBMOV, GKBDRW)

        RETURN
        END


        SUBROUTINE GKBMOV(X, Y)
C*  move
        REAL X, Y

        EXTERNAL GKBVEC

C*  switch to the move utility
        CALL GMOVE(X, Y, GKBVEC)

        RETURN
        END


        SUBROUTINE GKBDRW(X, Y)
C*  draw
        REAL X, Y

        EXTERNAL GKBVEC, GKBCOD

C*  switch to the dashed line generator
        CALL GDASH(X, Y, GKBVEC, GKBCOD)

        RETURN
        END


        SUBROUTINE GKBSPL(N, PX, PY, LTYPE, TNR)
C*  solid polyline

        REAL PX(N), PY(N)
        INTEGER N, LTYPE, TNR

        EXTERNAL GKBVEC, GKBCOD

C*  switch to the polyline utility
        CALL GPOLIN(N, PX, PY, LTYPE, TNR, GKBVEC, GKBCOD)

        RETURN
        END


        SUBROUTINE GKBFA(N, PX, PY, TNR)
C*  fill area

        INTEGER N
        REAL PX(N), PY(N)
        INTEGER TNR

        REAL YRES

        EXTERNAL GKBPL

C*  switch to the fill area utility
        YRES = 1.0 / 600.0
        CALL GFILLA(N, PX, PY, TNR, GKBPL, YRES)

        RETURN
        END


        SUBROUTINE GKBSFA(N, PX, PY, TNR)
C*  solid fill area

        INTEGER N
        REAL PX(N), PY(N)
        INTEGER TNR

        REAL YRES

        EXTERNAL GKBPL

C*  switch to the solid fill area utility
        YRES = 1.0 / 600.0
        CALL GKSFA(N, PX, PY, TNR, GKBPL, YRES)

        RETURN
        END


        SUBROUTINE GKBSDT(WN, VP, DPI)
C*  set up device transformation

        REAL WN(4), VP(4)
        INTEGER DPI, COLI
        REAL X, Y

        REAL A, B, C, D, E, F, G, H
        INTEGER VALUE, WIDTH, HEIGHT, X1, Y1, X2, Y2, TMP
        INTEGER DX, XINC, DY, YINC, XPLOT, YPLOT, RUNC
        LOGICAL LANDSC

        SAVE

        INCLUDE 'gksdefs.i'

        E = (VP(2) - VP(1)) / (WN(2) - WN(1))
        F = (DPI * 10 - 1) / 0.254
        G = (VP(4) - VP(3)) / (WN(4) - WN(3))
        H = (DPI * 8 - 1) / 0.2032

        A = E * F
        B = F * (VP(1) - WN(1) * E)
        C = G * H
        D = H * (VP(3) - WN(3) * G)

        WIDTH  = MAX(MIN(NINT(A * (WN(2) - WN(1))) + 1, DPI * 10), 100)
        HEIGHT = MAX(MIN(NINT(C * (WN(4) - WN(3))) + 1, DPI * 8), 100)

        LANDSC = (WIDTH .GT. DPI * 8)
        IF (LANDSC) THEN
          TMP = WIDTH
          WIDTH = HEIGHT
          HEIGHT = TMP
        ELSE
          C = -C
          D = HEIGHT - D - 1
        END IF

        CALL GKBSBM(WIDTH, HEIGHT)

        RETURN


        ENTRY GKBSCI(COLI)
C*  set color index

        IF (COLI .GE. 1) THEN
          VALUE = 1
        ELSE
          VALUE = 0
        END IF

        RETURN


        ENTRY GKBVEC(X, Y)
C*  re-initialize vector drawing sequence

        X1 = INT(A * X + B)
        Y1 = INT(C * Y + D)
        IF (LANDSC) THEN
          TMP = X1
          X1 = Y1
          Y1 = TMP
        END IF

        RETURN


        ENTRY GKBCOD(X, Y)

C*  device transformation
        X2 = INT(A * X + B)
        Y2 = INT(C * Y + D)
        IF (LANDSC) THEN
          TMP = X2
          X2 = Y2
          Y2 = TMP
        END IF

        RUNC = 0
        DX = ABS(X1 - X2)
        IF (X2 .GT. X1) XINC =  1
        IF (X2 .EQ. X1) XINC =  0
        IF (X2 .LT. X1) XINC = -1
        DY = ABS(Y1 - Y2)
        IF (Y2 .GT. Y1) YINC =  1
        IF (Y2 .EQ. Y1) YINC =  0
        IF (Y2 .LT. Y1) YINC = -1
        XPLOT = X1
        YPLOT = Y1

        IF (DX .GT. DY) THEN
C*  iterate x
          CALL GKBSP(XPLOT, YPLOT, VALUE)
    1     IF (XPLOT .NE. X2) THEN
            XPLOT = XPLOT + XINC
            RUNC = RUNC + DY
            IF (RUNC .GE. DX - RUNC) THEN
              YPLOT = YPLOT + YINC
              RUNC = RUNC - DX
            END IF
            CALL GKBSP(XPLOT, YPLOT, VALUE)
            GOTO 1
          END IF
        ELSE
C*  iterate y
          CALL GKBSP(XPLOT, YPLOT, VALUE)
    2     IF (YPLOT .NE. Y2) THEN
            YPLOT = YPLOT + YINC
            RUNC = RUNC + DX
            IF (RUNC .GE. DY - RUNC) THEN
              XPLOT = XPLOT +XINC
              RUNC = RUNC - DY
            END IF
            CALL GKBSP(XPLOT, YPLOT, VALUE)
            GOTO 2
          END IF
        END IF

        X1 = X2
        Y1 = Y2

        RETURN
        END


        SUBROUTINE GKBSCO(ID)
C*  set connection identifier

        INTEGER ID
        INTEGER X, Y, VALUE

        INTEGER SIZE
        PARAMETER (SIZE = 750 * 600 / 8)

        INTEGER CONID, WIDTH, HEIGHT, COLS, NBYTES
        CHARACTER*1 BITMAP(-10:SIZE)
        INTEGER J, K, MASK, BITS(0:7), LW, LH
        CHARACTER LF, BACK, W*8, H*8, HEADER*11

        SAVE

        DATA BITS /128, 64, 32, 16, 8, 4, 2, 1/

        CONID = ID

        LF = CHAR(10)
        BACK = CHAR(0)

        RETURN


        ENTRY GKBSBM(X, Y)
C*  set bitmap

        WIDTH = X
        HEIGHT = Y
        X = 8 * INT(X / 8.0 + 0.9)
        Y = 8 * INT(Y / 8.0 + 0.9)
        COLS = X
        NBYTES = MIN(X * Y / 8, SIZE)

        RETURN


        ENTRY GKBCBM
C*  clear bitmap

        DO 1, K = 1, NBYTES
          BITMAP(K) = BACK
    1   CONTINUE

        RETURN


        ENTRY GKBSP(X, Y, VALUE)
C*  set pixel

        J = Y * COLS + X
        K = J / 8 + 1
        IF (K .GE. 1 .AND. K .LE. NBYTES) THEN
          MASK = BITS(MOD(J, 8))
          IF (VALUE .NE. 0) THEN
            BITMAP(K) = CHAR(IOR(ICHAR(BITMAP(K)), MASK))
          ELSE
            BITMAP(K) = CHAR(IAND(ICHAR(BITMAP(K)), 255 - MASK))
          END IF
        END IF

        RETURN


        ENTRY GKBWBM
C*  write bitmap

        CALL GDEC(WIDTH, LW, W)
        CALL GDEC(HEIGHT, LH, H)

        HEADER = 'P4'//LF//W(1:LW)//' '//H(1:LH)//LF
        DO 2, K = 1, 11
          BITMAP(K - 11) = HEADER(K:K)
    2   CONTINUE

        CALL BINOUT(CONID, 11 + NBYTES, BITMAP)

        RETURN
        END
