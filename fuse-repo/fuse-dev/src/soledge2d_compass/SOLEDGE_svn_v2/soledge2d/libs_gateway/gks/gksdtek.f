C*
C* Copyright @ 1984 - 1993   Josef Heinen
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

        SUBROUTINE GKDTK1 (FCTID,DX,DY,DIMX,IA,LR1,R1,LR2,R2,LC,CHARS)
C*  GKS logical device driver for TEKtronix 401x series terminals

        INTEGER LC, LR1, LR2
        INTEGER FCTID,DX,DY,DIMX,IA(*)
        REAL R1(*),R2(*)
        CHARACTER*(*) CHARS

        EXTERNAL GK1PL,GK1SFA

        INTEGER I, J, BORDER, ICOLOR, ISTAT, ITERM, LTYPE, NP, NCHARS
        REAL X, Y, YRES, PX(2), PY(2)
        INTEGER ERRIND,TNR
        
        SAVE

C*  include GKS symbol definitons
        INCLUDE 'gksdefs.i'

C*  connection and type
        INTEGER CONID,WTYPE
C*  workstation state
        INTEGER STATE
C*  workstation transformation
        REAL WINDOW(4),VIEWPT(4)

C*  terminal control sequences
        CHARACTER ACSEQ*6,DASEQ*6
        INTEGER LACSEQ,LDASEQ

        GOTO (999,999,  2,  3,  4,  5,  6,999,  8,999,
     *         10,999, 12, 13, 14, 15,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999, 54, 55,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999, 81, 82,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999) FCTID+1

        IF (DX .LT. 0 .OR. DY .LT. 0 .OR. DIMX .LT. 0 .OR.
     *      LR1 .LT. 0 .OR. LR2 .LT. 0 .OR. LC .LT. 0)
     *      STOP 'GKS: internal inconsistency'
        GOTO 999

C*  open workstation
    2   CONTINUE
        CONID = IA(2)
        WTYPE = IA(3)

        IF (WTYPE .EQ. 72) THEN
C*   72: TEK4014
          ACSEQ(1:1) = CHAR(29)
          LACSEQ = 1
          DASEQ(1:1) = CHAR(31)
          LDASEQ = 1
        ELSE IF (WTYPE .EQ. 201) THEN
C*  201: TAB 132/15-G
C*  activate workstation - enter vector mode
          ACSEQ(1:1) = CHAR(29)
          LACSEQ = 1
C*  deactivate workstation - enter transparent mode
          ACSEQ(1:5) = CHAR(29)//CHAR(27)//'"0g'
          LACSEQ = 5
        ELSE IF (WTYPE .EQ. 204) THEN
C*  204: Monterey MG200
C*  activate workstation - enter graphic mode
          ACSEQ(1:4) = CHAR(27)//CHAR(92)//'5'//CHAR(29)
          LACSEQ = 4
C*  deactivate workstation - enter native mode
          DASEQ(1:1) = CHAR(24)
          LDASEQ = 1
        ELSE
C*   38: LN03 PLUS
C*  207: IBM Personal Computer
C*  activate workstation - modus TEK 4010/4014
          ACSEQ(1:6) = CHAR(27)//'[?38h'
          LACSEQ = 6
C*  deactivate workstation - modus VT220
          DASEQ(1:6) = CHAR(27)//'[?38l'
          LDASEQ = 6
        END IF

C*  set up connection
        CALL GK1SCO (CONID)
        IF (WTYPE .EQ. 38) THEN
          CALL GK1PB (ACSEQ(1:LACSEQ))
        END IF

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
        CALL GK1SDT (WINDOW,VIEWPT)
        GOTO 999

C*  close workstation
    3   CONTINUE
        IF (WTYPE .NE. 38) CALL GK1PB (ACSEQ(1:LACSEQ))
        CALL GK1PB (CHAR(27)//CHAR(12))
        CALL GK1PB (DASEQ(1:LDASEQ))
        CALL GK1U
        GOTO 999

C*  activate workstation
    4   CONTINUE
        IF (WTYPE .NE. 38) CALL GK1PB (ACSEQ(1:LACSEQ))
        STATE = GACTIV
        GOTO 999

C*  deactivate workstation
    5   CONTINUE
        IF (WTYPE .NE. 38) THEN
          CALL GK1PB (DASEQ(1:LDASEQ))
          CALL GK1U
        END IF
        STATE = GINACT
        GOTO 999

C*  clear workstation
    6   CONTINUE
        IF (WTYPE .NE. 38) THEN
          IF (STATE .NE. GACTIV) CALL GK1PB (ACSEQ(1:LACSEQ))
          CALL GK1PB (CHAR(27)//CHAR(12))
          CALL GK1PB (DASEQ(1:LDASEQ))
          CALL GK1U
          IF (STATE .EQ. GACTIV) CALL GK1PB (ACSEQ(1:LACSEQ))
        ELSE IF (WTYPE .EQ. 201) THEN
          CALL GK1PB (CHAR(29)//CHAR(27)//'"6g!ERA G'//CHAR(13)//
     *      CHAR(27)//'"0g')
          CALL GK1U
          IF (STATE .EQ. GACTIV) CALL GK1PB (ACSEQ(1:LACSEQ))
        ELSE
          CALL GK1PB (CHAR(27)//CHAR(12))
        END IF
        GOTO 999

C*  update workstation
    8   CONTINUE
        IF (WTYPE .NE. 38) THEN
          CALL GK1PB (DASEQ(1:LDASEQ))
          CALL GK1U
          CALL GK1PB (ACSEQ(1:LACSEQ))
        ELSE
          CALL GK1U
        END IF
        GOTO 999

C*  message
   10   CONTINUE
        IF (WTYPE .NE. 38) THEN
          IF (STATE .EQ. GACTIV) CALL GK1PB (DASEQ(1:LDASEQ))
          CALL GK1PB (CHARS)
          CALL GK1PB (CHAR(13)//CHAR(10))
          CALL GK1U
          IF (STATE .EQ. GACTIV) CALL GK1PB (ACSEQ(1:LACSEQ))
        ELSE
          CALL GK1PB (CHARS)
          CALL GK1PB (CHAR(13)//CHAR(10))
          CALL GK1U
        END IF
        GOTO 999

C*  polyline
   12   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire current line type
          CALL GQLN (ERRIND,LTYPE)
          CALL GQPLCI (ERRIND,ICOLOR)
          IF (ICOLOR .EQ. 0) LTYPE = 0
          IF (LTYPE .NE. GLSOLI) CALL GK1SLN (LTYPE)
          CALL GK1PL (IA(1),R1,R2,LTYPE,TNR)
          IF (LTYPE .NE. GLSOLI) CALL GK1SLN (GLSOLI)
        END IF
        GOTO 999

C*  polymarker
   13   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire polymarker color index
          CALL GQPMCI (ERRIND,ICOLOR)
          IF (ICOLOR .EQ. 0) CALL GK1SLN (ICOLOR)
          CALL GPOLMK (IA(1),R1,R2,GK1PL)
          IF (ICOLOR .EQ. 0) CALL GK1SLN (GLSOLI)
        END IF
        GOTO 999

C*  text
   14   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
          CALL GQTXCI (ERRIND,ICOLOR)
          IF (ICOLOR .EQ. 0) CALL GK1SLN (ICOLOR)
          NCHARS = LEN(CHARS)
          CALL GSIMTX (R1,R2,NCHARS,CHARS,GK1PL,GK1SFA)
          IF (ICOLOR .EQ. 0) CALL GK1SLN (GLSOLI)
        END IF
        GOTO 999

C*  fill area
   15   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
          CALL GQFACI (ERRIND,ICOLOR)
          IF (ICOLOR .EQ. 0) CALL GK1SLN (ICOLOR)
          IF (ICOLOR .LE. 1) THEN
            YRES = 1.0/768.0
            CALL GFILLA (IA(1),R1,R2,TNR,GK1PL,YRES)
          ELSE
            BORDER = 0
            CALL GK1PL (IA(1),R1,R2,BORDER,TNR)
          END IF
          IF (ICOLOR .EQ. 0) CALL GK1SLN (GLSOLI)
        END IF
        GOTO 999

C*  set workstation window
   54   CONTINUE
        WINDOW(1) = R1(1)
        WINDOW(2) = R1(2)
        WINDOW(3) = R2(1)
        WINDOW(4) = R2(2)
        CALL GK1SDT (WINDOW,VIEWPT)
        GOTO 999

C*  set workstation viewport
   55   CONTINUE
        VIEWPT(1) = R1(1)
        VIEWPT(2) = R1(2)
        VIEWPT(3) = R2(1)
        VIEWPT(4) = R2(2)
        CALL GKFVP (VIEWPT, 0.256, 0.192)
        CALL GK1SDT (WINDOW,VIEWPT)
        GOTO 999

C*  request locator
   81   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          CALL GK1GIN (CONID,ISTAT,ITERM,X,Y)
          IA(1) = ISTAT
          IA(4) = ITERM
          R1(1) = X
          R2(1) = Y
        END IF
        GOTO 999

C*  request stroke
   82   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
          NP = IA(3)
          DO 8201 I = 1,NP
            CALL GK1GIN (CONID,ISTAT,ITERM,X,Y)
            IF (ISTAT .NE. GOK) GOTO 8202
            R1(I) = X
            R2(I) = Y
            J = I-1
            IF (I .EQ. 1) THEN
              J = 1
              R1(2) = X
              R2(2) = Y
            END IF
            PX(1) = R1(J)
            PY(1) = R2(J)
            PX(2) = R1(J+1)
            PY(2) = R2(J+1)
            CALL GK1PL (2,PX,PY,GLSOLI,0)
            CALL GK1U
 8201     CONTINUE
 8202     IF (I.GT.2) ISTAT = GOK
          IA(1) = ISTAT
          IA(3) = I-1
          IA(4) = ITERM
        END IF

  999   RETURN
        END


        SUBROUTINE GK1PL (N,PX,PY,LTYPE,TNR)
C*  polyline

        REAL PX(N),PY(N)
        INTEGER TNR
        INTEGER LTYPE, N

        EXTERNAL GK1VEC,GK1COD

C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,LTYPE,TNR,GK1VEC,GK1COD)

        RETURN
        END


        SUBROUTINE GK1FA (N,PX,PY,TNR)
C*  fill area

        INTEGER N
        REAL PX(N),PY(N)
        INTEGER TNR

        REAL YRES

        EXTERNAL GK1PL

C*  switch to the fill area utility
        YRES = 1.0/768.0
        CALL GFILLA (N,PX,PY,TNR,GK1PL,YRES)

        RETURN
        END


        SUBROUTINE GK1SFA (N,PX,PY,TNR)
C*  solid fill area

        INTEGER N
        REAL PX(N),PY(N)
        INTEGER TNR

        REAL YRES

        EXTERNAL GK1PL

C*  switch to the solid fill area utility
        YRES = 1.0/768.0
        CALL GKSFA (N,PX,PY,TNR,GK1PL,YRES)

        RETURN
        END


        SUBROUTINE GK1SDT (WN,VP)
C*  set up device transformation

        REAL WN(4), VP(4), XCUR, YCUR
        INTEGER CONID
        INTEGER I, ISTAT, ITERM, IX, IY, JX, JY, NCHARS
        REAL X, Y, A, B, C, D, E, F, G, H

C*  terminal control sequences
        CHARACTER CHR*6, LCSEQ*3
        INTEGER IN(4), FREE

        SAVE

        INCLUDE 'gksdefs.i'

C*  locator - enter graphics input mode
        LCSEQ(1:3) = CHAR(29)//CHAR(27)//CHAR(26)

        E = (VP(2)-VP(1))/(WN(2)-WN(1))
        F = (1024-1)/0.256
        G = (VP(4)-VP(3))/(WN(4)-WN(3))
        H = (768-1)/0.192

        A = E*F
        B = F*(VP(1)-WN(1)*E)
        C = G*H
        D = H*(VP(3)-WN(3)*G)

        RETURN


        ENTRY GK1DT (X,Y,JX,JY)

C*  device transformation
        JX = NINT(A*X+B)
        JY = NINT(C*Y+D)

        RETURN


        ENTRY GK1VEC (X,Y)
C*  re-initialize vector drawing sequence

        CALL GK1BUF (FREE)
        IF (FREE .LE. 5) CALL GK1U

        CALL GK1PB (CHAR(29))


        ENTRY GK1COD (X,Y)

C*  device transformation
        IX = NINT(A*X+B)
        IY = NINT(C*Y+D)

        CALL GK1PCH (IX,IY)

        CALL GK1BUF (FREE)
        IF (FREE .LE. 4) THEN
            CALL GK1U
            CALL GK1PB (CHAR(29))
            CALL GK1PCH (IX,IY)
        END IF

        RETURN


        ENTRY GK1GIN (CONID,ISTAT,ITERM,XCUR,YCUR)
C*  graphic input

        CALL GK1U
        CALL BUFOUT (CONID,3,LCSEQ)

        NCHARS = 6
        CALL BUFIN (CONID,ISTAT,0,NCHARS,CHR)
        IF (ISTAT .EQ. GOK) THEN
          ITERM = ICHAR(CHR(1:1))
          DO 1 I = 2,5
            IN(I-1) = ICHAR(CHR(I:I))
   1      CONTINUE
C*  decode the input buffer
          IX = MOD(IN(1),32)*32+MOD(IN(2),32)
          IY = MOD(IN(3),32)*32+MOD(IN(4),32)
          XCUR = (IX-B)/A
          YCUR = (IY-D)/C
        ELSE
          ITERM = 0
        END IF

        RETURN
        END


        SUBROUTINE GK1PCH (IX,IY)
C*  calculate plot characters to arrive at IX,IY

        INTEGER IX, IY 

        CHARACTER CHR*4

C*  order is HIY,LOY,HIX,LOX
        CHR(1:1) = CHAR(MOD(IY/32,32)+ICHAR(' '))
        CHR(2:2) = CHAR(MOD(IY,32)+ICHAR('`'))
        CHR(3:3) = CHAR(MOD(IX/32,32)+ICHAR(' '))
        CHR(4:4) = CHAR(MOD(IX,32)+ICHAR('@'))

C*  output plot characters
        CALL GK1PB (CHR)

        RETURN
        END


        SUBROUTINE GK1SLN (LTYPE)
C*  set linetype

        INTEGER LTYPE
        CHARACTER LNTBL(-8:4)

C*  line type table
        DATA LNTBL /'`','`','`','`','g','f','e','d','x',
     &    '`','c','a','b'/

C*  there is no support for linetypes -9..-30
        IF (LTYPE .GE. 8) THEN
          CALL GK1PB (CHAR(27)//LNTBL(LTYPE))
        ELSE
          CALL GK1PB (CHAR(27)//LNTBL(2))
        END IF

        RETURN
        END


        SUBROUTINE GK1SCO (ID)
C*  set up connection

        INTEGER I, ID, IPNTR
        INTEGER CONID, FREE

        CHARACTER*(*) BUFF
        CHARACTER*500 IOBUFF

        SAVE

        DATA IPNTR /0/

        CONID = ID

        RETURN


        ENTRY GK1PB (BUFF)
C*  pack buffer

        IF (LEN(BUFF) .GT. 500-IPNTR-1) THEN
C*  xfer buffer
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = CHAR(10)

          CALL BUFOUT (CONID,IPNTR,IOBUFF)
          IPNTR = 0
        END IF

        DO 1 I = 1,LEN(BUFF)
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = BUFF(I:I)
   1    CONTINUE

        RETURN


        ENTRY GK1BUF (FREE)
C*  return free buffer space

        FREE = 500-IPNTR-1

        RETURN


        ENTRY GK1U
C*  update

        IPNTR = IPNTR+1
        IOBUFF(IPNTR:IPNTR) = CHAR(10)

        CALL BUFOUT (CONID,IPNTR,IOBUFF)
        IPNTR = 0

        RETURN
        END
