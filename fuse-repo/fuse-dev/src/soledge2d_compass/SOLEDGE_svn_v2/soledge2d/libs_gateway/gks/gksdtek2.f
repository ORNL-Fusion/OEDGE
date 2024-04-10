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

        SUBROUTINE GKDTK2 (FCTID,DX,DY,DIMX,IA,LR1,R1,LR2,R2,LC,CHARS)
C*  GKS logical device driver for TEKtronix 42xx series terminals

        INTEGER LC, LR1, LR2
        INTEGER FCTID,DX,DY,DIMX,IA(*)
        REAL R1(*),R2(*)
        CHARACTER*(*) CHARS

        EXTERNAL GK2PL,GK2SCI,GK2FA

        INTEGER I, J, LTYPE, NP, NCHARS
        INTEGER ICOLOR, ISTAT, ITERM
        INTEGER ERRIND,STYLE,TNR
        REAL B, G, R, PX(2), PY(2), X, Y, YRES
        
        SAVE

C*  include GKS symbol definitions
        INCLUDE 'gksdefs.i'

C*  connection and type
        INTEGER CONID
C*  workstation state
        INTEGER STATE
C*  workstation transformation
        REAL WINDOW(4),VIEWPT(4)

C*  terminal control sequences
        CHARACTER OPSEQ(25),CLSEQ(11),ACSEQ(14),DASEQ(5),CLRSEQ(2)
        CHARACTER IC(1),CRLF(2)

C*  open workstation
C*      select text precision "STRING"
C*      set color mode
C*      set communication region color index
C*      disable communication region buffer
        DATA OPSEQ /'?','M','Q','1','?','T','M','1','1','1',
     &    '?','L','I','1','0','0','?','T','D','1','0','?','K','A','0'/
C*  close workstation - set communication region color index
        DATA CLSEQ /'?','L','I','3','0','0','?','T','D','3','0'/
C*  activate workstation - select TEK mode
        DATA ACSEQ /'?','?','?','?','?','?','?','?','?','?',
     &  '?','%','!','0'/
C*  deactivate workstation - select ANSI mode
        DATA DASEQ /'?','%','!','1','?'/

        GOTO (999,999,  2,  3,  4,  5,  6,999,  8,999,
     *         10,999, 12, 13, 14, 15, 16,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999, 48,999,
     *        999,999,999,999, 54, 55,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999, 81, 82,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999) FCTID+1

        IF (LR1 .LT. 0 .OR. LR2 .LT. 0 .OR. LC .LT. 0)
     *      STOP 'GKS: internal inconsistency'
        GOTO 999

C*  open workstation
    2   CONTINUE
        CONID = IA(2)

C*  open workstation
        OPSEQ(1) = CHAR(27)
        OPSEQ(5) = CHAR(27)
        OPSEQ(11)= CHAR(27)
        OPSEQ(17)= CHAR(27)
        OPSEQ(22)= CHAR(27)
C*  close workstation
        CLSEQ(1) = CHAR(27)
        CLSEQ(7) = CHAR(27)
C*  activate workstation
        DO 2001, I=1,10
          ACSEQ(I) = CHAR(22)
 2001   CONTINUE
        ACSEQ(11)= CHAR(27)
C*  deactivate workstation
        DASEQ(1) = CHAR(27)
        DASEQ(5) = CHAR(13)
C*  clear workstation
        CLRSEQ(1)= CHAR(27)
        CLRSEQ(2)= CHAR(12)
C*  newline
        CRLF(1)= CHAR(13)
        CRLF(2)= CHAR(10)

C*  initialize colormap
        CALL GK2ICM
C*  set up connection
        CALL GK2SCO (CONID)
        CALL GK2PB (14,ACSEQ)
        CALL GK2PB (25,OPSEQ)
C*  change foreground/background color
        CALL GK2SCR (1,0.0,0.0,0.0)
        CALL GK2SCR (0,1.0,1.0,1.0)
        CALL GK2PB (5,DASEQ)
        CALL GK2U

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
        CALL GK2SDT (WINDOW,VIEWPT)
        GOTO 999

C*  close workstation
    3   CONTINUE
        CALL GK2PB (14,ACSEQ)
        CALL GK2PB (2,CLRSEQ)
        CALL GK2PB (11,CLSEQ)
C*  change foreground/background color
        CALL GK2SCR (1,1.0,1.0,1.0)
        CALL GK2SCR (0,0.0,0.0,0.0)
        CALL GK2PB (5,DASEQ)
        CALL GK2U
        GOTO 999

C*  activate workstation
    4   CONTINUE
        CALL GK2PB (14,ACSEQ)
        STATE = GACTIV
        GOTO 999

C*  deactivate workstation
    5   CONTINUE
        CALL GK2PB (5,DASEQ)
        CALL GK2U
        STATE = GINACT
        GOTO 999

C*  clear workstation
    6   CONTINUE
        IF (STATE .NE. GACTIV) CALL GK2PB (14,ACSEQ)
        CALL GK2PB (2,CLRSEQ)
        CALL GK2PB (5,DASEQ)
        CALL GK2U
        IF (STATE .EQ. GACTIV) CALL GK2PB (14,ACSEQ)
        GOTO 999

C*  update workstation
    8   CONTINUE
        CALL GK2PB (5,DASEQ)
        CALL GK2U
        CALL GK2PB (14,ACSEQ)
        GOTO 999

C*  message
   10   CONTINUE
        IF (STATE .EQ. GACTIV) CALL GK2PB (5,DASEQ)
        DO 100 I=1,LEN(CHARS)
          IC(1) = CHARS(I:I)
          CALL GK2PB (1,IC)
  100   CONTINUE
        CALL GK2PB (2,CRLF)
        CALL GK2U
        IF (STATE .EQ. GACTIV) CALL GK2PB (14,ACSEQ)
        GOTO 999

C*  polyline
   12   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire current line type, linewidth, and color index
          CALL GQLN (ERRIND,LTYPE)
          IF (LTYPE .NE. GLSOLI) CALL GK2SLN (LTYPE)
          CALL GQPLCI (ERRIND,ICOLOR)
          CALL GK2SCI (ICOLOR)
          CALL GK2PL (IA(1),R1,R2,LTYPE,TNR)
          IF (LTYPE .NE. GLSOLI) CALL GK2SLN (GLSOLI)
        END IF
        GOTO 999

C*  polymarker
   13   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire polymarker color index
          CALL GQPMCI (ERRIND,ICOLOR)
          CALL GK2SCI (ICOLOR)
          CALL GPOLMK (IA(1),R1,R2,GK2PL)
        END IF
        GOTO 999

C*  text
   14   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
          CALL GQTXCI (ERRIND,ICOLOR)
          CALL GK2SCI (ICOLOR)
          NCHARS = LEN(CHARS)
          CALL GSIMTX (R1,R2,NCHARS,CHARS,GK2PL,GK2FA)
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
          CALL GK2SCI (ICOLOR)
C*  inquire current fill area interior style
          CALL GQFAIS (ERRIND,STYLE)
          IF (STYLE .NE. GSOLID) THEN
            YRES = 1.0/768.0
            CALL GFILLA (IA(1),R1,R2,TNR,GK2PL,YRES)
          ELSE
            CALL GK2FA (IA(1),R1,R2,TNR)
          END IF
        END IF
        GOTO 999

C*  cell array
   16   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
          CALL GCELLA (R1(1),R1(2),R2(1),R2(2),DX,DY,
     *      DIMX,IA,GK2SCI,GK2FA)
        END IF
        GOTO 999

C*  set color representation
   48   CONTINUE
        I = IA(2)
        IF (I.GE.0 .AND. I.LE.8) THEN
          R = R1(1)
          G = R1(2)
          B = R1(3)
          IF (STATE .EQ. GINACT) CALL GK2PB (14,ACSEQ)
          CALL GK2SCR (I,R,G,B)
          CALL GK2PB (5,DASEQ)
          CALL GK2U
          IF (STATE .EQ. GACTIV) CALL GK2PB (14,ACSEQ)
        END IF
        GOTO 999

C*  set workstation window
   54   CONTINUE
        WINDOW(1) = R1(1)
        WINDOW(2) = R1(2)
        WINDOW(3) = R2(1)
        WINDOW(4) = R2(2)
        CALL GK2SDT (WINDOW,VIEWPT)
        GOTO 999

C*  set workstation viewport
   55   CONTINUE
        VIEWPT(1) = R1(1)
        VIEWPT(2) = R1(2)
        VIEWPT(3) = R2(1)
        VIEWPT(4) = R2(2)
        CALL GKFVP (VIEWPT, 0.256, 0.192)
        CALL GK2SDT (WINDOW,VIEWPT)
        GOTO 999

C*  request locator
   81   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          CALL GK2GIN (CONID,ISTAT,ITERM,X,Y)
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
            CALL GK2GIN (CONID,ISTAT,ITERM,X,Y)
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
            CALL GK2PL (2,PX,PY,GLSOLI,0)
            CALL GK2U
 8201     CONTINUE
 8202     IF (I.GT.2) ISTAT = GOK
          IA(1) = ISTAT
          IA(3) = I-1
          IA(4) = ITERM
        END IF

  999   RETURN
        END


        SUBROUTINE GK2PL (N,PX,PY,LTYPE,TNR)
C*  polyline

        REAL PX(N),PY(N)
        INTEGER N,LTYPE,TNR

        EXTERNAL GK2VEC,GK2COD

C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,LTYPE,TNR,GK2VEC,GK2COD)

        RETURN
        END


        SUBROUTINE GK2FA (N,PX,PY,TNR)
C*  fill area

        REAL PX(N),PY(N)
        INTEGER N
        INTEGER TNR

        EXTERNAL GK2V1,GK2C1

C*  begin panel
        CALL GK2BP
C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,0,TNR,GK2V1,GK2C1)
C*  end panel
        CALL GK2EP

        RETURN
        END


        SUBROUTINE GK2SDT (WN,VP)
C*  set up device transformation

        INTEGER I, ISTAT, ITERM, NCHARS, IX, IY, JX, JY
        REAL A, B, C, D, E, F, G, H, X, Y
        REAL WN(4), VP(4), XCUR, YCUR
        INTEGER CONID

C*  terminal control sequences
        CHARACTER PLSEQ(1),LF(1)
        CHARACTER PSEL(3),PBEG(3),PEND(3)
        CHARACTER LCSEQ*2 

        INTEGER COLI,IN(4)
        LOGICAL PFLAG

        CHARACTER CHR*6

        SAVE

        INCLUDE 'gksdefs.i'

C*  select fill pattern
        DATA PSEL /'?','M','P'/
C*  begin panel
        DATA PBEG /'?','L','P'/
C*  end panel
        DATA PEND /'?','L','E'/

        DATA PFLAG /.FALSE./
        DATA COLI /0/

C*  polyline - enter vector mode
        PLSEQ(1) = CHAR(29)
C*  locator - enter graphics input mode
        LCSEQ(1:1) = CHAR(27)
        LCSEQ(2:2) = CHAR(26)
C*  bypass cancelation character
        LF(1) = CHAR(10)

        PSEL(1)  = CHAR(27)
        PBEG(1)  = CHAR(27)
        PEND(1)  = CHAR(27)

        E = (VP(2)-VP(1))/(WN(2)-WN(1))
        F = (1024-1)/0.256
        G = (VP(4)-VP(3))/(WN(4)-WN(3))
        H = (768-1)/0.192

        A = E*F
        B = F*(VP(1)-WN(1)*E)
        C = G*H
        D = H*(VP(3)-WN(3)*G)

        RETURN


        ENTRY GK2DT (X,Y,JX,JY)

C*  device transformation
        JX = NINT(A*X+B)
        JY = NINT(C*Y+D)

        RETURN


        ENTRY GK2VEC (X,Y)
C*  re-initialize vector drawing sequence

        CALL GK2PB (1,PLSEQ)


        ENTRY GK2COD (X,Y)

C*  device transformation
        IX = NINT(A*X+B)
        IY = NINT(C*Y+D)
        CALL GK2PCH (IX,IY)

        RETURN


        ENTRY GK2GIN (CONID,ISTAT,ITERM,XCUR,YCUR)
C*  graphic input

        CALL GK2U
        CALL BUFOUT (CONID,2,LCSEQ)

        NCHARS = 6
        CALL BUFIN (CONID,ISTAT,0,NCHARS,CHR)
C*  send bypass cancelation character
        CALL GK2PB (1,LF)
        CALL GK2U

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


        ENTRY GK2BP
C*  begin panel

C*  inquire current colour index
        CALL GK2QCI (COLI)
C*  select fill pattern
        CALL GK2PB (3,PSEL)
        CALL GK2PI (-COLI)

        CALL GK2PB (3,PBEG)
        PFLAG = .TRUE.

        RETURN


        ENTRY GK2EP
C*  end panel

        CALL GK2PB (3,PEND)

        RETURN


        ENTRY GK2V1 (X,Y)
C*  re-initialize vector drawing sequence

        IF (.NOT. PFLAG) CALL GK2PB (1,PLSEQ)


        ENTRY GK2C1 (X,Y)
C*  device transformation

        IX = NINT(A*X+B)
        IY = NINT(C*Y+D)
        CALL GK2PCH (IX,IY)

        IF (PFLAG) THEN
           CALL GK2PB (1,PLSEQ)
           CALL GK2PCH (IX,IY)
           PFLAG = .FALSE.
        END IF

        RETURN
        END


        SUBROUTINE GK2PCH (IX,IY)
C*  calculate plot characters to arrive at IX,IY

        INTEGER IX, IY
        CHARACTER CHR(4)

C*  order is HIY,LOY,HIX,LOX
        CHR(1) = CHAR(MOD(IY/32,32)+ICHAR(' '))
        CHR(2) = CHAR(MOD(IY,32)+ICHAR('`'))
        CHR(3) = CHAR(MOD(IX/32,32)+ICHAR(' '))
        CHR(4) = CHAR(MOD(IX,32)+ICHAR('@'))

C*  output plot characters
        CALL GK2PB (4,CHR)

        RETURN
        END


        SUBROUTINE GK2SLN (LTYPE)
C*  set linetype

        INTEGER LTYPE
        CHARACTER SLNSEQ(4),LNTBL(-8:4)

C*  set line style
        DATA SLNSEQ /'?','M','V','0'/
C*  line style table
        DATA LNTBL /'0','0','0','0','7','6','5','4','0',
     &    '0','3','1','2'/

        SLNSEQ(1) = CHAR(27)
C*  there is no support for linetypes -9..-30
        IF (LTYPE .GE. -8) THEN
          SLNSEQ(4) = LNTBL(LTYPE)
        ELSE
          SLNSEQ(4) = LNTBL(2)
        END IF
        CALL GK2PB (4,SLNSEQ)

        RETURN
        END


        SUBROUTINE GK2SCR (INDEX,R,G,B)
C*  set color representation

        INTEGER IB, IG, IR
        INTEGER INDEX
        REAL R,G,B

        CHARACTER CSEQ(6),DACSEQ(5)

C*  assign colors
        DATA CSEQ /'?','T','G','1','4','0'/
        DATA DACSEQ /'?','T','F','4','0'/

        CSEQ(1) = CHAR(27)
        DACSEQ(1)= CHAR(27)

        IF (INDEX .LE. 8) THEN
          IR = NINT(R*100)
          IG = NINT(G*100)
          IB = NINT(B*100)

          CSEQ(6) = CHAR(ICHAR('0')+INDEX)
          CALL GK2PB (6,CSEQ)
          CALL GK2PI (IR)
          CALL GK2PI (IG)
          CALL GK2PI (IB)

          IF (INDEX .LE. 7) THEN
            DACSEQ(5) = CHAR(ICHAR('0')+INDEX)
            CALL GK2PB (5,DACSEQ)
            CALL GK2PI (IR)
            CALL GK2PI (IG)
            CALL GK2PI (IB)
          END IF
        END IF

        RETURN
        END


        SUBROUTINE GK2PI (I)
C*  pack integer

        INTEGER I,J
        CHARACTER CHR(1)

        J = ABS(I)
        IF (J.GT.15) THEN
          CHR(1) = CHAR(ICHAR('@')+MOD(J/16,64))
          CALL GK2PB (1,CHR)
        END IF
        IF (I.GE.0) THEN
          CHR(1) = CHAR(ICHAR('0')+MOD(J,16))
        ELSE
          CHR(1) = CHAR(ICHAR(' ')+MOD(J,16))
        END IF
        CALL GK2PB (1,CHR)

        RETURN
        END


        SUBROUTINE GK2ICM
C*  initialize colormap

        INTEGER ICOLOR

        CHARACTER SCSEQ(4)

        INTEGER CMAP(0:15)
        INTEGER ML(0:11), COLOR

        INTEGER CI

C*  set color index
        DATA SCSEQ /'?','M','L','1'/
C*  colormap       W  B  R  G  B  C  Y  M  R* G* B* C* Y* M*
        DATA CMAP /0, 1, 2, 3, 4, 5, 7, 6, 8,10,12,11, 9,13,14,15/
C*  map location
        DATA ML /4,10,7,13,2,8,6,12,3,9,5,11/

        DATA COLOR /1/

        SCSEQ(1) = CHAR(27)
        COLOR = -1

        RETURN


        ENTRY GK2SCI (ICOLOR)
C*  set color index

        IF (ICOLOR .GE. 8) THEN
            CI = ML(MOD(ICOLOR-8,72)/10)
        ELSE
            CI = ICOLOR
        END IF

        IF (CI .NE. COLOR) THEN
            SCSEQ(4) = CHAR(ICHAR('0')+CMAP(CI))
            CALL GK2PB (4,SCSEQ)
            COLOR = CI
        END IF


        ENTRY GK2QCI (ICOLOR)
C*  inquire color index

        ICOLOR = COLOR

        RETURN
        END


        SUBROUTINE GK2SCO (ID)
C*  set up connection

        INTEGER I, ID, IPNTR, LEN
        CHARACTER IOBUFF*80
        CHARACTER BUFF(LEN)

        INTEGER CONID

        SAVE

        DATA IPNTR /0/

        CONID = ID

        RETURN


        ENTRY GK2PB (LEN,BUFF)
C*  pack buffer

        IF (IPNTR+LEN .GT. 80) THEN
C*  xfer buffer
          CALL BUFOUT (CONID,IPNTR,IOBUFF)
          IPNTR = 0
        END IF

        DO 1 I = 1,LEN
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = BUFF(I)
   1    CONTINUE

        RETURN


        ENTRY GK2U
C*  update

        CALL BUFOUT (CONID,IPNTR,IOBUFF)
        IPNTR = 0

        RETURN
        END
