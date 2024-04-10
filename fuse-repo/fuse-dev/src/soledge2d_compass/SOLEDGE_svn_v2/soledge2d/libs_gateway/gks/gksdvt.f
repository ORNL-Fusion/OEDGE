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

        SUBROUTINE GKDVT3 (FCTID,DX,DY,DIMX,IA,LR1,R1,LR2,R2,LC,CHARS)
C*  GKS logical device driver for VT3xx series terminals

        INTEGER LC, LR1, LR2
        INTEGER FCTID,DX,DY,DIMX,IA(*)
        REAL R1(*),R2(*)
        CHARACTER*(*) CHARS

        CHARACTER*1 ESC,BS

        EXTERNAL GK3PL,GK3SCI,GK3FA

        INTEGER I, ICOLOR, ISTAT, ITERM, J, LTYPE, NP, NCHARS
        REAL X, Y, YRES
        INTEGER ERRIND, TNR, STYLE
        REAL R, G, B, PX(2), PY(2)
        
        SAVE

C*  include GKS symbol definitions
        INCLUDE 'gksdefs.i'

C*  connection and type
        INTEGER CONID,WSTYPE
C*  workstation state
        INTEGER STATE
C*  workstation transformation
        REAL WINDOW(4),VIEWPT(4)


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
        WSTYPE = IA(3)

        ESC = CHAR(27)
        BS = CHAR(92)

C*  initialize colormap
        CALL GK3ICM
C*  set up connection
        CALL GK3SCO (CONID)

C*  select reverse video, enter ReGIS mode
        CALL GK3PB (ESC//'[?5h'//ESC//'P1pW(V)')
C*  set colormap
        IF (WSTYPE .EQ. 17) THEN
          CALL GK3PB ('S(M1(AH0L50S100)2(AH120L50S100)')
          CALL GK3PB ('3(AH240L50S100)4(AH60L50S100)')
          CALL GK3PB ('5(AH300L50S100)6(AH180L50S100)')
          CALL GK3PB ('9(AH30L50S100)10(AH150L50S100)')
          CALL GK3PB ('11(AH270L50S100)12(AH90L50S100)')
          CALL GK3PB ('13(AH330L50S100)14(AH210L50S100))')
        END IF
        CALL GK3PB (ESC//BS)
        CALL GK3U

C*  set default workstation window
        WINDOW(1) = 0.
        WINDOW(2) = 1.
        WINDOW(3) = 0.
        WINDOW(4) = 1.

C*  set default workstation viewport
        VIEWPT(1) = 0.048
        VIEWPT(2) = 0.192
        VIEWPT(3) = 0.
        VIEWPT(4) = 0.144

C*  set up device transformation
        CALL GK3SDT (WINDOW,VIEWPT)

        GOTO 999

C*  close workstation
    3   CONTINUE
        CALL GK3U
        GOTO 999

C*  activate workstation
    4   CONTINUE
C*  enter ReGIS mode
        CALL GK3PB (ESC//'Pp')
        STATE = GACTIV
        GOTO 999

C*  deactivate workstation
    5   CONTINUE
C*  exit ReGIS mode
        CALL GK3PB (ESC//BS)
        CALL GK3U
        STATE = GINACT
        GOTO 999

C*  clear workstation
    6   CONTINUE
        IF (STATE .EQ. GACTIV) CALL GK3PB (ESC//BS)
        CALL GK3PB (ESC//'[H'//ESC//'[J')
        CALL GK3U
        IF (STATE .EQ. GACTIV) CALL GK3PB (ESC//'Pp')
        GOTO 999

C*  update workstation
    8   CONTINUE
        CALL GK3PB (ESC//BS)
        CALL GK3U
        CALL GK3PB (ESC//'Pp')
        GOTO 999

C*  message
   10   CONTINUE
        IF (STATE .EQ. GACTIV) CALL GK3PB (ESC//BS)
        CALL GK3PB (ESC//'[H')
        CALL GK3PB (CHARS)
        CALL GK3PB (CHAR(13)//CHAR(10))
        CALL GK3U
        IF (STATE .EQ. GACTIV) CALL GK3PB (ESC//'Pp')
        GOTO 999

C*  polyline
   12   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire current line type and color index
          CALL GQLN (ERRIND,LTYPE)
          CALL GQPLCI (ERRIND,ICOLOR)
          CALL GK3SCI (ICOLOR)
          IF (LTYPE .NE. GLSOLI) CALL GK3SLN (LTYPE)
          CALL GK3PL (IA(1),R1,R2,LTYPE,TNR)
          IF (LTYPE .NE. GLSOLI) CALL GK3SLN (GLSOLI)
        END IF
        GOTO 999

C*  polymarker
   13   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire polymarker color index
          CALL GQPMCI (ERRIND,ICOLOR)
          CALL GK3SCI (ICOLOR)
          CALL GPOLMK (IA(1),R1,R2,GK3PL)
        END IF
        GOTO 999

C*  text
   14   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
          CALL GQTXCI (ERRIND,ICOLOR)
          CALL GK3SCI (ICOLOR)
          NCHARS = LEN(CHARS)
          CALL GSIMTX (R1,R2,NCHARS,CHARS,GK3PL,GK3FA)
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
          CALL GK3SCI (ICOLOR)
C*  inquire current fill area interior style
          CALL GQFAIS (ERRIND,STYLE)
          IF (STYLE .NE. GSOLID) THEN
            YRES = 1.0/480.0
            CALL GFILLA (IA(1),R1,R2,TNR,GK3PL,YRES)
          ELSE
            CALL GK3FA (IA(1),R1,R2,TNR)
          END IF
        END IF
        GOTO 999

C*  cell array
   16   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
          CALL GCELLA (R1(1),R1(2),R2(1),R2(2),DX,DY,
     *      DIMX,IA,GK3SCI,GK3FA)
        END IF
        GOTO 999

C*  set color representation
   48   CONTINUE
        I = IA(2)
        IF (I.GE.0 .AND. I.LE.15) THEN
          R = R1(1)
          G = R1(2)
          B = R1(3)
          CALL GK3SCR (I,R,G,B)
        END IF
        GOTO 999

C*  set workstation window
   54   CONTINUE
        WINDOW(1) = R1(1)
        WINDOW(2) = R1(2)
        WINDOW(3) = R2(1)
        WINDOW(4) = R2(2)
        CALL GK3SDT (WINDOW,VIEWPT)
        GOTO 999

C*  set workstation viewport
   55   CONTINUE
        VIEWPT(1) = R1(1)
        VIEWPT(2) = R1(2)
        VIEWPT(3) = R2(1)
        VIEWPT(4) = R2(2)
        CALL GKFVP (VIEWPT, 0.24, 0.144)
        CALL GK3SDT (WINDOW,VIEWPT)
        GOTO 999

C*  request locator
   81   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
          CALL GK3GIN (CONID,ISTAT,ITERM,X,Y)
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
            CALL GK3GIN (CONID,ISTAT,ITERM,X,Y)
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
            CALL GK3PL (2,PX,PY,GLSOLI,0)
            CALL GK3U
 8201     CONTINUE
 8202     IF (I.GT.2) ISTAT = GOK
          IA(1) = ISTAT
          IA(3) = I-1
          IA(4) = ITERM
        END IF

  999   RETURN
        END


        SUBROUTINE GK3PL (N,PX,PY,LTYPE,TNR)
C*  polyline

        INTEGER LTYPE, N
        REAL PX(N),PY(N)
        INTEGER TNR

        EXTERNAL GK3MOV,GK3DRW

C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,LTYPE,TNR,GK3MOV,GK3DRW)

        RETURN
        END


        SUBROUTINE GK3MOV (X,Y)
C*  move

        EXTERNAL GK3VEC

        REAL X, Y
C*  switch to the move utility
        CALL GMOVE (X,Y,GK3VEC)

        RETURN
        END


        SUBROUTINE GK3DRW (X,Y)
C*  draw

        EXTERNAL GK3VEC,GK3COD
        REAL X, Y

C*  switch to the dashed line generator
        CALL GDASH (X,Y,GK3VEC,GK3COD)

        RETURN
        END


        SUBROUTINE GK3FA (N,PX,PY,TNR)
C*  fill area

        INTEGER N
        REAL PX(N),PY(N)
        INTEGER TNR

        EXTERNAL GK3MOV,GK3DRW

C*  enable polygon fill
        CALL GK3PB ('F(')
C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,0,TNR,GK3MOV,GK3DRW)
        CALL GK3PB (')')

        RETURN
        END


        SUBROUTINE GK3SDT (WN,VP)
C*  set up device transformation

        INTEGER ISTAT, ITERM, I, J, JJ, IX, IY, JX, JY, NCHARS
        REAL A, B, C, D, E, F, G, H, X, Y
        REAL WN(4), VP(4), XCUR, YCUR
        INTEGER CONID, CR

        CHARACTER CHR*16

        SAVE

C*  include GKS symbol definitions
        INCLUDE 'gksdefs.i'

        DATA CR /8192/

        E = (VP(2)-VP(1))/(WN(2)-WN(1))
        F = (800-1)/0.24
        G = (VP(4)-VP(3))/(WN(4)-WN(3))
        H = (480-1)/0.144

        A = E*F
        B = F*(VP(1)-WN(1)*E)
        C = G*H
        D = H*(VP(3)-WN(3)*G)

        RETURN


        ENTRY GK3DT (X,Y,JX,JY)

        JX =     NINT(A*X+B)
        JY = 479-NINT(C*Y+D)

        RETURN


        ENTRY GK3VEC (X,Y)
C*  re-initialize vector drawing sequence

        CALL GK3PB ('P')
        IX =     NINT(A*X+B)
        IY = 479-NINT(C*Y+D)

        CALL GK3PCH (IX,IY)

        RETURN


        ENTRY GK3COD (X,Y)
C*  device transformation

        CALL GK3PB ('V')
        IX =     NINT(A*X+B)
        IY = 479-NINT(C*Y+D)
        CALL GK3PCH (IX,IY)

        RETURN


        ENTRY GK3GIN (CONID,ISTAT,ITERM,XCUR,YCUR)
C*  graphic input

        CALL GK3U
C*  enter one-shot graphics input mode
        CALL BUFOUT (CONID,12,'R(I0)R(P(I))')

        NCHARS = 16
        CALL BUFIN (CONID,ISTAT,CR,NCHARS,CHR)
        IF (ISTAT .EQ. GOK) THEN
          IF (CHR(1:1) .LT. ' ')
     *      CALL BUFIN (CONID,ISTAT,CR,NCHARS,CHR)
          ITERM = ICHAR(CHR(1:1))
C*  decode the input buffer
          I = INDEX(CHR,'[')+1
          J = INDEX(CHR,',')-1
          JJ = INDEX(CHR(I:J),'.')-1
          IF (JJ .GT. 0) J = JJ
          CALL GSTR(J-I+1,CHR(I:J),IX)

          I = INDEX(CHR,',')+1
          J = INDEX(CHR,']')-1
          JJ = INDEX(CHR(I:J),'.')-1
          IF (JJ .GT. 0) J = JJ
          CALL GSTR(J-I+1,CHR(I:J),IY)

          XCUR =     (IX-B)/A
          YCUR = (479-IY-D)/C
        ELSE
          ITERM = 0
        END IF

        RETURN
        END


        SUBROUTINE GK3PCH (IX,IY)
C*  encode coordinates

        INTEGER IX, IY
        CHARACTER XSTR*5,YSTR*5
        INTEGER XLEN,YLEN

        CALL GDEC (IX,XLEN,XSTR)
        CALL GDEC (IY,YLEN,YSTR)

        CALL GK3PB ('['//XSTR(1:XLEN)//','//YSTR(1:YLEN)//']')

        RETURN
        END


        SUBROUTINE GK3SLN (LTYPE)
C*  set linetype

        INTEGER LTYPE
        CHARACTER*1 LNTBL(-8:4)

C*  line style table
        DATA LNTBL /'1','1','1','1','8','7','6','5','1',
     *    '1','2','4','3'/

C*  there is no support for linetypes -9..-30
        IF (LTYPE .GE. 8) THEN
          CALL GK3PB ('W(P'//LNTBL(LTYPE)//')')
        ELSE
          CALL GK3PB ('W(P'//LNTBL(1)//')')
        END IF

        RETURN
        END


        SUBROUTINE GK3ICM
C*  initialize colormap

        INTEGER ICOLOR
        REAL R, G, B

        INTEGER ML(0:11), PALETT(0:15), CMAP(0:15)
        INTEGER CI, COLOR

        INTEGER ILEN
        CHARACTER ISTR*3

        SAVE

C*  colormap         W  D  R  G  B  C  Y  M  R* G* B* C* Y* M*
        DATA PALETT /7, 0, 2, 3, 1, 5, 6, 4,10,11, 9,13,14,12, 8,15/
C*  map location
        DATA ML /4,10,7,13,2,8,6,12,3,9,5,11/

        DATA COLOR /1/

        DO 1, I = 0, 15
            CMAP(I) = PALETT(I)
   1    CONTINUE

        COLOR = -1

        RETURN


        ENTRY GK3SCI (ICOLOR)
C*  set color index

        IF (ICOLOR .GE. 8) THEN
            CI = ML(MOD(ICOLOR-8,72)/6)
        ELSE
            CI = ICOLOR
        END IF

        IF (CI .NE. COLOR) THEN
            CALL GDEC (CMAP(CI),ILEN,ISTR)
            CALL GK3PB ('W(I'//ISTR(1:ILEN)//')')
            COLOR = CI
        END IF

        RETURN


        ENTRY GK3SCR (ICOLOR, R, G, B)
C*  set color representation

        IF (ICOLOR .LT. 8)
     *      CALL GKRMAP (ICOLOR, R, G, B, 8, PALETT, CMAP)

        RETURN
        END


        SUBROUTINE GK3SCO (ID)
C*  set up connection

        INTEGER I, ID, IPNTR
        CHARACTER*(*) BUFF

        CHARACTER IOBUFF*80
        INTEGER CONID

        SAVE

        DATA IPNTR /0/

        CONID = ID

        RETURN


        ENTRY GK3PB (BUFF)
C*  pack buffer

        IF (IPNTR+LEN(BUFF) .GT. 80) THEN
C*  xfer buffer
          CALL BUFOUT (CONID,IPNTR,IOBUFF)
          IPNTR = 0
        END IF

        DO 1 I = 1,LEN(BUFF)
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = BUFF(I:I)
   1    CONTINUE

        RETURN


        ENTRY GK3U
C*  update

        CALL BUFOUT (CONID,IPNTR,IOBUFF)
        IPNTR = 0

        RETURN
        END
