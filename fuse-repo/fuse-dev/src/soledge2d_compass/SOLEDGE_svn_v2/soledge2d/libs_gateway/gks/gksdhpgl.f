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

        SUBROUTINE GKDHP (FCTID,DX,DY,DIMX,IA,LR1,R1,LR2,R2,LC,CHARS)
C*  GKS logical device handler for HP-GL devices

        INTEGER LC, LR1, LR2
        INTEGER FCTID,DX,DY,DIMX,IA(*)
        REAL R1(*),R2(*)
        CHARACTER*(*) CHARS

        EXTERNAL GKHPL,GKHSPL,GKHSFA

        REAL YRES
        INTEGER ERRIND,TNR,COLI,STYLE,LTYPE,BORDER,NCHARS
        LOGICAL EMPTY

        SAVE

C*  include GKS symbol definitions
        INCLUDE 'gksdefs.i'

C *     Workstation State List
C               connection and type
        INTEGER CONID
C               workstation state
        INTEGER STATE
C               workstation transformation
        REAL WINDOW(4),VIEWPT(4)


        GOTO (999,999,  2,  3,  4,  5,  6,999,  8,999,
     *        999,999, 12, 13, 14, 15,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999, 54, 55,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999) FCTID+1

        IF (LR1 .LT. 0 .OR. LR2 .LT. 0 .OR. LC .LT. 0 .OR.
     *      DX .LT. 0 .OR. DY .LT. 0 .OR. DIMX .LT. 0)
     *      STOP 'GKS: internal inconsistency'
        GOTO 999

C*  open workstation
    2   CONTINUE
        CONID = IA(2)

C*  set up connection
        CALL GKHSCO (CONID)
        CALL GKHPB (
     *    CHAR(27)//'.Y'//CHAR(27)//'.N;19:'//CHAR(27)//'.I;;17:PS4;'
     *    //'IN;')
C*  select default pen
        CALL GKHPEN (-1)

C*  set default workstation window
        WINDOW(1) = 0.
        WINDOW(2) = 1.
        WINDOW(3) = 0.
        WINDOW(4) = 1.

C*  set default workstation viewport
        VIEWPT(1) = 0.046
        VIEWPT(2) = 0.226
        VIEWPT(3) = 0.005
        VIEWPT(4) = 0.185

C*  set up device transformation
        CALL GKHSDT (WINDOW,VIEWPT)

        EMPTY = .TRUE.
        GOTO 999

C*  close workstation
    3   CONTINUE
CCC        CALL GKHPAC
        CALL GKHPU (1)
        CALL GKHPB ('AF;SP;')
        CALL GKHU
        GOTO 999

C*  activate workstation
    4   CONTINUE
        STATE = GACTIV
        GOTO 999

C*  deactivate workstation
    5   CONTINUE
        STATE = GINACT
        CALL GKHPU (1)
        CALL GKHU
        GOTO 999

C*  clear workstation
    6   CONTINUE
        IF (.NOT.EMPTY) THEN
          CALL GKHPU (1)
          CALL GKHPB ('AF;')
          CALL GKHU
          EMPTY = .TRUE.
        END IF
        GOTO 999

C*  update workstation
    8   CONTINUE
        CALL GKHPU (1)
        CALL GKHU
        GOTO 999

C*  polyline
   12   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire line type, linewidth, and color index
          CALL GQLN (ERRIND,LTYPE)
          CALL GQPLCI (ERRIND,COLI)
          CALL GKHPEN (COLI)
C*  polyline
          CALL GKHPL (IA(1),R1,R2,LTYPE,TNR)
          CALL GKHPU (1)
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  polymarker
   13   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire current polymarker color index
          CALL GQPMCI (ERRIND,COLI)
          CALL GKHPEN (COLI)
C*  polymarker
          CALL GPOLMK (IA(1),R1,R2,GKHPL)
          CALL GKHPU (1)
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  text
   14   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire current text color index
          CALL GQTXCI (ERRIND,COLI)
          CALL GKHPEN (COLI)
C*  text
          NCHARS = LEN(CHARS)
          CALL GSIMTX (R1,R2,NCHARS,CHARS,GKHSPL,GKHSFA)
          CALL GKHPU (1)
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  fill area
   15   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
C*  inquire current normalization transformation
          CALL GQCNTN (ERRIND,TNR)
C*  set up device transformation
          CALL GSDT (WINDOW,VIEWPT)
C*  inquire current fill area color index
          CALL GQFACI (ERRIND,COLI)
          CALL GKHPEN (COLI)
C*  inquire current fill area interior style
          CALL GQFAIS (ERRIND,STYLE)
          IF (STYLE .NE. GSOLID) THEN
            YRES = 10.0/7600.0
            CALL GFILLA (IA(1),R1,R2,TNR,GKHPL,YRES)
          ELSE
            BORDER = 0
            CALL GKHPL (IA(1),R1,R2,BORDER,TNR)
          END IF
          CALL GKHPU (1)
          EMPTY = .FALSE.
        END IF
        GOTO 999

C*  set workstation window
   54   CONTINUE
        WINDOW(1) = R1(1)
        WINDOW(2) = R1(2)
        WINDOW(3) = R2(1)
        WINDOW(4) = R2(2)
        CALL GKHSDT (WINDOW,VIEWPT)
        GOTO 999

C*  set workstation viewport
   55   CONTINUE
        VIEWPT(1) = R1(1)
        VIEWPT(2) = R1(2)
        VIEWPT(3) = R2(1)
        VIEWPT(4) = R2(2)
        CALL GKHSDT (WINDOW,VIEWPT)
        GOTO 999

  999   RETURN
        END


        SUBROUTINE GKHPL (N,PX,PY,LTYPE,TNR)
C*  polyline

        REAL PX(N),PY(N)
        INTEGER TNR
        INTEGER LTYPE, N

        EXTERNAL GKHMOV,GKHDRW

C*  switch to the polyline utility
        CALL GPOLIN (N,PX,PY,LTYPE,TNR,GKHMOV,GKHDRW)

        RETURN
        END


        SUBROUTINE GKHFA (N,PX,PY,TNR)
C*  fill area

        INTEGER N
        REAL PX(N),PY(N)
        INTEGER TNR

        EXTERNAL GKHPL

        REAL YRES

C*  switch to the fill area utility
        YRES = 10.0/7600.0
        CALL GFILLA (N,PX,PY,TNR,GKHPL,YRES)

        RETURN
        END


        SUBROUTINE GKHSPL (N,PX,PY,LTYPE,TNR)
C*  solid polyline

        INTEGER N
        REAL PX(N),PY(N)
        INTEGER LTYPE,TNR

        EXTERNAL GKHVEC,GKHCOD

C*  switch to the solid polyline utility
        CALL GPOLIN (N,PX,PY,LTYPE,TNR,GKHVEC,GKHCOD)

        RETURN
        END


        SUBROUTINE GKHSFA (N,PX,PY,TNR)
C*  solid fill area

        INTEGER N
        REAL PX(N),PY(N)
        INTEGER TNR

        EXTERNAL GKHSPL

        REAL YRES

C*  switch to the solid fill area utility
        YRES = 10.0/7600.0
        CALL GKSFA (N,PX,PY,TNR,GKHSPL,YRES)

        RETURN
        END


        SUBROUTINE GKHMOV (X,Y)
C*  move routine

        EXTERNAL GKHVEC
        REAL X, Y

C*  switch to the move utility
        CALL GMOVE (X,Y,GKHVEC)

        RETURN
        END


        SUBROUTINE GKHDRW (X,Y)
C*  draw routine

        EXTERNAL GKHVEC,GKHCOD
        REAL X, Y
 
C*  switch to the dashed line generator
        CALL GDASH (X,Y,GKHVEC,GKHCOD)

        RETURN
        END


        SUBROUTINE GKHSDT (WN,VP)
C*  set up device transformation

        INTEGER IX, IY, JX, JY, L 
        REAL A, B, C, D, E, F, G, H, X, Y
        REAL WN(4), VP(4)

        SAVE

        INCLUDE 'gksdefs.i'

        E = (VP(2)-VP(1))/(WN(2)-WN(1))
        F = 10870/0.272
        G = (VP(4)-VP(3))/(WN(4)-WN(3))
        H = 7600/0.19

        A = E*F
        B = F*(VP(1)-WN(1)*E)
        C = G*H
        D = H*(VP(3)-WN(3)*G)

        RETURN


        ENTRY GKHDT (X,Y,JX,JY)

C*  device transformation
        JX = NINT(A*X+B)
        JY = NINT(C*Y+D)

        RETURN


        ENTRY GKHVEC (X,Y)
C*  re-initialize vector drawing sequence

        CALL GKHPU (1)
        L = 0


        ENTRY GKHCOD (X,Y)

C*  device transformation
        IX = NINT(A*X+B)
        IY = NINT(C*Y+D)

        IF (L .EQ. 0) THEN
          L = 1
        ELSE
          IF (L .EQ. 1) CALL GKHPU (0)
          L = 2
        END IF
        CALL GKHPCH (IX,IY)
        
        RETURN
        END


        SUBROUTINE GKHPCH (IX,IY)
C*  encode coordinates

        INTEGER IX, IY
        CHARACTER XSTR*16,YSTR*16
        INTEGER XLEN,YLEN

        CALL GDEC (IX,XLEN,XSTR)
        CALL GDEC (IY,YLEN,YSTR)

        CALL GKHPB ('PA'//XSTR(1:XLEN)//','//YSTR(1:YLEN)//';')

        RETURN
        END


        SUBROUTINE GKHPEN (COLI)
C*  select pen

        INTEGER COLI

        CHARACTER STR*1
        INTEGER PEN
        DATA PEN /0/

        IF (COLI .NE. PEN) THEN
          IF (COLI .LE. 8) THEN
C*  select pen
            PEN = ABS(COLI)
            STR = CHAR(ICHAR('0')+PEN)
            CALL GKHPU (-1)
            CALL GKHPB ('SP'//STR//';')
          END IF
        END IF

        RETURN
        END


        SUBROUTINE GKHPU (STATE)
C*  pen control

        INTEGER STATE

        INTEGER PENUP
        DATA PENUP /0/

        IF (STATE .NE. PENUP) THEN
          IF (STATE .NE. 0) THEN
            CALL GKHPB ('PU;')
            PENUP = 1
          ELSE
            CALL GKHPB ('PD;')
            PENUP = 0
          END IF
        END IF

        RETURN
        END


CCC        SUBROUTINE GKHPAC
CCC
CCC        INTEGER NCHARS
CCC        CHARACTER BUFFER*100
CCC
CCC        CALL GKINFO (NCHARS, BUFFER)
CCC
CCC        IF (NCHARS .GT. 0) THEN
CCC          CALL GKHPB ('SP1;PA0,0;LO11;SR0.5,1;SS;')
CCC          CALL GKHPB ('LB'//BUFFER(1:24)//'  BY  '//
CCC     *        BUFFER(36:NCHARS)//CHAR(3))
CCC          CALL GKHPB ('PU;')
CCC        END IF
CCC
CCC        RETURN
CCC        END


        SUBROUTINE GKHSCO (ID)
C*  set up connection

        INTEGER I, ID, IPNTR
        CHARACTER*(*) BUFF

        CHARACTER IOBUFF*500
        INTEGER CONID

        SAVE

        DATA IPNTR /0/

        CONID = ID

        RETURN


        ENTRY GKHPB (BUFF)
C*  pack buffer

        IF (LEN(BUFF)+2 .GT. 500-IPNTR) THEN
C*  xfer buffer
          CALL BUFOUT (CONID,IPNTR,IOBUFF)
          IPNTR = 0
        END IF

        DO 1 I = 1,LEN(BUFF)
          IPNTR = IPNTR+1
          IOBUFF(IPNTR:IPNTR) = BUFF(I:I)
   1    CONTINUE

        IPNTR = IPNTR+1
        IOBUFF(IPNTR:IPNTR) = CHAR(10)

        RETURN


        ENTRY GKHU
C*  update

        CALL BUFOUT (CONID,IPNTR,IOBUFF)
        IPNTR = 0

        RETURN
        END
