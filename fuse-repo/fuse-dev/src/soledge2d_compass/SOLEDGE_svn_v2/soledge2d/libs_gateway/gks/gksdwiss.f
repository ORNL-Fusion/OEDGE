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

        SUBROUTINE GKDWIS (FCTID,DX,DY,DIMX,IA,LR1,R1,LR2,R2,LC,CHARS)
C*  GKS Workstation Independent Segment Storage

        INTEGER LC, LR1, LR2
        INTEGER FCTID,DX,DY,DIMX,IA(*)
        REAL R1(*),R2(*)
        CHARACTER*(*) CHARS(*)

        SAVE

C*  include GKS symbol definitions
        INCLUDE 'gksdefs.i'

C *     Workstation State List
        INTEGER STATE
C               connection identifier
        INTEGER CONID
C               segment number
        INTEGER SGNUM

        CHARACTER*100 TMP
        INTEGER LTMP
        CHARACTER*132 STR(1)
        INTEGER L, MIA, MLR1, MLR2, MLC

        GOTO (999,999,  2,  3,  4,  5,999,999,999,333,
     *        333,333,333,333,333,333,333,333,333,333,
     *        333,333,333,333,333,333,333,333,333,333,
     *        333,333,333,333,333,333,333,333,333,333,
     *        333,999,999,999,999,999,999,999,333,333,
     *        333,333,333,333,333,999, 56, 57,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999,999,999,999,
     *        999,999,999,999,999,999,999) FCTID+1
        GOTO 999

C*  open workstation
    2   CONTINUE
        CALL GKTMP (LTMP, TMP)

        CONID = IA(2)
        CALL GKDSC (CONID, LTMP, TMP)

        SGNUM = 0
        GOTO 999

C*  close workstation
    3   CONTINUE
        OPEN (CONID, FILE=TMP(1:LTMP), FORM='UNFORMATTED',
     *      STATUS='UNKNOWN')
        CLOSE (CONID, STATUS='DELETE')
        GOTO 999

C*  activate workstation
    4   CONTINUE
        STATE = GACTIV
        GOTO 999

C*  deactivate workstation
    5   CONTINUE
        STATE = GINACT
        GOTO 999

C*  create segment
   56   CONTINUE
        SGNUM = IA(1)

        OPEN (CONID, FILE=TMP(1:LTMP), FORM='UNFORMATTED',
     *      STATUS='UNKNOWN')

        MIA = 1
        MLR1 = 1
        MLR2 = 1
        MLC = 1
        GOTO 999

C*  close segment
   57   CONTINUE
        SGNUM = 0

        CLOSE (CONID)

        DIMX = MIA
        LR1 = MLR1
        LR2 = MLR2
        LC = MLC
        GOTO 999

  333   CONTINUE
        IF (STATE .EQ. GACTIV) THEN
            IF (SGNUM .NE. 0) THEN
                DO 444 L = 1,LC
                    STR(L) = CHARS(L)
  444           CONTINUE
                MIA = MAX(MIA, DX*DY)
                MLR1 = MAX(MLR1, LR1)
                MLR2 = MAX(MLR2, LR2)
                MLC = MAX(MLC, LC)
                WRITE (CONID) SGNUM, FCTID,
     *              DX, DY, DIMX, (IA(L), L=1,DX*DY),
     *              LR1, (R1(L), L=1,LR1), LR2, (R2(L), L=1,LR2),
     *              LC, (STR(L), L=1,LC)
            END IF
        END IF

  999   CONTINUE

        RETURN
        END



        SUBROUTINE GKDCSG (WKID, SEGN, IA, R1, R2, CHARS)
C           copy segment to workstation

        INTEGER WKID, SEGN, IA(*)
        REAL R1(*), R2(*)
        CHARACTER*132 CHARS(*)

        INTEGER ICONID, ILTMP
        CHARACTER*(*) ITMP

        INTEGER L, LC, LR1, LR2
        INTEGER FCTID, DX, DY, DIMX

        INTEGER I, CONID, SGNUM
        CHARACTER TMP*100
        INTEGER LTMP

        DATA CONID /88/
        DATA TMP /'fort.88'/
        DATA LTMP /7/

        CALL GKDSWK (WKID)

        OPEN (CONID, FILE=TMP(1:LTMP), FORM='UNFORMATTED',
     *      STATUS='UNKNOWN')
        REWIND CONID

  333   CONTINUE
        READ (CONID, END=999, ERR=999) SGNUM, FCTID,
     *      DX, DY, DIMX, (IA(L), L=1,DX*DY),
     *      LR1, (R1(L), L=1,LR1), LR2, (R2(L), L=1,LR2),
     *      LC, (CHARS(L), L=1,LC)

        IF (SEGN .EQ. 0 .OR. SGNUM .EQ. SEGN) THEN

            GOTO (333,333,333,333,333,333,333,333,333,333,
     *            333,333, 12, 13, 14, 15, 16,333,333, 19,
     *             20, 21,333, 23, 24, 25,333, 27, 28, 29,
     *             30, 31, 32, 33, 34,333, 36, 37, 38,333,
     *            333,333,333,333,333,333,333,333,333, 49,
     *             50,333, 52, 53,333,333,333,333,333,333,
     *            333,333,333,333,333,333,333,333,333,333,
     *            333,333,333,333,333,333,333,333,333,333,
     *            333,333,333,333,333,333,333,333,333,333,
     *            333,333,333,333,333,333,333,333,333,333,
     *            333,333,333,333,333,333,333) FCTID+1
            GOTO 333
   12       CALL GPL (IA,R1,R2)
            GOTO 333
   13       CALL GPM (IA,R1,R2)
            GOTO 333
   14       DO 140, I = 132, 1, -1
               IF (CHARS(1)(I:I) .NE. ' ') THEN
                   CALL GTX (R1,R2,CHARS(1)(1:I))
                   GOTO 333
               END IF
  140       CONTINUE
            GOTO 333
   15       CALL GFA (IA,R1,R2)
            GOTO 333
   16       CALL GCA (R1(1),R2(1),R1(2),R2(2),DX,DY,1,1,DIMX,DY,IA)
            GOTO 333
   19       CALL GSLN (IA)
            GOTO 333
   20       CALL GSLWSC (R1)
            GOTO 333
   21       CALL GSPLCI (IA)
            GOTO 333
   23       CALL GSMK (IA)
            GOTO 333
   24       CALL GSMKSC (R1)
            GOTO 333
   25       CALL GSPMCI (IA)
            GOTO 333
   27       CALL GSTXFP (IA(1),IA(2))
            GOTO 333
   28       CALL GSCHXP (R1)
            GOTO 333
   29       CALL GSCHSP (R1)
            GOTO 333
   30       CALL GSTXCI (IA)
            GOTO 333
   31       CALL GSCHH (R1)
            GOTO 333
   32       CALL GSCHUP (R1,R2)
            GOTO 333
   33       CALL GSTXP (IA)
            GOTO 333
   34       CALL GSTXAL (IA(1),IA(2))
            GOTO 333
   36       CALL GSFAIS (IA)
            GOTO 333
   37       CALL GSFASI (IA)
            GOTO 333
   38       CALL GSFACI (IA)
            GOTO 333
   49       CALL GSWN (IA,R1(1),R1(2),R2(1),R2(2))
            GOTO 333
   50       CALL GSVP (IA,R1(1),R1(2),R2(1),R2(2))
            GOTO 333
   52       CALL GSELNT (IA)
            GOTO 333
   53       CALL GSCLIP (IA)
            GOTO 333

        ELSE
            GOTO 333
        END IF

  999   CONTINUE

        CLOSE (CONID)

        CALL GKDSWK (0)

        RETURN


        ENTRY GKDSC (ICONID, ILTMP, ITMP)
C               set connection identifier

        CONID = ICONID
        LTMP = ILTMP
        TMP = ITMP

        RETURN
        END
