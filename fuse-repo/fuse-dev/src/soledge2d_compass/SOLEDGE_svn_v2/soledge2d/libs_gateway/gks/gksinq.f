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
C       GKS inquiry functions
C

        SUBROUTINE GQOPS (OPSTA)
C               inquire operating state value

        INTEGER NBSP
        PARAMETER (NBSP = 160)

        INTEGER MAXSTR
        PARAMETER (MAXSTR = 255)

        INTEGER WKID,CONID,WTYPE,WKCAT,WKSTAT,N
        REAL QX,QY,RX,RY
        CHARACTER*(*) CHARS
        INTEGER NCHARS

        INTEGER FLAG(13),INDEX,TYPE,FONT,PREC,PATH,ALH,ALV,STYLE
        REAL FACTOR,SPACE,HEIGHT,UX,UY

        INTEGER TNR,CLSW,I
        REAL CLRT(4)

        REAL CPX,CPY,TRX(4),TRY(4)
        INTEGER SEGN

        INTEGER OFF,ON
        PARAMETER (OFF=0,ON=1)

        INTEGER OPSTA,LEV,ERRIND,NUMBER,MAXTNR
        INTEGER NCOLI,COLA,NPCI
        INTEGER WKIND,DCUNIT,OL,LX,LY
        REAL WN(4),VP(4)

        CHARACTER TEXT*255
        INTEGER NTEXT

C               external functions
        LOGICAL GANY
        INTEGER GORD

C               Predefined text fonts, text precisions, fill area
C               interior styles and fill area style indices
        INTEGER PFONT(6), PPREC(6), PINTS(5), PSTYLI(5)

C               GKS symbol definitions
        INCLUDE 'gksdefs.i'
C               GKS description table
        INCLUDE 'gksdescr.i'
C               GKS state list
        INCLUDE 'gksstate.i'

        SAVE

        DATA PFONT /1,1,1,-2,-3,-4/
        DATA PPREC /0,1,2,2,2,2/
        DATA PINTS /0,1,3,3,3/
        DATA PSTYLI /1,1,1,2,3/


        OPSTA = STATE

        RETURN



        ENTRY GQLVKS (ERRIND,LEV)
C               inquire level of GKS

        LEV = LEVEL
        ERRIND = OFF

        RETURN



        ENTRY GQEWK (N,ERRIND,NUMBER,WTYPE)
C               inquire list element of available workstation types

        IF (N .LT. 1 .OR. N .GT. NWSTY) THEN
          ERRIND = ON
        ELSE
          NUMBER = NWSTY
          WTYPE = LWSTY(N)
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQMNTN (ERRIND,MAXTNR)
C               inquire maximum normalization transformation number

        MAXTNR = MXNTNR
        ERRIND = OFF

        RETURN



        ENTRY GQOPWK (N,ERRIND,OL,WKID)
C               inquire set member of open workstations

        IF (N .LT. 1 .OR. N .GT. MNOPWS) THEN
          ERRIND = ON
        ELSE
          CALL GNUM (SOPWS,MNOPWS,OL)
          WKID = SOPWS(N,1)
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQACWK (N,ERRIND,OL,WKID)
C               inquire set member of active workstations

        IF (N .LT. 1 .OR. N .GT. MNACWS) THEN
          ERRIND = ON
        ELSE
          CALL GNUM (SACWS,MNACWS,OL)
          WKID = SACWS(N)
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQCHH (ERRIND,HEIGHT)
C               inquire character height

        HEIGHT = CHH
        ERRIND = OFF

        RETURN



        ENTRY GQCHUP (ERRIND,UX,UY)
C               inquire character up vector

        UX = CHUP(1)
        UY = CHUP(2)
        ERRIND = OFF

        RETURN



        ENTRY GQTXP (ERRIND,PATH)
C               inquire text path

        PATH = TXP
        ERRIND = OFF

        RETURN



        ENTRY GQTXAL (ERRIND,ALH,ALV)
C               inquire text alignment

        ALH = TXAL(1)
        ALV = TXAL(2)
        ERRIND = OFF

        RETURN



        ENTRY GQASF (ERRIND,FLAG)
C               inquire aspect source flags

        DO 1 I = 1,13
          FLAG(I) = ASF(I)
   1    CONTINUE
        ERRIND = OFF

        RETURN



        ENTRY GQPLI (ERRIND,INDEX)
C               inquire polyline index

        INDEX = LINDEX
        ERRIND = OFF

        RETURN



        ENTRY GQLN (ERRIND,TYPE)
C               inquire line type

        IF (KERNEL) THEN
          IF (ASF(1) .EQ. GINDIV) THEN
            TYPE = LTYPE
          ELSE
            TYPE = LINDEX
          END IF
        ELSE
          TYPE = LTYPE
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQLWSC (ERRIND,FACTOR)
C               inquire linewidth scale factor

        IF (KERNEL) THEN
          IF (ASF(2) .EQ. GINDIV) THEN
            FACTOR = LWIDTH
          ELSE
            FACTOR = 1.0
          END IF
        ELSE
          FACTOR = LWIDTH
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQPLCI (ERRIND,INDEX)
C               inquire polyline color index

        IF (KERNEL) THEN
          IF (ASF(3) .EQ. GINDIV) THEN
            INDEX = PLCOLI
          ELSE
            INDEX = 1
          END IF
        ELSE
          INDEX = PLCOLI
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQPMI (ERRIND,INDEX)
C               inquire polymarker index

        INDEX = MINDEX
        ERRIND = OFF

        RETURN



        ENTRY GQMK (ERRIND,TYPE)
C               inquire marker type

        IF (KERNEL) THEN
          IF (ASF(4) .EQ. GINDIV) THEN
            TYPE = MTYPE
          ELSE
            TYPE = MINDEX
          END IF
        ELSE
          TYPE = MTYPE
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQMKSC (ERRIND,FACTOR)
C               inquire marker size scale factor

        IF (KERNEL) THEN
          IF (ASF(5) .EQ. GINDIV) THEN
            FACTOR = MSZSC
          ELSE
            FACTOR = 1.0
          END IF
        ELSE
          FACTOR = MSZSC
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQPMCI (ERRIND,INDEX)
C               inquire polymarker color index

        IF (KERNEL) THEN
          IF (ASF(6) .EQ. GINDIV) THEN
            INDEX = PMCOLI
          ELSE
            INDEX = 1
          END IF
        ELSE
          INDEX = PMCOLI
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQTXI (ERRIND,INDEX)
C               inquire text index

        INDEX = TINDEX
        ERRIND = OFF

        RETURN



        ENTRY GQTXFP (ERRIND,FONT,PREC)
C               inquire text font and precision

        IF (KERNEL) THEN
          IF (ASF(7) .EQ. GINDIV) THEN
            FONT = TXFONT
            PREC = TXPREC
          ELSE
            FONT = PFONT(TINDEX)
            PREC = PPREC(TINDEX)
          END IF
        ELSE
          FONT = TXFONT
          PREC = TXPREC
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQCHXP (ERRIND,FACTOR)
C               inquire character expansion factor

        IF (KERNEL) THEN
          IF (ASF(8) .EQ. GINDIV) THEN
            FACTOR = CHXP
          ELSE
            FACTOR = 1.0
          END IF
        ELSE
          FACTOR = CHXP
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQCHSP (ERRIND,SPACE)
C               inquire character spacing

        IF (KERNEL) THEN
          IF (ASF(9) .EQ. GINDIV) THEN
            SPACE = CHSP
          ELSE
            SPACE = 0.0
          END IF
        ELSE
          SPACE = CHSP
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQTXCI (ERRIND,INDEX)
C               inquire text color index

        IF (KERNEL) THEN
          IF (ASF(10) .EQ. GINDIV) THEN
            INDEX = TXCOLI
          ELSE
            INDEX = 1
          END IF
        ELSE
          INDEX = TXCOLI
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQFAI (ERRIND,INDEX)
C               inquire fill area index

        INDEX = FINDEX
        ERRIND = OFF

        RETURN



        ENTRY GQFAIS (ERRIND,STYLE)
C               inquire fill area interior style

        IF (KERNEL) THEN
          IF (ASF(11) .EQ. GINDIV) THEN
            STYLE = INTS
          ELSE
            STYLE = PINTS(FINDEX)
          END IF
        ELSE
          STYLE = INTS
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQFASI (ERRIND,INDEX)
C               inquire fill area style index

        IF (KERNEL) THEN
          IF (ASF(12) .EQ. GINDIV) THEN
            INDEX = STYLI
          ELSE
            INDEX = PSTYLI(FINDEX)
          END IF
        ELSE
          INDEX = STYLI
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQFACI (ERRIND,INDEX)
C               inquire fill area color index

        IF (KERNEL) THEN
          IF (ASF(13) .EQ. GINDIV) THEN
            INDEX = FACOLI
          ELSE
            INDEX = 1
          END IF
        ELSE
          INDEX = FACOLI
        END IF
        ERRIND = OFF

        RETURN



        ENTRY GQCNTN (ERRIND,TNR)
C               inquire current normalization transformation number

        TNR = CNTNR
        ERRIND = OFF

        RETURN



        ENTRY GQNT (TNR,ERRIND,WN,VP)
C               inquire normalization transformation

        IF (TNR .LT. 0 .OR. TNR .GT. MXNTNR) THEN
          ERRIND = ON
        ELSE
          DO 2 I = 1,4
            WN(I) = WINDOW(I,TNR)
            VP(I) = VIEWPT(I,TNR)
   2      CONTINUE
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQCLIP (ERRIND,CLSW,CLRT)
C               inquire clipping indicator

        CLSW = CLIP

        DO 3 I = 1,4
          IF (CLSW .EQ. GCLIP) THEN
            CLRT(I) = VIEWPT(I,CNTNR)
          ELSE
            CLRT(I) = VIEWPT(I,0)
          END IF
   3    CONTINUE

        ERRIND = OFF

        RETURN



        ENTRY GQWKC (WKID,ERRIND,CONID,WTYPE)
C               inquire workstation connection and type

        IF (.NOT. GANY(SOPWS,MNOPWS,WKID)) THEN
          ERRIND = ON
        ELSE
          WKIND = GORD(SOPWS,MNOPWS,WKID)
          CONID = SOPWS(WKIND,2)
          WTYPE = SOPWS(WKIND,3)
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQWKS (WKID,ERRIND,WKSTAT)
C               inquire workstation state

        IF (.NOT. GANY(SOPWS,MNOPWS,WKID)) THEN
          ERRIND = ON
        ELSE
          IF (GANY(SACWS,MNACWS,WKID)) THEN
            WKSTAT = GACTIV
          ELSE
            WKSTAT = GINACT
          END IF
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQWKCA (WTYPE,ERRIND,WKCAT)
C               inquire workstation category

        IF (.NOT. GANY(LWSTY,NWSTY,WTYPE)) THEN
          ERRIND = ON
        ELSE
          WKIND = GORD(LWSTY,NWSTY,WTYPE)
          WKCAT = WSCAT(WKIND)
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQCF (WTYPE,ERRIND,NCOLI,COLA,NPCI)
C               inquire workstation category

        IF (.NOT. GANY(LWSTY,NWSTY,WTYPE)) THEN
          ERRIND = ON
        ELSE
          NCOLI = 980
          COLA = GCOLOR
          NPCI = 8
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQDSP (WTYPE,ERRIND,DCUNIT,RX,RY,LX,LY)
C               inquire display space size
        ENTRY GQMDS (WTYPE,ERRIND,DCUNIT,RX,RY,LX,LY)
C               inquire maximum display surface

        IF (.NOT. GANY(LWSTY,NWSTY,WTYPE)) THEN
          ERRIND = ON
        ELSE
          WKIND = GORD(LWSTY,NWSTY,WTYPE)
          DCUNIT = GMETRE
          RX = GDC(1,WKIND)
          RY = GDC(2,WKIND)
          LX = GRU(1,WKIND)
          LY = GRU(2,WKIND)
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQTXX (WKID,QX,QY,CHARS,ERRIND,CPX,CPY,TRX,TRY)
C               inquire text extend

        IF (.NOT. GANY(SOPWS,MNOPWS,WKID)) THEN
          ERRIND = ON
        ELSE
          NTEXT = MIN(LEN(CHARS),MAXSTR-1)
          DO 4 I = 1,NTEXT
            IF (CHARS(I:I) .NE. CHAR(NBSP)) THEN
              TEXT(I:I) = CHARS(I:I)
            ELSE
              TEXT(I:I) = ' '
            END IF
   4      CONTINUE
          I = NTEXT + 1
          TEXT(I:I) = CHAR(0)
          CALL GQTEXT (QX,QY,TEXT(1:NTEXT),CPX,CPY,TRX,TRY)
          ERRIND = OFF
        END IF

        RETURN



        ENTRY GQTXXS (WKID,QX,QY,NCHARS,CHARS,ERRIND,CPX,CPY,TRX,TRY)
C               inquire text extend (FORTRAN 77 subset)

        IF (.NOT. GANY(SOPWS,MNOPWS,WKID)) THEN
          ERRIND = ON
        ELSE
          NTEXT = MIN(NCHARS,MAXSTR-1)
          DO 5 I = 1,NTEXT
            IF (CHARS(I:I) .NE. CHAR(NBSP)) THEN
              TEXT(I:I) = CHARS(I:I)
            ELSE
              TEXT(I:I) = ' '
            END IF
   5      CONTINUE
          I = NTEXT + 1
          TEXT(I:I) = CHAR(0)
          CALL GQTEXT (QX,QY,TEXT(1:NTEXT),CPX,CPY,TRX,TRY)
          ERRIND = OFF
        END IF

        RETURN


        ENTRY GQOPSG (ERRIND,SEGN)
C           inquire name of open segment

        IF (STATE .EQ. GSGOP) THEN
          SEGN = OPSG
          ERRIND = OFF
        ELSE
          ERRIND = ON
        END IF

        RETURN


        ENTRY GQSGWK (WKID,N,ERRIND,OL,SEGN)

        ERRIND = ON

        RETURN
        END
