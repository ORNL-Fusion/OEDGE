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
C*  J.Heinen@FZ-juelich.de.
C*
C*
C       GLIGKS V4.5
C
C       Graphic Kernel System
C               ISO/DIS 7942

C       Reentrant FORTRAN 77 Function Interface for
C       GKS Version 7.4 Level 0b


C       written by      J.Heinen, FZ Juelich
C       date:           April 1984


        SUBROUTINE GOPKS (FILE, BUFFER)
C               Open GKS

        INTEGER NBSP
        PARAMETER (NBSP = 160)

        REAL FEPS
        PARAMETER (FEPS = 1.0E-06)

        INTEGER MAXSTR
        PARAMETER (MAXSTR = 255)

        INTEGER FILE, BUFFER, WKID, CONID, WTYPE, COFL, REFL, LDR
        CHARACTER*(*) DATREC(LDR)

        INTEGER N, DX, DY, SCOL, SROW, NCOL, NROW, COLIA(DX,DY)
        REAL PX(N), PY(N), QX, QY, RX, RY
        CHARACTER*(*) CHARS
        INTEGER NCHARS

        INTEGER INDEX, TYPE, FONT, PREC, PATH, ALH, ALV, STYLE
        REAL WIDTH, FACTOR, SPACE, HEIGHT, UX, UY

        INTEGER FLAG(13), CI
        REAL CR, CG, CB

        INTEGER TNR, CLSW
        REAL XMIN, XMAX, YMIN, YMAX

        INTEGER LCDNR, SKDNR, STDNR, CHDNR
        INTEGER LOSTR, STAT, NP, CHNR, DEFMOD, IRGMOD
        CHARACTER*(*) STR
        REAL SX(NP), SY(NP)

        INTEGER IDX, IA(4), LR1, LR2
        REAL R1(5), R2(5)
        INTEGER LC, I
        PARAMETER (LC=0)
        CHARACTER C*1

        INTEGER FUNID, DIMIDR, MAXODR, LENODR
        CHARACTER*80 IDR(DIMIDR), ODR(MAXODR)

        INTEGER SEGN

        REAL X0, Y0, TX, TY, PHI, FX, FY, TRAN(2,3)
        INTEGER ISW

        INTEGER WKIND, IPTR(2)
        INTEGER MIA, MR1, MR2, MC
        DOUBLE PRECISION DALIGN

        INTEGER ICONID, IWTYPE, IMASK, ITYPE, ISTYLE
        LOGICAL WISS
        REAL XORG, YORG, XPOINT, YPOINT, XSHIFT, YSHIFT, SINF, COSF

        LOGICAL OPENED

        CHARACTER TEXT*255
        INTEGER NTEXT

C               external functions
        LOGICAL GANY, GALL
        INTEGER GORD

C *     GKS Error State List
C               error file
        INTEGER ERRFIL
C               buffers
        INTEGER BUFF
C               identification of the GKS procedure which caused the
C               error detection
        INTEGER FCTID
C               error number
        INTEGER ERRNO

C               GKSGRAL marker types
        INTEGER GRALMK(-114:-101)

C               GKSGRAL hatch styles
        INTEGER GRALHA(-106:-101)
C               GDDM hatch styles
        INTEGER GDDMHA(-6:-1)

C               segment attributes, GKS attributes
        INTEGER SGATTR(20), GKATTR(20)

C               GKS symbol definitions
        INCLUDE 'gksdefs.i'
C               GKS description table
        INCLUDE 'gksdescr.i'
C               GKS state list
        INCLUDE 'gksstate.i'

        EQUIVALENCE (DALIGN, IPTR(1))

C               force the LINKER to search the object module for the
C               block data subprogram (required for VAX FORTRAN)
        EXTERNAL GKSDAT
 
C               backward compatibility flag for GR/GR3-Software
        COMMON /GKSOPT/ GR_BC
CDEC$ PSECT /GKSOPT/ NOSHR

        SAVE

        DATA GRALMK /-20,-15,-18,-5,-17,-3,-7,
     *               -1,-16,-19,-14,-12,-2,-6/

        DATA GRALHA /7,7,1,8,8,2/
        DATA GDDMHA /4,10,3,9,2,1/

        DATA IA(4) /0/


        FCTID = EOPKS

        INQUIRE (UNIT=FILE, OPENED=OPENED)
        IF (OPENED) THEN
          ERRFIL = FILE
        ELSE
          ERRFIL = 2 + 100
        END IF
        BUFF = BUFFER

        IF (STATE .NE. GGKCL) THEN
C               GKS not in proper state. GKS must be in the state GKCL
          ERRNO = 1
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE
C               call the devive driver link routine
          IDX = 1
          IA(1) = ERRFIL
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)


C               open file which contains character fonts
          CALL OPFONT(GR_BC)

          STATE = GGKOP

C               delete all workstation identifiers from the set of open
C               and the set of active workstations
          CALL GINI (SOPWS,MNOPWS)
          NOPWS = 0
          CALL GINI (SACWS,MNACWS)
          NACWS = 0

          X11 = .FALSE.
          WISS = .FALSE.

C               initialize aspect source flags
          DO 1 I = 1,13
            ASF(I) = GBUNDL
   1      CONTINUE

C               set default bundle table indices
          LINDEX = 1
          MINDEX = 1
          TINDEX = 1
          FINDEX = 1

C               set default global output attributes
          CHH = 0.01
          CHUP(1) = 0.
          CHUP(2) = 1.
          TXP = GRIGHT
          TXAL(1) = GAHNOR
          TXAL(2) = GAVNOR

C               set default polyline attributes
          LTYPE = GLSOLI
          LWIDTH = 1.
          PLCOLI = 1
C               set default polymarker attributes
          MTYPE = GPOINT
          MSZSC = 1.
          PMCOLI = 1
C               set default text attributes
          TXFONT = 1
          TXPREC = GSTRP
          CHXP = 1.
          CHSP = 0.
          TXCOLI = 1
C               set default fill area attributes
          INTS = GHOLLO
          STYLI = 1
          FACOLI = 1

C               initialize normalization transformations
          DO 2 CNTNR = 0,MXNTNR
            LNTNR(CNTNR) = CNTNR
            WINDOW(1,CNTNR) = 0.
            WINDOW(2,CNTNR) = 1.
            WINDOW(3,CNTNR) = 0.
            WINDOW(4,CNTNR) = 1.
            VIEWPT(1,CNTNR) = 0.
            VIEWPT(2,CNTNR) = 1.
            VIEWPT(3,CNTNR) = 0.
            VIEWPT(4,CNTNR) = 1.
C               set up normalization transformation
            CALL GSNT (CNTNR,WINDOW(1,CNTNR),VIEWPT(1,CNTNR))
   2      CONTINUE
          CNTNR = 0

C               set cliping indicator
          CLIP = GCLIP

C               set segment transformation
          MAT(1,1) = 1.
          MAT(1,2) = 0.
          MAT(2,1) = 0.
          MAT(2,2) = 1.
          MAT(1,3) = 0.
          MAT(2,3) = 0.
          CALL GSST (MAT)

          IA(4) = 0

        END IF

        RETURN



        ENTRY GCLKS
C               close GKS

        FCTID = ECLKS
        IF (STATE .NE. GGKOP) THEN
C               GKS not in proper state. GKS must be in the state GKOP
          ERRNO = 2
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

C               call the devive driver link routine
          IDX = 0
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,0,0,IA,LR1,R1,LR2,R2,LC,C)

          CALL CLFONT
          STATE = GGKCL

        END IF

        RETURN



        ENTRY GOPWK (WKID,CONID,WTYPE)
C               open workstation

        FCTID = EOPWK
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IMASK = WTYPE/65536 * 65536
            IWTYPE = MOD(WTYPE,65536)
            IF (IWTYPE .EQ. 0) THEN
C               get default workstation type
                CALL GTWSTY(IWTYPE)
            END IF

            IF (.NOT. GANY (LWSTY,NWSTY, IWTYPE)) THEN
C               specified workstation type is invalid
              ERRNO = 22
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is open
                ERRNO = 24
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

                IF (.NOT. GANY (SOPWS,MNOPWS, 0)) THEN
C               specified workstation cannot be opened
                  ERRNO = 26
                  CALL GERHND(ERRNO,FCTID,ERRFIL)

                ELSE

                  IF (IWTYPE .EQ. 5 .AND. WISS) THEN
C               WISS is already open
                    ERRNO = 28
                    CALL GERHND(ERRNO,FCTID,ERRFIL)

                  ELSE

                    ICONID = CONID
C               check file descriptor
                    CALL GKFD (ICONID)

C               add workstation identifier to the set of open
C               workstations
                    CALL GADD (WKID,SOPWS,MNOPWS)
                    NOPWS = NOPWS + 1

                    IPTR(1) = CONID
                    IPTR(2) = 0
                    CALL GKLTOI (CONID, IPTR(1))
                    WKIND = GORD (SOPWS,MNOPWS,WKID)
                    SOPWS(WKIND,2) = CONID
                    SOPWS(WKIND,3) = IWTYPE
                    SOPWS(WKIND,4) = GORD (LWSTY,NWSTY,IWTYPE)
                    SOPWS(WKIND,5) = IPTR(1)
                    SOPWS(WKIND,6) = IPTR(2)

                    IF (STATE .EQ. GGKOP) STATE = GWSOP

C               call the devive driver link routine
                    IDX = 3
                    IIA(1) = WKID
                    IIA(2) = ICONID
                    IIA(3) = IMASK + IWTYPE
                    LR1 = 0
                    LR2 = 0
                    CALL GKDDLK(FCTID,IDX,1,IDX,IIA,LR1,R1,LR2,R2,LC,C)
                    IF (IIA(1) .NE. 0 .OR. IIA(2) .NE. 0) THEN

                      SOPWS(WKIND,5) = IIA(1)
                      SOPWS(WKIND,6) = IIA(2)

                      IF ((IWTYPE.GE.210 .AND. IWTYPE.LE.215) .OR.
     *                     IWTYPE.EQ.218 .OR.
     *                    (IWTYPE.GE.230 .AND. IWTYPE.LE.233)) THEN
                        X11 = .TRUE.
                        WKIND = GORD(LWSTY,NWSTY, IWTYPE)
                        GDC(1,WKIND) = R1(1)
                        GDC(2,WKIND) = R2(1)
                        GRU(1,WKIND) = INT(R1(2))
                        GRU(2,WKIND) = INT(R2(2))
                      ELSE IF (IWTYPE.EQ.5) THEN
                        WISS = .TRUE.
                      END IF

                    ELSE

C               delete workstation identifier from the
C               set of open workstations
                      CALL GDEL (WKID,SOPWS,MNOPWS)
                      NOPWS = NOPWS - 1

                      IF (GALL (SOPWS,MNOPWS, 0)) STATE = GGKOP

C               open failed
                      ERRNO = 905
                      CALL GERHND(ERRNO,FCTID,ERRFIL)

                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GCLWK (WKID)
C               close workstation

        FCTID = ECLWK
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (GANY (SACWS,MNACWS, WKID)) THEN
C               specified workstation is active
                ERRNO = 29
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

C                 call the devive driver link routine
                IDX = 1
                IA(1) = WKID
                LR1 = 0
                LR2 = 0
                CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

                WKIND = GORD (SOPWS,MNOPWS,WKID)
                IWTYPE = SOPWS(WKIND,3)

                IF ((IWTYPE.GE.210 .AND. IWTYPE.LE.215) .OR.
     *               IWTYPE.EQ.218 .OR.
     *              (IWTYPE.GE.230 .AND. IWTYPE.LE.233)) THEN
                  X11 = .FALSE.
                ELSE IF (IWTYPE.EQ.5) THEN
                  WISS = .FALSE.
                END IF

C               delete workstation identifier from the set of open
C               workstations
                CALL GDEL (WKID,SOPWS,MNOPWS)
                NOPWS = NOPWS - 1

C               bump data structure
                IF (WKIND .LT. MNOPWS) THEN
                  DO 8 I = WKIND+1,MNOPWS
                    SOPWS(I-1,2) = SOPWS(I,2)
                    SOPWS(I-1,3) = SOPWS(I,3)
                    SOPWS(I-1,4) = SOPWS(I,4)
                    SOPWS(I-1,5) = SOPWS(I,5)
                    SOPWS(I-1,6) = SOPWS(I,6)
   8              CONTINUE
                END IF

                IF (GALL (SOPWS,MNOPWS, 0)) STATE = GGKOP

              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GACWK (WKID)
C               activate workstation

        FCTID = EACWK
        IF (STATE .NE. GWSOP .AND. STATE .NE. GWSAC) THEN
C               GKS not in proper state. GKS must be either in the state
C               WSOP or in the state WSAC
          ERRNO = 6
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (GANY (SACWS,MNACWS, WKID)) THEN
C               specified workstation is active
                ERRNO = 29
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

C               add workstation identifier to the set of active
C               workstations
                CALL GADD (WKID,SACWS,MNACWS)
                NACWS = NACWS + 1

                STATE = GWSAC
C               call the devive driver link routine
                IDX = 1
                IA(1) = WKID
                LR1 = 0
                LR2 = 0
                CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GDAWK (WKID)
C               deactivate workstation

        FCTID = EDAWK
        IF (STATE .NE. GWSAC) THEN
C               GKS not in proper state. GKS must be in the state WSAC
          ERRNO = 3
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SACWS,MNACWS, WKID)) THEN
C               specified workstation is not active
              ERRNO = 30
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

C               call the devive driver link routine
              IDX = 1
              IA(1) = WKID
              LR1 = 0
              LR2 = 0
              CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

C               delete workstation identifier from the set of active
C               workstations
              CALL GDEL (WKID,SACWS,MNACWS)

              NACWS = NACWS - 1
              IF (NACWS .EQ. 0) STATE = GWSOP

            END IF
          END IF
        END IF

        RETURN



        ENTRY GCLRWK (WKID,COFL)
C               clear workstation

        FCTID = ECLRWK
        IF (STATE .NE. GWSOP .AND. STATE .NE. GWSAC) THEN
C               GKS not in proper state. GKS must be either in the state
C               WSOP or in the state WSAC
          ERRNO = 6
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

C               call the devive driver link routine
              IDX = 2
              IA(1) = WKID
              IA(2) = COFL
              LR1 = 0
              LR2 = 0
              CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

            END IF
          END IF
        END IF

        RETURN



        ENTRY GUWK (WKID,REFL)
C               update workstation

        FCTID = EUWK
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

C               call the devive driver link routine
              IDX = 2
              IA(1) = WKID
              IA(2) = REFL
              LR1 = 0
              LR2 = 0
              CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

            END IF
          END IF
        END IF

        RETURN



        ENTRY GMSG (WKID,CHARS)
C               message

        FCTID = EMSG
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          NTEXT = MIN(LEN(CHARS),MAXSTR-1)
          DO 3 I = 1,NTEXT
            IF (CHARS(I:I) .NE. CHAR(NBSP)) THEN
              TEXT(I:I) = CHARS(I:I)
            ELSE
              TEXT(I:I) = ' '
            END IF
   3      CONTINUE
          I = NTEXT + 1
          TEXT(I:I) = CHAR(0)

C               call the devive driver link routine
          IDX = 1
          IA(1) = WKID
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,1,TEXT(1:NTEXT))

        END IF

        RETURN



        ENTRY GMSGS (WKID,NCHARS,CHARS)
C               message (FORTRAN 77 subset)

        FCTID = EMSG
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          NTEXT = MIN(NCHARS,MAXSTR-1)
          DO 4 I = 1,NTEXT
            IF (CHARS(I:I) .NE. CHAR(NBSP)) THEN
              TEXT(I:I) = CHARS(I:I)
            ELSE
              TEXT(I:I) = ' '
            END IF
   4      CONTINUE
          I = NTEXT + 1
          TEXT(I:I) = CHAR(0)

C               call the devive driver link routine
          IDX = 1
          IA(1) = WKID
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,1,TEXT(1:NTEXT))

        END IF

        RETURN



        ENTRY GPL (N,PX,PY)
C               polyline

        FCTID = EPL
        IF (STATE .LT. GWSAC) THEN
C               GKS not in proper state. GKS must be either in the state
C               WSAC or in the state SGOP
          ERRNO = 5
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE
          IF (N .LT. 2) THEN
C               number of points is invalid
            ERRNO = 100
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IDX = 1
            IA(1) = N
C               call the devive driver link routine
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,N,PX,N,PY,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GPM (N,PX,PY)
C               polymarker

        FCTID = EPM
        IF (STATE .LT. GWSAC) THEN
C               GKS not in proper state. GKS must be either in the state
C               WSAC or in the state SGOP
          ERRNO = 5
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (N .LT. 1) THEN
C               number of points is invalid
            ERRNO = 100
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IDX = 1
            IA(1) = N
C               call the devive driver link routine
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,N,PX,N,PY,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GTX (QX,QY,CHARS)
C               text

        FCTID = ETX
        IF (STATE .LT. GWSAC) THEN
C               GKS not in proper state. GKS must be either in the state
C               WSAC or in the state SGOP
          ERRNO = 5
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          NTEXT = MIN(LEN(CHARS),MAXSTR-1)
          DO 9 I = 1,NTEXT
            IF (CHARS(I:I) .NE. CHAR(NBSP)) THEN
              TEXT(I:I) = CHARS(I:I)
            ELSE
              TEXT(I:I) = ' '
            END IF
   9      CONTINUE
          I = NTEXT + 1
          TEXT(I:I) = CHAR(0)

          IDX = 0
          R1(1) = QX
          R2(1) = QY
C               call the devive driver link routine
          CALL GKDDLK(FCTID,IDX,0,0,IA,1,R1,1,R2,1,TEXT(1:NTEXT))

        END IF

        RETURN



        ENTRY GTXS (QX,QY,NCHARS,CHARS)
C               text (FORTRAN 77 subset)

        FCTID = ETX
        IF (STATE .LT. GWSAC) THEN
C               GKS not in proper state. GKS must be either in the state
C               WSAC or in the state SGOP
          ERRNO = 5
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          NTEXT = MIN(NCHARS,MAXSTR-1)
          DO 10 I = 1,NTEXT
            IF (CHARS(I:I) .NE. CHAR(NBSP)) THEN
              TEXT(I:I) = CHARS(I:I)
            ELSE
              TEXT(I:I) = ' '
            END IF
  10      CONTINUE
          I = NTEXT + 1
          TEXT(I:I) = CHAR(0)

          IDX = 0
          R1(1) = QX
          R2(1) = QY
C               call the devive driver link routine
          CALL GKDDLK(FCTID,IDX,0,0,IA,1,R1,1,R2,1,TEXT(1:NTEXT))

        END IF

        RETURN



        ENTRY GFA (N,PX,PY)
C               fill area

        FCTID = EFA
        IF (STATE .LT. GWSAC) THEN
C               GKS not in proper state. GKS must be either in the state
C               WSAC or in the state SGOP
          ERRNO = 5
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (N .LT. 3) THEN
C               number of points is invalid
            ERRNO = 100
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IDX = 1
            IA(1) = N
C               call the devive driver link routine
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,N,PX,N,PY,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GCA (QX,QY,RX,RY,DX,DY,SCOL,SROW,NCOL,NROW,COLIA)
C               cell array

        FCTID = ECA
        IF (STATE .LT. GWSAC) THEN
C               GKS not in proper state. GKS must be either in the state
C               WSAC or in the state SGOP
          ERRNO = 5
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (DX .LT. 1 .OR. DY .LT. 1) THEN
C               dimensions of color index array are invalid
            ERRNO = 84
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

C               call the devive driver link routine
            R1(1) = QX
            R2(1) = QY
            R1(2) = RX
            R2(2) = RY
            CALL GKDDLK(FCTID,NCOL,NROW,DX,COLIA(SCOL,SROW),
     *        2,R1,2,R2,0,C)

          END IF
        END IF

        RETURN



        ENTRY GSASF (FLAG)
C               set aspect source flags

        FCTID = ESASF
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          DO 13 I = 1,13
            ASF(I) = FLAG(I)
  13      CONTINUE

        END IF

        RETURN



        ENTRY GSPLI (INDEX)
C               set polyline index

        FCTID = ESPLI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .LT. 1 .OR. INDEX .GT. 5) THEN
C               Polyline index is invalid
            ERRNO = 60
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            LINDEX = INDEX
C               call the devive driver link routine
            IDX = 1
            IA(1) = INDEX
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSLN (TYPE)
C               set linetype

        FCTID = ESLN
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

C               this implementation supports eight linetypes
C               not defined by GKS Version 7.4
          IF (TYPE .LT. -30 .OR. TYPE .EQ. 0 .OR. TYPE .GT. GLDASD) THEN
C               Linetype is invalid
            ERRNO = 62
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (TYPE .NE. LTYPE) THEN

            LTYPE = TYPE
C               call the devive driver link routine
            IDX = 1
            IA(1) = TYPE
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSLWSC (WIDTH)
C               set linewidth scale factor

        FCTID = ESLWSC
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE IF (WIDTH .NE. LWIDTH) THEN

          LWIDTH = WIDTH
C               call the devive driver link routine
          IDX = 0
          LR1 = 1
          R1(1) = WIDTH
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,0,0,IA,LR1,R1,LR2,R2,LC,C)

        END IF

        RETURN



        ENTRY GSPLCI (INDEX)
C               set polyline color index

        FCTID = ESPLCI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .LT. 0) THEN
C               color index is invalid
            ERRNO = 65
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (INDEX .NE. PLCOLI) THEN

            PLCOLI = INDEX
C               call the devive driver link routine
            IDX = 1
            IA(1) = INDEX
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSPMI (INDEX)
C               set polymarker index

        FCTID = ESPMI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .LT. 1 .OR. INDEX .GT. 5) THEN
C               Polymarker index is invalid
            ERRNO = 64
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            MINDEX = INDEX
C               call the devive driver link routine
            IDX = 1
            IA(1) = INDEX
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSMK (TYPE)
C               set marker TYPE

        FCTID = ESMK
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (TYPE .GE. -114 .AND. TYPE .LE. -101) THEN
            ITYPE = GRALMK(TYPE)
          ELSE
            ITYPE = TYPE
          END IF
C               this implementation supports twenty marker types
C               not defined by GKS Version 7.4
          IF (ITYPE.LT.-20 .OR. ITYPE.EQ.0 .OR. ITYPE.GT.GXMARK) THEN
C               marker type is invalid
            ERRNO = 66
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (ITYPE .NE. MTYPE) THEN

            MTYPE = ITYPE
C               call the devive driver link routine
            IDX = 1
            IA(1) = ITYPE
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSMKSC (FACTOR)
C               set marker size scale factor

        FCTID = ESMKSC
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE IF (FACTOR .NE. MSZSC) THEN

          MSZSC = FACTOR
C               call the devive driver link routine
          IDX = 0
          LR1 = 1
          R1(1) = FACTOR
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,0,0,IA,LR1,R1,LR2,R2,LC,C)

        END IF

        RETURN



        ENTRY GSPMCI (INDEX)
C               set polymarker color index

        FCTID = ESPMCI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .LT. 0) THEN
C               color index is invalid
            ERRNO = 85
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (INDEX .NE. PMCOLI) THEN

            PMCOLI = INDEX
C               call the devive driver link routine
            IDX = 1
            IA(1) = INDEX
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSTXI (INDEX)
C               set text index

        FCTID = ESTXI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .LT. 1 .OR. INDEX .GT. 6) THEN
C               Text index is invalid
            ERRNO = 68
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            TINDEX = INDEX
C               call the devive driver link routine
            IDX = 1
            IA(1) = INDEX
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSTXFP (FONT,PREC)
C               set text font and precision

        FCTID = ESTXFP
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (FONT .EQ. 0) THEN
C               text font is invalid
            ERRNO = 70
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (FONT .NE. TXFONT .OR. PREC .NE. TXPREC) THEN

            TXFONT = FONT
            TXPREC = PREC
C               call the devive driver link routine
            IDX = 2
            IA(1) = FONT
            IA(2) = PREC
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSCHXP (FACTOR)
C               set character expansion factor

        FCTID = ESCHXP
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (FACTOR .LE. 0.) THEN
C               character expansion factor is invalid
            ERRNO = 72
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (FACTOR .NE. CHXP) THEN

            CHXP = FACTOR
C               call the devive driver link routine
            IDX = 0
            LR1 = 1
            R1(1) = FACTOR
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,0,0,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSCHSP (SPACE)
C               set character spacing

        FCTID = ESCHSP
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE IF (SPACE .NE. CHSP) THEN

          CHSP = SPACE
C               call the devive driver link routine
          IDX = 0
          LR1 = 1
          R1(1) = SPACE
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,0,0,IA,LR1,R1,LR2,R2,LC,C)

        END IF

        RETURN



        ENTRY GSTXCI (INDEX)
C               set text color index

        FCTID = ESTXCI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .LT. 0) THEN
C               color index is invalid
            ERRNO = 85
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (INDEX .NE. TXCOLI) THEN

            TXCOLI = INDEX
C               call the devive driver link routine
            IDX = 1
            IA(1) = INDEX
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSCHH (HEIGHT)
C               set character height

        FCTID = ESCHH
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (HEIGHT .LE. 0.) THEN
C               character height is invalid
            ERRNO = 73
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (HEIGHT .NE. CHH) THEN

            CHH = HEIGHT
C               call the devive driver link routine
            IDX = 0
            LR1 = 1
            R1(1) = HEIGHT
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,0,0,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSCHUP (UX,UY)
C               set character up vector

        FCTID = ESCHUP
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

            IF (ABS(UX) .LE. FEPS .AND. ABS(UY) .LE. FEPS) THEN 
C               character up vector is invalid
            ERRNO = 74
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (UX .NE. CHUP(1) .OR. UY .NE. CHUP(2)) THEN

            CHUP(1) = UX
            CHUP(2) = UY
C               call the devive driver link routine
            IDX = 0
            LR1 = 1
            R1(1) = UX
            LR2 = 1
            R2(1) = UY
            CALL GKDDLK(FCTID,IDX,0,0,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSTXP (PATH)
C               set text path

        FCTID = ESTXP
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE IF (PATH .NE. TXP) THEN

          TXP = PATH
C               call the devive driver link routine
          IDX = 1
          IA(1) = PATH
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

        END IF

        RETURN



        ENTRY GSTXAL (ALH,ALV)
C               set text alignment

        FCTID = ESTXAL
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE IF (ALH .NE. TXAL(1) .OR. ALV .NE. TXAL(2)) THEN

          TXAL(1) = ALH
          TXAL(2) = ALV
C               call the devive driver link routine
          IDX = 2
          IA(1) = ALH
          IA(2) = ALV
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

        END IF

        RETURN



        ENTRY GSFAI (INDEX)
C               set fill area index

        FCTID = ESFAI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .LT. 1 .OR. INDEX .GT. 5) THEN
C               Fill area index is invalid
            ERRNO = 75
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            FINDEX = INDEX
C               call the devive driver link routine
            IDX = 1
            IA(1) = INDEX
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSFAIS (STYLE)
C               set fill area interior style

        FCTID = ESFAIS
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE IF (STYLE .NE. INTS) THEN

          INTS = STYLE
C               call the devive driver link routine
          IDX = 1
          IA(1) = STYLE
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

        END IF

        RETURN



        ENTRY GSFASI (INDEX)
C               set fill area style index

        FCTID = ESFASI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .GE. -106 .AND. INDEX .LE. -101) THEN
            ISTYLE = GRALHA(INDEX)
          ELSE IF (INDEX .GE. -6 .AND. INDEX .LE. -1) THEN
            ISTYLE = GDDMHA(INDEX)
          ELSE
            ISTYLE = INDEX
          END IF

          IF (ISTYLE .LT. 0) THEN
C               style index is invalid
            ERRNO = 78
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (ISTYLE .NE. STYLI) THEN

            STYLI = ISTYLE
C               call the devive driver link routine
            IDX = 1
            IA(1) = ISTYLE
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSFACI (INDEX)
C               set fill area color index

        FCTID = ESFACI
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the states
C               GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (INDEX .LT. 0) THEN
C               color index is invalid
            ERRNO = 85
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (INDEX .NE. FACOLI) THEN

            FACOLI = INDEX
C               call the devive driver link routine
            IDX = 1
            IA(1) = INDEX
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSCR (WKID,CI,CR,CG,CB)
C               set color representation

        FCTID = ESCR
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states GKOP,WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (CI .LT. 0) THEN
C               color index is invalid
                ERRNO = 85
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

                IF (CR .LT. 0. OR. CR .GT. 1. .OR.
     *              CG .LT. 0. OR. CG .GT. 1. .OR.
     *              CB .LT. 0. OR. CB .GT. 1.) THEN
C               color is invalid
                  ERRNO = 88
                  CALL GERHND(ERRNO,FCTID,ERRFIL)

                ELSE

C               call the devive driver link routine
                  IDX = 2
                  IA(1) = WKID
                  IA(2) = CI
                  LR1 = 3
                  R1(1) = CR
                  R1(2) = CG
                  R1(3) = CB
                  LR2 = 0
                  CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

                END IF
              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GSWN (TNR,XMIN,XMAX,YMIN,YMAX)
C               set window

        FCTID = ESWN
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (TNR .LT. 1 .OR. TNR .GT. MXNTNR) THEN
C               transformation number is invalid
            ERRNO = 50
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (XMIN .GE. XMAX .OR. YMIN .GE. YMAX) THEN
C               rectangle definition is invalid
              ERRNO = 51
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              WINDOW(1,TNR) = XMIN
              WINDOW(2,TNR) = XMAX
              WINDOW(3,TNR) = YMIN
              WINDOW(4,TNR) = YMAX
C               set up normalization transformation
              CALL GSNT (TNR,WINDOW(1,TNR),VIEWPT(1,TNR))

C               call the devive driver link routine
              IDX = 1
              IA(1) = TNR
              LR1 = 2
              R1(1) = XMIN
              R1(2) = XMAX
              LR2 = 2
              R2(1) = YMIN
              R2(2) = YMAX
              CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

            END IF
          END IF
        END IF

        RETURN



        ENTRY GSVP (TNR,XMIN,XMAX,YMIN,YMAX)
C               set viewport

        FCTID = ESVP
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (TNR .LT. 1 .OR. TNR .GT. MXNTNR) THEN
C               transformation number is invalid
            ERRNO = 50
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (XMIN .GE. XMAX .OR. YMIN .GE. YMAX) THEN
C               rectangle definition is invalid
              ERRNO = 51
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (XMIN .LT. 0. .OR. XMAX .GT. 1. .OR.
     *            YMIN .LT. 0. .OR. YMAX .GT. 1.) THEN
C               viewport is not within the NDC unit square
                ERRNO = 52
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

                VIEWPT(1,TNR) = XMIN
                VIEWPT(2,TNR) = XMAX
                VIEWPT(3,TNR) = YMIN
                VIEWPT(4,TNR) = YMAX
C               set up normalization transformation
                CALL GSNT (TNR,WINDOW(1,TNR),VIEWPT(1,TNR))

C               call the devive driver link routine
                IDX = 1
                IA(1) = TNR
                LR1 = 2
                R1(1) = XMIN
                R1(2) = XMAX
                LR2 = 2
                R2(1) = YMIN
                R2(2) = YMAX
                CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GSELNT (TNR)
C               select normalization transformation

        FCTID = ESELNT
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (TNR .LT. 0 .OR. TNR .GT. MXNTNR) THEN
C               transformation number is invalid
            ERRNO = 50
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE IF (TNR .NE. CNTNR) THEN

            CNTNR = TNR
C               call the devive driver link routine
            IDX = 1
            IA(1) = TNR
            LR1 = 0
            LR2 = 0
            CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          END IF
        END IF

        RETURN



        ENTRY GSCLIP (CLSW)
C               set clipping indicator

        FCTID = ESCLIP
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE IF (CLSW .NE. CLIP) THEN

          CLIP = CLSW
C               call the devive driver link routine
          IDX = 1
          IA(1) = CLSW
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

        END IF

        RETURN



        ENTRY GSWKWN (WKID,XMIN,XMAX,YMIN,YMAX)
C               set workstation window

        FCTID = ESWKWN
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (XMIN .GE. XMAX .OR. YMIN .GE. YMAX) THEN
C               rectangle definition is invalid
                ERRNO = 51
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

                IF (XMIN .LT. 0. .OR. XMAX .GT. 1. .OR.
     *              YMIN .LT. 0. .OR. YMAX .GT. 1.) THEN
C               workstation window is not within the NDC unit square
                  ERRNO = 53
                  CALL GERHND(ERRNO,FCTID,ERRFIL)

                ELSE

C               call the devive driver link routine
                  IDX = 1
                  IA(1) = WKID
                  LR1 = 2
                  R1(1) = XMIN
                  R1(2) = XMAX
                  LR2 = 2
                  R2(1) = YMIN
                  R2(2) = YMAX
                  CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

                END IF
              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GSWKVP (WKID,XMIN,XMAX,YMIN,YMAX)
C               set workstation viewport

        FCTID = ESWKVP
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (XMIN .GE. XMAX .OR. YMIN .GE. YMAX) THEN
C               rectangle definition is invalid
                ERRNO = 51
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

C               call the devive driver link routine
                IDX = 1
                IA(1) = WKID
                LR1 = 2
                R1(1) = XMIN
                R1(2) = XMAX
                LR2 = 2
                R2(1) = YMIN
                R2(2) = YMAX
                CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GINLC (WKID,LCDNR,TNR,QX,QY,TYPE,XMIN,XMAX,YMIN,YMAX,
     *    LDR,DATREC)
C               initialize locator

        FCTID = EINLC
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

C               call the devive driver link routine
              IDX = 4
              IA(1) = WKID
              IA(2) = LCDNR
              IA(3) = TNR
              IA(4) = TYPE
              LR1 = 3
              R1(1) = QX
              R1(2) = XMIN
              R1(3) = XMAX
              LR2 = 3
              R2(1) = QY
              R2(2) = YMIN
              R2(3) = YMAX
              CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,1,DATREC(1))

            END IF
          END IF
        END IF

        RETURN



        ENTRY GRQLC (WKID,LCDNR,STAT,TNR,QX,QY)
C               request locator

        FCTID = ERQLC
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

C               call the devive driver link routine
              IDX = 2
              IA(1) = WKID
              IA(2) = LCDNR
              LR1 = 1
              LR2 = 1
              R1(1) = QX
              R2(1) = QY
              CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)
              STAT = IA(1)
              TNR = 0
              QX = R1(1)
              QY = R2(1)
            END IF
          END IF
        END IF

        RETURN



        ENTRY GRQSK (WKID,SKDNR,N,STAT,TNR,NP,SX,SY)
C               request stroke

        FCTID = ERQSK
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

C               call the devive driver link routine
              IDX = 3
              IA(1) = WKID
              IA(2) = SKDNR
              IA(3) = N
              LR1 = N
              LR2 = N
              CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,SX,LR2,SY,LC,C)
              STAT = IA(1)
              TNR = 0
              NP = IA(3)
            END IF
          END IF
        END IF

        RETURN



        ENTRY GRQST (WKID,STDNR,STAT,LOSTR,STR)
C               request string

        FCTID = ERQST
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

C               call the devive driver link routine
              IDX = 2
              IA(1) = WKID
              IA(2) = STDNR
              LR1 = 0
              LR2 = 0
              CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,1,STR)
              STAT = IA(1)
              LOSTR = IA(2)
            END IF
          END IF
        END IF

        RETURN



        ENTRY GCRSG (SEGN)
C               create segment

        FCTID = ECRSG
        IF (STATE .NE. GWSAC) THEN
C               GKS not in proper state. GKS must be in the state WSAC
          ERRNO = 3
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

C               call the devive driver link routine
          IDX = 1
          IA(1) = SEGN
          LR1 = 0
          LR2 = 0
          CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

          STATE = GSGOP
          OPSG = SEGN

C               save segment attributes
          DO 20, I = 1,20
            SGATTR(I) = GKSL(I)
   20     CONTINUE

        END IF

        RETURN



        ENTRY GASGWK (WKID, SEGN)
C               associate segment with workstation

        ENTRY GCSGWK (WKID, SEGN)
C               copy segment to workstation

        FCTID = ECSGWK
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. WISS) THEN
C               WISS is not open
              ERRNO = 27
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (.NOT. GANY (SACWS,MNACWS, WKID)) THEN
C               specified workstation is not active
                ERRNO = 30
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

C               save GKS attributes, restore segment attributes
                DO 30, I = 1,20
                  GKATTR(I) = GKSL(I)
                  GKSL(I) = SGATTR(I)
   30           CONTINUE

C               copy segment
                CALL GKCSG(WKID, SEGN, MIA, MR1, MR2, MC, 0)

C               restore saved GKS attributes
                DO 32, I = 1,20
                  GKSL(I) = GKATTR(I)
   32           CONTINUE

              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GRSGWK (WKID)
C               redraw all segments on workstation

        FCTID = ERSGWK
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (WISS) THEN

              IF (.NOT. GANY (SACWS,MNACWS, WKID)) THEN
C               specified workstation is not active
                ERRNO = 30
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

C               save GKS attributes, restore segment attributes
                DO 40, I = 1,20
                  GKATTR(I) = GKSL(I)
                  GKSL(I) = SGATTR(I)
   40           CONTINUE

C               redraw segments
                CALL GKCSG(WKID, 0, MIA, MR1, MR2, MC, 1)

C               restore saved GKS attributes
                DO 42, I = 1,20
                  GKSL(I) = GKATTR(I)
   42           CONTINUE

              END IF
            END IF
          END IF
        END IF

        RETURN



        ENTRY GCLSG
C               close segment

        FCTID = ECLSG
        IF (STATE .NE. GSGOP) THEN
C               GKS not in proper state. GKS must be in the state SGOP
          ERRNO = 4
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

C               call the devive driver link routine
          IDX = 0
          MIA = 0
          MR1 = 0
          MR2 = 0
          MC = 0
          CALL GKDDLK(FCTID,IDX,1,MIA,IA,MR1,R1,MR2,R2,MC,C)

          STATE = GWSAC

        END IF

        RETURN



        ENTRY GEVTM (X0,Y0,TX,TY,PHI,FX,FY,ISW,TRAN)
C               evaluate transformation matrix

        FCTID = EEVTM
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (ISW .EQ. GWC) THEN

            XORG = 0
            YORG = 0
            CALL GNT(XORG,YORG,CNTNR)

            XPOINT = X0
            YPOINT = Y0
            CALL GNT(XPOINT,YPOINT,CNTNR)

            XSHIFT = TX
            YSHIFT = TY
            CALL GNT(XSHIFT,YSHIFT,CNTNR)

            XSHIFT = XSHIFT - XORG
            YSHIFT = YSHIFT - YORG

          ELSE

            XPOINT = X0
            YPOINT = Y0
            XSHIFT = TX
            YSHIFT = TY

          END IF

          COSF = COS(PHI)
          SINF = SIN(PHI)

          TRAN(1,1) = FX * COSF
          TRAN(1,2) = FX * SINF
          TRAN(2,1) = -FY * SINF
          TRAN(2,2) = FY * COSF
          TRAN(1,3) = XPOINT + XSHIFT - XPOINT * TRAN(1,1) -
     *      YPOINT * TRAN(2,1)
          TRAN(2,3) = YPOINT + YSHIFT - XPOINT * TRAN(1,2) -
     *      YPOINT * TRAN(2,2)

        END IF

        RETURN



        ENTRY GSSGT (SEGN,TRAN)
C               set segment transformation

        FCTID = ESSGT
        IF (STATE .LT. GGKOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states GKOP,WSOP,WSAC,SGOP
          ERRNO = 8
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE
C               set up segment transformation
          MAT(1,1) = TRAN(1,1)
          MAT(1,2) = TRAN(1,2)
          MAT(2,1) = TRAN(2,1)
          MAT(2,2) = TRAN(2,2)
          MAT(1,3) = TRAN(1,3)
          MAT(2,3) = TRAN(2,3)
          CALL GSST (MAT)

        END IF

        RETURN


        ENTRY GESC (FUNID, DIMIDR, IDR, MAXODR, LENODR, ODR)
C               escape

        TEXT = IDR(1)

        IF (FUNID .EQ. -400) THEN
C               inquire input terminator
            LENODR = 1
            ODR(1) = CHAR(IA(4))
        ELSE
            LENODR = 0
        END IF

        RETURN


        ENTRY GRQCH (WKID, CHDNR, STAT, CHNR)

        IA(1) = WKID
        IA(2) = CHDNR
        CHNR = 0

        RETURN


        ENTRY GDSG (SEGN)

        IA(1) = SEGN

        RETURN


        ENTRY GGTITM (WKID,TYPE,LENODR)

        FCTID = EGTITM
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              WKIND = GORD (SOPWS,MNOPWS,WKID)
              IF (SOPWS(WKIND,3) .NE. 3) THEN
C               specified workstation is not of category MI
                ERRNO = 34
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

C                 call the devive driver link routine
                IDX = 1
                IA(1) = WKID
                LR1 = 0
                LR2 = 0
                CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,LC,C)

                TYPE = IA(1)
                LENODR = IA(2)

              END IF
            END IF
          END IF
        END IF

        RETURN


        ENTRY GRDITM (WKID,LENIDR,MAXODR,ODR)

        FCTID = ERDITM
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (WKID .LE. 0) THEN
C               specified workstation identifier is invalid
            ERRNO = 20
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (.NOT. GANY (SOPWS,MNOPWS, WKID)) THEN
C               specified workstation is not open
              ERRNO = 25
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              WKIND = GORD (SOPWS,MNOPWS,WKID)
              IF (SOPWS(WKIND,3) .NE. 3) THEN
C               specified workstation is not of category MI
                ERRNO = 34
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

C                 call the devive driver link routine
                IDX = 3
                IA(1) = WKID
                IA(2) = LENIDR
                IA(3) = MAXODR
                LR1 = 0
                LR2 = 0
                CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,MAXODR,ODR)

              END IF
            END IF
          END IF
        END IF

        RETURN


        ENTRY GIITM (TYPE,LENIDR,DIMIDR,IDR)

        FCTID = EIITM
        IF (STATE .LT. GWSOP) THEN
C               GKS not in proper state. GKS must be in one of the
C               states WSOP,WSAC,SGOP
          ERRNO = 7
          CALL GERHND(ERRNO,FCTID,ERRFIL)

        ELSE

          IF (TYPE .LE. 0) THEN
C               item type is not a valid GKS item
            ERRNO = 164
            CALL GERHND(ERRNO,FCTID,ERRFIL)

          ELSE

            IF (LENIDR .LT. 8) THEN
C               item length is invalid
              ERRNO = 161
              CALL GERHND(ERRNO,FCTID,ERRFIL)

            ELSE

              IF (DIMIDR .LT. 1) THEN
C               metafile item is invalid
                ERRNO = 163
                CALL GERHND(ERRNO,FCTID,ERRFIL)

              ELSE

C                 call the devive driver link routine
                IDX = 3
                IA(1) = TYPE
                IA(2) = LENIDR
                IA(3) = DIMIDR
                LR1 = 0
                LR2 = 0
                CALL GKDDLK(FCTID,IDX,1,IDX,IA,LR1,R1,LR2,R2,DIMIDR,IDR)

              END IF
            END IF
          END IF
        END IF

        RETURN


        ENTRY GSDS (WKID, DEFMOD, IRGMOD)

        IA(1) = WKID
        IA(2) = DEFMOD
        IA(3) = IRGMOD

        RETURN


        ENTRY GGDP

        RETURN
        END


        BLOCK DATA GKSDAT

C               GKS symbol definitions
        INCLUDE 'gksdefs.i'
C               GKS description table
        INCLUDE 'gksdescr.i'
C               GKS state list
        INCLUDE 'gksstate.i'

        SAVE
C               operating state value
        DATA STATE /0/
C               list of available workstation types
        DATA LWSTY /  200,   201,
     *  204,   207,    82,    51,
     *   53,    72,    16,    17,
     *   61,    62,    63,    64,
     *  210,   211,   212,   213,
     *  214,   215,   216,   217,
     *  218,
     *  230,   231,   232,   233,
     *    7,     8,     5,    41,
     *   38,   103,   104,    92,
     *    2,     3,   101,   102/
C               workstation categories
        DATA WSCAT /    0,     2,
     *    2,     2,     2,     0,
     *    0,     2,     2,     2,
     *    0,     0,     0,     0,
     *    0,     2,     2,     2,
     *    0,     0,     0,     0,
     *    0,
     *    0,     2,     2,     2,
     *    4,     4,     3,     2,
     *    0,     0,     0,     0,
     *    4,     5,     0,     0/
C               maximum display surface size (DC)
        DATA GDC /      1.0,1.0,        0.256,0.192,
     *  0.256,0.192,    0.256,0.192,    0.256,0.192,    0.272,0.19,
     *  0.272,0.19,     0.256,0.192,    0.24,0.144,     0.24,0.144,
     *  0.28575,0.19685,0.28575,0.19685,0.28575,0.19685,0.28575,0.19685,
     *  0.333,0.281,    0.333,0.281,    0.333,0.281,    0.333,0.281,
     *  0.333,0.281,    0.333,0.281,    0.333,0.281,    0.333,0.281,
     *  0.333,0.281,
     *  0.333,0.281,    0.333,0.281,    0.333,0.281,    0.333,0.281,
     *  1.0,1.0,        1.0,1.0,        1.0,1.0,        0.333,0.281,
     *  0.256,0.192,    0.254,0.2032,   0.254,0.2032,   0.2794,0.2032,
     *  1.0,1.0,        1.0,1.0,        0.288,0.1984,   0.288,0.1984/
C               maximum display surface size (raster units)
        DATA GRU /      1012,835,       1024,768,
     *  1024,768,       1024,768,       1024,768,       10870,7600,
     *  10870,7600,     1024,768,       800,480,        800,480,
     *  6750,4650,      6750,4650,      6750,4650,      6750,4650,
     *  1024,864,       1024,864,       1024,864,       1024,864,
     *  1024,864,       1024,864,       1024,864,       1024,864,
     *  1024,864,
     *  1024,864,       1024,864,       1024,864,       1024,864,
     *  65536,65536,    65536,65536,    32767,32767,    1024,864,
     *  1024,768,       720,576,        750,600,        1980,1440,
     *  65536,65536,    65536,65536,    810,558,        810,558/

        END
