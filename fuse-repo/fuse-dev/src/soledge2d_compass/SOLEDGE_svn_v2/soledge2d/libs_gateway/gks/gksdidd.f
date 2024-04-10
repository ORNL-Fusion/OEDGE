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

        SUBROUTINE GKDDLK (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,LC,CHARS)
C               devive driver link routine

        INTEGER I, LC, LRX, LRY
        INTEGER FCTID, DX, DY, DIMX, IA(*), ID
        REAL RX(*), RY(*)
        CHARACTER*(*) CHARS

        INTEGER OK
        PARAMETER (OK = 0)

        SAVE
C               GKS symbol definitions
        INCLUDE 'gksdefs.i'
C               GKS description table
        INCLUDE 'gksdescr.i'
C               GKS state list
        INCLUDE 'gksstate.i'

        INTEGER WKID, IPTR(2)
        LOGICAL ALL
        DOUBLE PRECISION DALIGN

        EQUIVALENCE (DALIGN, IPTR(1))

        DATA WKID /0/

        CALL GERSET (OK, -1)

        KERNEL = .TRUE.

        GOTO (1000,1000,2000,2000,2000,2000,2000,1000,2000,1000,
     *        2000,2000,1000,1000,1000,1000,1000,1000,1000,1000,
     *        1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     *        1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     *        1000,1000,1000,1000,1000,1000,1000,1000,2000,1000,
     *        1000,1000,1000,1000,2000,2000,1000,1000,1000,1000,
     *        1000,1000,1000,1000,1000,1000,1000,1000,1000,2000,
     *        1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     *        1000,2000,2000,1000,1000,1000,2000,1000,1000,1000,
     *        1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     *        1000,1000,1000,1000,1000,1000,1000 ) FCTID+1
        GOTO 5000

 1000   CONTINUE
        ALL = .TRUE.
        GOTO 3000

 2000   CONTINUE
        ALL = .FALSE.

 3000   CONTINUE

        IF (NOPWS .GT. 0) THEN

          DO 4000 I = 1,NOPWS

            IF (ALL .OR. IA(1) .EQ. SOPWS(I,1)) THEN

              IF (WKID .NE. 0) THEN
                IF (WKID .NE. SOPWS(I,1)) GOTO 4000
              END IF
                
              GOTO (4000, 201,
     *               204, 207,  82,  51,
     *                53,  72,  16,  17,
     *                61,  62,  63,  64,
     *               210, 211, 212, 213,
     *               214, 215, 216, 217,
     *               218,
     *               230, 231, 232, 233,
     *                 7,   8,   5,  41,
     *                38, 103, 104,  92,
     *                 2,   3, 101, 102), SOPWS(I,4)
              GOTO 4000

C*  TEKtronix 401x Series
C*  TAB 132/15-G
C*  Monterey MG200
C*  IBM PC
C*  DIGITAL LN03 PLUS
   72         CONTINUE
  201         CONTINUE
  204         CONTINUE
  207         CONTINUE
   38         CALL GKDTK1 (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,LC,CHARS)
              GOTO 4000

C*  TEKtronix 41xx Series
   82         CALL GKDTK2 (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,LC,CHARS)
              GOTO 4000

C*  HP-GL
   51         CONTINUE
   53         CALL GKDHP (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,LC,CHARS)
              GOTO 4000

C*  VT3xx
   16         CONTINUE
   17         CALL GKDVT3 (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,LC,CHARS)
              GOTO 4000

C*  PostScript, Display PostScript, LJ250
   61         CONTINUE
   62         CONTINUE
   63         CONTINUE
   64         CONTINUE
   92         CALL GKDPS (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,LC,CHARS)
              GOTO 4000

C*  Adobe's Portable Document Format (PDF)
  101         CONTINUE
  102         IPTR(1) = SOPWS(I,5)
              IPTR(2) = SOPWS(I,6)
              CALL GKDPDF (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,
     *          LC,CHARS,IPTR(1))
              GOTO 4000

C*  X Windows (Xlib)
  210         CONTINUE
  211         CONTINUE
  212         CONTINUE
  213         CONTINUE
  214         CONTINUE
  215         CONTINUE
  216         CONTINUE
  217         CONTINUE
  218         CONTINUE
  230         CONTINUE
  231         CONTINUE
  232         CONTINUE
  233         IPTR(1) = SOPWS(I,5)
              IPTR(2) = SOPWS(I,6)
              CALL GKDXW (FCTID,DX,DY,DIMX,IA(1),LRX,RX(1),LRY,RY(1),
     *          LC,CHARS,IPTR(1))
              GOTO 4000

C*  GKSM Output Metafile
    2         IPTR(1) = SOPWS(I,5)
              IPTR(2) = SOPWS(I,6)
              CALL GKDMFO (FCTID,DX,DY,DIMX,IA(1),LRX,RX(1),LRY,RY(1),
     *          LC,CHARS,IPTR(1))
              GOTO 4000

C*  GKSM Input Metafile
    3         IPTR(1) = SOPWS(I,5)
              IPTR(2) = SOPWS(I,6)
              CALL GKDMFI (FCTID,DX,DY,DIMX,IA(1),LRX,RX(1),LRY,RY(1),
     *          LC,CHARS,IPTR(1))
              GOTO 4000

C*  CGM (Computer Graphics Metafile)
    7         CONTINUE
    8         IPTR(1) = SOPWS(I,5)
              IPTR(2) = SOPWS(I,6)
              CALL GKDCGM (FCTID,DX,DY,DIMX,IA(1),LRX,RX(1),LRY,RY(1),
     *          LC,CHARS,IPTR(1))
              GOTO 4000

C*  WISS (Workstation Independent Segment Storage)
    5         CALL GKDWIS (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,LC,CHARS)
              GOTO 4000

C*  IBM PC w\ VGA
C*  VAX UIS
   41         CALL GKDGA (FCTID,DX,DY,DIMX,IA(1),LRX,RX(1),LRY,RY(1),
     *          LC,CHARS)
              GOTO 4000

C*  PBM (Portable BitMap)
  103         CONTINUE
  104         CALL GKDBM (FCTID,DX,DY,DIMX,IA,LRX,RX,LRY,RY,LC,CHARS)
              GOTO 4000

            END IF

 4000     CONTINUE

        END IF

 5000   CONTINUE

        KERNEL = .FALSE.

        RETURN


        ENTRY GKDSWK (ID)
C               select workstation

        WKID = ID

        RETURN
        END
