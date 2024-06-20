C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C ----------------------------------------------------------------------
C  GRAXS - ACHSENPROGRAMM
C  AUTOR: G. EGERER
C  DATUM:  7.11.83
C  AENDERUNGEN:
C      7.12.83
C     18. 9.89   MOEGLICHKEIT ZUM UNTERDRUECKEN DER ACHSENBESCHRIFTUNG
C                MITTELS DER OPTIONEN "XLABEL" UND "YLABEL" HINZUGEFUEGT
C     02.11.89   TEXTX und TEXTXY muessen jetzt CHARACTER sein.
C     18.01.90   auch OPTION muss jetzt CHARACTER sein.
C     12.07.91   Netzlinien werden vom letzten GRDSH bestimmt (Groten).
C     12.07.91   Wenn Pixel innen: Label naeher an den Achsen.
C     20.08.91   Netzlinien nicht auf Achsen.
C     10.09.91   UPDATE: GROTEN
C     29.01.92   Groten: Beschriftungen besser zentriert mit GSTXAL
C     20.09.93   Groten: Neue Option: statt F | FORMAT auch G | GENAU
C                Das Format wird nicht verkuerzt, so wie geschrieben.
C                Wenn man in OPTION statt F=( ... ) bzw. FORMAT=( ... )
C                                    nun  G=( ... ) bzw.  GENAU=( ... )
C                schreibt, sollte das Format mit fuehrenden und
C                beendenden Nullen geschrieben werden.
C     20.09.93   Groten: Neue Option: R | RAND
C                Wenn man die Beschriftung der Achsen ueber den
C                Achsenrand hinaus zulassen will, muss man R oder RAND
C                in die Option schreiben (wie I topr INNEN).
C
C ----------------------------------------------------------------------
      SUBROUTINE GRAXS(LOPT,OPTION,LTXTX,TEXTX,LTXTXY,TEXTXY)

C                           * EINGABEPARAMETER:
C                           * LOPT   - LAENGE DER Zeichenkette
C                           *          "OPTION"
C                           * OPTION - Zeichenkette, DIE ALLE
C                           *          INFORMATIONEN BZGL. AUSSEHEN DER
C                           *          ACHSE ENTHAELT
C                           * LTXTX  - LAENGE DER Zeichenkette
C                           *          "TEXTX"
C                           * TEXTX  - TEXT, DER AN DIE X-ACHSE
C                           *          GESCHRIEBEN WERDEN SOLL
C                           * LTXTXY - LAENGE DER Zeichenkette
C                           *          "TEXTXY"
C                           * TEXTXY - TEXT, DER AN DIE Y-ACHSE
C                           *          GESCHRIEBEN WERDEN SOLL
      INTEGER     LOPT, LTXTX,  LTXTXY
      character(len=*) textx,textxy
      CHARACTER(len=*) option

C                           * COMMON-VARIABLEN:
      COMMON  /GRPP/ SCLCM(2,2), SCLV(2,2), FCTR, DSH(3), INTLIN,
     $               CHRHEI, CHRANG, INTSYM, ISTEP, DUMMY
CDEC$ PSECT /GRPP/ NOSHR
      INTEGER INTLIN, INTSYM, ISTEP,lop,ltxx,ltxxy
      REAL    SCLCM, SCLV, FCTR, DSH, CHRHEI, CHRANG, DUMMY
      COMMON      /GRAXCB/ CPRCTL,GENAU,RANDB
CDEC$ PSECT /GRAXCB/ NOSHR
      CHARACTER(len=1) CPRCTL,GENAU,RANDB
      COMMON  /GRAXDS/ SVDSH1, SVDSH2, SVDSH3
CDEC$ PSECT /GRAXDS/ NOSHR
      SAVE /GRPP/,/GRAXCB/, /GRAXDS/
      REAL             SVDSH1, SVDSH2, SVDSH3

C                           * PROGRAMMVARIABLEN:
      LOGICAL LEDGE(2), LINNEN, LXACHS
      INTEGER I, IACHSE(2), IANZFU, IBSCHR, ICOMPL, IFAK, IFORMT(2),
     $        IGROBU, ILABEL(2), INTS, IP(2), IPRECS(2), ISTFU,
     $        ITBLFU(0:3,3), ITEXT, J, LTEXT, MANKBZ, N
      REAL    AMINA, AXLEN(2), CHRAND, CHRSZC, DIFFGU, E0N, FAKCMV(2),
     $        FRSTGU, HELP, PKTLU(2), PKTRO(2), RANGLE, SLENGU

C                           * BENAMTE PROGRAMMKONSTANTEN:
      REAL      CPI, CRNIL
      PARAMETER (CPI = 3.1415926536, CRNIL = -1.)

C                           * INITIALISIERUNGEN:
C                           *    TABELLE, DIE FUER JEDE STUFE DER FEIN-
C                           *    UNTERTEILUNG (U) JEWEILS DIE ANZAHL DER
C                           *    FEINUNTERTEILUNGSSTRICHE ZWISCHEN 2
C                           *    GROBUNTERTEILUNGSSTRICHEN FUER JEDE
C                           *    GROBUNTERTEILUNG ENTHAELT:
C                           *       GROBUNTERTEILUNG = 1, 2, 5 *10**N
C                           *               U = 0      0  0  0
C                           *               U = 1      1  1  4
C                           *               U = 2      4  3  4
C                           *               U = 3      9  9  9
      DATA ((ITBLFU(I,J), J = 1,3), I = 0,3) /0, 0, 0,
     $                                        1, 1, 4,
     $                                        4, 3, 4,
     $                                        9, 9, 9/
C                           * ------------------------------------------

      CPRCTL = '1'

      CHRSZC = CHRHEI
      CHRAND = CHRANG
      INTS = INTSYM

C---- DIE NAECHSTEN ZEILEN NEU 18.01.90 GROTEN
      IF ( LOPT.LT.0 .OR. LOPT.GT.LEN(OPTION) ) THEN
         lop=LEN(option)
      ELSE
         lop=lopt
      ENDIF
      IF ( LTXTX.LT.0 .OR. LTXTX.GT.LEN(TEXTX) ) THEN
         ltxx=LEN(textx)
      ELSE
         ltxx=ltxtx
      ENDIF
      IF ( LTXTXY.LT.0 .OR. LTXTXY.GT.LEN(TEXTXY) ) THEN
         ltxxy=LEN(textxy)
      ELSE
         ltxxy=ltxtxy
      ENDIF
C                           * PRUEFEN UND BEREITSTELLEN DER INFORMA-
C                           * TIONEN DER Zeichenkette "OPTION"
      CALL GRAXOP(lop,OPTION,IACHSE,LINNEN,IFORMT,AMINA,ILABEL,ISTFU,
     $   ITEXT)
      IF ( AMINA .EQ. CRNIL ) AMINA = 2.*CHRSZC

      SVDSH1 = DSH(1)
      SVDSH2 = DSH(2)
      SVDSH3 = DSH(3)
      CALL GRDSH(1.,0.,1.)

      DO 1, I = 1,2
C                           * LAENGE DER X- (I = 1) BZW. Y-ACHSE (I = 2)
C                           * IN CM
         AXLEN(I) = ABS(SCLCM(I,2)-SCLCM(I,1))
C                           * FAKTOR ZUR UMRECHNUNG VON CM-ANGABEN IN
C                           * BENUTZERWERTE FUER X- BZW. Y-ACHSE
         FAKCMV(I) = ABS(SCLV(I,2)-SCLV(I,1))/AXLEN(I)
    1 CONTINUE
C                           * WINKEL EINES ZEICHENS IN BOGENMASS
      RANGLE = CPI/180.*CHRAND
      LXACHS = .TRUE.
C                           * ERZEUGEN DER ACHSEN (I = 1: X-ACHSE,
C                           * I = 2: Y-ACHSE):
      DO 2, I = 1,2
C                           * SOLL X- BZW. Y-ACHSE ERZEUGT WERDEN?
         IF ( IACHSE(I) .NE. 0 ) THEN
C                           * GROBUNTERTEILUNG FUER DIE ACHSE ERMITTELN
            CALL GRAXGU(SCLV(I,1),SCLV(I,2),AXLEN(I),CHRSZC,RANGLE,
     $         AMINA,IFORMT,MANKBZ,IFAK,N,IPRECS)
            E0N = 10.**N
C                           * ABSZISSE DES 1. GROBUNTERTEILUNGSSTRICHES
            FRSTGU = REAL(MANKBZ)*E0N
C                           * ABSTAND DER GROBUNTERTEILUNGSSTRICHE
            DIFFGU = REAL(IFAK)*E0N
C                           * ANZAHL DER FEINUNTERTEILUNGSSTRICHE
C                           * ZWISCHEN 2 GROBUNTERTEILUNGSSTRICHEN
            IF ( IFAK .EQ. 5 ) THEN
               IANZFU = ITBLFU(ISTFU,3)
            ELSE
               IANZFU = ITBLFU(ISTFU,IFAK)
            ENDIF
C                           * SOLLEN AN DEN GROBUNTERTEILUNGSSTRICHEN
C                           * NETZLINIEN GEZEICHNET WERDEN?
            IF ( IACHSE(I) .LT. 0 ) THEN
               IGROBU = 2
               IACHSE(I) = -IACHSE(I)
            ELSE
C                           * ES SOLLEN KEINE NETZLINIEN GEZEICHNET
C                           * WERDEN
               IGROBU = 1
            ENDIF
            ICOMPL = 3-I
C                           * LAENGE EINES GROBUNTERTEILUNGSSTRICHES IN
C                           * BENUTZEREINHEITEN
            SLENGU = CHRSZC*FAKCMV(ICOMPL)
C                           * LEDGE(1) = .TRUE., FALLS INNENLIEGENDE
C                           * UNTERTEILUNGSSTRICHE AUF DEM LINKEN RAND
C                           * GEZEICHNET WERDEN SOLLEN
            LEDGE(1) = ( ABS(IACHSE(ICOMPL)) .EQ. 2 ) .OR.
     $         ( IACHSE(ICOMPL) .EQ. 0 )
C                           * LEDGE(2) = .TRUE., FALLS INNENLIEGENDE
C                           * UNTERTEILUNGSSTRICHE AUF DEM RECHTEN RAND
C                           * GEZEICHNET WERDEN SOLLEN
            LEDGE(2) = ( ABS(IACHSE(ICOMPL)) .LE. 1 )
C                           * ABSZISSE DES LINKEN RANDES DER ACHSE
            PKTLU(1) = SCLV(I,1)
C                           * ABSZISSE DES RECHTEN RANDES DER ACHSE
            PKTRO(1) = SCLV(I,2)
C                           * SOLL NUR  E I N E  X- BZW. Y-ACHSE
C                           * GEZEICHNET WERDEN?
            IF ( IACHSE(I) .NE. 3 ) THEN
C                           * ORDINATE DER ACHSE
               PKTLU(2) = SCLV(ICOMPL,IACHSE(I))
C                           * ORDINATE DER ENDPUNKTE VON EVTL. ZU
C                           * ZEICHNENDEN NETZLINIEN
               PKTRO(2)= SCLV(ICOMPL,3-IACHSE(I))
C                           * ACHSE ZEICHNEN
               CALL GRAXDR(LXACHS,PKTLU,PKTRO,LEDGE,LINNEN,IGROBU,
     $            FRSTGU,DIFFGU,IANZFU,SLENGU,FAKCMV(I))
C                           * SOLL DIE ACHSENBESCHRIFTUNG FUER DIE ACHSE
C                           * UNTERDRUECKT WERDEN?
               IF (  ( ILABEL(I) .EQ. -3 ) .OR.
     $            ( ILABEL(I) .EQ. -IACHSE(I) )  ) THEN
                     IP(1) = 0
                     IP(2) = 0
               ELSE
C                           * DIE ACHSE SOLL BESCHRIFTET WERDEN
                  IP(1) = IPRECS(1)
                  IP(2) = IPRECS(2)
               ENDIF
            ELSE
C                           * ES SOLLEN  Z W E I  X- BZW. Y-ACHSEN
C                           * GEZEICHNET WERDEN
C                           * ORDINATE DER UNTEREN BZW. LINKEN ACHSE
               PKTLU(2) = SCLV(ICOMPL,1)
C                           * ORDINATE DER ENDPUNKTE VON EVTL. ZU
C                           * ZEICHNENDEN NETZLINIEN
               PKTRO(2) = SCLV(ICOMPL,2)
C                           * UNTERE BZW. LINKE ACHSE ZEICHNEN
               CALL GRAXDR(LXACHS,PKTLU,PKTRO,LEDGE,LINNEN,IGROBU,
     $            FRSTGU,DIFFGU,IANZFU,SLENGU,FAKCMV(I))
C                           * VERTAUSCHEN DER ORDINATEN
               HELP = PKTLU(2)
               PKTLU(2) = PKTRO(2)
               PKTRO(2) = HELP
C                           * OBERE BZW. RECHTE ACHSE IST MIT STRICHEN,
C                           * JEDOCH NICHT MIT NETZLINIEN, AN DER GROB-
C                           * UNTERTEILUNG ZU ZEICHNEN
               IGROBU = 1
C                           * OBERE BZW. RECHTE ACHSE ZEICHNEN
               CALL GRAXDR(LXACHS,PKTLU,PKTRO,LEDGE,LINNEN,IGROBU,
     $            FRSTGU,DIFFGU,IANZFU,SLENGU,FAKCMV(I))
C                           * SOLL DIE ACHSENBESCHRIFTUNG FUER BEIDE
C                           * ACHSEN UNTERDRUECKT WERDEN?
               IF ( ILABEL(I) .EQ. -3 ) THEN
                  IP(1) = 0
                  IP(2) = 0
               ELSE
C                           * MINDESTENS EINE DER BEIDEN ACHSEN SOLL
C                           * BESCHRIFTET WERDEN
                  IP(1) = IPRECS(1)
                  IP(2) = IPRECS(2)
C                           * BESTIMME DIE ACHSEN, DIE BESCHRIFTET
C                           * WERDEN SOLLEN
                  IACHSE(I) = 3+ILABEL(I)
               ENDIF
            ENDIF
C                           * ORDINATE DES UNTEREN RANDES DES ACHSEN-
C                           * BEREICHES
            PKTLU(2) = SCLV(ICOMPL,I)
C                           * ORDINATE DES OBEREN RANDES DES ACHSEN-
C                           * BEREICHES
            PKTRO(2) = SCLV(ICOMPL,ICOMPL)
C                           * ACHSE(N) BESCHRIFTEN
            IF ( I .EQ. 1 ) THEN
               IF ( ITEXT .GE. 1 ) THEN
                  LTEXT = LTXX
               ELSE
                  LTEXT = 0
               ENDIF
               CALL GRAXLB(LXACHS,PKTLU,PKTRO,FAKCMV,CHRSZC,IACHSE(1),
     $            IP,CHRAND,FRSTGU,DIFFGU,N,LTEXT,TEXTX,INTS,LINNEN)
            ELSE
               HELP = FAKCMV(1)
               FAKCMV(1) = FAKCMV(2)
               FAKCMV(2) = HELP
               IF ( IACHSE(2) .NE. 3 ) THEN
                  IBSCHR = 3-IACHSE(2)
               ELSE
                  IBSCHR = 3
               ENDIF
               CHRAND = CHRAND-90.
               IF ( ITEXT .EQ. 2 ) THEN
                  LTEXT = LTXXY
               ELSE
                  LTEXT = 0
               ENDIF
               CALL GRAXLB(LXACHS,PKTLU,PKTRO,FAKCMV,CHRSZC,IBSCHR,
     $            IP,CHRAND,FRSTGU,DIFFGU,N,LTEXT,TEXTXY,INTS,LINNEN)
            ENDIF
         ENDIF
C                           * WINKEL EINES ZEICHENS UM PI/2 VERMINDERN,
C                           * DA Y-ACHSE UNTER EINEM WINKEL VON -90 GRAD
C                           * GEZEICHNET WIRD
         RANGLE = RANGLE-CPI/2.
         LXACHS = .FALSE.
    2 CONTINUE

      CALL GRDSH(SVDSH1,SVDSH2,SVDSH3)

C                           * WURDEN MELDUNGEN AUSGEGEBEN?
      IF ( CPRCTL .NE. '1' ) THEN
C                           * ERZEUGE 12 LEERZEILEN
         DO 3, I = 1,6
            WRITE ( *, '( ''0'' )' )
    3    CONTINUE
      ENDIF
      END
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C ----------------------------------------------------------------------
C  GRAXOP - SCANNEN DER Zeichenkette "OPTION" UND BEREITSTELLEN
C           ALLER INFORMATIONEN BZGL. AUSSEHEN DER ACHSE
C           WIRD AUFGERUFEN VON: GRAXS
C  AUTOR: G. EGERER
C  DATUM:  1.12.83
C  AENDERUNGEN:
C     13. 4.84
C     18. 9.89   KEYWORTE "XLABEL" UND "YLABEL" HINZUGEFUEGT
C ----------------------------------------------------------------------
      SUBROUTINE GRAXOP(LOPT,OPTION,
     $   IACHSE,LINNEN,IFORMT,AMINA,ILABEL,ISTFU,ITEXT)

C                           * EINGABEPARAMETER:
C                           * LOPT   - LAENGE DER Zeichenkette
C                           *          "OPTION"
      INTEGER     LOPT
C                           * OPTION - Zeichenkette, DIE ALLE
C                           *          INFORMATIONEN BZGL. AUSSEHEN DER
C                           *          ACHSE IN FORM VON KEYWORTEN
C                           *          ENTHAELT
      CHARACTER(len=*) OPTION

C                           * AUSGABEPARAMETER:
C                           * IACHSE - WERTE DER KEYWORTE "XACHSE" UND
C                           *          "YACHSE":
C                           *          IACHSE(1)=WERT "XACHSE"
C                           *          IACHSE(2)=WERT "YACHSE"
C                           * LINNEN - GIBT AN, OB UNTERTEILUNGSSTRICHE
C                           *          INNEN LIEGEN SOLLEN (LINNEN =
C                           *          $TRUE.) ODER AUSSEN (LINNEN =
C                           *          $FALSE.)
C                           *          (KEYWORTE: "INNEN", "AUSSEN")
C                           * IFORMT - WERT DES KEYWORTES "FORMAT":
C                           *          IFORMT(1)=MAXIMALE ANZAHL DER
C                           *             ZIFFERN VOR DEM DEZIMALPUNKT
C                           *          IFORMT(2)=MAXIMALE ANZAHL DER
C                           *             ZIFFERN HINTER DEZIMALPUNKT
C                           * AMINA  - WERT DES KEYWORTES "MINABS"
C                           * ILABEL - WERTE DER KEYWORTE "XLABEL" UND
C                           *          "YLABEL":
C                           *          ILABEL(1)=WERT "XLABEL"
C                           *          ILABEL(2)=WERT "YLABEL"
C                           * ISTFU  - WERT DES KEYWORTES "U"
C                           * ITEXT  - GIBT AN, OB TEXT AN DIE ACHSE
C                           *          GESCHRIEBEN WERDEN SOLL (ITEXT =
C                           *          1 (TEXT AN X-ACHSE) ODER ITEXT =
C                           *          2 (TEXT AN X- UND Y-ACHSE)) ODER
C                           *          NICHT (ITEXT = 0)
C                           *          (KEYWORTE: "TEXT", "TEXTXY")
      INTEGER IACHSE(2), IFORMT(2), ILABEL(2), ISTFU, ITEXT
      LOGICAL LINNEN
      REAL    AMINA

C                           * COMMON-VARIABLEN:
      COMMON      /GRAXCB/ CPRCTL,GENAU,RANDB
CDEC$ PSECT /GRAXCB/ NOSHR
      SAVE /GRAXCB/
      CHARACTER(len=1) CPRCTL,GENAU,RANDB

C                           * BENAMTE PROGRAMMKONSTANTEN:
      INTEGER   CKWANZ
      PARAMETER (CKWANZ = 13)
      INTEGER   CNIL
      REAL      CRNIL
      PARAMETER (CNIL = -32766, CRNIL = -1.)

C                           * PROGRAMMVARIABLEN:
      CHARACTER (len=120) COPT
      CHARACTER (len=120) CSYMBL
      CHARACTER (len=122) CHELP
      CHARACTER (len=17)  CFMT
      CHARACTER (len=6)   CKWLTB(CKWANZ)
      CHARACTER (len=3)   CKWSTB(CKWANZ)
      CHARACTER (len=1)   CURCHR
      LOGICAL       LDLM
      INTEGER I, IAPOS, ICNT, ICUR, IEPOS, IFMP, IFMQ, IHELP, ILEN,
     $        INEXT, IOPTLE, ISCNTB(25,0:15), ISTATE, ISYMBL,
     $        ISYMTB(CKWANZ), IVLPOS, J, KLAPOS, KLLEVL, KLZPOS, KWIND
      REAL    HELP

C                           * INITIALISIERUNGEN:
C TABELLE ZUR SYNTAXDEFINITION:
C
C SYMBOL S0  "," "=" "(" ")" "-" "+" "." " " "E" ZIFF S1 S2  S3  S4  S5
C -----+----------------------------------------------------------------
C ZU-  |
C STAND|
C    1 |         -4                       1
C    2 |         -5                       2
C    3 |         -6                       3
C    4 |                     18  18       4     -23
C    5 |                      7   7   8   5      24
C    6 |             -9      18  18       6     -23
C    7 |                              8          24
C    8 |                                         25
C    9 |     13              10  10       9     -11
C   10 |                                        -11
C   11 |     13         -21              12      11
C   12 |     13         -21              12
C   13 |                     14  14      13     -15
C   14 |                                        -15
C   15 |                -21              16      15
C   16 |                -21              16
C   17 |                     18  18             -23
C   18 |                                        -23
C   19 |                                 19          21  22   1   2   3
C END- |
C ZUST.|
C   20 |                                 20          21  22   1   2   3
C   21 |    -19                          21
C   22 |    -19  -4                      22
C   23 |    -19                          21      23
C   24 |    -19                      25  21  17  24
C   25 |    -19                          21  17  25
C
C S0 (SYMBOLKLASSE 0): ENTHAELT ALLE SYMBOLE, DIE NICHT GESONDERT AUFGE-
C                      FUEHRT UND IN KEINER ANDEREN KLASSE ENTHALTEN
C                      SIND
C S1 (SYMBOLKLASSE 1): ENTHAELT DIE SYMBOLE (KEYWORTE) "I", "INNEN",
C                      "A", "AUSSEN", "T", "TEXT", "TXY" UND "TEXTXY"
C S2 (SYMBOLKLASSE 2): ENTHAELT DIE SYMBOLE (KEYWORTE) "X", "XACHSE",
C                      "Y" UND "YACHSE"
C S3 (SYMBOLKLASSE 3): ENTHAELT DIE SYMBOLE (KEYWORTE) "XLB", "XLABEL",
C                      "YLB", "YLABEL" UND "U"
C S4 (SYMBOLKLASSE 4): ENTHAELT DIE SYMBOLE (KEYWORTE) "M" UND "MINABS"
C S5 (SYMBOLKLASSE 5): ENTHAELT DIE SYMBOLE (KEYWORTE) "F" UND "FORMAT"
C ZIFF               : ENTHAELT ALLE ZIFFERN VON "0" - "9"
C
      DATA ((ISCNTB(I,J), J = 0,15), I = 1,19)
     $   / 2*CNIL, -4, 5*CNIL, 1, 7*CNIL,
     $   2*CNIL, -5, 5*CNIL, 2, 7*CNIL,
     $   2*CNIL, -6, 5*CNIL, 3, 7*CNIL,
     $   5*CNIL, 18, 18, CNIL, 4, CNIL, -23, 5*CNIL,
     $   5*CNIL, 7, 7, 8, 5, CNIL, 24, 5*CNIL,
     $   3*CNIL, -9, CNIL, 18, 18, CNIL, 6, CNIL, -23, 5*CNIL,
     $   7*CNIL, 8, 2*CNIL, 24, 5*CNIL,
     $   10*CNIL, 25, 5*CNIL,
     $   CNIL, 13, 3*CNIL, 10, 10, CNIL, 9, CNIL, -11, 5*CNIL,
     $   10*CNIL, -11, 5*CNIL,
     $   CNIL, 13, 2*CNIL, -21, 3*CNIL, 12, CNIL, 11, 5*CNIL,
     $   CNIL, 13, 2*CNIL, -21, 3*CNIL, 12, 7*CNIL,
     $   5*CNIL, 14, 14, CNIL, 13, CNIL, -15, 5*CNIL,
     $   10*CNIL, -15, 5*CNIL,
     $   4*CNIL, -21, 3*CNIL, 16, CNIL, 15, 5*CNIL,
     $   4*CNIL, -21, 3*CNIL, 16, 7*CNIL,
     $   5*CNIL, 18, 18, 3*CNIL, -23, 5*CNIL,
     $   10*CNIL, -23, 5*CNIL,
     $   8*CNIL, 19, 2*CNIL, 21, 22, 1, 2, 3 /
      DATA ((ISCNTB(I,J), J = 0,15), I = 20,25)
     $   / 8*CNIL, 20, 2*CNIL, 21, 22, 1, 2, 3,
     $   CNIL, -19, 6*CNIL, 21, 7*CNIL,
     $   CNIL, -19, -4, 5*CNIL, 22, 7*CNIL,
     $   CNIL, -19, 6*CNIL, 21, CNIL, 23, 5*CNIL,
     $   CNIL, -19, 5*CNIL, 25, 21, 17, 24, 5*CNIL,
     $   CNIL, -19, 6*CNIL, 21, 17, 25, 5*CNIL /
C
      DATA ISYMTB /12, 12, 11, 11, 15, 14, 13, 13, 13, 11, 11, 15, 11/
      DATA CKWLTB /'XACHSE',
     $             'YACHSE',
     $             'INNEN ',
     $             'AUSSEN',
     $             'FORMAT',
     $             'MINABS',
     $             'XLABEL',
     $             'YLABEL',
     $             '      ',
     $             'TEXT  ',
     $             'TEXTXY',
     $             'GENAU ',
     $             'RAND  '/
      DATA CKWSTB /'X  ',
     $             'Y  ',
     $             'I  ',
     $             'A  ',
     $             'F  ',
     $             'M  ',
     $             'XLB',
     $             'YLB',
     $             'U  ',
     $             'T  ',
     $             'TXY',
     $             'G  ',
     $             'R  '/
C                           * ------------------------------------------

C                           * VORBESETZUNG DER AUSGABEPARAMETER MIT
C                           * DEFAULTWERTEN:
       IACHSE(1) = 0
       IACHSE(2) = 0
       LINNEN = .FALSE.
       IFORMT(1) = 3
       IFORMT(2) = 2
       AMINA = CRNIL
       ILABEL(1) = 0
       ILABEL(2) = 0
       ISTFU = 2
       ITEXT = 2
       GENAU = '0'
       RANDB = '0'

      IF ( LOPT .GT. 120 ) THEN
         WRITE ( *, '( '''//CPRCTL//''',A,/,1X,A )' )
     $      '**GRAXS: Zeichenkette "OPTION" ZU LANG. '//
     $         'SIE WIRD RECHTS ABGESCHNITTEN',
     $      '         AUF DIE LAENGE VON 120 ZEICHEN.'
         CPRCTL = '0'
         IOPTLE = 120
      ELSE
         IOPTLE = LOPT
      ENDIF

      IF ( ioptle.GT.0) THEN
         copt=option(:ioptle)
      ELSE
         copt=' '
      ENDIF

C                           * KLAMMERLEVEL
      KLLEVL = 0
C                           * LDLM = .FALSE. ZEIGT AN, DAS "," NICHT ALS
C                           * DELIMITER ZWISCHEN KEYWORTEN GEBRAUCHT
C                           * WIRD UND IM FEHLERFALLE IGNORIERT WERDEN
C                           * SOLL (GILT Z.B. FUER "," IN "X=,I")
      LDLM = .FALSE.
      IAPOS = 1
C                           * KEYWORT-INDEX
      KWIND = CNIL
      KLAPOS = CNIL
      IVLPOS = CNIL
C                           * ANFANGSZUSTAND
      ISTATE = 20

      ICUR = 1
  100    CONTINUE
         IF ( ICUR .GT. IOPTLE ) THEN
            IF ( ISTATE .GE. 20 ) THEN
C                           * ISTATE IST EIN ENDZUSTAND
               ISTATE = -19
            ELSE
C                           * ISTATE IST KEIN ENDZUSTAND
               ISTATE = CNIL
            ENDIF
         ELSE
C                           * ERMITTLE INDEX DES AKTUELLEN SYMBOLS
            CURCHR = COPT(ICUR:ICUR)
            IF (  ( CURCHR .GE. 'A' ) .AND. ( CURCHR .LE. 'Z' )  ) THEN
               DO 118, I = ICUR+1, IOPTLE
                  IF (  ( COPT(I:I) .LT. 'A' ) .OR.
     $               ( COPT(I:I) .GT. 'Z' )  ) GO TO 119
  118          CONTINUE
  119          IF ( I-1 .EQ. ICUR ) THEN
C                           * EINZELNER BUCHSTABE (KEIN KEYWORT)?
                  ISYMBL = INDEX('E',CURCHR)
                  IF ( ISYMBL .NE. 0 ) ISYMBL = ISYMBL+8
               ELSE
                  ISYMBL = 0
               ENDIF
               IF ( ISYMBL .EQ. 0 ) THEN
C                           * KEYWORT?
                  CSYMBL = COPT(ICUR:I-1)
                  DO 128, KWIND = 1,CKWANZ
                     IF (  ( CSYMBL .EQ. CKWSTB(KWIND) ) .OR.
     $                  ( CSYMBL .EQ. CKWLTB(KWIND) )  ) THEN
                           ISYMBL = ISYMTB(KWIND)
                           GO TO 129
                     ENDIF
  128             CONTINUE
  129             CONTINUE
                  IF (KWIND.EQ.12) GENAU = '1'
                  IF (KWIND.EQ.13) RANDB = '1'
               ENDIF
               INEXT = I
            ELSE IF (  ( CURCHR .GE. '0' ) .AND. ( CURCHR .LE. '9' )  )
     $         THEN
                  IF (  ( ISTATE .EQ. 11 ) .OR. ( ISTATE .EQ. 15 ) .OR.
     $               ( ISTATE .EQ. 23 )  ) THEN
                        IF ( ICNT .LT. 8 ) THEN
                           ISYMBL = 10
                           ICNT = ICNT+1
                        ELSE
                           ISYMBL = 0
                        ENDIF
                  ELSE
                     ISYMBL = 10
                  ENDIF
                  INEXT = ICUR+1
            ELSE
               ISYMBL = INDEX(',=()-+. ',CURCHR)
               INEXT = ICUR+1
            ENDIF

            ISTATE = ISCNTB(ISTATE,ISYMBL)
         ENDIF

         IF ( ISTATE .LT. 0 ) THEN
            IF ( ISTATE .EQ. -19 ) THEN
C                           * AKTUELLES SYMBOL = "," AUF KLAMMERLEVEL 0
               IF (  ( KWIND .EQ. 1 ) .OR. ( KWIND .EQ. 2 )  ) THEN
C                           * LETZTES KEYWORT WAR "XACHSE" ODER "YACHSE"
                  I = KWIND
                  IF ( IVLPOS .NE. CNIL ) THEN
C                           * FUER "XACHSE" BZW. "YACHSE" WURDE EIN WERT
C                           * ANGEGEBEN
                     WRITE ( CFMT, '( ''(BN,I'',I3,'')'' )' )
     $                  ICUR-IVLPOS
                     READ ( COPT(IVLPOS:), CFMT ) IHELP
                     IF (  ( ABS(IHELP) .GT. 3 ) .OR.
     $                  ( IHELP .EQ. 0 )  ) THEN
                           WRITE ( *, '( '''//CPRCTL//''',A,/,1X,A )' )
     $                        '**GRAXS: UNZULAESSIGER WERT FUER "'//
     $                           CKWLTB(KWIND)//'" ANGEGEBEN.',
     $                        '         "'//CKWLTB(KWIND)//
     $                           '" WIRD IGNORIERT.'
                           CPRCTL = '0'
                     ELSE
                        IACHSE(I) = IHELP
                     ENDIF
                  ELSE
                     IACHSE(I) = 1
                  ENDIF
               ELSE IF ( KWIND .EQ. 3 ) THEN
C                           * LETZTES KEYWORT WAR "INNEN"
                  LINNEN = .TRUE.
               ELSE IF ( KWIND .EQ. 4 ) THEN
C                           * LETZTES KEYWORT WAR "AUSSEN"
                  LINNEN = .FALSE.
               ELSE IF ( KWIND .EQ. 5 .OR. KWIND .EQ. 12 ) THEN
C                           * LETZTES KEYWORT WAR "FORMAT"
                  I = INDEX(COPT(IVLPOS:ICUR-1),',')
                  IF ( I .NE. 0 ) THEN
                     WRITE ( CFMT,
     $                  '( ''(BN,I'',I3,'',1X,I'',I3,'')'' )' )
     $                     (IVLPOS+I-2)-KLAPOS, KLZPOS-(IVLPOS+I)
                     READ ( COPT(KLAPOS+1:), CFMT ) IFMP, IFMQ
                     IF (  ( IFMP .LT. 0 ) .OR. ( IFMQ .LT. 0 ) .OR.
     $                  ( IFMP+IFMQ .GT. 6 )  ) THEN
                           WRITE ( *, '( '''//CPRCTL//''',A )' )
     $                        '**GRAXS: UNZULAESSIGES FORMAT '//
     $                           'ANGEGEBEN. FORMAT WIRD IGNORIERT.'
                           CPRCTL = '0'
                     ELSE
                        IFORMT(1) = IFMP
                        IFORMT(2) = IFMQ
                     ENDIF
                  ELSE
                     IF ( KLAPOS .NE. CNIL ) THEN
                        IVLPOS = KLAPOS+1
                        IHELP = KLZPOS
                     ELSE
                        IHELP = ICUR
                     ENDIF
                     WRITE ( CFMT, '( ''(BN,I'',I3,'')'' )' )
     $                  IHELP-IVLPOS
                     READ ( COPT(IVLPOS:), CFMT ) IFMP
                     IF (  ( IFMP .LT. 0 ) .OR. ( IFMP .GT. 6 )  ) THEN
                        WRITE ( *, '( '''//CPRCTL//''',A )' )
     $                     '**GRAXS: UNZULAESSIGES FORMAT ANGEGEBEN. '//
     $                        'FORMAT WIRD IGNORIERT.'
                        CPRCTL = '0'
                     ELSE
                        IFORMT(1) = IFMP
                        IFORMT(2) = 0
                     ENDIF
                  ENDIF
               ELSE IF ( KWIND .EQ. 6 ) THEN
C                           * LETZTES KEYWORT WAR "MINABS"
                  I = INDEX(COPT(IVLPOS:ICUR-1),'E')
                  IF ( I .NE. 0 ) THEN
                     WRITE ( CFMT, '( ''(BN,F'',I3,''.0)'' )' ) I-1
                     READ ( COPT(IVLPOS:), CFMT ) HELP
                     IF ( HELP .NE. 0. ) THEN
                        WRITE ( CFMT, '( ''(BN,I'',I3,'')'' )' )
     $                     ICUR-(IVLPOS+I)
                        READ ( COPT(IVLPOS+I:), CFMT ) IHELP
                        IHELP = IHELP+INT(ALOG10(ABS(HELP)))
                        IF ( IABS(IHELP) .GT. 74 ) THEN
C                           * WEISE AUF HELP UNZULAESSIGEN WERT ZU
                           HELP = -1.
                        ELSE
                           WRITE ( CFMT, '( ''(BN,F'',I3,''.0)'' )' )
     $                        ICUR-IVLPOS
                           READ ( COPT(IVLPOS:), CFMT ) HELP
                        ENDIF
                     ELSE
C                           * HELP = 0
                     ENDIF
                  ELSE
                     WRITE ( CFMT, '( ''(BN,F'',I3,''.0)'' )' )
     $                  ICUR-IVLPOS
                     READ ( COPT(IVLPOS:), CFMT ) HELP
                  ENDIF
                  IF ( HELP .LT. 0. ) THEN
                     WRITE ( *, '( '''//CPRCTL//''',A )' )
     $                  '**GRAXS: UNZULAESSIGER WERT FUER "MINABS" '//
     $                     'WIRD IGNORIERT.'
                     CPRCTL = '0'
                  ELSE
                     AMINA = HELP
                  ENDIF
               ELSE IF (  ( KWIND .EQ. 7 ) .OR. ( KWIND .EQ. 8 )  ) THEN
C                           * LETZTES KEYWORT WAR "XLABEL" ODER "YLABEL"
                  WRITE ( CFMT, '( ''(BN,I'',I3,'')'' )' )
     $               ICUR-IVLPOS
                  READ ( COPT(IVLPOS:), CFMT ) IHELP
                  IF (  ( IHELP .LT. -3 ) .OR. ( IHELP .GT. -1 )  ) THEN
                     WRITE ( *, '( '''//CPRCTL//''',A )' )
     $                  '**GRAXS: UNZULAESSIGER WERT FUER "'//
     $                     CKWLTB(KWIND)//'" WIRD IGNORIERT.'
                     CPRCTL = '0'
                  ELSE
                     ILABEL(KWIND-6) = IHELP
                  ENDIF
               ELSE IF ( KWIND .EQ. 9 ) THEN
C                           * LETZTES KEYWORT WAR "U"
                  WRITE ( CFMT, '( ''(BN,I'',I3,'')'' )' )
     $               ICUR-IVLPOS
                  READ ( COPT(IVLPOS:), CFMT ) IHELP
                  IF (  ( IHELP .LT. 0 ) .OR. ( IHELP .GT. 3 )  ) THEN
                     WRITE ( *, '( '''//CPRCTL//''',A )' )
     $                  '**GRAXS: UNZULAESSIGER WERT FUER "U" WIRD '//
     $                     'IGNORIERT.'
                     CPRCTL = '0'
                  ELSE
                     ISTFU = IHELP
                  ENDIF
               ELSE IF ( KWIND .EQ. 10 ) THEN
C                           * LETZTES KEYWORT WAR "TEXT"
C                 IF ( ITEXT .EQ. 0 ) ITEXT = 1
               ELSE IF ( KWIND .EQ. 11 ) THEN
C                           * LETZTES KEYWORT WAR "TEXTXY"
C                 ITEXT = 2
               ELSE
C                           * LEERE EINGABE (KWIND = CNIL)
               ENDIF
               LDLM = .TRUE.
               IAPOS = ICUR
               KLAPOS = CNIL
               IVLPOS = CNIL
               ISTATE = -ISTATE
            ELSE IF (  ( ISTATE .GE. -6 ) .AND. ( ISTATE .LE. -4 )  )
     $         THEN
C                           * AKTUELLES SYMBOL = "="
                  IVLPOS = ICUR+1
                  ISTATE = -ISTATE
            ELSE IF (  ( ISTATE .EQ. -11 ) .OR. ( ISTATE .EQ. -15 ) .OR.
     $         ( ISTATE .EQ. -23 )  ) THEN
C                           * AKTUELLES SYMBOL = 1. ZIFFER EINER
C                           * INTEGERZAHL
                  ICNT = 0
                  ISTATE = -ISTATE
            ELSE IF ( ISTATE .EQ. -9 ) THEN
C                           * AKTUELLES SYMBOL = "("
               KLLEVL = KLLEVL+1
               KLAPOS = ICUR
               ISTATE = -ISTATE
            ELSE IF ( ISTATE .EQ. -21 ) THEN
C                           * AKTUELLES SYMBOL = ")"
               KLLEVL = KLLEVL-1
               KLZPOS = ICUR
               ISTATE = -ISTATE
            ELSE
C                           * ISTATE = CNIL

               I = INEXT-1
  130             CONTINUE
                  IF ( I .GT. IOPTLE ) THEN
                     IEPOS = IOPTLE
                     ICUR = I
                     GO TO 139
                  ENDIF
                  IF ( COPT(I:I) .EQ. ',' ) THEN
                     IF ( KLLEVL .EQ. 0 ) THEN
                        IF ( LDLM ) THEN
                           IEPOS = I-1
                           ISTATE = 19
                        ELSE
                           IEPOS = I
                           ISTATE = 20
                           KWIND = CNIL
                        ENDIF
                        KLAPOS = CNIL
                        IVLPOS = CNIL
                        INEXT = I+1
                        GO TO 139
                     ENDIF
                  ELSE IF ( COPT(I:I) .EQ. '(' ) THEN
                     KLLEVL = KLLEVL+1
                  ELSE IF ( COPT(I:I) .EQ. ')' ) THEN
                     KLLEVL = KLLEVL-1
                  ELSE
                  ENDIF
                  I = I+1
                  GO TO 130
  139          CONTINUE

               WRITE ( *, '( '''//CPRCTL//''',A,/,1X,A )' )
     $            '**GRAXS: Zeichenkette "OPTION" ENTHAELT '//
     $               'UNZULAESSIGE SYNTAX ODER',
     $            '         UNGUELTIGES KEYWORT:'
               CHELP = '"'//COPT(IAPOS:IEPOS)//'"'
               ILEN = IEPOS-IAPOS+3
               IF ( ILEN .GT. 63 ) THEN
                  WRITE ( *, '( 1X,2A )' ) '         ', CHELP(:63)
                  ILEN = ILEN-63
                  I = 64
               ELSE
                  I = 1
               ENDIF
               IF ( ILEN+16 .LE. 63 ) THEN
                  WRITE ( *, '( 1X,3A )' ) '         ',
     $               CHELP(I:I+ILEN-1), ' WIRD IGNORIERT.'
               ELSE
                  WRITE ( *, '( 1X,2A )' ) '         ',
     $               CHELP(I:I+ILEN-1)
                  WRITE ( *, '( 1X,A )' ) '         WIRD IGNORIERT.'
               ENDIF
               CPRCTL = '0'

               IAPOS = IEPOS+1
            ENDIF
         ENDIF

         IF ( ICUR .GT. IOPTLE ) GO TO 199
         ICUR = INEXT
         GO TO 100
  199 CONTINUE

      END
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C ----------------------------------------------------------------------
C  GRAXGU - ERMITTELT DIE GROBUNTERTEILUNG FUER EINE ACHSE
C           WIRD AUFGERUFEN VON: GRAXS
C  AUTOR: G. EGERER, DATUM: 10. 8.83, LETZTE AENDERUNG: 11.12.85
C ----------------------------------------------------------------------
      SUBROUTINE GRAXGU(XAVAL,XBVAL,AXLEN,CHRHEI,CHRANG,AMINA,IFORMT,
     $   MANKBZ,IFAK,N,IPRECS)

C                           * EINGABEPARAMETER:
C                           * XAVAL  - ANFANGSWERT DES WERTEBEREICHES
C                           *          DER ACHSE
C                           * XBVAL  - ENDWERT DES WERTEBEREICHES
C                           *          DER ACHSE
C                           * AXLEN  - LAENGE DER ACHSE IN CM
C                           * CHRHEI - GROESSE EINES ZEICHENS DER
C                           *          ACHSENBESCHRIFTUNG IN CM
C                           * CHRANG - WINKEL EINES ZEICHENS IN
C                           *          BOGENMASS
C                           * AMINA  - MINIMALABSTAND ZWISCHEN 2 ZAHLEN
C                           *          DER BESCHRIFTUNG IN CM
C                           * IFORMT - FORMATANGABE DES BENUTZERS FUER
C                           *          DIE BESCHRIFTUNG
C                           *          ( FORMAT=(IFORMT(1),IFORMT(2)) )
      REAL    XAVAL, XBVAL, AXLEN, CHRHEI, CHRANG, AMINA
      INTEGER IFORMT(2)

C                           * AUSGABEPARAMETER:
C                           * MANKBZ,
C                           * IFAK,
C                           * N      - KLEINSTE BESCHRIFTUNGSZAHL =
C                           *          MANKBZ*10**N,
C                           *          ERMITTELTE GROBUNTERTEILUNG =
C                           *          IFAK*10**N
C                           * IPRECS - MAXIMALE GENAUIGKEIT DER
C                           *          BESCHRIFTUNGSZAHLEN:
C                           *          IPRECS(1)=ANZAHL DER STELLEN VOR
C                           *             DEM DEZIMALPUNKT
C                           *          IPRECS(2)=ANZAHL DER STELLEN
C                           *             HINTER DEM DEZIMALPUNKT
C                           *          BEI BESCHRIFTUNG IM E-FORMAT:
C                           *          IPRECS(1)=ANZAHL DER SIGNIFIKAN-
C                           *             TEN STELLEN
C                           *          IPRECS(2)=-32766
      INTEGER MANKBZ, IFAK, N, IPRECS(2)
C        COMMON - Block
      COMMON      /GRAXCB/ CPRCTL,GENAU,RANDB
CDEC$ PSECT /GRAXCB/ NOSHR
      SAVE /GRAXCB/
      CHARACTER(len=1) CPRCTL,GENAU,RANDB
C                           * PROGRAMMVARIABLEN:
      CHARACTER (len=12) CHELP
      LOGICAL          LCEIL(2)
      INTEGER I, IEXP, IFAKA, IFMP, IFMQ, IHELP, ILOGBZ(2), IND, INDFAK,
     $        IPRENA, ISIGNW(2), IWIDTH(2), J, MANTIS, MAXLOG, MAXWID,
     $        NA
      REAL    DIFFCM, DIFFV, E0NA, FAK(3), HELP, RANGEV, WIDCM, X(2),
     $        XA, XB

C                           * INITIALISIERUNGEN:
      DATA FAK /1., 2., 5./
C                           * ------------------------------------------

      IF ( XAVAL .GT. XBVAL ) THEN
         XA = XBVAL
         XB = XAVAL
      ELSE
         XA = XAVAL
         XB = XBVAL
      ENDIF

      IF (  ( IFORMT(1) .EQ. 0 ) .AND. ( IFORMT(2) .EQ. 0 )  ) THEN
         IFMP = 3
         IFMQ = 2
      ELSE
         IFMP = IFORMT(1)
         IFMQ = IFORMT(2)
      ENDIF

      RANGEV = XB-XA
C                           * ERMITTLUNG EINER MINIMALEN DIFFERENZ
C                           * DER BESCHRIFTUNGSZAHLEN (ANZAHL DER GROB-
C                           * UNTERTEILUNGEN < 150 ANGENOMMEN):
      DIFFV = RANGEV/150.
C                           * NA = FLOOR(ALOG10(DIFFV))
C                           *    = (ANZAHL DER STELLEN VON DIFFV)-1
      HELP = ALOG10(DIFFV)
      IF ( HELP .GE. 0. ) THEN
         NA = INT(HELP)
      ELSE
         NA = INT(HELP+79.)-79
      ENDIF
      E0NA = 10.**NA
      HELP = DIFFV/E0NA
      IF ( HELP .GE. 5. ) THEN
         INDFAK = 1
         E0NA = E0NA*10.
         NA = NA+1
      ELSE IF ( HELP .GE. 2. ) THEN
         INDFAK = 3
      ELSE
         INDFAK = 2
      ENDIF

      DO 1888, I = 1,4
         DO 1188, IND = INDFAK,3
            DIFFV = FAK(IND)*E0NA
            IF ( DIFFV .GT. RANGEV ) THEN
               IFAK = IFAKA
               N = IPRENA
               IPRECS(1) = 0
               IPRECS(2) = 0
               GO TO 1999
            ENDIF
            IFAKA = INT(FAK(IND))
            DIFFCM = AXLEN/RANGEV*DIFFV
C                           * ERMITTLUNG DER GROESSTEN UND KLEINSTEN
C                           * BESCHRIFTUNGSZAHL:
C                           * GROESSTE BESCHRIFTUNGSZAHL =
C                           *    FLOOR(XB/DIFFV)*DIFFV
C                           * KLEINSTE BESCHRIFTUNGSZAHL =
C                           *    CEIL(XA/DIFFV)*DIFFV
C                           * (BEIDE ZAHLEN WERDEN IN DER FORM
C                           *    MANTIS*10**NA
C                           * ERMITTELT;
C                           * ILOGBZ(J) = ALOG10(BESCHRIFTUNGSZAHL)
C                           *           = (ANZAHL DER STELLEN DER
C                           *                BESCHRIFTUNGSZAHL)-1   )
            X(1) = XB
            X(2) = XA
            LCEIL(1) = .FALSE.
            LCEIL(2) = .TRUE.
            DO 1118, J = 1,2
               IF ( X(J) .EQ. 0. ) THEN
                  MANTIS = 0
                  IWIDTH(J) = -32766
               ELSE
                  IF ( X(J) .LT. 0. ) THEN
                     ISIGNW(J) = 1
                     X(J) = -X(J)
                     LCEIL(J) = .NOT. LCEIL(J)
                  ELSE
                     ISIGNW(J) = 0
                  ENDIF
                  WRITE ( CHELP, '( 6P,E12.5E3 )' ) X(J)
                  READ ( CHELP(:6), '( I6 )' ) MANTIS
                  MANTIS = MANTIS*10/IFAKA
                  READ ( CHELP(9:), '( I4 )' ) IEXP
                  IEXP = IEXP-1-NA
C**                         * MANTIS*10**IEXP = (XB BZW. XA)/DIFFV
                  IF ( IEXP .LE. -7 ) THEN
C                           * ABS( (XA BZW. XB)/DIFFV ) < 1
                     IF ( LCEIL(J) ) THEN
                        MANTIS = IFAKA
                        ILOGBZ(J) = NA
                        IWIDTH(J) = ISIGNW(J)+NA
                     ELSE
                        MANTIS = 0
                        IWIDTH(J) = -32766
                     ENDIF
                  ELSE IF ( IEXP .LT. 0 ) THEN
                     IHELP = MANTIS/10**(-IEXP)
                     IF ( LCEIL(J) ) THEN
                        IF ( MANTIS-IHELP*10**(-IEXP) .NE. 0 ) THEN
C                           * ((XB BZW. XA)/DIFFV) IST KEINE GANZE ZAHL
                           MANTIS = (IHELP+1)*IFAKA
                        ELSE
C                           * ((XB BZW. XA)/DIFFV) IST EINE GANZE ZAHL
                           MANTIS = IHELP*IFAKA
                        ENDIF
                        ILOGBZ(J) = IDINT(DLOG10(DBLE(MANTIS))+1.D-8)+
     $                     NA
                        IWIDTH(J) = ISIGNW(J)+ILOGBZ(J)
                     ELSE
                        IF ( IHELP .EQ. 0 ) THEN
                           MANTIS = 0
                           IWIDTH(J) = -32766
                        ELSE
                           MANTIS = IHELP*IFAKA
                           ILOGBZ(J) = IDINT(DLOG10(DBLE(MANTIS))+
     $                        1.D-8)+NA
                           IWIDTH(J) = ISIGNW(J)+ILOGBZ(J)
                        ENDIF
                     ENDIF
                  ELSE
                     MANTIS = MANTIS*IFAKA
                     ILOGBZ(J) = IDINT(DLOG10(DBLE(MANTIS))+1.D-8)+
     $                  IEXP+NA
                     IWIDTH(J) = ISIGNW(J)+ILOGBZ(J)
                     MANTIS = MANTIS*10**IEXP
                  ENDIF
               ENDIF
 1118       CONTINUE
            IF ( ISIGNW(2) .EQ. 0 ) THEN
               MANKBZ = MANTIS
            ELSE
               MANKBZ = -MANTIS
            ENDIF

C                           * ERMITTLUNG DER GROESSTEN LAENGE EINER
C                           * BESCHRIFTUNGSZAHL (HABEN KLEINSTE UND
C                           * GROESSTE BESCHRIFTUNGSZAHL DIESELBE LAENGE
C                           * MUESSEN DIE GROESSEN MAXLOG UND MAXWID
C                           * BEZOGEN AUF DIE GROESSTE BESCHRIFTUNGSZAHL
C                           * GESETZT WERDEN, DENN DIE GROESSTE BE-
C                           * SCHRIFTUNGSZAHL KANN MEHR STELLEN HABEN,
C                           * WENN SIE POSITIV IST UND DIE KLEINSTE
C                           * NEGATIV (Z.B. -4 ... 10).):
            IF ( IWIDTH(1) .GE. IWIDTH(2) ) THEN
               IF ( IWIDTH(1) .EQ. -32766 ) THEN
C                           * KLEINSTE BESCHRIFTUNGSZAHL = GROESSTE BE-
C                           * SCHRIFTUNGSZAHL = 0
                  MAXLOG = -32766
               ELSE
                  MAXLOG = ILOGBZ(1)+1
                  MAXWID = ISIGNW(1)
               ENDIF
            ELSE
               MAXLOG = ILOGBZ(2)+1
               MAXWID = ISIGNW(2)
            ENDIF
C                           * LAENGSTE BESCHRIFTUNGSZAHL UNGLEICH 0?
            IF ( MAXLOG .NE. -32766 ) THEN
C                           * PASST LAENGSTE BESCHRIFTUNGSZAHL IN DAS
C                           * VORGEGEBENE FORMAT?
               IF (  ( MAXLOG .LE. IFMP ) .AND. ( -NA .LE. IFMQ )  )
     $            THEN
                     IF ( NA .GE. 0 ) THEN
C                           * KEINE STELLEN HINTER DEZIMALPUNKT
                        IPRECS(1) = MAXLOG
                        IPRECS(2) = 0
C                           * LAENGE DER ZAHL = ANZAHL DER STELLEN VOR
C                           * DEZIMALPUNKT (EVTL. +1 FUER VORZEICHEN):
                        MAXWID = MAXWID+MAXLOG
                     ELSE
C                           * NA < 0; -NA STELLEN HINTER DEZIMALPUNKT
                        IF (GENAU.EQ.'0') THEN
                           IPRECS(1) = MAX(MAXLOG,0)
                           IPRECS(2) = -NA
                        ELSE
                           IPRECS(1) = IFORMT(1)
                           IPRECS(2) = IFORMT(2)
                        ENDIF
C                           * LAENGE DER ZAHL = ANZAHL ALLER STELLEN
C                           * +1 FUER DEZIMALPUNKT (EVTL. +1 FUER
C                           * VORZEICHEN):
                        MAXWID = MAXWID+IPRECS(1)-NA+1
                     ENDIF
                  ELSE
C                           * LAENGSTE BESCHRIFTUNGSZAHL PASST NICHT
C                           * IN DAS VORGEGEBENE FORMAT
C                           * ZAHL IM E-FORMAT MIT ALLEN STELLEN
C                           * DARSTELLBAR?
                     IHELP = MAXLOG-NA
                     IF ( IHELP .LE. IFMP+IFMQ ) THEN
                        IPRECS(1) = IHELP
                        IPRECS(2) = -32766
                        MAXWID = IHELP
C                           * ANWEISUNG IN VORAUSGEHENDER ZEILE FALSCH
C                           * UND DURCH FOLGENDE ZU ERSETZEN?
C                       MAXWID = MAXWID+IHELP
                     ELSE
                        MAXWID = -32766
                     ENDIF
                  ENDIF
            ELSE
C                           * LAENGSTE BESCHRIFTUNGSZAHL = 0
               IPRECS(1) = 1
               IPRECS(2) = 0
               MAXWID = 1
            ENDIF
            IF ( MAXWID .NE. -32766 ) THEN
               WIDCM = ( ABS(COS(CHRANG))*REAL(MAXWID)+
     $            ABS(SIN(CHRANG)) )*CHRHEI+AMINA
               IF ( WIDCM .LE. DIFFCM ) THEN
                  IFAK = IFAKA
                  N = NA
                  GO TO 1999
               ENDIF
            ENDIF
            IPRENA = NA
 1188    CONTINUE
         INDFAK = 1
         E0NA = E0NA*10.
         NA = NA+1
 1888 CONTINUE
 1999 CONTINUE

      IF (  ( IFORMT(1) .EQ. 0 ) .AND. ( IFORMT(2) .EQ. 0 )  ) THEN
         IPRECS(1) = 0
         IPRECS(2) = 0
      ENDIF

      END
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C ----------------------------------------------------------------------
C  GRAXDR - ZEICHNET EINE ACHSE
C           WIRD AUFGERUFEN VON: GRAXS
C  AUTOR: G. EGERER, DATUM: 30. 8.83, LETZTE AENDERUNG: 16. 9.83
C                   Groten:                  Aenderung  12.07.91
C ----------------------------------------------------------------------
      SUBROUTINE GRAXDR(LXACHS,PKTLU,PKTRO,LEDGE,LINNEN,IGROBU,
     $   FRSTGU,DIFFGU,IANZFU,SLENGU,FAKCMV)

C                           * EINGABEPARAMETER:
C                           * LXACHS - GIBT AN, OB EINE X-ACHSE
C                           *          (LXACHS = .TRUE.) ODER EINE
C                           *          Y-ACHSE (LXACHS = .FALSE.)
C                           *          GEZEICHNET WERDEN SOLL
C                           * PKTLU  - LINKE UNTERE ECKE DES ZEICHEN-
C                           *          FELDES FUER DIE ACHSE
C                           *          (XLU = PKTLU(1); YLU = PKTLU(2) )
C                           * PKTRO  - RECHTE OBERE ECKE DES ZEICHEN-
C                           *          FELDES
C                           * LEDGE  - GIBT AN, OB INNENLIEGENDE UNTER-
C                           *          TEILUNGSSTRICHE AUF DEN RANDBE-
C                           *          GRENZUNGEN GEZEICHNET WERDEN
C                           *          SOLLEN (LEDGE(I) = .TRUE., I =
C                           *          1 (LINKER RAND), 2 (RECHTER
C                           *          RAND)) ODER NICHT (LEDGE(I) =
C                           *          $FALSE.)
C                           * LINNEN - GIBT AN, OB UNTERTEILUNGSSTRICHE
C                           *          INNEN LIEGEN SOLLEN (LINNEN =
C                           *          $TRUE.) ODER AUSSEN (LINNEN =
C                           *          $FALSE.)
C                           * IGROBU - GIBT AN, OB GROBUNTERTEILUNGS-
C                           *          STRICHE GEZEICHNET (IGROBU = 1),
C                           *          NICHT GEZEICHNET (IGROBU = 0)
C                           *          ODER NETZLINIEN (IGROBU = 2)
C                           *          GEZEICHNET WERDEN SOLLEN
C                           * FRSTGU - ABSZISSE DES 1. GROBUNTER-
C                           *          TEILUNGSSTRICHES
C                           * DIFFGU - ABSTAND DER GROBUNTERTEILUNGS-
C                           *          STRICHE
C                           * IANZFU - ANZAHL DER FEINUNTERTEILUNGS-
C                           *          STRICHE ZWISCHEN 2 GROBUNTER-
C                           *          TEILUNGSSTRICHEN
C                           * SLENGU - LAENGE EINES GROBUNTERTEILUNGS-
C                           *          STRICHES
C                           * FAKCMV - FAKTOR ZUR UMRECHNUNG VON CM-
C                           *          ANGABEN IN BENUTZERWERTE
      LOGICAL LXACHS, LEDGE(2), LINNEN
      REAL    PKTLU(2), PKTRO(2), FRSTGU, DIFFGU, SLENGU, FAKCMV
      INTEGER IGROBU, IANZFU

C                           * PROGRAMMVARIABLEN:
      INTEGER I, IANZ, IEDGEA, IEDGEB,IGU,NGU
      REAL    DIFFU, EPSLN, SLGU, TOLRNZ, XA, XB, XFU, XGU, XGUA,
     $        YDRWFU, YDRWGU, YJMPFU, YJMPGU, YLU, YJMPHU, YDRWHU
C                           * COMMON-Block
      COMMON  /GRAXDS/ SVDSH1, SVDSH2, SVDSH3
CDEC$ PSECT /GRAXDS/ NOSHR
      SAVE /GRAXDS/
      REAL             SVDSH1, SVDSH2, SVDSH3
C                           * BENAMTE PROGRAMMKONSTANTEN:
      REAL      CEPSCM, CTOLCM
      PARAMETER (CEPSCM = 2.E-4, CTOLCM = 0.04)
C                           * ------------------------------------------

      EPSLN = CEPSCM*FAKCMV
      TOLRNZ = CTOLCM*FAKCMV

C                           * ZEICHNEN DER ACHSE:
      YLU = PKTLU(2)
      IF ( LXACHS ) THEN
         CALL GRJMP(PKTLU(1),YLU)
         CALL GRDRW(PKTRO(1),YLU)
      ELSE
         CALL GRJMP(YLU,PKTLU(1))
         CALL GRDRW(YLU,PKTRO(1))
      ENDIF

C                           * ERMITTLUNG DER ORDINATEN FUER DIE
C                           * UNTERTEILUNGSSTRICHE:
      IF ( YLU .GT. PKTRO(2) ) THEN
         SLGU = -SLENGU
      ELSE
         SLGU = SLENGU
      ENDIF
      IF ( LINNEN ) THEN
         YJMPHU = PKTRO(2)-SLGU
         YDRWHU = YLU+SLGU
         YJMPGU = YDRWHU
         YDRWGU = YLU
         YJMPFU = YLU
         YDRWFU = YLU+SLGU/2.
      ELSE
         YJMPGU = YLU
         YJMPHU = PKTRO(2)
         YDRWGU = YLU-SLGU
         YDRWHU = YLU
         YJMPFU = YLU
         YDRWFU = YLU-SLGU/2.
      ENDIF

C                           * ERMITTLUNG DER RICHTUNG, IN DER DIE
C                           * UNTERTEILUNGSSTRICHE AN DIE ACHSE GE-
C                           * ZEICHNET WERDEN SOLLEN:
      IF ( PKTLU(1) .GT. PKTRO(1) ) THEN
         XA = PKTRO(1)
         XB = PKTLU(1)
         IEDGEA = 2
         IEDGEB = 1
      ELSE
         XA = PKTLU(1)
         XB = PKTRO(1)
         IEDGEA = 1
         IEDGEB = 2
      ENDIF

C                           * ERMITTLUNG DES ABSTANDES DER
C                           * FEINUNTERTEILUNGSSTRICHE:
C                           * FEINUNTERTEILUNG FUER DIE ACHSE
C                           * GEWUENSCHT?
      IF ( IANZFU .GT. 0 ) THEN
         DIFFU = DIFFGU/REAL(IANZFU+1)
      ELSE
C                           * FEINUNTERTEILUNG NICHT GEWUENSCHT
C                           * (ABSTAND DER FEINUNTERTEILUNGSSTRICHE
C                           * WIRD GROESSER ALS DER ABSTAND DER
C                           * GROBUNTERTEILUNGSSTRICHE GEWAEHLT)
         DIFFU = DIFFGU*1.1
      ENDIF

C                           * LIEGT 1. GROBUNTERTEILUNGSSTRICH NICHT AUF
C                           * RANDBEGRENZUNG?
      IF ( FRSTGU .GE. XA+TOLRNZ ) THEN
C                           * ANZAHL DER FEINUNTERTEILUNGSSTRICHE VOR
C                           * 1. GROBUNTERTEILUNGSSTRICH
         IANZ = INT( (FRSTGU-(XA-EPSLN))/DIFFU )
C                           * ABSZISSE DES 1. FEINUNTERTEILUNGSSTRICHES
         XFU = FRSTGU-REAL(IANZ)*DIFFU
C                           * SOLLEN INNENLIEGENDE UNTERTEILUNGSSTRICHE
C                           * AUF RANDBEGRENZUNG NICHT GEZEICHNET WER-
C                           * DEN, SOLLEN UNTERTEILUNGSSTRICHE INNEN
C                           * LIEGEN UND LIEGT 1. FEINUNTERTEILUNGS-
C                           * STRICH AUF RANDBEGRENZUNG?
         IF (  ( .NOT. LEDGE(IEDGEA) ) .AND. LINNEN .AND.
     $      ( XFU .LT. XA+TOLRNZ )  ) THEN
               IANZ = IANZ-1
               XFU = XFU+DIFFU
            ENDIF
         XGUA = FRSTGU
      ELSE
C                           * 1. GROBUNTERTEILUNGSSTRICH LIEGT AUF
C                           * RANDBEGRENZUNG
C                           * SOLLEN INNENLIEGENDE UNTERTEILUNGSSTRICHE
C                           * AUF RANDBEGRENZUNG GEZEICHNET WERDEN?
         IF ( LEDGE(IEDGEA) ) THEN
C                           * KEINE FEINUNTERTEILUNGSSTRICHE VOR
C                           * 1. GROBUNTERTEILUNGSSTRICH
            IANZ = 0
            XGUA = FRSTGU
         ELSE
C                           * INNENLIEGENDE UNTERTEILUNGSSTRICHE AUF
C                           * RANDBEGRENZUNG SOLLEN NICHT GEZEICHNET
C                           * WERDEN
C                           * LIEGEN UNTERTEILUNGSSTRICHE AUSSEN?
            IF ( .NOT. LINNEN ) THEN
C                           * ZEICHNE 1. GROBUNTERTEILUNGSSTRICH, JEDOCH
C                           * KEINE NETZLINIE
               IF ( LXACHS ) THEN
                  CALL GRJMP(FRSTGU,YLU)
                  CALL GRDRW(FRSTGU,YDRWGU)
               ELSE
                  CALL GRJMP(YLU,FRSTGU)
                  CALL GRDRW(YDRWGU,FRSTGU)
               ENDIF
            ENDIF
C                           * ANZAHL DER FEINUNTERTEILUNGSSTRICHE VOR
C                           * 2. GROBUNTERTEILUNGSSTRICH
            IANZ = IANZFU
            XFU = FRSTGU+DIFFU
C                           * ABSZISSE DES 2. GROBUNTERTEILUNGSSTRICHES
            XGUA = FRSTGU+DIFFGU
         ENDIF
      ENDIF
C                           * ZEICHNEN DER UNTERTEILUNGSSTRICHE:

      NGU = (XB-TOLRNZ-XGUA)/DIFFGU
      XGU = XGUA
      DO 10 IGU = 0,NGU
C                           * FEINUNTERTEILUNGSSTRICHE
         DO 2, I = 1,IANZ
            IF ( LXACHS ) THEN
               CALL GRJMP(XFU,YJMPFU)
               CALL GRDRW(XFU,YDRWFU)
            ELSE
               CALL GRJMP(YJMPFU,XFU)
               CALL GRDRW(YDRWFU,XFU)
            ENDIF
            XFU = XFU+DIFFU
    2    CONTINUE
         IANZ = IANZFU
         XFU = XGU+DIFFU
C                           * GROBUNTERTEILUNGSSTRICHE
         IF ( IGROBU .NE. 0 ) THEN
            IF ( LXACHS ) THEN
               IF ( IGROBU.EQ.2 .AND. XGU.NE.XA ) THEN
                  CALL GRDSH(SVDSH1,SVDSH2,SVDSH3)
                  CALL GRJMP(XGU,YJMPHU)
                  CALL GRDRW(XGU,YDRWHU)
                  CALL GRDSH(1.,0.,1.)
                  CALL GRJMP(XGU,YDRWHU)
               ELSE
                  CALL GRJMP(XGU,YJMPGU)
               ENDIF
               CALL GRDRW(XGU,YDRWGU)
            ELSE
               IF ( IGROBU.EQ.2 .AND. XGU.NE.XA   ) THEN
                  CALL GRDSH(SVDSH1,SVDSH2,SVDSH3)
                  CALL GRJMP(YJMPHU,XGU)
                  CALL GRDRW(YDRWHU,XGU)
                  CALL GRDSH(1.,0.,1.)
                  CALL GRJMP(YDRWHU,XGU)
               ELSE
                  CALL GRJMP(YJMPGU,XGU)
               ENDIF
               CALL GRDRW(YDRWGU,XGU)
            ENDIF
         ENDIF
         XGU = XGU + DIFFGU
   10 CONTINUE

C                           * LIEGT LETZTER GROBUNTERTEILUNGSSTRICH AUF
C                           * RANDBEGRENZUNG?
      IF ( XGU .LE. XB+EPSLN ) THEN
C                           * FEINUNTERTEILUNGSSTRICHE VOR LETZTEM
C                           * GROBUNTERTEILUNGSSTRICH
         DO 3, I = 1,IANZFU
            IF ( LXACHS ) THEN
               CALL GRJMP(XFU,YJMPFU)
               CALL GRDRW(XFU,YDRWFU)
            ELSE
               CALL GRJMP(YJMPFU,XFU)
               CALL GRDRW(YDRWFU,XFU)
            ENDIF
            XFU = XFU+DIFFU
    3    CONTINUE
C                           * LETZTER GROBUNTERTEILUNGSSTRICH
         IF ( IGROBU .NE. 0 ) THEN
C                           * SOLLEN INNENLIEGENDE UNTERTEILUNGSSTRICHE
C                           * AUF RANDBEGRENZUNG GEZEICHNET WERDEN?
            IF ( LEDGE(IEDGEB) ) THEN
               IF ( LXACHS ) THEN
                  CALL GRJMP(XGU,YJMPGU)
                  CALL GRDRW(XGU,YDRWGU)
               ELSE
                  CALL GRJMP(YJMPGU,XGU)
                  CALL GRDRW(YDRWGU,XGU)
               ENDIF
            ELSE
C                           * INNENLIEGENDE UNTERTEILUNGSSTRICHE AUF
C                           * RANDBEGRENZUNG SOLLEN NICHT GEZEICHNET
C                           * WERDEN
C                           * LIEGEN UNTERTEILUNGSSTRICHE AUSSEN?
               IF ( .NOT. LINNEN ) THEN
C                           * ZEICHNE LETZTEN GROBUNTERTEILUNGSSTRICH,
C                           * JEDOCH KEINE NETZLINIE
                  IF ( LXACHS ) THEN
                     CALL GRJMP(XGU,YLU)
                     CALL GRDRW(XGU,YDRWGU)
                  ELSE
                     CALL GRJMP(YLU,XGU)
                     CALL GRDRW(YDRWGU,XGU)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ELSE
C                           * LETZTER GROBUNTERTEILUNGSSTRICH LIEGT
C                           * NICHT AUF RANDBEGRENZUNG
C                           * (ER WURDE BEREITS GEZEICHNET)
C                           * ANZAHL DER FEINUNTERTEILUNGSSTRICHE HINTER
C                           * LETZTEM GROBUNTERTEILUNGSSTRICH
         IANZ = INT( ((XB+EPSLN)-(XGU-DIFFGU))/DIFFU )
C                           * SOLLEN INNENLIEGENDE UNTERTEILUNGSSTRICHE
C                           * AUF RANDBEGRENZUNG NICHT GEZEICHNET WER-
C                           * DEN, SOLLEN UNTERTEILUNGSSTRICHE INNEN
C                           * LIEGEN UND LIEGT LETZTER FEINUNTER-
C                           * TEILUNGSSTRICH AUF RANDBEGRENZUNG?
         IF (  ( .NOT. LEDGE(IEDGEB) ) .AND. LINNEN .AND.
     $      ( XFU+REAL(IANZ-1)*DIFFU .GT. XB-TOLRNZ )  )
     $         IANZ = IANZ-1
C                           * FEINUNTERTEILUNGSSTRICHE HINTER LETZTEM
C                           * GROBUNTERTEILUNGSSTRICH
         DO 4, I = 1,IANZ
            IF ( LXACHS ) THEN
               CALL GRJMP(XFU,YJMPFU)
               CALL GRDRW(XFU,YDRWFU)
            ELSE
               CALL GRJMP(YJMPFU,XFU)
               CALL GRDRW(YDRWFU,XFU)
            ENDIF
            XFU = XFU+DIFFU
    4    CONTINUE
      ENDIF
      END
C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C ----------------------------------------------------------------------
C  GRAXLB - ZEICHNET DIE LABELS FUER EINE ACHSE
C           WIRD AUFGERUFEN VON: GRAXS
C  AUTOR: G. EGERER, DATUM:  3.11.83, LETZTE AENDERUNG: 13. 4.84
C ----------------------------------------------------------------------
      SUBROUTINE GRAXLB(LXACHS,PKTLU,PKTRO,FAKCMV,CHRSZC,IBSCHR,IPRECS,
     $   CHRAND,FRSTBZ,DIFFBZ,N,LTEXT,TEXT,INTSYM,LINNEN)

C                           * EINGABEPARAMETER:
C                           * LXACHS - GIBT AN, OB EINE X-ACHSE
C                           *          (LXACHS = .TRUE.) ODER EINE
C                           *          Y-ACHSE (LXACHS = .FALSE.)
C                           *          BESCHRIFTET WERDEN SOLL
C                           * PKTLU  - LINKE UNTERE ECKE DES ACHSEN-
C                           *          BEREICHES
C                           *          (XLU = PKTLU(1); YLU = PKTLU(2) )
C                           * PKTRO  - RECHTE OBERE ECKE DES ACHSEN-
C                           *          BEREICHES
C                           * FAKCMV - FAKTOREN ZUR UMRECHNUNG VON CM-
C                           *          ANGABEN IN BENUTZERWERTE FUER
C                           *          X-ACHSE (FAKCMV(1)) UND Y-ACHSE
C                           *          (FAKCMV(2))
C                           * CHRSZC - GROESSE EINES ZEICHENS IN CM
C                           * IBSCHR - GIBT AN, OB LABELS NUR UNTERHALB
C                           *          (IBSCHR = 1), NUR OBERHALB
C                           *          (IBSCHR = 2) ODER SOWOHL UNTER-
C                           *          ALS AUCH OBERHALB DES ACHSEN-
C                           *          BEREICHES (IBSCHR = 3) ERZEUGT
C                           *          WERDEN SOLLEN
C                           * IPRECS - MAXIMALE GENAUIGKEIT DER
C                           *          BESCHRIFTUNGSZAHLEN:
C                           *          IPRECS(1)=ANZAHL DER STELLEN VOR
C                           *             DEM DEZIMALPUNKT
C                           *          IPRECS(2)=ANZAHL DER STELLEN
C                           *             HINTER DEM DEZIMALPUNKT
C                           *          BEI BESCHRIFTUNG IM E-FORMAT:
C                           *          IPRECS(1)=ANZAHL DER SIGNIFIKAN-
C                           *             TEN STELLEN
C                           *          IPRECS(2)=-32766
C                           * CHRAND - WINKEL EINES ZEICHENS DER ACHSEN-
C                           *          BESCHRIFTUNG IN GRAD
C                           * FRSTBZ - KLEINSTE BESCHRIFTUNGSZAHL
C                           * DIFFBZ - POSITIVE DIFFERENZ ZWEIER
C                           *          BESCHRIFTUNGSZAHLEN
C                           * N      - EXPONENT BEI BESCHRIFTUNG IM
C                           *          E-FORMAT
C                           * LTEXT  - LAENGE DER Zeichenkette
C                           *          "TEXT"
C                           * TEXT   - TEXT, DER AN DIE ACHSE
C                           *          GESCHRIEBEN WERDEN SOLL
C                           * INTSYM - INTENSITAET EINES ZEICHENS
C                           * LINNEN - ob Pixel innen gezeichnet werden
      LOGICAL LXACHS, LINNEN
      REAL    PKTLU(2), PKTRO(2), FAKCMV(2), CHRSZC, CHRAND, FRSTBZ,
     $        DIFFBZ
      INTEGER IBSCHR, IPRECS(2), N, LTEXT, INTSYM
      CHARACTER(len=*) TEXT

C                           * COMMON-VARIABLEN:
      COMMON      /GRAXCB/ CPRCTL,GENAU,RANDB
CDEC$ PSECT /GRAXCB/ NOSHR
      SAVE /GRAXCB/
      CHARACTER (len=1) CPRCTL,GENAU,RANDB

C                           * PROGRAMMVARIABLEN:
      CHARACTER (len=11) CFMT
      CHARACTER (len=8)  CLABEL
      CHARACTER (len=1)  CHELP
      LOGICAL      LVARFM
      INTEGER IAPOS, ICHRNO, IEPOS, MAXLEN,NBZ,IBZ
      REAL    ABZPV, ABZSV, ALBWID, AMAXWD, BZ, CHRANR, CHSZVX, CHSZVY,
     $        COSANC, COSCSZ, ERRLIM, FX, FXDFBZ, FY, HELP, ORDPV,
     $        ORDSV, SINANC, SINCSZ, XDISTV, XEPOS, XLEFT, XMAX, XPOS1,
     $        XPOS2, XRIGHT, XSWIDV, XTEST, YBOT, YDISTV, YEPOS1,
     $        YEPOS2, YINVA1, YINVA2, YPOS1, YPOS2, YTOP, YUP, FDIST

C                           * BENAMTE PROGRAMMKONSTANTEN:
      REAL      CFDSTA, CFDSTI, CFEXSZ, CYUP
      PARAMETER (CFDSTA=2.,CFEXSZ=2./3.,CYUP=1.-2./3.*CFEXSZ,CFDSTI=1.)
C                           * ZWISCHENRAUMBREITE FUER ZEICHEN DER
C                           * GROESSE 1 = 1-1/SQRT(3)
      REAL      CSPACE
      PARAMETER (CSPACE = 1.-1./1.7320508076)
      REAL      CEPS1, CTLRNZ
      PARAMETER (CEPS1 = 2.E-4, CTLRNZ = 5.E-5)
      REAL      CPI
      PARAMETER (CPI = 3.1415926536)
C                           * ------------------------------------------

      AMAXWD = 0.

      CHSZVX = CHRSZC*FAKCMV(1)
      CHSZVY = CHRSZC*FAKCMV(2)

      FDIST = CFDSTA
      IF ( LINNEN ) FDIST = CFDSTI
      XDISTV = CHSZVX*FDIST
      YDISTV = CHSZVY*FDIST

      IF ( PKTLU(1) .GT. PKTRO(1) ) THEN
C                           * BENUTZERWERTE VON LINKS NACH RECHTS
C                           * ABFALLEND
C                           * GROESSTE ABSZISSE DES ACHSENBEREICHES
         XMAX = PKTLU(1)
C                           * RICHTUNGSAENDERUNGSFAKTOR FUER DEN WECHSEL
C                           * VOM CM- ZUM BENUTZERWERTE(VALUE)-SYSTEM
         FX = -1.
         XDISTV = -XDISTV
         CHSZVX = -CHSZVX
      ELSE
C                           * BENUTZERWERTE VON LINKS NACH RECHTS
C                           * WACHSEND
         XMAX = PKTRO(1)
         FX = 1.
      ENDIF

C                           * ABSZISSE ZUM LINKEN RAND DES
C                           * BESCHRIFTUNGSBEREICHES
      XLEFT = PKTLU(1)-XDISTV
C                           * ABSZISSE ZUM RECHTEN RAND
      XRIGHT = PKTRO(1)+XDISTV

      IF ( PKTLU(2) .GT. PKTRO(2) ) THEN
C                           * BENUTZERWERTE VON UNTEN NACH OBEN
C                           * ABFALLEND
         FY = -1.
         YDISTV = -YDISTV
         CHSZVY = -CHSZVY
      ELSE
C                           * BENUTZERWERTE VON UNTEN NACH OBEN WACHSEND
         FY = 1.
      ENDIF
C                           * SOLL DIE ACHSE BESCHRIFTET WERDEN?
      IF (  ( IPRECS(1) .NE. 0 ) .OR. ( IPRECS(2) .NE. 0 )  ) THEN

C                           * UMRECHNUNG VON CHRAND IN BOGENMASS
         CHRANR = CPI/180.*CHRAND

         SINANC = SIN(CHRANR)
         COSANC = COS(CHRANR)
         IF ( ABS(SINANC) .LE. CTLRNZ ) THEN
            SINANC = 0.
         ELSE
            IF ( ABS(COSANC) .LE. CTLRNZ ) COSANC = 0.
         ENDIF

C                           * ORDINATE DES ABSTANDSPUNKTES S BEZUEGLICH
C                           * DES VALUE-SYSTEMS (WINKEL ZWISCHEN X-ACHSE
C                           * UND DER VERBINDUNG DES NULLPUNKTES MIT S
C                           * IM CM-SYSTEM = CHRAND, ENTFERNUNG DES
C                           * PUNKTES S VOM NULLPUNKT IM CM-SYSTEM =
C                           * CHRSZC)
         SINCSZ = SINANC*CHRSZC
         ORDSV = SINCSZ*FAKCMV(2)*FY
C                           * ABSZISSE VON S BZGL. DES VALUE-SYSTEMS
         COSCSZ =COSANC*CHRSZC
         ABZSV = COSCSZ*FAKCMV(1)*FX

C                           * ORDINATE DES EBENENPUNKTES P BZGL. DES
C                           * VALUE-SYSTEMS (WINKEL ZWISCHEN X-ACHSE UND
C                           * DER VERBINDUNG DES NULLPUNKTES MIT P IM
C                           * CM-SYSTEM = CHRAND+90, ENTFERNUNG DES
C                           * PUNKTES P VOM NULLPUNKT IM CM-SYSTEM =
C                           * CHRSZC)
         ORDPV = COSCSZ*FAKCMV(2)*FY
C                           * ABSZISSE VON P BZGL. DES VALUE-SYSTEMS
         ABZPV = SINCSZ*FAKCMV(1)*FX

C                           * ERMITTLUNG DER ORDINATEN (YPOS1 (ODER
C                           * YINVA1) UND EVTL. YPOS2 (ODER YINVA2)) DER
C                           * LINKEN UNTEREN ECKPUNKTE DER ZU ERZEUGEN-
C                           * DEN LABEL, DER ORDINATE ZUM UNTEREN RAND
C                           * (YBOT) UND/ODER ZUM OBEREN RAND DES BE-
C                           * SCHRIFTUNGSBEREICHES (YTOP)
C                           * (VON DER BREITE EINES LABELS ABHAENGIGE
C                           * ORDINATEN WERDEN UNVOLLSTAENDIG ERMITTELT;
C                           * ES FEHLT JEWEILS DAS VON DER BREITE AB-
C                           * HAENGIGE TEILSTUECK; YINVA1, YINVA2 SIND
C                           * UNVOLLSTAENDIGE ORDINATEN):
         IF ( SINANC .GE. 0. ) THEN
            IF ( COSANC .GE. 0. ) THEN
C                           * 0 <= CHRAND <= 90
               IF ( IBSCHR .NE. 2 ) THEN
                  IF ( SINANC .EQ. 0. ) THEN
C                           * CHRAND = 0
                     YPOS1 = PKTLU(2)-YDISTV-ORDPV
                     YBOT = YPOS1
                  ELSE
                     YINVA1 = PKTLU(2)-YDISTV-ORDPV
                     YBOT = YINVA1
                  ENDIF
                  IF ( IBSCHR .EQ. 3 ) THEN
                     YPOS2 = PKTRO(2)+YDISTV
                     YTOP = YPOS2+ORDPV
                  ENDIF
               ELSE
                  YPOS1 = PKTRO(2)+YDISTV
                  YTOP = YPOS1+ORDPV
               ENDIF
            ELSE
C                           * 90 < CHRAND <= 180
               IF ( IBSCHR .NE. 2 ) THEN
                  IF ( SINANC .EQ. 0. ) THEN
C                           * CHRAND = 180
                     YPOS1 = PKTLU(2)-YDISTV
                     YBOT = YPOS1+ORDPV
                  ELSE
                     YINVA1 = PKTLU(2)-YDISTV
                     YBOT = YINVA1+ORDPV
                  ENDIF
                  IF ( IBSCHR .EQ. 3 ) THEN
                     YPOS2 = PKTRO(2)+YDISTV-ORDPV
                     YTOP = YPOS2
                  ENDIF
               ELSE
                  YPOS1 = PKTRO(2)+YDISTV-ORDPV
                  YTOP = YPOS1
               ENDIF
            ENDIF
         ELSE
            IF ( COSANC .LE. 0. ) THEN
C                           * 180 < CHRAND <= 270
               IF ( IBSCHR .NE. 2 ) THEN
                  YPOS1 = PKTLU(2)-YDISTV
                  YBOT = YPOS1+ORDPV
                  IF ( IBSCHR .EQ. 3 ) THEN
                     YINVA2 = PKTRO(2)+YDISTV-ORDPV
                     YTOP = YINVA2
                  ENDIF
               ELSE
                  YINVA1 = PKTRO(2)+YDISTV-ORDPV
                  YTOP = YINVA1
               ENDIF
            ELSE
C                           * 270 < CHRAND < 360
               IF ( IBSCHR .NE. 2 ) THEN
                  YPOS1 = PKTLU(2)-YDISTV-ORDPV
                  YBOT = YPOS1
                  IF ( IBSCHR .EQ. 3 ) THEN
                     YINVA2 = PKTRO(2)+YDISTV
                     YTOP = YINVA2+ORDPV
                  ENDIF
               ELSE
                  YINVA1 = PKTRO(2)+YDISTV
                  YTOP = YINVA1+ORDPV
               ENDIF
            ENDIF
         ENDIF

C                           * ERZEUGUNG DES FORMATS FUER DIE UMWANDLUNG
C                           * DER BESCHRIFTUNGSZAHLEN IN ZEICHENKETTEN
C                           * MITTELS INTERNEM I/O (LVARFM = .TRUE.
C                           * ZEIGT AN, DASS DAS FORMAT EVTL. NICHT FUER
C                           * ALLE UMZUWANDELNDE ZAHLEN BENUTZT WERDEN
C                           * KANN):
         LVARFM = .FALSE.
         IF ( IPRECS(2) .GT. 0 ) THEN
            IF ( IPRECS(1) .EQ. 0 ) THEN
               IF ( FRSTBZ .LT. 0. ) THEN
C                           * DATENFELDLAENGE W FUER FORMAT = ANZAHL DER
C                           * STELLEN HINTER DEM DEZIMALPUNKT +1 (FUER
C                           * ".") +1 (FUER "-")
                  IEPOS = IPRECS(2)+2
                  WRITE ( CFMT, '( ''(F'',I1,''.'',I1,'')'' )' ) IEPOS,
     $               IPRECS(2)
C                           * FUER 0, SOWIE POSITIVE ZAHLEN IST DIESES
C                           * FORMAT NICHT ZU BENUTZEN, DAMIT DIE
C                           * ZEICHENKETTENFORM VON 0 KEINEN DEZIMAL-
C                           * PUNKT UND KEINE NACHPUNKTSTELLEN, VON
C                           * POSITIVEN ZAHLEN KEINE "0" VOR DEM DEZI-
C                           * MALPUNKT ENTHALTEN
                  LVARFM = .TRUE.
                  ERRLIM = 0.5*10.**(-IPRECS(2))
               ELSE IF ( FRSTBZ .EQ. 0. ) THEN
C                           * W = 2 (ERFORDERLICH, DA REAL-ZAHL NICHT
C                           * OHNE "." MITTELS INTERNEM I/O IN ZEICHEN-
C                           * KETTE UMGEWANDELT WERDEN KANN; "." WIRD
C                           * HIER JEDOCH NICHT BESTANDTEIL DER ZEICHEN-
C                           * KETTENFORM DER ZAHL)
                  IEPOS = 1
                  CFMT = '(F2.0)'
                  LVARFM = .TRUE.
                  ERRLIM = 0.5*10.**(-IPRECS(2))
               ELSE
C                           * W = ANZAHL DER STELLEN HINTER DEM DEZIMAL-
C                           * PUNKT +1 (FUER ".")
                  IEPOS = IPRECS(2)+1
                  WRITE ( CFMT, '( ''(F'',I1,''.'',I1,'')'' )' ) IEPOS,
     $               IPRECS(2)
               ENDIF
            ELSE
C                           * W = MAXIMALE ANZAHL DER STELLEN VOR DEM
C                           * DEZIMALPUNKT + ANZAHL DER STELLEN HINTER
C                           * DEM DEZIMALPUNKT +1 (FUER ".") +1 (FUER
C                           * MOEGLICHES "-")
               IEPOS = IPRECS(1)+IPRECS(2)+2
               WRITE ( CFMT, '( ''(F'',I1,''.'',I1,'')'' )' ) IEPOS,
     $            IPRECS(2)
            ENDIF
         ELSE IF ( IPRECS(2) .EQ. -32766 ) THEN
C                           * BESCHRIFTUNG IM E-FORMAT (LABELS SIND DIE
C                           * MANTISSEN DER BESCHRIFTUNGSZAHLEN):
C                           * W = MAXIMALE ANZAHL DER SIGNIFIKANTEN
C                           * STELLEN +1 (FUER MOEGLICHES "-") +1 (FUER
C                           * ".", DER JEDOCH NICHT BESTANDTEIL DER MIT
C                           * DIESEM FORMAT ZU ERZEUGENDEN ZEICHENKETTEN
C                           * WIRD)
            IEPOS = IPRECS(1)+1
            WRITE ( CFMT, '( ''('',I3,''P,F'',I1,''.0)'' )' ) -N,
     $         IEPOS+1
         ELSE
C                           * IPRECS(2) = 0
C                           * W = MAXIMALE ANZAHL DER STELLEN VOR DEM
C                           * DEZIMALPUNKT +1 (FUER MOEGLICHES "-")
C                           * +1 (FUER ".", DER JEDOCH NICHT BESTANDTEIL
C                           * DER MIT DIESEM FORMAT ZU ERZEUGENDEN ZEI-
C                           * CHENKETTEN WIRD)
            IEPOS = IPRECS(1)+1
            WRITE ( CFMT, '( ''(F'',I1,''.0)'' )' ) IEPOS+1
         ENDIF

         FXDFBZ = FX*DIFFBZ
C                           * SCHLEIFE FUER LABELERZEUGUNG:
         NBZ = (XMAX+CEPS1*FAKCMV(1)-FRSTBZ)/DIFFBZ
         DO 10 IBZ = 0,NBZ
            BZ = FRSTBZ + DIFFBZ*IBZ

C                           * UMWANDLUNG DER BESCHRIFTUNGSZAHL BZ IN
C                           * ZEICHENKETTE:
            IF ( LVARFM ) THEN
               IF ( BZ .LE. -ERRLIM ) THEN
C                           * BZ < 0
               ELSE IF ( BZ .LT. ERRLIM ) THEN
C                           * BZ = 0
C                           * ENTHAELT CFMT FORMAT FUER NEGATIVE ZAHLEN?
                  IF ( IEPOS .NE. 1 ) THEN
C                           * ERZEUGE DAS FUER 0 ZU BENUTZENDE FORMAT
                     IEPOS = 1
                     CFMT = '(F2.0)'
                  ENDIF
               ELSE
C                           * BZ > 0
C                           * CFMT ENTHAELT FORMAT FUER NEGATIVE ZAHLEN
C                           * ODER FUER 0
C                           * ERZEUGE FORMAT FUER POSITIVE ZAHLEN
                  IEPOS = IPRECS(2)+1
                  WRITE ( CFMT, '( ''(F'',I1,''.'',I1,'')'' )' ) IEPOS,
     $               IPRECS(2)
                  LVARFM = .FALSE.
               ENDIF
            ENDIF

            WRITE ( CLABEL, CFMT ) BZ

            IF (GENAU .EQ. '1') THEN
               IAPOS = INDEX(CLABEL,' .')
               IF (IAPOS .GT. 0) THEN
                  CLABEL(IAPOS:IAPOS+1) = '0.'
               ELSE
                  IAPOS = INDEX(CLABEL,' -.')
                  IF (IAPOS .GT. 0) CLABEL(IAPOS:IAPOS+2) = '-0.'
               ENDIF
            ENDIF

C                           * SUCHEN DES 1. VON BLANK VERSCHIEDENEN
C                           * ZEICHENS IN DER ERZEUGTEN ZEICHENKETTE
            DO 2, IAPOS = 1, IEPOS-1
               IF ( CLABEL(IAPOS:IAPOS) .NE. ' ' ) GO TO 3
    2       CONTINUE

    3       ICHRNO = IEPOS-IAPOS+1
C                           * BREITE DES LABELS FUER BZ BEZOGEN AUF
C                           * ZEICHEN DER GROESSE 1
            ALBWID = REAL(ICHRNO)-CSPACE
C                           * ABSZISSE DES PUNKTES S' BZGL. DES VALUE-
C                           * SYSTEMS (WINKEL ZWISCHEN X-ACHSE UND DER
C                           * VERBINDUNG DES NULLPUNKTES MIT S' IM CM-
C                           * SYSTEM = CHRAND, ENTFERNUNG DES PUNKTES
C                           * S' VOM NULLPUNKT IM CM-SYSTEM =
C                           * ALBWID*CHRSZC)
            XSWIDV = ALBWID*ABZSV

C                           * ERMITTLUNG DER ABSZISSEN (XPOS1 UND EVTL.
C                           * XPOS2), SOWIE DER BISHER NOCH UNVOLLSTAEN-
C                           * DIG ERMITTELTEN ORDINATEN DER LINKEN
C                           * UNTEREN ECKPUNKTE DER ZU ZEICHNENDEN LABEL
C                           * FUER BZ
            IF ( SINANC .EQ. 0. ) THEN
C                           * CHRAND = 0 | CHRAND = 180
               XPOS1 = BZ-XSWIDV/2.
               IF ( IBSCHR .EQ. 3 ) XPOS2 = XPOS1
            ELSE
               IF (  ( IBSCHR .NE. 2 ) .EQV. ( SINANC .GT. 0. )  ) THEN
                  XPOS1 = BZ-XSWIDV+ABZPV/2.
                  IF ( IBSCHR .EQ. 3 ) XPOS2 = XPOS1+XSWIDV
                  YPOS1 = YINVA1-ALBWID*ORDSV
C                           * ALBWID*ORDSV = ORDINATE VON S'
               ELSE
                  XPOS1 = BZ+ABZPV/2.
                  IF ( IBSCHR .EQ. 3 ) THEN
                     XPOS2 = XPOS1-XSWIDV
                     YPOS2 = YINVA2-ALBWID*ORDSV
C                           * ALBWID*ORDSV = ORDINATE VON S'
                  ENDIF
               ENDIF
            ENDIF

            IF ( IBSCHR .NE. 3 ) THEN
C                           * E I N E  ACHSE IST MIT LABELS ZU VERSEHEN

C                           * LIEGT BZ-FXDFBZ LINKS VOM BESCHRIFTUNGS-
C                           * BEREICH?
C                           * (BZ-FXDFBZ LIEGT LINKS VON ALLEN ABSZISSEN
C                           * DES ZEICHENBEREICHES DES LABELS FUER BZ)
               IF ( ((BZ-FXDFBZ)-XLEFT)*FX .LT. 0. ) THEN
C                           * TESTE, OB DER AM WEITESTEN LINKS LIEGENDE
C                           * ECKPUNKT DES LABELS FUER BZ LINKS VOM
C                           * BESCHRIFTUNGSBEREICH LIEGT:
                  IF ( SINANC .GE. 0. ) THEN
                     IF ( COSANC .GE. 0. ) THEN
C                           * 0 <= CHRAND <= 90
                        XTEST = XPOS1-ABZPV
                     ELSE
C                           * 90 < CHRAND <= 180
                        XTEST = XPOS1+XSWIDV-ABZPV
                     ENDIF
                  ELSE
                     IF ( COSANC .LE. 0. ) THEN
C                           * 180 < CHRAND <= 270
                        XTEST = XPOS1+XSWIDV
                     ELSE
C                           * 270 < CHRAND < 360
                        XTEST = XPOS1
                     ENDIF
                  ENDIF

                  IF ((XTEST-XLEFT)*FX.LT.0. .AND. RANDB.EQ.'0')GOTO 9
               ENDIF

C                           * LIEGT BZ+FXDFBZ (LIEGT RECHTS VON ALLEN
C                           * ABSZISSEN DES LABELBEREICHES) RECHTS VOM
C                           * BESCHRIFTUNGSBEREICH?
               IF ( ((BZ+FXDFBZ)-XRIGHT)*FX .GT. 0. ) THEN
C                           * TESTE, OB DER AM WEITESTEN RECHTS LIEGENDE
C                           * LABELECKPUNKT RECHTS VOM BESCHRIFTUNGS-
C                           * BEREICH LIEGT:
                  IF ( SINANC .GE. 0. ) THEN
                     IF ( COSANC .GE. 0. ) THEN
C                           * 0 <= CHRAND <= 90
                        XTEST = XPOS1+XSWIDV
                     ELSE
C                           * 90 < CHRAND <= 180
                        XTEST = XPOS1
                     ENDIF
                  ELSE
                     IF ( COSANC .LE. 0. ) THEN
C                           * 180 < CHRAND <= 270
                        XTEST = XPOS1-ABZPV
                     ELSE
C                           * 270 < CHRAND < 360
                        XTEST = XPOS1+XSWIDV-ABZPV
                     ENDIF
                  ENDIF

                  IF((XTEST-XRIGHT)*FX.GT.0. .AND. RANDB.EQ.'0' )GOTO 9
               ENDIF

C                           * ERMITTLE GROESSTE BREITE EINES LABELS
               IF ( ALBWID .GT. AMAXWD ) AMAXWD = ALBWID

C                           * ZEICHNEN DES LABELS
               IF ( LXACHS ) THEN
                  CALL GRTXT(XPOS1,YPOS1,ICHRNO,CLABEL(IAPOS:IEPOS))
               ELSE
                  CALL GRTXT(YPOS1,XPOS1,ICHRNO,CLABEL(IAPOS:IEPOS))
               ENDIF

            ELSE
C                           * Z W E I  ACHSEN SIND MIT LABELS ZU
C                           * VERSEHEN

               IF ( ((BZ-FXDFBZ)-XLEFT)*FX .LT. 0. ) THEN
C                           * TESTE, OB DER ECKPUNKT, DER VON ALLEN
C                           * ECKPUNKTEN DER 2 LABEL FUER BZ AM WEITE-
C                           * STEN LINKS LIEGT, LINKS VOM BESCHRIFTUNGS-
C                           * BEREICH LIEGT:
                  IF ( SINANC .GE. 0. ) THEN
                     IF ( COSANC .GE. 0. ) THEN
C                           * 0 <= CHRAND <= 90
                        XTEST = XPOS1-ABZPV
                     ELSE
C                           * 90 < CHRAND <= 180
                        XTEST = XPOS2+XSWIDV-ABZPV
                     ENDIF
                  ELSE
                     IF ( COSANC .LE. 0. ) THEN
C                           * 180 < CHRAND <= 270
                        XTEST = XPOS1+XSWIDV
                     ELSE
C                           * 270 < CHRAND < 360
                        XTEST = XPOS2
                     ENDIF
                  ENDIF
                  IF ( (XTEST-XLEFT)*FX .LT. 0. ) GO TO 9
               ENDIF

               IF ( ((BZ+FXDFBZ)-XRIGHT)*FX .GT. 0. ) THEN
C                           * TESTE, OB DER ECKPUNKT, DER VON ALLEN
C                           * ECKPUNKTEN DER 2 LABEL AM WEITESTEN
C                           * RECHTS LIEGT, RECHTS VOM BESCHRIFTUNGS-
C                           * BEREICH LIEGT:
                  IF ( SINANC .GE. 0. ) THEN
                     IF ( COSANC .GE. 0. ) THEN
C                           * 0 <= CHRAND <= 90
                        XTEST = XPOS2+XSWIDV
                     ELSE
C                           * 90 < CHRAND <= 180
                        XTEST = XPOS1
                     ENDIF
                  ELSE
                     IF ( COSANC .LE. 0. ) THEN
C                           * 180 < CHRAND <= 270
                        XTEST = XPOS2-ABZPV
                     ELSE
C                           * 270 < CHRAND < 360
                        XTEST = XPOS1+XSWIDV-ABZPV
                     ENDIF
                  ENDIF
                  IF ( (XTEST-XRIGHT)*FX .GT. 0. ) GO TO 9
               ENDIF

C                           * ERMITTLE GROESSTE BREITE EINES LABELS
               IF ( ALBWID .GT. AMAXWD ) AMAXWD = ALBWID

C                           * ZEICHNEN DER LABEL
               IF ( LXACHS ) THEN
                  CALL GRTXT(XPOS1,YPOS1,ICHRNO,CLABEL(IAPOS:IEPOS))
                  CALL GRTXT(XPOS2,YPOS2,ICHRNO,CLABEL(IAPOS:IEPOS))
               ELSE
                  CALL GRTXT(YPOS1,XPOS1,ICHRNO,CLABEL(IAPOS:IEPOS))
                  CALL GRTXT(YPOS2,XPOS2,ICHRNO,CLABEL(IAPOS:IEPOS))
               ENDIF

            ENDIF

    9       CONTINUE
   10    CONTINUE

      ENDIF

C                           * WURDE ACHSE BESCHRIFTET?
      IF ( AMAXWD .NE. 0. ) THEN

C                           * HILFSGROESSE ZUR ERMITTLUNG DER BISHER
C                           * NOCH UNVOLLSTAENDIG ERMITTELTEN ORDINATEN
C                           * YBOT UND YTOP:
C                           *    YBOT = YBOT-HELP
C                           *    YTOP = YTOP+HELP
C                           * (HELP = BETRAGSMAESSIG GROESSTE ORDINATE
C                           * ALLER PUNKTE S')
         HELP = FY*ABS(AMAXWD*ORDSV)

C                           * EXPONENTENANGABE ZU ZEICHNEN?
         IF (  ( IPRECS(2) .EQ. -32766 ) .AND. ( N .NE. 0 )  ) THEN

C                           * UMWANDLUNG DES EXPONENTEN N IN ZEICHEN-
C                           * KETTE:
            WRITE ( CLABEL, '( I3 )' ) N
C                           * SUCHEN DES 1. VON BLANK VERSCHIEDENEN
C                           * ZEICHENS IN DER ERZEUGTEN ZEICHENKETTE
            DO 4, IAPOS = 1,2
               IF ( CLABEL(IAPOS:IAPOS) .NE. ' ' ) GO TO 5
    4       CONTINUE

    5       ICHRNO = 4-IAPOS
C                           * ABSZISSE FUER DEN EXPONENTEN
            XEPOS = XRIGHT-(REAL(ICHRNO)-CSPACE)*CFEXSZ*CHSZVX
C                           * ABSZISSE FUER EXPONENTENANGABE
            XPOS1 = XEPOS-3.*CHSZVX

            YDISTV = YDISTV/2.

C                           * ERMITTLUNG DER ORDINATEN FUER EXPONENTEN-
C                           * ANGABE (YPOS1 UND/ODER YPOS2; ANGABE IST
C                           * FUER "IBSCHR = 3" 2-MAL ZU ZEICHNEN), FUER
C                           * DEN EXPONENTEN (YEPOS1 UND/ODER YEPOS2),
C                           * SOWIE DER ORDINATE ZUM UNTEREN RAND (YBOT)
C                           * UND/ODER ZUM OBEREN RAND DES BESCHRIF-
C                           * TUNGSBEREICHES (YTOP):
            YUP = CYUP*CHSZVY
            IF ( IBSCHR .NE. 2 ) THEN
               YPOS1 = YBOT-HELP-YDISTV-CHSZVY
               YBOT = YPOS1
               YEPOS1 = YPOS1+YUP
            ENDIF
            IF ( IBSCHR .NE. 1 ) THEN
               YPOS2 = YTOP+HELP+YDISTV
               YTOP = YPOS2+CHSZVY
               YEPOS2 = YPOS2+YUP
            ENDIF

C                           * ZEICHNEN DER EXPONENTENANGABE
            IF ( LXACHS ) THEN
               CALL GRCHRC(CHRSZC,0.,INTSYM)
               IF ( IBSCHR .NE. 2 ) CALL GRTXT(XPOS1,YPOS1,3,'*10')
               IF ( IBSCHR .NE. 1 ) CALL GRTXT(XPOS1,YPOS2,3,'*10')
               CALL GRCHRC(CFEXSZ*CHRSZC,0.,INTSYM)
               IF ( IBSCHR .NE. 2 )
     $            CALL GRTXT(XEPOS,YEPOS1,ICHRNO,CLABEL(IAPOS:3))
               IF ( IBSCHR .NE. 1 )
     $            CALL GRTXT(XEPOS,YEPOS2,ICHRNO,CLABEL(IAPOS:3))
            ELSE
               CALL GRCHRC(CHRSZC,90.,INTSYM)
               IF ( IBSCHR .NE. 2 ) CALL GRTXT(YPOS1,XPOS1,3,'*10')
               IF ( IBSCHR .NE. 1 ) CALL GRTXT(YPOS2,XPOS1,3,'*10')
               CALL GRCHRC(CFEXSZ*CHRSZC,90.,INTSYM)
               IF ( IBSCHR .NE. 2 )
     $            CALL GRTXT(YEPOS1,XEPOS,ICHRNO,CLABEL(IAPOS:3))
               IF ( IBSCHR .NE. 1 )
     $            CALL GRTXT(YEPOS2,XEPOS,ICHRNO,CLABEL(IAPOS:3))
            ENDIF

         ELSE
C                           * ACHSE WURDE BESCHRIFTET, JEDOCH IST KEINE
C                           * EXPONENTENANGABE ZU ZEICHNEN
            IF ( IBSCHR .NE. 2 ) YBOT = YBOT-HELP
            IF ( IBSCHR .NE. 1 ) YTOP = YTOP+HELP

         ENDIF

      ELSE
C                           * ACHSE WURDE NICHT BESCHRIFTET
         IF ( IBSCHR .NE. 2 ) YBOT = PKTLU(2)
         IF ( IBSCHR .NE. 1 ) YTOP = PKTRO(2)
         YDISTV = YDISTV+CHSZVY

      ENDIF

C                           * SOLL AN DIE ACHSE TEXT GESCHRIEBEN WERDEN?
      IF ( LTEXT .GT. 0 ) THEN

C                           * MAXIMALE LAENGE FUER EINEN AN DIE ACHSE ZU
C                           * SCHREIBENDEN TEXT
         MAXLEN = INT((XRIGHT-XLEFT)/CHSZVX+CSPACE)
         IF ( MAXLEN .GT. 0 ) THEN
C                           * ZU SCHREIBENDER TEXT ZU LANG?
            IF ( LTEXT .GT. MAXLEN ) THEN
               IF ( LXACHS ) THEN
                  CHELP = 'X'
               ELSE
                  CHELP = 'Y'
               ENDIF
               WRITE ( *, '( '''//CPRCTL//''',A )' )
     $            '**GRAXS: TEXT FUER '//CHELP//
     $               '-ACHSE ZU LANG. ER WIRD RECHTS ABGESCHNITTEN.'
               CPRCTL = '0'
               ICHRNO = MAXLEN
            ELSE
               ICHRNO = LTEXT
            ENDIF

C                           * ABSZISSE FUER TEXT (TEXT WIRD ZENTRIERT)
            XPOS1 = (XLEFT+XRIGHT)*.5

C                           * ERMITTLUNG DER ORDINATE UND ZEICHNEN DES
C                           * TEXTES:
            CALL GSTXAL(2,4)
            IF ( LXACHS ) THEN
               IF ( IBSCHR .NE. 2 ) THEN
                  YPOS1 = YBOT-YDISTV-CHSZVY
               ELSE
                  YPOS1 = YTOP+YDISTV
               ENDIF
               CALL GRCHRC(CHRSZC,0.,INTSYM)
               CALL GRTXT(XPOS1,YPOS1,ICHRNO,TEXT)
            ELSE
               IF ( IBSCHR .NE. 1 ) THEN
                  YPOS1 = YTOP+YDISTV
               ELSE
                  YPOS1 = YBOT-YDISTV-CHSZVY
               ENDIF
               CALL GRCHRC(CHRSZC,90.,INTSYM)
               CALL GRTXT(YPOS1,XPOS1,ICHRNO,TEXT)
            ENDIF
            CALL GSTXAL(0,4)

         ELSE
C                           * MAXLEN = 0
            IF ( LXACHS ) THEN
               CHELP = 'X'
            ELSE
               CHELP = 'Y'
            ENDIF
            WRITE ( *, '( '''//CPRCTL//''',A,/,1X,A )' )
     $         '**GRAXS: TEXT FUER '//CHELP//
     $            '-ACHSE WIRD WEGGELASSEN, DA DER VERFUEGBARE PLATZ',
     $         '         ZUM SCHREIBEN EINES TEXTES NICHT AUSREICHT.'
            CPRCTL = '0'

         ENDIF

      ENDIF

      IF ( LXACHS ) THEN
         HELP = CHRAND
      ELSE
         HELP = CHRAND+90.
      ENDIF
      CALL GRCHRC(CHRSZC,HELP,INTSYM)

      END
