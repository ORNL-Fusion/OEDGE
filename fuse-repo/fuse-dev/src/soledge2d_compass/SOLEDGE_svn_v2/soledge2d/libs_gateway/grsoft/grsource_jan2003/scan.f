C@PROCESS OPT(2) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C ----------------------------------------------------------------------
C  SCAN - SCANNEN DES TEXTES "CTEXT"
C  AUTOR: G. EGERER, DATUM: 26. 9.85, LETZTE AENDERUNG: 17. 1.88
C ----------------------------------------------------------------------
      SUBROUTINE SCAN(CTEXT,CELEMT,IETAB)

C                           * EINGABEPARAMETER:
C                           * CTEXT  - TEXTZEILE
C                           * CELEMT - FORMALE BESCHREIBUNG DER IN DER
C                           *          TEXTZEILE ERLAUBTEN ELEMENTE
C                           *          DURCH EINE ZEICHENKETTE. DIESE
C                           *          DARF NUR DIE FOLGENDEN ZEICHEN
C                           *          ENTHALTEN:
C                           *                     | STEHT FUER
C                           *             ZEICHEN | DAS ELEMENT
C                           *            ---------+--------------------
C                           *                C    | TEXT ODER IN ' EIN-
C                           *                     | GESCHLOSSENE
C                           *                     | ZEICHENKETTE
C                           *                I    | INTEGER-ZAHL MIT
C                           *                     | GROSSEM WERTE-
C                           *                     | BEREICH (INTEGER*4)
C                           *                S    | INTEGER-ZAHL MIT
C                           *                     | KLEINEM WERTE-
C                           *                     | BEREICH (INTEGER*2)
C                           *                R    | EINFACH GENAUE
C                           *                     | REAL-ZAHL
C                           *                D    | DOPPELT GENAUE
C                           *                     | REAL-ZAHL
C                           *          DIE ELEMENTE IN DER TEXTZEILE
C                           *          MUESSEN DURCH LEERZEICHEN
C                           *          VONEINANDER GETRENNT SEIN.

      CHARACTER(len=*) CTEXT, CELEMT

C                           * AUSGABEPARAMETER:
C                           * IETAB  - TABELLE, DIE JEWEILS DIE ANFANGS-
C                           *          UND DIE ENDPOSITION IN DER TEXT-
C                           *          ZEILE VON ALLEN GUELTIGEN ELE-
C                           *          MENTEN ENTHAELT.
C                           *          FUER UNGUELTIGE ELEMENTE ENT-
C                           *          HAELT DIESE TABELLE DIE FOLGENDE
C                           *          INFORMATION:
C                           *                      | ANF.- | END-
C                           *           FEHLER     | POS.  | POSITION
C                           *          ------------+-------+-----------
C                           *           SYNTAX-    |   -1  | POSITION
C                           *           FEHLER     |       | DES 1. UN-
C                           *                      |       | GUELTIGEN
C                           *                      |       | ZEICHENS
C                           *          ------------+-------+-----------
C                           *           WERT       |   -2  | KORREKTE
C                           *           AUSSERHALB |       | ENDPOSI-
C                           *           DES WERTE- |       | TION DER
C                           *           BEREICHES  |       | ZAHL
C                           *          ------------+-------+-----------
C                           *           NICHT ZU-  |   -3  |     0
C                           *           LAESSIGE   |       |
C                           *           FORMALE    |       |
C                           *           ELEMENTBE- |       |
C                           *           SCHREIBUNG |       |
      INTEGER IETAB(2,*)

C                           * BENAMTE PROGRAMMKONSTANTEN:

C                           * KONSTANTEN ZUR FESTLEGUNG DES WERTE-
C                           * BEREICHES FUER STANDARD-INTEGER-ZAHLEN
C                           * (INTEGER*4)
      INTEGER            IMAXLI
      PARAMETER          (IMAXLI = 10)
CIBMVSI*4      PARAMETER          (IMAXLI = 10)
CIBMPCI*4      PARAMETER          (IMAXLI = 10)
      CHARACTER(len=IMAXLI) CIMINA, CIMAXA
      PARAMETER          (CIMINA = '2147483647',
     .                    CIMAXA = '2147483647')
CIBMVSI*4      PARAMETER          (CIMINA = '2147483647',
CIBMVSI*4     .                    CIMAXA = '2147483647')
CIBMPCI*4      PARAMETER          (CIMINA = '2147483648',
CIBMPCI*4     .                    CIMAXA = '2147483647')

C                           * KONSTANTEN ZUR FESTLEGUNG DES WERTE-
C                           * BEREICHES FUER SHORT-INTEGER-ZAHLEN
C                           * (INTEGER*2)
      INTEGER            IMAXLS
      PARAMETER          (IMAXLS = 5)
CIBMVSI*2      PARAMETER          (IMAXLS = 5)
CIBMPCI*2      PARAMETER          (IMAXLS = 5)
      CHARACTER(len=IMAXLS) CSMINA, CSMAXA
      PARAMETER          (CSMINA = '32767',
     .                    CSMAXA = '32767')
CIBMVSI*2      PARAMETER          (CSMINA = '32767',
CIBMVSI*2     .                    CSMAXA = '32767')
CIBMPCI*2      PARAMETER          (CSMINA = '32768',
CIBMPCI*2     .                    CSMAXA = '32767')

C                           * KONSTANTEN ZUR FESTLEGUNG DES WERTE-
C                           * BEREICHES FUER EINFACH GENAUE REAL-ZAHLEN
      INTEGER   IRMINE, IRMAXE, IRMINM, IRMAXM
      PARAMETER (IRMINE =  -79,
     .           IRMAXE =   75,
     .           IRMINM = 5398,
     .           IRMAXM = 7236)
CIBMVSR*4      PARAMETER (IRMINE =  -79,
CIBMVSR*4     .           IRMAXE =   75,
CIBMVSR*4     .           IRMINM = 5398,
CIBMVSR*4     .           IRMAXM = 7236)
CIBMPCR*4      PARAMETER (IRMINE =  -38,
CIBMPCR*4     .           IRMAXE =   38,
CIBMPCR*4     .           IRMINM = 1176,
CIBMPCR*4     .           IRMAXM = 3401)

C                           * KONSTANTEN ZUR FESTLEGUNG DES WERTE-
C                           * BEREICHES FUER DOPPELT GENAUE REAL-ZAHLEN
      INTEGER   IDMINE, IDMAXE, IDMINM, IDMAXM
      PARAMETER (IDMINE =  -79,
     .           IDMAXE =   75,
     .           IDMINM = 5398,
     .           IDMAXM = 7236)
CIBMVSR*8      PARAMETER (IDMINE =  -79,
CIBMVSR*8     .           IDMAXE =   75,
CIBMVSR*8     .           IDMINM = 5398,
CIBMVSR*8     .           IDMAXM = 7236)
CIBMPCR*8      PARAMETER (IDMINE = -308,
CIBMPCR*8     .           IDMAXE =  308,
CIBMPCR*8     .           IDMINM = 2226,
CIBMPCR*8     .           IDMAXM = 1796)

C                           * SONSTIGE KONSTANTEN
      INTEGER   S01, S02, S03, S04, S05, S06, S07, S08, S09, S10,
     .          S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,
     .          S21, S22, S23, S24, S25
      PARAMETER (S01 =  1, S02 =  2, S03 =  3, S04 =  4, S05 =  5,
     .           S06 =  6, S07 =  7, S08 =  8, S09 =  9, S10 = 10,
     .           S11 = 11, S12 = 12, S13 = 13, S14 = 14, S15 = 15,
     .           S16 = 16, S17 = 17, S18 = 18, S19 = 19, S20 = 20,
     .           S21 = 21, S22 = 22, S23 = 23, S24 = 24, S25 = 25)
      INTEGER   S05M, S08M, S09M, S10M, S11M, S12M, S15M, S17M, S18M,
     .          S20M, S25M
      PARAMETER (S05M = -S05, S08M = -S08, S09M = -S09, S10M = -S10,
     .           S11M = -S11, S12M = -S12, S15M = -S15, S17M = -S17,
     .           S18M = -S18, S20M = -S20, S25M = -S25)
      INTEGER   SNIL
      PARAMETER (SNIL = -32766)

C                           * PROGRAMMVARIABLEN:
      INTEGER I, ICUR, IEANZ, IELEMT, ISCNTB(24,0:8), ISTATE, ISYMBL,
     .        ITXTLE, J

      SAVE ISCNTB

      LOGICAL LVALID

C                           * INITIALISIERUNGEN:
C TABELLE ZUR SYNTAXDEFINITION:
C
C SYMBOL       S0  " "   "-"   "+"   "."   "D"   "E"   "'"   ZIFF
C                                          "d"   "e"
C ---------+--------------------------------------------------------
C ANFANGS- |
C ZUSTAENDE|
C C  S01   | S15M  S01  S15M  S15M  S15M  S15M  S15M  S05M   S15M
C IS S02   |   -   S02  S08M  S08M    -     -     -     -    S17M
C R  S03   |   -   S03  S09M  S09M  S10M    -     -     -    S18M
C D  S04   |   -   S04  S11M  S11M  S12M    -     -     -    S20M
C FOLGE-   |
C ZUSTAENDE|
C C  S05   |  S07  S07   S07   S07   S07   S07   S07   S06    S07
C C  S06   |   -    -     -     -     -     -     -    S07     -
C C  S07   |  S07  S07   S07   S07   S07   S07   S07   S16    S07
C IS S08   |   -    -     -     -     -     -     -     -     S17
C R  S09   |   -    -     -     -    S10    -     -     -     S18
C R  S10   |   -    -     -     -     -     -     -     -     S19
C D  S11   |   -    -     -     -    S12    -     -     -     S20
C D  S12   |   -    -     -     -     -     -     -     -     S21
C RD S13   |   -    -    S14   S14    -     -     -     -     S22
C RD S14   |   -    -     -     -     -     -     -     -     S22
C END-     |
C ZUSTAENDE|
C C  S15   |  S15 S25M   S15   S15   S15   S15   S15   S15    S15
C C  S16   |   -  S25M    -     -     -     -     -    S07     -
C IS S17   |   -  S25M    -     -     -     -     -     -     S17
C R  S18   |   -  S25M    -     -    S19    -    S13    -     S18
C R  S19   |   -  S25M    -     -     -     -    S13    -     S19
C D  S20   |   -  S25M    -     -    S21   S13   S13    -     S20
C D  S21   |   -  S25M    -     -     -    S13   S13    -     S21
C RD S22   |   -  S25M    -     -     -     -     -     -     S23
C RD S23   |   -  S25M    -     -     -     -     -     -     S24
C RD S24   |   -  S25M    -     -     -     -     -     -      -
C
C S0 (SYMBOLKLASSE 0): ENTHAELT ALLE SYMBOLE, DIE NICHT GESONDERT AUFGE-
C                      FUEHRT SIND
C ZIFF               : ENTHAELT ALLE ZIFFERN VON "0" - "9"
C
      DATA ((ISCNTB(I,J), J = 0,8), I = 1,19)
     1   / S15M,  S01, S15M, S15M, S15M, S15M, S15M, S05M, S15M,
     2     SNIL,  S02, S08M, S08M, SNIL, SNIL, SNIL, SNIL, S17M,
     3     SNIL,  S03, S09M, S09M, S10M, SNIL, SNIL, SNIL, S18M,
     4     SNIL,  S04, S11M, S11M, S12M, SNIL, SNIL, SNIL, S20M,
     5      S07,  S07,  S07,  S07,  S07,  S07,  S07,  S06,  S07,
     6     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S07, SNIL,
     7      S07,  S07,  S07,  S07,  S07,  S07,  S07,  S16,  S07,
     8     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S17,
     9     SNIL, SNIL, SNIL, SNIL,  S10, SNIL, SNIL, SNIL,  S18,
     +     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S19,
     1     SNIL, SNIL, SNIL, SNIL,  S12, SNIL, SNIL, SNIL,  S20,
     2     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S21,
     3     SNIL, SNIL,  S14,  S14, SNIL, SNIL, SNIL, SNIL,  S22,
     4     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S22,
     5      S15, S25M,  S15,  S15,  S15,  S15,  S15,  S15,  S15,
     6     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL,  S07, SNIL,
     7     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S17,
     8     SNIL, S25M, SNIL, SNIL,  S19, SNIL,  S13, SNIL,  S18,
     9     SNIL, S25M, SNIL, SNIL, SNIL, SNIL,  S13, SNIL,  S19 /
      DATA ((ISCNTB(I,J), J = 0,8), I = 20,24)
     +   / SNIL, S25M, SNIL, SNIL,  S21,  S13,  S13, SNIL,  S20,
     1     SNIL, S25M, SNIL, SNIL, SNIL,  S13,  S13, SNIL,  S21,
     2     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S23,
     3     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S24,
     4     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL /
C                           * ------------------------------------------
      ICUR = 0
      ITXTLE = LEN(CTEXT)
      IEANZ = LEN(CELEMT)
      DO 1888, I = 1, IEANZ
         IELEMT = INDEX('CIRDS',CELEMT(I:I))
         IF ( IELEMT .EQ. 0 ) IELEMT = INDEX('cirds',CELEMT(I:I))
         IF ( IELEMT .EQ. 5 ) THEN
            ISTATE = 2
         ELSE
            ISTATE = IELEMT
         ENDIF
         IF ( ISTATE .EQ. 0 ) THEN
            IETAB(1,I) = -3
            IETAB(2,I) = 0
            WRITE(*,*) 'SCAN (W): Invalid descriptor "',CELEMT(I:I),
     .         '" ignored.'
            GOTO 1888
         ENDIF
 1100       CONTINUE
            ICUR = ICUR+1
            IF ( ICUR .GT. ITXTLE ) THEN
               IF ( ISTATE .GE. 15 ) THEN
C                           * ISTATE IST EIN ENDZUSTAND
                  ISTATE = S25M
               ELSE
C                           * ISTATE IST KEIN ENDZUSTAND
                  ISTATE = SNIL
               ENDIF
            ELSE
C                           * ERMITTLE INDEX DES AKTUELLEN SYMBOLS
               ISYMBL = INDEX(' -+.DEde''0123456789',CTEXT(ICUR:ICUR))
               IF ( ISYMBL .GE. 7 ) THEN
                  ISYMBL = ISYMBL - 2
                  ISYMBL = MIN(ISYMBL,8)
               ENDIF
               ISTATE = ISCNTB(ISTATE,ISYMBL)
            ENDIF
            IF ( ISTATE .LT. 0 ) THEN
               IF ( ISTATE .EQ. S25M ) THEN
                  IETAB(2,I) = ICUR-1
                  LVALID = .TRUE.
                  IF ( IELEMT .EQ. 2 ) THEN
                     CALL SCANIR(CIMINA,CIMAXA,
     .                  CTEXT(IETAB(1,I):IETAB(2,I)),LVALID)
                  ELSE IF ( IELEMT .EQ. 3 ) THEN
                     CALL SCANRR(IRMINM,IRMINE,IRMAXM,IRMAXE,
     .                  CTEXT(IETAB(1,I):IETAB(2,I)),LVALID)
                  ELSE IF ( IELEMT .EQ. 4 ) THEN
                     CALL SCANRR(IDMINM,IDMINE,IDMAXM,IDMAXE,
     .                  CTEXT(IETAB(1,I):IETAB(2,I)),LVALID)
                  ELSE IF ( IELEMT .EQ. 5 ) THEN
                     CALL SCANIR(CSMINA,CSMAXA,
     .                  CTEXT(IETAB(1,I):IETAB(2,I)),LVALID)
                  ENDIF
                  IF ( .NOT. LVALID ) IETAB(1,I) = -2
               ELSE IF ( ISTATE .NE. SNIL ) THEN
                  IETAB(1,I) = ICUR
                  ISTATE = -ISTATE
               ELSE
                  IETAB(1,I) = -1
                  IETAB(2,I) = ICUR
 1110                CONTINUE
                     IF ( ICUR .GT. ITXTLE ) GOTO 1119
                     IF ( CTEXT(ICUR:ICUR) .EQ. ' ' ) GOTO 1119
                     ICUR = ICUR+1
                     GOTO 1110
 1119             CONTINUE
               ENDIF
            ENDIF
            IF ( ICUR .GT. ITXTLE ) GOTO 1199
            IF ( ISTATE .GE. 0 ) GOTO 1100
 1199    CONTINUE
 1888 CONTINUE
      RETURN
      END
C@PROCESS OPT(2) NOSDUMP NOGOSTMT
C ----------------------------------------------------------------------
C  SCANIR - UEBERPRUEFEN DES WERTES EINER INTEGER-ZAHL AUF EINER
C           ZEICHENKETTE
C           AUFRUF:  IM UNTERPROGRAMM "SCAN"
C  AUTOR: G. EGERER, DATUM: 10.12.87, LETZTE AENDERUNG: 13. 1.88
C ----------------------------------------------------------------------
      SUBROUTINE SCANIR(CIMINA,CIMAXA,CINT,LVALID)

C                           * EINGABEPARAMETER:
C                           * CIMINA - ZEICHENKETTE MIT DEM INTEGER-WERT
C                           *          FUER DIE UNTERE GRENZE DES WERTE-
C                           *          BEREICHES
C                           *          (UNTERE GRENZE = -CIMINA)
C                           * CIMAXA - ZEICHENKETTE MIT DEM INTEGER-WERT
C                           *          FUER DIE OBERE GRENZE DES WERTE-
C                           *          BEREICHES
C                           *          (OBERE GRENZE = +CIMAXA)
C                           * CINT   - ZEICHENKETTE MIT DER ZU PRUEFEN-
C                           *          DEN INTEGER-ZAHL
      CHARACTER(len=*) CIMINA, CIMAXA, CINT

C                           * AUSGABEPARAMETER:
C                           * LVALID - HAT DEN WERT ".TRUE.", WENN DIE
C                           *          INTEGER-ZAHL AUF "CINT" EINEN
C                           *          WERT INNERHALB DES WERTE-
C                           *          BEREICHES HAT, SONST ".FALSE."
      LOGICAL LVALID

C                           * PROGRAMMVARIABLEN:
      INTEGER   CMAXLE
      PARAMETER (CMAXLE = 10)
      CHARACTER(len=CMAXLE) CMAXVL, CINPVL

      INTEGER I, IENDS, ILENT, ISTRTS, ISTRTT

C                           * ------------------------------------------
      IENDS = LEN(CINT)
      IF ( CINT(1:1) .EQ. '-' ) THEN
         CMAXVL = CIMINA
         ILENT = LEN(CIMINA)
         ISTRTS = 2
      ELSE
         CMAXVL = CIMAXA
         ILENT = LEN(CIMAXA)
         IF ( CINT(1:1) .EQ. '+' ) THEN
            ISTRTS = 2
         ELSE
            ISTRTS = 1
         ENDIF
      ENDIF
      IF ( ILENT .GT. CMAXLE ) THEN
         WRITE(*,*) 'Fehler in Routine: SCANIR'
         WRITE(*,*) 'Konstante "CMAXLE" zu klein.'
         WRITE(*,*) 'Erforderlicher Minimalwert fuer "CMAXLE":',ILENT
         WRITE(*,*) 'Wenden Sie sich bitte an die Programmberatung, '//
     .      'Tel. 6658'
         STOP 16
      ENDIF
  100    CONTINUE
         IF (  ( INDEX('123456789',CINT(ISTRTS:ISTRTS)) .GT. 0 ) .OR.
     .      ( ISTRTS .EQ. IENDS )  ) GOTO 199
         ISTRTS = ISTRTS + 1
         GOTO 100
  199 CONTINUE
      ISTRTT = ILENT-(IENDS-ISTRTS)
      IF ( ISTRTT .GE. 1 ) THEN
         CINPVL(ISTRTT:ILENT) = CINT(ISTRTS:IENDS)
         DO 288, I = 1, ISTRTT-1
            CINPVL(I:I) = '0'
  288    CONTINUE
         I = 1
  300       CONTINUE
            IF ( ( CINPVL(I:I) .NE. CMAXVL(I:I) ) .OR.
     .         ( I .EQ. ILENT )  ) GOTO 399
            I = I + 1
            GOTO 300
  399    CONTINUE
         LVALID = ( CINPVL(I:I) .LE. CMAXVL(I:I) )
      ELSE
         LVALID = .FALSE.
      ENDIF
      RETURN
      END
C@PROCESS OPT(2) NOSDUMP NOGOSTMT
C ----------------------------------------------------------------------
C  SCANRR - UEBERPRUEFEN DES WERTES EINER REAL-ZAHL AUF EINER ZEICHEN-
C           KETTE
C           AUFRUF:  IM UNTERPROGRAMM "SCAN"
C  AUTOR: G. EGERER, DATUM: 15.12.87, LETZTE AENDERUNG: 18. 1.88
C ----------------------------------------------------------------------
      SUBROUTINE SCANRR(IMIMAN,IMIEXP,IMAMAN,IMAEXP,CREAL,LVALID)

C                           * EINGABEPARAMETER:
C                           * IMIMAN - ZAHL ZWISCHEN 1000 UND 9999, DIE
C                           *          ZUR BESTIMMUNG DER KLEINSTEN
C                           *          REAL-ZAHL GROESSER 0 DIENT (S.U.)
C                           * IMIEXP - EXPONENT DER KLEINSTEN REAL-ZAHL
C                           *          GROESSER 0
C                           *
C                           *          KLEINSTE REAL-ZAHL GROSSER 0 =
C                           *          IMIMAN/1000 * 10**IMIEXP
C                           *
C                           * IMAMAN - ZAHL ZWISCHEN 1000 UND 9999, DIE
C                           *          ZUR BESTIMMUNG DER GROESSTEN
C                           *          REAL-ZAHL DIENT (S.U.)
C                           * IMAEXP - EXPONENT DER GROESSTEN REAL-ZAHL
C                           *
C                           *          GROESSTE REAL-ZAHL =
C                           *          (IMAMAN+1)/1000 * 10**IMAEXP
C                           *
C                           * CREAL  - ZEICHENKETTE MIT DER ZU PRUEFEN-
C                           *          DEN REAL-ZAHL
      INTEGER       IMIMAN, IMIEXP, IMAMAN, IMAEXP
      CHARACTER(len=*) CREAL

C                           * AUSGABEPARAMETER:
C                           * LVALID - HAT DEN WERT ".TRUE.", WENN DIE
C                           *          REAL-ZAHL AUF "CREAL" EINEN WERT
C                           *          INNERHALB DES WERTEBEREICHES HAT,
C                           *          SONST ".FALSE."
      LOGICAL LVALID

C                           * PROGRAMMVARIABLEN:
      INTEGER I, ICNT, IEND, IEXP, IMANTI, IPOS
      LOGICAL LFOUND
      CHARACTER(len=4) C1234

      SAVE C1234

C                           * INITIALISIERUNGEN:
      DATA C1234 /'1234'/
C                           * ------------------------------------------
      IEND = LEN(CREAL)
      IPOS = MAX(1,IEND-4)
      LFOUND = .FALSE.
  100    CONTINUE
         IF (  ( LFOUND ) .OR. ( IPOS .GT. IEND )  ) GOTO 199
         IF ( INDEX('EeDd',CREAL(IPOS:IPOS)) .GT. 0 ) THEN
            LFOUND = .TRUE.
         ELSE
            IPOS = IPOS + 1
         ENDIF
         GOTO 100
  199 CONTINUE
      IF ( LFOUND ) THEN
         I = IEND - IPOS
          READ(CREAL(IPOS+1:IEND),'(I'//C1234(I:I)//')') IEXP
      ELSE
         IEXP = 0
      ENDIF

      IEND = IPOS - 1
      IPOS = 1
      LFOUND = .FALSE.
  200    CONTINUE
         IF (  ( LFOUND ) .OR. ( IPOS .GT. IEND )  ) GOTO 299
         IF ( INDEX('123456789',CREAL(IPOS:IPOS)) .GT. 0 ) THEN
            LFOUND = .TRUE.
         ELSE
            IPOS = IPOS + 1
         ENDIF
         GOTO 200
  299 CONTINUE
      IF ( LFOUND ) THEN
         I = INDEX(CREAL(1:IEND),'.')
         IF ( I .EQ. 0 ) I = IEND + 1
         IF ( I .GT. IPOS ) THEN
            IEXP = IEXP + I - IPOS - 1
         ELSE
            IEXP = IEXP + I - IPOS
         ENDIF
         IF (  ( IEXP .LT. IMAEXP ) .AND. ( IEXP .GT. IMIEXP )  ) THEN
            LVALID = .TRUE.
         ELSE IF ( ( IEXP .GT. IMAEXP ) .OR. ( IEXP .LT. IMIEXP ) ) THEN
            LVALID = .FALSE.
         ELSE
            IMANTI = 0
            ICNT = 0
  300          CONTINUE
               I = INDEX('0123456789',CREAL(IPOS:IPOS))
               IF ( I .GT. 0 ) THEN
                  IMANTI = IMANTI*10 + (I-1)
                  ICNT = ICNT + 1
               ENDIF
               IPOS = IPOS + 1
               IF (  ( IPOS .GT. IEND ) .OR. ( ICNT .EQ. 4 )  ) GOTO 399
               GOTO 300
  399       CONTINUE
            IMANTI = IMANTI * 10**(4-ICNT)
            IF ( IEXP .EQ. IMAEXP ) THEN
               LVALID = ( IMANTI .LE. IMAMAN )
            ELSE
               LVALID = ( IMANTI .GE. IMIMAN )
            ENDIF
         ENDIF
      ELSE
         LVALID = .TRUE.
      ENDIF
      RETURN
      END
