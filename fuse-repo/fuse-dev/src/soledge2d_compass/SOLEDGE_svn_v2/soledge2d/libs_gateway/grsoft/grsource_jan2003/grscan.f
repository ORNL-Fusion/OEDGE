C@PROCESS OPT(2) NOSDUMP NOGOSTMT                                        SCA0001
C ----------------------------------------------------------------------SCA00020
C  SCAN - SCANNEN DES TEXTES "CTEXT"                                    SCA00030
C  AUTOR: G. EGERER, DATUM: 26. 9.85, LETZTE AENDERUNG: 17. 1.88        SCA00040
C  Update : Busch 9.12.91 (Version von 17.1.88 als Source und umbenannt
C                          als GRSCAN)
C ----------------------------------------------------------------------SCA00050
      SUBROUTINE GRSCAN(CTEXT,CELEMT,IETAB)                             SCA00060
                                                                        SCA00070
C                           * EINGABEPARAMETER:                         SCA00080
C                           * CTEXT  - TEXTZEILE                        SCA00090
C                           * CELEMT - FORMALE BESCHREIBUNG DER IN DER  SCA00100
C                           *          TEXTZEILE ERLAUBTEN ELEMENTE     SCA00110
C                           *          DURCH EINE ZEICHENKETTE. DIESE   SCA00120
C                           *          DARF NUR DIE FOLGENDEN ZEICHEN   SCA00130
C                           *          ENTHALTEN:                       SCA00140
C                           *                     | STEHT FUER          SCA00150
C                           *             ZEICHEN | DAS ELEMENT         SCA00160
C                           *            ---------+-------------------- SCA00170
C                           *                C    | TEXT ODER IN ' EIN- SCA00180
C                           *                     | GESCHLOSSENE        SCA00190
C                           *                     | ZEICHENKETTE        SCA00200
C                           *                I    | INTEGER-ZAHL MIT    SCA00210
C                           *                     | GROSSEM WERTE-      SCA00220
C                           *                     | BEREICH (INTEGER*4) SCA00230
C                           *                S    | INTEGER-ZAHL MIT    SCA00240
C                           *                     | KLEINEM WERTE-      SCA00250
C                           *                     | BEREICH (INTEGER*2) SCA00260
C                           *                R    | EINFACH GENAUE      SCA00270
C                           *                     | REAL-ZAHL           SCA00280
C                           *                D    | DOPPELT GENAUE      SCA00290
C                           *                     | REAL-ZAHL           SCA00300
C                           *          DIE ELEMENTE IN DER TEXTZEILE    SCA00310
C                           *          MUESSEN DURCH LEERZEICHEN        SCA00320
C                           *          VONEINANDER GETRENNT SEIN.       SCA00330
                                                                        SCA00340
      CHARACTER(len=*) CTEXT, CELEMT                                       SCA00350
                                                                        SCA00360
C                           * AUSGABEPARAMETER:                         SCA00370
C                           * IETAB  - TABELLE, DIE JEWEILS DIE ANFANGS-SCA00380
C                           *          UND DIE ENDPOSITION IN DER TEXT- SCA00390
C                           *          ZEILE VON ALLEN GUELTIGEN ELE-   SCA00400
C                           *          MENTEN ENTHAELT.                 SCA00410
C                           *          FUER UNGUELTIGE ELEMENTE ENT-    SCA00420
C                           *          HAELT DIESE TABELLE DIE FOLGENDE SCA00430
C                           *          INFORMATION:                     SCA00440
C                           *                      | ANF.- | END-       SCA00450
C                           *           FEHLER     | POS.  | POSITION   SCA00460
C                           *          ------------+-------+----------- SCA00470
C                           *           SYNTAX-    |   -1  | POSITION   SCA00480
C                           *           FEHLER     |       | DES 1. UN- SCA00490
C                           *                      |       | GUELTIGEN  SCA00500
C                           *                      |       | ZEICHENS   SCA00510
C                           *          ------------+-------+----------- SCA00520
C                           *           WERT       |   -2  | KORREKTE   SCA00530
C                           *           AUSSERHALB |       | ENDPOSI-   SCA00540
C                           *           DES WERTE- |       | TION DER   SCA00550
C                           *           BEREICHES  |       | ZAHL       SCA00560
C                           *          ------------+-------+----------- SCA00570
C                           *           NICHT ZU-  |   -3  |     0      SCA00580
C                           *           LAESSIGE   |       |            SCA00590
C                           *           FORMALE    |       |            SCA00600
C                           *           ELEMENTBE- |       |            SCA00610
C                           *           SCHREIBUNG |       |            SCA00620
      INTEGER IETAB(2,*)                                                SCA00630
                                                                        SCA00640
C                           * BENAMTE PROGRAMMKONSTANTEN:               SCA00650
                                                                        SCA00660
C                           * KONSTANTEN ZUR FESTLEGUNG DES WERTE-      SCA00670
C                           * BEREICHES FUER STANDARD-INTEGER-ZAHLEN    SCA00680
C                           * (INTEGER*4)                               SCA00690
      INTEGER            IMAXLI                                         SCA00700
      PARAMETER          (IMAXLI = 10)                                  SCA00710
CIBMVSI*4      PARAMETER          (IMAXLI = 10)                         SCA00720
CIBMPCI*4      PARAMETER          (IMAXLI = 10)                         SCA00730
      CHARACTER(len=IMAXLI) CIMINA, CIMAXA                                 SCA00740
      PARAMETER          (CIMINA = '2147483647',                        SCA00750
     .                    CIMAXA = '2147483647')                        SCA00760
CIBMVSI*4      PARAMETER          (CIMINA = '2147483647',               SCA00770
CIBMVSI*4     .                    CIMAXA = '2147483647')               SCA00780
CIBMPCI*4      PARAMETER          (CIMINA = '2147483648',               SCA00790
CIBMPCI*4     .                    CIMAXA = '2147483647')               SCA00800
                                                                        SCA00810
C                           * KONSTANTEN ZUR FESTLEGUNG DES WERTE-      SCA00820
C                           * BEREICHES FUER SHORT-INTEGER-ZAHLEN       SCA00830
C                           * (INTEGER*2)                               SCA00840
      INTEGER            IMAXLS                                         SCA00850
      PARAMETER          (IMAXLS = 5)                                   SCA00860
CIBMVSI*2      PARAMETER          (IMAXLS = 5)                          SCA00870
CIBMPCI*2      PARAMETER          (IMAXLS = 5)                          SCA00880
      CHARACTER(len=IMAXLS) CSMINA, CSMAXA                                 SCA00890
      PARAMETER          (CSMINA = '32767',                             SCA00900
     .                    CSMAXA = '32767')                             SCA00910
CIBMVSI*2      PARAMETER          (CSMINA = '32767',                    SCA00920
CIBMVSI*2     .                    CSMAXA = '32767')                    SCA00930
CIBMPCI*2      PARAMETER          (CSMINA = '32768',                    SCA00940
CIBMPCI*2     .                    CSMAXA = '32767')                    SCA00950
                                                                        SCA00960
C                           * KONSTANTEN ZUR FESTLEGUNG DES WERTE-      SCA00970
C                           * BEREICHES FUER EINFACH GENAUE REAL-ZAHLEN SCA00980
      INTEGER   IRMINE, IRMAXE, IRMINM, IRMAXM                          SCA00990
      PARAMETER (IRMINE =  -79,                                         SCA01000
     .           IRMAXE =   75,                                         SCA01010
     .           IRMINM = 5398,                                         SCA01020
     .           IRMAXM = 7236)                                         SCA01030
CIBMVSR*4      PARAMETER (IRMINE =  -79,                                SCA01040
CIBMVSR*4     .           IRMAXE =   75,                                SCA01050
CIBMVSR*4     .           IRMINM = 5398,                                SCA01060
CIBMVSR*4     .           IRMAXM = 7236)                                SCA01070
CIBMPCR*4      PARAMETER (IRMINE =  -38,                                SCA01080
CIBMPCR*4     .           IRMAXE =   38,                                SCA01090
CIBMPCR*4     .           IRMINM = 1176,                                SCA01100
CIBMPCR*4     .           IRMAXM = 3401)                                SCA01110
                                                                        SCA01120
C                           * KONSTANTEN ZUR FESTLEGUNG DES WERTE-      SCA01130
C                           * BEREICHES FUER DOPPELT GENAUE REAL-ZAHLEN SCA01140
      INTEGER   IDMINE, IDMAXE, IDMINM, IDMAXM                          SCA01150
      PARAMETER (IDMINE =  -79,                                         SCA01160
     .           IDMAXE =   75,                                         SCA01170
     .           IDMINM = 5398,                                         SCA01180
     .           IDMAXM = 7236)                                         SCA01190
CIBMVSR*8      PARAMETER (IDMINE =  -79,                                SCA01200
CIBMVSR*8     .           IDMAXE =   75,                                SCA01210
CIBMVSR*8     .           IDMINM = 5398,                                SCA01220
CIBMVSR*8     .           IDMAXM = 7236)                                SCA01230
CIBMPCR*8      PARAMETER (IDMINE = -308,                                SCA01240
CIBMPCR*8     .           IDMAXE =  308,                                SCA01250
CIBMPCR*8     .           IDMINM = 2226,                                SCA01260
CIBMPCR*8     .           IDMAXM = 1796)                                SCA01270
                                                                        SCA01280
C                           * SONSTIGE KONSTANTEN                       SCA01290
      INTEGER   S01, S02, S03, S04, S05, S06, S07, S08, S09, S10,       SCA01300
     .          S11, S12, S13, S14, S15, S16, S17, S18, S19, S20,       SCA01310
     .          S21, S22, S23, S24, S25                                 SCA01320
      PARAMETER (S01 =  1, S02 =  2, S03 =  3, S04 =  4, S05 =  5,      SCA01330
     .           S06 =  6, S07 =  7, S08 =  8, S09 =  9, S10 = 10,      SCA01340
     .           S11 = 11, S12 = 12, S13 = 13, S14 = 14, S15 = 15,      SCA01350
     .           S16 = 16, S17 = 17, S18 = 18, S19 = 19, S20 = 20,      SCA01360
     .           S21 = 21, S22 = 22, S23 = 23, S24 = 24, S25 = 25)      SCA01370
      INTEGER   S05M, S08M, S09M, S10M, S11M, S12M, S15M, S17M, S18M,   SCA01380
     .          S20M, S25M                                              SCA01390
      PARAMETER (S05M = -S05, S08M = -S08, S09M = -S09, S10M = -S10,    SCA01400
     .           S11M = -S11, S12M = -S12, S15M = -S15, S17M = -S17,    SCA01410
     .           S18M = -S18, S20M = -S20, S25M = -S25)                 SCA01420
      INTEGER   SNIL                                                    SCA01430
      PARAMETER (SNIL = -32766)                                         SCA01440
                                                                        SCA01450
C                           * PROGRAMMVARIABLEN:                        SCA01460
      INTEGER I, ICUR, IEANZ, IELEMT, ISCNTB(24,0:8), ISTATE, ISYMBL,   SCA01470
     .        ITXTLE, J                                                 SCA01480
                                                                        SCA01490
      SAVE ISCNTB                                                       SCA01500
                                                                        SCA01510
      LOGICAL LVALID                                                    SCA01520
                                                                        SCA01530
C                           * INITIALISIERUNGEN:                        SCA01540
C TABELLE ZUR SYNTAXDEFINITION:                                         SCA01550
C                                                                       SCA01560
C SYMBOL       S0  " "   "-"   "+"   "."   "D"   "E"   "'"   ZIFF       SCA01570
C                                          "d"   "e"                    SCA01580
C ---------+--------------------------------------------------------    SCA01590
C ANFANGS- |                                                            SCA01600
C ZUSTAENDE|                                                            SCA01610
C C  S01   | S15M  S01  S15M  S15M  S15M  S15M  S15M  S05M   S15M       SCA01620
C IS S02   |   -   S02  S08M  S08M    -     -     -     -    S17M       SCA01630
C R  S03   |   -   S03  S09M  S09M  S10M    -     -     -    S18M       SCA01640
C D  S04   |   -   S04  S11M  S11M  S12M    -     -     -    S20M       SCA01650
C FOLGE-   |                                                            SCA01660
C ZUSTAENDE|                                                            SCA01670
C C  S05   |  S07  S07   S07   S07   S07   S07   S07   S06    S07       SCA01680
C C  S06   |   -    -     -     -     -     -     -    S07     -        SCA01690
C C  S07   |  S07  S07   S07   S07   S07   S07   S07   S16    S07       SCA01700
C IS S08   |   -    -     -     -     -     -     -     -     S17       SCA01710
C R  S09   |   -    -     -     -    S10    -     -     -     S18       SCA01720
C R  S10   |   -    -     -     -     -     -     -     -     S19       SCA01730
C D  S11   |   -    -     -     -    S12    -     -     -     S20       SCA01740
C D  S12   |   -    -     -     -     -     -     -     -     S21       SCA01750
C RD S13   |   -    -    S14   S14    -     -     -     -     S22       SCA01760
C RD S14   |   -    -     -     -     -     -     -     -     S22       SCA01770
C END-     |                                                            SCA01780
C ZUSTAENDE|                                                            SCA01790
C C  S15   |  S15 S25M   S15   S15   S15   S15   S15   S15    S15       SCA01800
C C  S16   |   -  S25M    -     -     -     -     -    S07     -        SCA01810
C IS S17   |   -  S25M    -     -     -     -     -     -     S17       SCA01820
C R  S18   |   -  S25M    -     -    S19    -    S13    -     S18       SCA01830
C R  S19   |   -  S25M    -     -     -     -    S13    -     S19       SCA01840
C D  S20   |   -  S25M    -     -    S21   S13   S13    -     S20       SCA01850
C D  S21   |   -  S25M    -     -     -    S13   S13    -     S21       SCA01860
C RD S22   |   -  S25M    -     -     -     -     -     -     S23       SCA01870
C RD S23   |   -  S25M    -     -     -     -     -     -     S24       SCA01880
C RD S24   |   -  S25M    -     -     -     -     -     -      -        SCA01890
C                                                                       SCA01900
C S0 (SYMBOLKLASSE 0): ENTHAELT ALLE SYMBOLE, DIE NICHT GESONDERT AUFGE-SCA01910
C                      FUEHRT SIND                                      SCA01920
C ZIFF               : ENTHAELT ALLE ZIFFERN VON "0" - "9"              SCA01930
C                                                                       SCA01940
      DATA ((ISCNTB(I,J), J = 0,8), I = 1,19)                           SCA01950
     1   / S15M,  S01, S15M, S15M, S15M, S15M, S15M, S05M, S15M,        SCA01960
     2     SNIL,  S02, S08M, S08M, SNIL, SNIL, SNIL, SNIL, S17M,        SCA01970
     3     SNIL,  S03, S09M, S09M, S10M, SNIL, SNIL, SNIL, S18M,        SCA01980
     4     SNIL,  S04, S11M, S11M, S12M, SNIL, SNIL, SNIL, S20M,        SCA01990
     5      S07,  S07,  S07,  S07,  S07,  S07,  S07,  S06,  S07,        SCA02000
     6     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S07, SNIL,        SCA02010
     7      S07,  S07,  S07,  S07,  S07,  S07,  S07,  S16,  S07,        SCA02020
     8     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S17,        SCA02030
     9     SNIL, SNIL, SNIL, SNIL,  S10, SNIL, SNIL, SNIL,  S18,        SCA02040
     +     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S19,        SCA02050
     1     SNIL, SNIL, SNIL, SNIL,  S12, SNIL, SNIL, SNIL,  S20,        SCA02060
     2     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S21,        SCA02070
     3     SNIL, SNIL,  S14,  S14, SNIL, SNIL, SNIL, SNIL,  S22,        SCA02080
     4     SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S22,        SCA02090
     5      S15, S25M,  S15,  S15,  S15,  S15,  S15,  S15,  S15,        SCA02100
     6     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL,  S07, SNIL,        SCA02110
     7     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S17,        SCA02120
     8     SNIL, S25M, SNIL, SNIL,  S19, SNIL,  S13, SNIL,  S18,        SCA02130
     9     SNIL, S25M, SNIL, SNIL, SNIL, SNIL,  S13, SNIL,  S19 /       SCA02140
      DATA ((ISCNTB(I,J), J = 0,8), I = 20,24)                          SCA02150
     +   / SNIL, S25M, SNIL, SNIL,  S21,  S13,  S13, SNIL,  S20,        SCA02160
     1     SNIL, S25M, SNIL, SNIL, SNIL,  S13,  S13, SNIL,  S21,        SCA02170
     2     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S23,        SCA02180
     3     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL,  S24,        SCA02190
     4     SNIL, S25M, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL, SNIL /       SCA02200
C                           * ------------------------------------------SCA02210
      ICUR = 0                                                          SCA02220
      ITXTLE = LEN(CTEXT)                                               SCA02230
      IEANZ = LEN(CELEMT)                                               SCA02240
      DO 1888, I = 1, IEANZ                                             SCA02250
         IELEMT = INDEX('CIRDS',CELEMT(I:I))                            SCA02260
         IF ( IELEMT .EQ. 0 ) IELEMT = INDEX('cirds',CELEMT(I:I))       SCA02270
         IF ( IELEMT .EQ. 5 ) THEN                                      SCA02280
            ISTATE = 2                                                  SCA02290
         ELSE                                                           SCA02300
            ISTATE = IELEMT                                             SCA02310
         ENDIF                                                          SCA02320
         IF ( ISTATE .EQ. 0 ) THEN                                      SCA02330
            IETAB(1,I) = -3                                             SCA02340
            IETAB(2,I) = 0                                              SCA02350
            WRITE(*,*) 'SCAN (W): Invalid descriptor "',CELEMT(I:I),    SCA02360
     .         '" ignored.'                                             SCA02370
            GOTO 1888                                                   SCA02380
         ENDIF                                                          SCA02390
 1100       CONTINUE                                                    SCA02400
            ICUR = ICUR+1                                               SCA02410
            IF ( ICUR .GT. ITXTLE ) THEN                                SCA02420
               IF ( ISTATE .GE. 15 ) THEN                               SCA02430
C                           * ISTATE IST EIN ENDZUSTAND                 SCA02440
                  ISTATE = S25M                                         SCA02450
               ELSE                                                     SCA02460
C                           * ISTATE IST KEIN ENDZUSTAND                SCA02470
                  ISTATE = SNIL                                         SCA02480
               ENDIF                                                    SCA02490
            ELSE                                                        SCA02500
C                           * ERMITTLE INDEX DES AKTUELLEN SYMBOLS      SCA02510
               ISYMBL = INDEX(' -+.DEde''0123456789',CTEXT(ICUR:ICUR))  SCA02520
               IF ( ISYMBL .GE. 7 ) THEN                                SCA02530
                  ISYMBL = ISYMBL - 2                                   SCA02540
                  ISYMBL = MIN(ISYMBL,8)                                SCA02550
               ENDIF                                                    SCA02560
               ISTATE = ISCNTB(ISTATE,ISYMBL)                           SCA02570
            ENDIF                                                       SCA02580
            IF ( ISTATE .LT. 0 ) THEN                                   SCA02590
               IF ( ISTATE .EQ. S25M ) THEN                             SCA02600
                  IETAB(2,I) = ICUR-1                                   SCA02610
                  LVALID = .TRUE.                                       SCA02620
                  IF ( IELEMT .EQ. 2 ) THEN                             SCA02630
                     CALL GRSCIR(CIMINA,CIMAXA,                         SCA02640
     .                  CTEXT(IETAB(1,I):IETAB(2,I)),LVALID)            SCA02650
                  ELSE IF ( IELEMT .EQ. 3 ) THEN                        SCA02660
                     CALL GRSCRR(IRMINM,IRMINE,IRMAXM,IRMAXE,           SCA02670
     .                  CTEXT(IETAB(1,I):IETAB(2,I)),LVALID)            SCA02680
                  ELSE IF ( IELEMT .EQ. 4 ) THEN                        SCA02690
                     CALL GRSCRR(IDMINM,IDMINE,IDMAXM,IDMAXE,           SCA02700
     .                  CTEXT(IETAB(1,I):IETAB(2,I)),LVALID)            SCA02710
                  ELSE IF ( IELEMT .EQ. 5 ) THEN                        SCA02720
                     CALL GRSCIR(CSMINA,CSMAXA,                         SCA02730
     .                  CTEXT(IETAB(1,I):IETAB(2,I)),LVALID)            SCA02740
                  ENDIF                                                 SCA02750
                  IF ( .NOT. LVALID ) IETAB(1,I) = -2                   SCA02760
               ELSE IF ( ISTATE .NE. SNIL ) THEN                        SCA02770
                  IETAB(1,I) = ICUR                                     SCA02780
                  ISTATE = -ISTATE                                      SCA02790
               ELSE                                                     SCA02800
                  IETAB(1,I) = -1                                       SCA02810
                  IETAB(2,I) = ICUR                                     SCA02820
 1110                CONTINUE                                           SCA02830
                     IF ( ICUR .GT. ITXTLE ) GOTO 1119                  SCA02840
                     IF ( CTEXT(ICUR:ICUR) .EQ. ' ' ) GOTO 1119         SCA02850
                     ICUR = ICUR+1                                      SCA02860
                     GOTO 1110                                          SCA02870
 1119             CONTINUE                                              SCA02880
               ENDIF                                                    SCA02890
            ENDIF                                                       SCA02900
            IF ( ICUR .GT. ITXTLE ) GOTO 1199                           SCA02910
            IF ( ISTATE .GE. 0 ) GOTO 1100                              SCA02920
 1199    CONTINUE                                                       SCA02930
 1888 CONTINUE                                                          SCA02940
      RETURN                                                            SCA02950
      END                                                               SCA02960
C@PROCESS OPT(2) NOSDUMP NOGOSTMT                                        SCA0297
C ----------------------------------------------------------------------SCA02980
C  GRSCIR - UEBERPRUEFEN DES WERTES EINER INTEGER-ZAHL AUF EINER        SCA02990
C           ZEICHENKETTE                                                SCA03000
C           AUFRUF:  IM UNTERPROGRAMM "SCAN"                            SCA03010
C  AUTOR: G. EGERER, DATUM: 10.12.87, LETZTE AENDERUNG: 13. 1.88        SCA03020
C ----------------------------------------------------------------------SCA03030
      SUBROUTINE GRSCIR(CIMINA,CIMAXA,CINT,LVALID)                      SCA03040
                                                                        SCA03050
C                           * EINGABEPARAMETER:                         SCA03060
C                           * CIMINA - ZEICHENKETTE MIT DEM INTEGER-WERTSCA03070
C                           *          FUER DIE UNTERE GRENZE DES WERTE-SCA03080
C                           *          BEREICHES                        SCA03090
C                           *          (UNTERE GRENZE = -CIMINA)        SCA03100
C                           * CIMAXA - ZEICHENKETTE MIT DEM INTEGER-WERTSCA03110
C                           *          FUER DIE OBERE GRENZE DES WERTE- SCA03120
C                           *          BEREICHES                        SCA03130
C                           *          (OBERE GRENZE = +CIMAXA)         SCA03140
C                           * CINT   - ZEICHENKETTE MIT DER ZU PRUEFEN- SCA03150
C                           *          DEN INTEGER-ZAHL                 SCA03160
      CHARACTER(len=*) CIMINA, CIMAXA, CINT                                SCA03170
                                                                        SCA03180
C                           * AUSGABEPARAMETER:                         SCA03190
C                           * LVALID - HAT DEN WERT ".TRUE.", WENN DIE  SCA03200
C                           *          INTEGER-ZAHL AUF "CINT" EINEN    SCA03210
C                           *          WERT INNERHALB DES WERTE-        SCA03220
C                           *          BEREICHES HAT, SONST ".FALSE."   SCA03230
      LOGICAL LVALID                                                    SCA03240
                                                                        SCA03250
C                           * PROGRAMMVARIABLEN:                        SCA03260
      INTEGER   CMAXLE                                                  SCA03270
      PARAMETER (CMAXLE = 10)                                           SCA03280
      CHARACTER(len=CMAXLE) CMAXVL, CINPVL                                 SCA03290
                                                                        SCA03300
      INTEGER I, IENDS, ILENT, ISTRTS, ISTRTT                           SCA03310
                                                                        SCA03320
C                           * ------------------------------------------SCA03330
      IENDS = LEN(CINT)                                                 SCA03340
      IF ( CINT(1:1) .EQ. '-' ) THEN                                    SCA03350
         CMAXVL = CIMINA                                                SCA03360
         ILENT = LEN(CIMINA)                                            SCA03370
         ISTRTS = 2                                                     SCA03380
      ELSE                                                              SCA03390
         CMAXVL = CIMAXA                                                SCA03400
         ILENT = LEN(CIMAXA)                                            SCA03410
         IF ( CINT(1:1) .EQ. '+' ) THEN                                 SCA03420
            ISTRTS = 2                                                  SCA03430
         ELSE                                                           SCA03440
            ISTRTS = 1                                                  SCA03450
         ENDIF                                                          SCA03460
      ENDIF                                                             SCA03470
      IF ( ILENT .GT. CMAXLE ) THEN                                     SCA03480
         WRITE(*,*) 'Fehler in Routine: GRSCIR'                         SCA03490
         WRITE(*,*) 'Konstante "CMAXLE" zu klein.'                      SCA03500
         WRITE(*,*) 'Erforderlicher Minimalwert fuer "CMAXLE":',ILENT   SCA03510
         WRITE(*,*) 'Wenden Sie sich bitte an die Programmberatung, '// SCA03520
     .      'Tel. 6658'                                                 SCA03530
         STOP 16                                                        SCA03540
      ENDIF                                                             SCA03550
  100    CONTINUE                                                       SCA03560
         IF (  ( INDEX('123456789',CINT(ISTRTS:ISTRTS)) .GT. 0 ) .OR.   SCA03570
     .      ( ISTRTS .EQ. IENDS )  ) GOTO 199                           SCA03580
         ISTRTS = ISTRTS + 1                                            SCA03590
         GOTO 100                                                       SCA03600
  199 CONTINUE                                                          SCA03610
      ISTRTT = ILENT-(IENDS-ISTRTS)                                     SCA03620
      IF ( ISTRTT .GE. 1 ) THEN                                         SCA03630
         CINPVL(ISTRTT:ILENT) = CINT(ISTRTS:IENDS)                      SCA03640
         DO 288, I = 1, ISTRTT-1                                        SCA03650
            CINPVL(I:I) = '0'                                           SCA03660
  288    CONTINUE                                                       SCA03670
         I = 1                                                          SCA03680
  300       CONTINUE                                                    SCA03690
            IF ( ( CINPVL(I:I) .NE. CMAXVL(I:I) ) .OR.                  SCA03700
     .         ( I .EQ. ILENT )  ) GOTO 399                             SCA03710
            I = I + 1                                                   SCA03720
            GOTO 300                                                    SCA03730
  399    CONTINUE                                                       SCA03740
         LVALID = ( CINPVL(I:I) .LE. CMAXVL(I:I) )                      SCA03750
      ELSE                                                              SCA03760
         LVALID = .FALSE.                                               SCA03770
      ENDIF                                                             SCA03780
      RETURN                                                            SCA03790
      END                                                               SCA03800
C@PROCESS OPT(2) NOSDUMP NOGOSTMT                                        SCA0381
C ----------------------------------------------------------------------SCA03820
C  GRSCRR - UEBERPRUEFEN DES WERTES EINER REAL-ZAHL AUF EINER ZEICHEN-  SCA03830
C           KETTE                                                       SCA03840
C           AUFRUF:  IM UNTERPROGRAMM "SCAN"                            SCA03850
C  AUTOR: G. EGERER, DATUM: 15.12.87, LETZTE AENDERUNG: 18. 1.88        SCA03860
C ----------------------------------------------------------------------SCA03870
      SUBROUTINE GRSCRR(IMIMAN,IMIEXP,IMAMAN,IMAEXP,CREAL,LVALID)       SCA03880
                                                                        SCA03890
C                           * EINGABEPARAMETER:                         SCA03900
C                           * IMIMAN - ZAHL ZWISCHEN 1000 UND 9999, DIE SCA03910
C                           *          ZUR BESTIMMUNG DER KLEINSTEN     SCA03920
C                           *          REAL-ZAHL GROESSER 0 DIENT (S.U.)SCA03930
C                           * IMIEXP - EXPONENT DER KLEINSTEN REAL-ZAHL SCA03940
C                           *          GROESSER 0                       SCA03950
C                           *                                           SCA03960
C                           *          KLEINSTE REAL-ZAHL GROSSER 0 =   SCA03970
C                           *          IMIMAN/1000 * 10**IMIEXP         SCA03980
C                           *                                           SCA03990
C                           * IMAMAN - ZAHL ZWISCHEN 1000 UND 9999, DIE SCA04000
C                           *          ZUR BESTIMMUNG DER GROESSTEN     SCA04010
C                           *          REAL-ZAHL DIENT (S.U.)           SCA04020
C                           * IMAEXP - EXPONENT DER GROESSTEN REAL-ZAHL SCA04030
C                           *                                           SCA04040
C                           *          GROESSTE REAL-ZAHL =             SCA04050
C                           *          (IMAMAN+1)/1000 * 10**IMAEXP     SCA04060
C                           *                                           SCA04070
C                           * CREAL  - ZEICHENKETTE MIT DER ZU PRUEFEN- SCA04080
C                           *          DEN REAL-ZAHL                    SCA04090
      INTEGER       IMIMAN, IMIEXP, IMAMAN, IMAEXP                      SCA04100
      CHARACTER(len=*) CREAL                                               SCA04110
                                                                        SCA04120
C                           * AUSGABEPARAMETER:                         SCA04130
C                           * LVALID - HAT DEN WERT ".TRUE.", WENN DIE  SCA04140
C                           *          REAL-ZAHL AUF "CREAL" EINEN WERT SCA04150
C                           *          INNERHALB DES WERTEBEREICHES HAT,SCA04160
C                           *          SONST ".FALSE."                  SCA04170
      LOGICAL LVALID                                                    SCA04180
                                                                        SCA04190
C                           * PROGRAMMVARIABLEN:                        SCA04200
      INTEGER I, ICNT, IEND, IEXP, IMANTI, IPOS                         SCA04210
      LOGICAL LFOUND                                                    SCA04220
      CHARACTER (len=4) C1234                                                 SCA04230
                                                                        SCA04240
      SAVE C1234                                                        SCA04250
                                                                        SCA04260
C                           * INITIALISIERUNGEN:                        SCA04270
      DATA C1234 /'1234'/                                               SCA04280
C                           * ------------------------------------------SCA04290
      IEND = LEN(CREAL)                                                 SCA04300
      IPOS = MAX(1,IEND-4)                                              SCA04310
      LFOUND = .FALSE.                                                  SCA04320
  100    CONTINUE                                                       SCA04330
         IF (  ( LFOUND ) .OR. ( IPOS .GT. IEND )  ) GOTO 199           SCA04340
         IF ( INDEX('EeDd',CREAL(IPOS:IPOS)) .GT. 0 ) THEN              SCA04350
            LFOUND = .TRUE.                                             SCA04360
         ELSE                                                           SCA04370
            IPOS = IPOS + 1                                             SCA04380
         ENDIF                                                          SCA04390
         GOTO 100                                                       SCA04400
  199 CONTINUE                                                          SCA04410
      IF ( LFOUND ) THEN                                                SCA04420
         I = IEND - IPOS                                                SCA04430
          READ(CREAL(IPOS+1:IEND),'(I'//C1234(I:I)//')') IEXP           SCA04440
      ELSE                                                              SCA04450
         IEXP = 0                                                       SCA04460
      ENDIF                                                             SCA04470
                                                                        SCA04480
      IEND = IPOS - 1                                                   SCA04490
      IPOS = 1                                                          SCA04500
      LFOUND = .FALSE.                                                  SCA04510
  200    CONTINUE                                                       SCA04520
         IF (  ( LFOUND ) .OR. ( IPOS .GT. IEND )  ) GOTO 299           SCA04530
         IF ( INDEX('123456789',CREAL(IPOS:IPOS)) .GT. 0 ) THEN         SCA04540
            LFOUND = .TRUE.                                             SCA04550
         ELSE                                                           SCA04560
            IPOS = IPOS + 1                                             SCA04570
         ENDIF                                                          SCA04580
         GOTO 200                                                       SCA04590
  299 CONTINUE                                                          SCA04600
      IF ( LFOUND ) THEN                                                SCA04610
         I = INDEX(CREAL(1:IEND),'.')                                   SCA04620
         IF ( I .EQ. 0 ) I = IEND + 1                                   SCA04630
         IF ( I .GT. IPOS ) THEN                                        SCA04640
            IEXP = IEXP + I - IPOS - 1                                  SCA04650
         ELSE                                                           SCA04660
            IEXP = IEXP + I - IPOS                                      SCA04670
         ENDIF                                                          SCA04680
         IF (  ( IEXP .LT. IMAEXP ) .AND. ( IEXP .GT. IMIEXP )  ) THEN  SCA04690
            LVALID = .TRUE.                                             SCA04700
         ELSE IF ( ( IEXP .GT. IMAEXP ) .OR. ( IEXP .LT. IMIEXP ) ) THENSCA04710
            LVALID = .FALSE.                                            SCA04720
         ELSE                                                           SCA04730
            IMANTI = 0                                                  SCA04740
            ICNT = 0                                                    SCA04750
  300          CONTINUE                                                 SCA04760
               I = INDEX('0123456789',CREAL(IPOS:IPOS))                 SCA04770
               IF ( I .GT. 0 ) THEN                                     SCA04780
                  IMANTI = IMANTI*10 + (I-1)                            SCA04790
                  ICNT = ICNT + 1                                       SCA04800
               ENDIF                                                    SCA04810
               IPOS = IPOS + 1                                          SCA04820
               IF (  ( IPOS .GT. IEND ) .OR. ( ICNT .EQ. 4 )  ) GOTO 399SCA04830
               GOTO 300                                                 SCA04840
  399       CONTINUE                                                    SCA04850
            IMANTI = IMANTI * 10**(4-ICNT)                              SCA04860
            IF ( IEXP .EQ. IMAEXP ) THEN                                SCA04870
               LVALID = ( IMANTI .LE. IMAMAN )                          SCA04880
            ELSE                                                        SCA04890
               LVALID = ( IMANTI .GE. IMIMAN )                          SCA04900
            ENDIF                                                       SCA04910
         ENDIF                                                          SCA04920
      ELSE                                                              SCA04930
         LVALID = .TRUE.                                                SCA04940
      ENDIF                                                             SCA04950
      RETURN                                                            SCA04960
      END                                                               SCA04970
