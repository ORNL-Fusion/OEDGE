C***********************************************************************
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRFONT(IFONT)
C     GKS FONT
C     Beim GKS von der TU Berlin gibt es die Fonts 1-18,20
C     mit proportionaler gerader Schrift
C     Bei GTSGRAL GKS gbit es dei Fonts -1 bis -11,-13,-51
C     mit proportionaler gerader Schrift und zusaetzlich
C     alle diese Fonts aequidistant sowohl gerade als auch kursiv
C
C     Umsetzungstabelle
C
C     TUB Berlin                GTSGRAL
C----------------------------------------------------------------------
C     Font 1 Standard            -1
C     Font 2 Standard thick      -2
C     Font 3 Roman               -3
C     Font 4 Roman thick         -4
C     Font 5 Roman outline       nicht vorhanden, wie Roman thick -4
C     Font 6 Grieche             nur      >
C     Font 7 Grieche thick       ein      = -13
C            neu - KFA erstellt                    jetzt neuer Font:-14
C     Font 8 Grieche outline     Grieche  >
C            neu - KFA erstellt                    jetzt neuer Font:-14
C     Font 9 Italic              -5
C     Font10 Italic thick        -6
C     Font11 Italic outline      nicht vorhanden, wie Italic thick -6
C     Font12 Script              -7
C     Font13 Script thick        -8
C     Font14 kyrilisch????       nicht vorhanden ==> Fuellfont -51
C     Font15 Gothic German       -10
C     Font16 Gothic English      -9
C     Font17 Gothic Italian      -11
C     Font18 kyrilisch????       nicht vorhanden ==> Fuellfont -51
C     Font20 Sonderzeichen       nicht vorhanden ==> Fuellfont -51
C
C Update : 26.7.91. ICOLOR --> PP(16)
C Update: 6.4 92 seit Febr. 92 gibt es den Sonderzeichenfont -14
C                die Zeichen koennen nicht automatisch umgesetzt
C                werden, die die Zuordnung der Ansteuerung
C                unterschiedlich ist
C                NEU: Font 14 GREEK complex-diese Zeichen werden
C                automatisch umgesetzt
C Update: 12.10.92
c     falls kein GTSGRAL-GKS und Font<0 dann Standard - Font=1 waehlen
C     wegen Probleme beim XGKS (BUS Error bei SUN und Strichcode
C     aehnlicher Output auf AIX mit XGKS)
C Update fuer GLIGKS 11.8.93 Busch
C
      INTEGER IFONT ,IFO  ,ICOLOR, INTLIN
      INTEGER PRUEF , NEWFNT(20)
      CHARACTER GKSTYP*7
      LOGICAL GRALGKS, FLGROT
      INTEGER ERRIND, FONT, PREC
      COMMON /GRPP/ pp(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRGKS/ GKSTYP
CDEC$ PSECT /GRGKS/ NOSHR
      SAVE /GRPP/ ,/GRGKS/
      EQUIVALENCE (PP(13),INTLIN),(PP(16),ICOLOR),(PP(17),SIZMRK)
      EQUIVALENCE (PP(9),IF),(PP(18),FLGROT)
cc    DATA NEWFNT/-1,-2,-3,-4,-4,-13,-13,-13,-5,-6,-6,-7,-8,-51,-10,
      DATA NEWFNT/-1,-2,-3,-4,-4,-13,-14,-14,-5,-6,-6,-7,-8,-51,-10,
     <            -9,-11,-51,-1,-51/
C-----------------------------------------------------------------------
c
c     Unterscheidung der GKS Versionen
C     Es gibt keine Query Funktion fuer den Namen der Implemantation
C     Abhilfe: GTSGRAL hat die MO Workstaion 321, 322, 323,
C              und die MI Workstations 324,325,326 auf der CRAY
C              und die MO WK 321,322,323 im CMS
C              PRUEF=3 CMS
C              PRUEF=6 CRAY
C     Falls diese vorhanden sind wird angenommen,
C     dass das GKS von GTSGRAL ist


      IFO=IFONT
      PRUEF =0
c     query Anzahl WK's
      CALL GQEWK(1,IERR,NUMBER,IWKT)
      DO 47 I=1,NUMBER
      CALL GQEWK(I,IERR,NUMBER,IWKT)
      IF (IWKT.EQ.321.OR.IWKT.EQ.322.OR.IWKT.EQ.323.OR.
     <    IWKT.EQ.324.OR.IWKT.EQ.325.OR.IWKT.EQ.326) PRUEF=PRUEF+1
47    CONTINUE

      IF ( PRUEF.EQ.6 .OR. PRUEF.EQ.3 )   GRALGKS=.TRUE.

C 11.8.93 Busch (wie GRALGKS behandeln , Fonts werden intern gemapped
      IF (GKSTYP.EQ.'GLIGKS' ) GRALGKS=.TRUE.

      IF ( GRALGKS .AND.IFONT.GT.0.and.IFONT.LE.20 ) IFO=NEWFNT(IFONT)

c     falls kein GTSGRAL-GKS und Font<0 dann Standard - Font=1 waehlen
C     wegen Probleme beim XGKS (BUS ERror bei SUN)
      IF ( .NOT.GRALGKS.AND.IFONT.LT.0) IFO=1
      IF ( GRALGKS.AND.IFONT.LT.0) IFO=IFONT

      CALL GQTXFP(ERRIND, FONT, PREC)
      CALL GSTXFP(IFO, PREC)
c Frabe --> COMMON GRPP
      IF = IFO
      RETURN
      END
C
      SUBROUTINE GRPREC(IPREC)
      INTEGER IPREC
      INTEGER ERRIND, FONT, PREC
      CALL GQTXFP(ERRIND, FONT, PREC)
      CALL GSTXFP(FONT, IPREC)
      RETURN
      END
