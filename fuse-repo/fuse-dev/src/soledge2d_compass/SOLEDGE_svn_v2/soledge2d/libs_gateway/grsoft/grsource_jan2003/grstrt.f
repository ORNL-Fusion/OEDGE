C@PROCESS OPT(3) NOSDUMP NOGOSTMT
***********************************************************************
** *** COPYRIGHT 1990 BY FORSCHUNGSZENTRUM (KFA) JUELICH GMBH    ***  *
**                                                                    *
** NAME:    GRSTRT                                                    *
**                                                                    *
** ZWECK:   STARTROUTINE GR-SOFTWARE auf UNIX Workstation             *
**                       GLIGKS, XGKS auf AIX und CRAY                *
**                                                                    *
** AUTOR :   M. BUSCH,        (ma.busch@kfa-juelich.de)               *
**                                                                    *
** DATUM :   18.10. 1990                                              *
** UPDATE:   22.11. 1990                                              *
** UPDATE:   20.12. 1990 beliebig bildgroesee grsclp                  *
** UPDATE:   18. 2. 1991 Eingabe des Augabedevices und Hintergrund    *
**                       GKS IMplementation und CGM Name              *
**                       WISS, TERM und CGM TYP wird automatisch      *
**                       Default: XGKS, weisser Hintergrund, Ausgabe  *
**                       Terminal
**                       erfragt                                      *
* Update: 26. 2.92 Busch COMMON GRREST eingefuegt
* Update: 18. 3.92 Busch COMMON GRTXPO eingefuegt, current position
*                        fuer Text neu (frueher: nur eine curr. pos.
** UPDATE:   13. 4. 1992 SUNGKS kann: Terminal, CGM, PS               *
**                       setenv GRGKS SUNGKS                          *
**                       setenv GRBG  SCHWARZ                         *
**                              Background  schwarz fuer TERM, CGM    *
**                              und weiss fuer PS                     *
**                       setenv GRCGM name.cgm |  name.ps             *
**                       setenv GRDEVICE TERM | CGm | PS              *
** UPDATE:   10. 8. 1992 Busch character spacing generell fuer xgks   *
** UPDATE:    8. 2. 1993 Busch COMMON GRCOM1 mit GKDEV1=5003, damit   *
**                             fuer interaktives gr3pan das X-terminal*
**                             wie sonst das 3270G Terminal behandelt *
**                             wird                                   *
**                             COMMON GRDEV    neu                    *
**                             COMMON GRCOM1   neu                    *
** Update    12. 8. 93 Busch   Anpassung fuer GLIGKS, GTSGRAL GKS     *
**                             XGKS und SUNGKS ueber Env. Variable    *
**                             GRGKS                                  *
***********************************************************************
***********************************************************************
C
      SUBROUTINE GRSTRT(CAMERA, DDNUMB)

      INTEGER CAMERA, DDNUMB,  LASF(13), ERRUN, CONID,
     1        WKTYP, WISCO, WISS, BACK,  DEV1, DEV2,
     2        GKDEV1, GKDEV2, IWK1, IWK2, IWKDUM, IGRSHOW, IGR3,
     3        SHOWPR, IGR3PL, FLPIC, RAHMEN, INTLIN, ICOLOR, IFONT

      REAL XLN(300), YLN(300), XMR(300), YMR(300)

      CHARACTER DEVICE*4, GRDEVVAL*4, GKSTYP*7

      LOGICAL EINZEL, PROMPT, FLGROT
CCC
CCC   COMMONBLOCK KURVEF,GRBLD, USW.
      COMMON /GRCIL/ CILBER(8)
CDEC$ PSECT /GRCIL/ NOSHR
      COMMON /GRDEV/ BACK, GRDEVVAL
CDEC$ PSECT /GRDEV/ NOSHR
      COMMON /GRWK/ IXGKS, IWK, IWKWIS
CDEC$ PSECT /GRWK/ NOSHR
      COMMON /GRGKS/ GKSTYP
CDEC$ PSECT /GRGKS/ NOSHR
      COMMON /SCALE/ XMAXDC, XUNITS, YUNITS
CDEC$ PSECT /SCALE/ NOSHR
      COMMON /MKRTP/ MKTYP
CDEC$ PSECT /MKRTP/ NOSHR
      COMMON /GRREST/ MAXPKT, NRLN, XCURR, YCURR, XLN, YLN, NRMR, XMR,
     1       YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /GRTXPO/ XTCURR, YTCURR
CDEC$ PSECT /GRTXPO/ NOSHR
C
CCC   COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRPIC/ FLPIC, NSCLC, NSCLV, NSCLP, RAHMEN, XMAXCM, YMAXCM,
     1       XDCPIC, YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON /GKSCM1/ DEV1, DEV2, GKDEV1, GKDEV2, IWK1, IWK2, IWKDUM,
     1       CONID, EINZEL, PROMPT
CDEC$ PSECT /GKSCM1/ NOSHR
      COMMON /GRCOM1/ IGRSHOW, ERRUN, IGR3, SHOWPR, IGR3PL
CDEC$ PSECT /GRCOM1/ NOSHR
      SAVE /GRCIL/ ,/GRDEV/ , /GRWK/, /GRGKS/,/SCALE/, /MKRTP/,/GRREST/,
     $     /GRTXPO/,/GRPP/,/GRPIC/,/GKSCM1/ ,/GRCOM1/

      EQUIVALENCE (PP(13),INTLIN), (PP(16),ICOLOR), (PP(17),SIZMRK)
      EQUIVALENCE (PP(9),IFONT), (PP(18),FLGROT)

      DATA LASF /13*1/


C     Indentifikation der UNIX Version
      IXGKS = 1234567890

      IWK=1

      ERRUN = 93

C     FLPIC=0   GRSCLP PICTURE SIZE WURDE NICHT AUFGERUFEN
C     NSCLC=0   ZAEHLER GRSCLC AUFRUFE - GRSCLP WIRD IGNORIERT WENN NIC
C               VOR ERSTEM GRSCLC AUFGERUFEN WURDE (BILD)
C     NSCLP=0   ZAEHLER GRSCLP AUFRUFE
C     PP(9)    HAT KEINE BEDEUTUNG MEHR (UEBERALL XMAXCM,YMAXCM
C              EINGESETZT
      NSCLC = 0
      NSCLV = 0
      NSCLP = 0
      RAHMEN = 0
      FLPIC = 0
C     Max. Anzahl Punkte pro GRLN, bzw Marker bei GPM
      MAXPKT = 300

      PP(1) = 0.
      PP(2) = 0.
      PP(3) = 39.5
      PP(4) = 28.7
      PP(5) = 0.
      PP(6) = 0.
      PP(7) = 39.5
      PP(8) = 28.7
c frueher: Rueckvergroesserungsfaktor
C jetzt Font
c     PP(9)=10.5
      IFONT = 1
      PP(10) = 1.
      PP(11) = 0.
      PP(12) = 1.
      INTLIN = 16
      PP(14) = 0.3
      PP(15) = 0.
c pp(16)
      ICOLOR = 1
c pp(17)
      SIZMRK = 0.2
c pp(18)
      FLGROT = .FALSE.

C Hintergrund weiss Standard (GLIGKS,XGKS,CGI(GRALGKS)
C
      BACK = 1

C--------------------------------------------------------------------
C GET Environment
C--------------------------------------------------------------------

      CALL GRGETE


C     FALLS PICTURE SIZE NICHT GESETZT - > STANDARD
      XMAXCM = 39.5
      YMAXCM = 28.7
      XDCPIC = 1.

c Voreinstellung testweise auf NT - da ohne WISS sonst nichts geht
      IF (GKSTYP .EQ. 'GLIGKS'.OR.GKSTYP .EQ. 'GRALGKS' ) THEN
        YDCPIC = 28.7/39.5
      ELSE
        YDCPIC = 1.
      END IF

CCC   INITIALISIERUNG DES CILBER-COMMON
CCC
      CILBER(1) = 0E0
      CILBER(2) = 0E0
      CILBER(3) = 0E0
      CILBER(4) = 0E0
      CILBER(5) = 0E0
      CILBER(6) = 0E0
      CILBER(7) = 0E0
      CILBER(8) = 0E0


CCC   INITIALISIEREN MARKERTYP
      MKTYP = 999
CCC
CCC   8. 2. 93
      SHOWPR = 1
CCC
CCC
CCC   OPEN GKS
CCC
      CONID = DDNUMB
      WKTYP = CAMERA
*UPDATE 7.11.90
      IBUF = 0


C schliesse GKS
      CALL GECLKS

      CALL GOPKS(ERRUN, IBUF)


CCC   SET ASPECT SOURCE FLAG
CCC   ALLE ENTRIES = 1 (INDIVIDUAL)
      CALL GSASF(LASF)

C UNIX 9.2.93
      GKDEV1 = WKTYP
      CALL GSCLIP(1)


      IF (GKSTYP .EQ. 'GLIGKS') THEN
        CALL GRCONF(WKTYP)

      ELSEIF ( GKSTYP.EQ.'GRALGKS') THEN
C OPEN WK
        CALL GOPWK(IWK, CONID, WKTYP)
        CALL GACWK(IWK)
        CALL GCRSG(1)


C    (SUNGKS, XGKS)
      ELSE

        CALL GRCNV(WKTYP, CONID, WISS, WISCO, DEVICE)
        IWKWIS = 2
        IWK = 1
C open AND ACTIVATE WISS

        CALL GOPWK(IWKWIS, WISCO, WISS)
        CALL GACWK(IWKWIS)

C OPEN WK
        CALL GOPWK(IWK, CONID, WKTYP)
C       CALL GACWK(IWK)



C     Definition der 8 GR Farben (SUNGKS,XGKS)
C     wenn kein CGM erstellt wird, werden in GRSCR die ersten 8 Farben
C     gesetzt

        IF (DEVICE .NE. 'CGM') CALL GRSCR(IWK)


        CALL GCRSG(1)

      END IF


      VWXP = 0.
      VWYP = 0.
      VWXQ = 1.
      IF (GKSTYP .EQ. 'GLIGKS'.OR.GKSTYP .EQ. 'GRALGKS') THEN
        VWYQ = 28.7/39.5

C     (SUNGKS, XGKS)
      ELSE
        VWYQ = 1.
      END IF

      INT = 1
      CALL GQWKCA(WKTYP, IERR,ICAT)
      IF (GKSTYP .EQ. 'XGKS') THEN
        IF (DEVICE .EQ. 'TERM') CALL GQDSP(WKTYP, IERR, IDUN, XDC, YDC,
     1     IX, IY)
        VWXQ = 1.
        VWYQ = YDC/XDC
      END IF

      CALL GSVP(INT, VWXP, VWXQ, VWYP, VWYQ)
C     DEFAULT WINDOW DIN A3
      CALL GSWN(INT, PP(5), PP(7), PP(6), PP(8))
CCC   SELECT TRANSFORMATIONS NUMMER
      CALL GSELNT(INT)
C     Workstation Window fuer Metafile


C     nicht fuer WISS
      IF (ICAT.NE.3) CALL GSWKWN( IWK,0.,VWXQ,0.,VWYQ)

CCC   INITIALISIEREN DER COMMON-WERTE
CCC
      XUNITS = 1E0
      YUNITS = 1E0
C     Initialisieren Current Position fuer GRLN(C),GRCHN(C)
      XCURR = 0.
      YCURR = 0.
C     Initialisieren Current Position fuer GRTXT(C)
      XTCURR = 0.
      YTCURR = 0.
      NRLN = 0
      NRMR = 0
CCC   SET TEXT FONT UND PRECISION
      CALL GSTXFP(1, 2)
      CALL GRCHRC(PP(14), PP(15), IDUMMY)
CCC   TEST STIFT 1 DEFAULT (WEISS)
      CALL GRNWPN(1)
CCC   char. spacing nach Test von Herrn Groten
CCC   10.8.92
      IF (GKSTYP.EQ.'XGKS') CALL GSCHSP(0.3)
      END

C--------------------------------------------------------------------
C GET Environment
C
C Default: GLIGKS
C          keine Env Variable erforderlich
C
C          XGKS
C          export GRGKS=XGKS
C
C          SUNGKS
C          export GRGKS=SUNGKS
C          export GRBG=SCHWARZ
C
C          CGI mit GRALGKS auf der CRAY
C          GRGKS wird erkannt wenn zuvor das Script gks_init aufge-
C          rufen wurde
C--------------------------------------------------------------------

      SUBROUTINE GRGETE
      INTEGER BACK

      CHARACTER GRDEVICE*8, GRDEVVAL*4, GRBG*4, GRBGV*7, GRGKS*5,
     1          GRGKSV*7, GRCGM*5, GRCGMV*20,GTSHV*14,GTSHOME*7
      COMMON /GRCGME/ GRCGMV
CDEC$ PSECT /GRCGME/ NOSHR
      COMMON /GRDEV/ BACK, GRDEVVAL
CDEC$ PSECT /GRDEV/ NOSHR
      COMMON /GRGKS/ GRGKSV
CDEC$ PSECT /GRGKS/ NOSHR

      SAVE /GRCGME/,/GRDEV/, /GRGKS/


      GRDEVICE = 'GRDEVICE'
      GRBG = 'GRBG'
      GRGKS = 'GRGKS'
      GRCGM = 'GRCGM'
      GTSHOME = 'GTSHOME'

C     GRDEVICE Type der Ausgabe: TERM = X Workstation, Terminal, INOUT
C                               CGM  = Ausgabe auf CGM
C             Ausgabe ist entweder auf Terminal oder auf CGM moeglich
C                               Default: Terminal
C     GRBG    Backgroundfarbe des GKS:  SCHWARZ=0, WEISS=1
C                                      Default: WEISS=1
C             SUNGKS: scharzer Hintergrund, XGKS: weisser Hintergrund
C     GRGKS   Default: XGKS
C             moeglich: XGKS, SUNGKS, GLIGKS, GRALGKS
C     GRCGM   Name des CGM, wird durch Extension .CGM erweitert
C             Default name: GR.CGM
C     GTSHOME CGI auf der CRAY , wenn initialisiert, dann dort cgi
C              directory
C             siehe gks_init auf CRAY - dort
C             wird gesetzt durch . /usr/local/cgi/bin/gralenv.sh

      grdevval = ' '
      grbgv = ' '
      grgksv = ' '
      grcgmv = ' '
      gtshv = ' '

Cmb 5.11.93
C Probleme bei VMS System  (GETENV nur fuer XGKS, SUNGKS notwendig)
C --> Kommentar
cjh
C      CALL GETENV(GRDEVICE, GRDEVVAL)
C      CALL GETENV(GRBG, GRBGV)
C      CALL GETENV(GRGKS, GRGKSV)
C      CALL GETENV(GRCGM, GRCGMV)
C      CALL GETENV(GTSHOME, GTSHV)
cjh
      grdevice = ' '
      grbg = ' '
      grgks = ' '
      grcgm = ' '


C Default GKS=GLIGKS
      IF (GRGKSV.EQ.' ') THEN
         BACK=1
         GRGKSV='GLIGKS'
         goto 999
      ENDIF

      IF (GRBGV .EQ. 'SCHWARZ') THEN
        BACK = 0

      ELSE
        BACK = 1
      END IF


999   END

C--------------------------------------------------------------------
C     nur aufgerufen wenn nicht GLIGKS, dh. XGKS, SUNGKS
C
      SUBROUTINE GRCNV(WKTYP, CONID, WISS, WISCO, DEVICE)
C--------------------------------------------------------------------
      INTEGER WKTYP, CONID, WISS, WISCO, CGMTYP, CGMCND, BACK, OUTTYP,
     1        OUTCND, IOCGM

      CHARACTER GRDEVVAL*4, GRGKSV*7, GRCGMV*20, CGMNAME*20, DEVICE*4
      COMMON /GRCGME/ GRCGMV
CDEC$ PSECT /GRCGME/ NOSHR
      COMMON /GRDEV/ BACK, GRDEVVAL
CDEC$ PSECT /GRDEV/ NOSHR
      COMMON /GRGKS/ GRGKSV
CDEC$ PSECT /GRGKS/ NOSHR
      SAVE /GRCGME/,/GRDEV/ ,/GRGKS/
      DATA iocgm/2/
c     number of available workstations
      CALL GQEWK(1, IERR, NWT, IWT)
      DO 1, I = 1, NWT
        CALL GQEWK(I, IERR, NWT, IWT)
C     workstation categorie
        CALL GQWKCA(IWT, IERR, ICAT)

c     icat=output, input,outin, wiss, mo,mi)
c         =0     , 1    ,2    , 3   , 4 ,5 )

        IF (ICAT .EQ. 2) THEN
          OUTTYP = IWT
          OUTCND = 0

        ELSE IF (ICAT .EQ. 3) THEN
          WISS = IWT
          WISCO = 90
C TEST EVTL AUCH MIT 0
        ELSE IF (ICAT .EQ. 4) THEN
          CGMTYP = IWT
          CGMCND = IOCGM
        END IF

    1 CONTINUE
c

      IF (GRDEVVAL .EQ. 'CGM') THEN
        WKTYP = CGMTYP
        CONID = CGMCND
        DEVICE = 'CGM'

      ELSE IF (GRDEVVAL .EQ. 'TERM') THEN
        WKTYP = OUTTYP
        CONID = OUTCND
        DEVICE = 'TERM'

      ELSE
        WKTYP = OUTTYP
        CONID = OUTCND
        DEVICE = 'TERM'
      END IF


C     CGM entweder vomBenutzer gesetzen Name oder Default
      IF (GRCGMV .NE. ' ') THEN
        CGMNAME = GRCGMV

      ELSE
        CGMNAME = 'gr.cgm'
      END IF

      IF (GRDEVVAL .EQ. 'CGM') OPEN (IOCGM,file=CGMNAME,
     *  status='unknown')

C     wegen Abfrage in GRSHOW
      IF (GRDEVVAL .EQ. ' ') GRDEVVAL = 'TERM'

      END


C----------------------------------------------------------------------
C      Setzen der ersten 8 Farben (nur SUNGKS, XGKS)
C
      SUBROUTINE GRSCR(IWK)
C--------------------------------------------------------------------
      INTEGER BACK
      CHARACTER GRDEVVAL*4
      COMMON /GRDEV/ BACK, GRDEVVAL
CDEC$ PSECT /GRDEV/ NOSHR
       SAVE /GRDEV/

C Farbe 0
      CALL GSCR(IWK, 0, 0., 0., 0.)

c Farbe 1 weiss
c xgks 1.Farbe schwarz wegen weissem Hintergrund
      IF (BACK .EQ. 1) THEN
        CALL GSCR(IWK, 1, 0., 0., 0.)
C XGKS Farbe 8: immer schwarz bei GR-Software
        CALL GSCR(IWK, 8, 0.,0.,0.)

C SunGKS 1.Farbe weiss wegen schwarzem Hintergrund
      ELSE
        CALL GSCR(IWK, 1, 1., 1., 1.)
C Default GKS=XGKS
C SunGKS Farbe 8: schwarz
        CALL GSCR(IWK, 8, 0., 0., 0.)
      END IF
C Farbe 2: rot
      CALL GSCR(IWK, 2, 1., 0., 0.)
C Farbe 3: blau
      CALL GSCR(IWK, 3, 0., 0., 1.)
C Frabe 4: gruen
      CALL GSCR(IWK, 4, 0., 1., 0.)
c Farbe 5: pink
      CALL GSCR(IWK, 5, 1., 0., 1.)
C Farbe 6: gelb
      CALL GSCR(IWK, 6, 1., 1., 0.)
C Farbe 7: tuerkis
      CALL GSCR(IWK, 7, 0., 1., 1.)

      END
