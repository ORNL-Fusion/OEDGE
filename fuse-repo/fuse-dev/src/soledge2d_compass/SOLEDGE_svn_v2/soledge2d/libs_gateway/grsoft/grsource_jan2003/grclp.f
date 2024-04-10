C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C-----------------------------------------------------------------------
C     Subroutine zum Clippen
C     an den gerade gueltigen Grenzen (GRSCLC und GRSCLV)
C
C     Parameter: ICLIP           Integer*4, INPUT
C
C     Bedeutung: ICLIP  = 0  --> Ausschalten des Clippen
C                ICLIP ^= 0  --> Anschalten des Clippen
C     Achtung: bei eingestelltem Clipping kann z.B. keine Achsen-
C              beschriftung gemacht werden
C     Sinnvolle Anwendung: z.B. Kurve, die eine Polstelle hat darstellen
C                          Vorgehen: Achse und Kurve getrennt zeichnen
C                          zuerst ohne Clipping die Achse und dann mit
C                          Clipping die Kurve
C     Author:  M. Busch (nach Vorlage U. Funk  IFF011)
C     changed: M. Busch 2. 4. 91  ( COMMON GRPIC )
C     CHANGED: M. BUSCH 18. 6. 91
C     CHANGED: M. BUSCH 21.10. 91 Klippen mit Drehen des Bildes
C                       (nach Vorlage von Hern Kleefeld, AVR, REW085)
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C-----------------------------------------------------------------------


      SUBROUTINE GRCLP (ICLIP)

c     IMPLICIT NONE
      INTEGER ICLIP, ERRIND ,NRLN
      INTEGER FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN
      REAL    XMAXCM,YMAXCM, XDCPIC, YDCPIC, XCURR, YCURR
      REAL    VXMIN, VXMAX, VYMIN, VYMAX, FAKX, FAKY
      real    WIND(4), VIEW(4)
      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      LOGICAL FLGROT
      COMMON /GRPP / PP(18)
CDEC$ PSECT /GRPP / NOSHR
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMAXCM,YMAXCM, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      EQUIVALENCE  (PP(18),FLGROT)
      SAVE /GRPP / ,  /GRREST/, /GRPIC/


C---- Ausgabe evt. vorhandener GRJMP- und GRDRW-Daten

      IF(NRLN.GT.1) THEN
        CALL GPL(NRLN,XLN,YLN)
        XCURR=XLN(NRLN)
        YCURR=YLN(NRLN)
        NRLN=0
      ENDIF

      IF(NRMR.GT.0) THEN
        CALL GPM(NRMR,XMR,YMR)
        NRMR=0
      ENDIF

C-----------------------------------------------------------------------
C     ICLIP = 0 --> Ausschalten des Clippen
C-----------------------------------------------------------------------
      IF (ICLIP.EQ.0) THEN
        CALL GSELNT (1)

      ELSE

C-----------------------------------------------------------------------
C     ICLIP ^= 0 --> Anschalten des Clippen
C-----------------------------------------------------------------------

c21.10.91. IF (FLGROT) THEN ... ENDIF neu eingefuegt

      IF ( FLGROT ) THEN

C        Koordinatentransformationen
C        ---------------------------
         X1CM = XMAXCM - PP(4)
         X2CM = XMAXCM - PP(2)
         Y1CM = PP(1)
         Y2CM = PP(3)

         x STEIG   =   ( PP(7) - PP(5) ) / ( PP(3) - PP(1) )
         X1VL = PP(5) + x STEIG * ( X1CM - PP(1) )
         X2VL = PP(5) + x STEIG * ( X2CM - PP(1) )

         y STEIG   =   ( PP(8) - PP(6) ) / ( PP(4) - PP(2) )
         Y1VL = PP(6) + y STEIG  * ( Y1CM - PP(2) )
         Y2VL = PP(6) + y STEIG  * ( Y2CM - PP(2) )

C       1. 'Scale cm': Lage des gedrehten Plotfensters in Einheiten des
C          n i c h t  gedrehten Koordinatensystems
C       ---------------------------------------------------------------
         CALL GRSCLC(X1CM, Y1CM, X2CM, Y2CM)

C       2. 'Scale Value' mit Werten, die dem gedrehten Plotfenster im
C          n i c h t  gedrehten Koordinatensystem zukaemen
C       -------------------------------------------------------------
         CALL GRSCLV(X1VL, Y1VL, X2VL, Y2VL)

        ENDIF


C       3. Dann erst Clipping einschalten (FLGROT=.TRUE.)
C       -------------------------------------------------
C       bzw. bei nicht gedreht: Nur Clipping einschalten
C       ------------------------------------------------


        CALL GQNT(1,ERRIND,WIND,VIEW)
        FAKX  =(VIEW(2)-VIEW(1))/XMAXCM
        VXMIN = PP(1)*FAKX
        VXMAX = PP(3)*FAKX
        FAKY  =(VIEW(4)-VIEW(3))/YMAXCM
        VYMIN = PP(2)*FAKY
        VYMAX = PP(4)*FAKY
        VXMIN = AMAX1(0.,VXMIN)
        VYMIN = AMAX1(0.,VYMIN)
        VXMAX = AMIN1(1.0,VXMAX)
        VYMAX = AMIN1(1.0,VYMAX)
        CALL GSVP   (2,VXMIN,VXMAX,VYMIN,VYMAX)
        CALL GSWN   (2,PP(5),PP(7),PP(6),PP(8))
        CALL GSELNT (2)

c21.10.91. IF (FLGROT) THEN ... ENDIF neu eingefuegt

      IF ( FLGROT) THEN
C       4. Jetzt 'Scale cm' und 'Scale Value' mit den urspruengl. Werten
C       ----------------------------------------------------------------
         CALL GRSCLC( PP(1), PP(2), PP(3), PP(4) )
         CALL GRSCLV( PP(5), PP(6), PP(7), PP(8) )
      ENDIF

      ENDIF
C
      RETURN
      END
