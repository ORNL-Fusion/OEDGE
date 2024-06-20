C@PROCESS NOSDUMP NOGOSTMT OPT(3)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRAXLOG(X10,Y10,X20,Y20,S10,S20,LINKS0,ACHSE)
************************************************************************
*     KOORDINATENACHSE (LOGARITHMISCH) IN BELIEBIGER RICHTUNG:         *
*     (X1,Y1)   SIND DER ANFANG DER ACHSE                              *
*     (X2,Y2)   SIND DER ENDPUNKT DER ACHSE                            *
*     S1,S2     SIND DIE WERTE DER ENDPUNKTE; 0<S1<=S2                 *
*     LINKS     GIBT INFORMATION UEBER DIE BESCHRIFTUNG MIT ZAHLEN:    *
*               .TRUE.  - BESCHRIFTUNG LINKS BZGL. DER ORIENTIERUNG    *
*               .FALSE. - BESCHRIFTUNG RECHTS BZGL. DER ORIENTIERUNG   *
*     ACHSE     1       - DIE GANZE ACHSE WIRD GEPLOTTET               *
*               0       - DIE (VORHANDENE) ACHSE WIRD NUR BESCHRIFTET  *
*               2       - DIE ACHSE WIRD gezeichnet und nur markiert   *
*             UPDATE: 29.9.89 GROTEN (ACHSENBEMASSUNG WEGGELASSEN)     *
*             UPDATE: 28.6.90 BUSCH: FEHLER BEI GSTXAL KORRIGIERT      *
*             UPDATE: 18.6.91 Groten: Abstaende der Schrift verbessert *
*             UPDATE: 17.7.91 Groten: COMMON /GRLOGY/ usw.             *
*             UPDATE:  5.5.92 Groten: Mehr Unterteilungen bei <= 2 Dec.*
*             UPDATE:  5.5.92 Groten: Gitterlinien auch zwischen Decad.*
************************************************************************
C---- UMRECHNUNG BOGENMASS-->GRAD; PI; LOG10(2), LOG10(5)
      PARAMETER (GRDFAK=57.29578, PI2TEL=1.570796, HAL=.5)

      CHARACTER(len=3) ZAHL
      LOGICAL LINKS,X2X5,LINKS0,SUBDEC,NETZ,LETZ
      integer achse
      LOGICAL YINNEN,YNETZ,YL
      REAL FL(8)
      CHARACTER(len=8) ZL
      COMMON /GRLOGY/ ISLOGY,YINNEN,YNETZ,YL,DSH1,DSH2,DSH3
CDEC$ PSECT /GRLOGY/ NOSHR
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR

      SAVE /GRLOGY/ , /GRPP/

      DATA FL/.30103004,.698970079,.47712125,.602059991,.77815125,
     $        .845098040,.90308999,.954242509/
      DATA ZL/'25346789'/


      IF ( S10*S20 .LE. 0. .OR. S10.EQ.S20 ) THEN
        WRITE (*,*) 'GRAXLOG: PARAMETER 5 AND/OR 6 IS WRONG:',S10,S20
        STOP
      ENDIF
      X1=X10
      Y1=Y10
      X2=X20
      Y2=Y20
      SUBDEC=S10.GT.0. .AND. S20.GT.0.
      S1=ABS(S10)
      S2=ABS(S20)
      LINKS=LINKS0
      IF (S1.GT.S2) THEN
         S1=ABS(S20)
         S2=ABS(S10)
         X1=X20
         X2=X10
         Y1=Y20
         Y2=Y10
         LINKS=.NOT.LINKS0
      ENDIF
C---- RETTUNG DER AKTUELLEN PLOTPARAMETER
      PP5=PP(5)
      PP6=PP(6)
      PP7=PP(7)
      PP8=PP(8)
      PP10=PP(10)
      PP11=PP(11)
      PP12=PP(12)
      PP15=PP(15)
      CALL GQTXAL(IER,IX0,IY0)

      CSIZE=PP(14)
      X1CM=PP(1)+(X1-PP5)*(PP(3)-PP(1))/(PP7-PP5)
      Y1CM=PP(2)+(Y1-PP6)*(PP(4)-PP(2))/(PP8-PP6)
      X2CM=PP(1)+(X2-PP5)*(PP(3)-PP(1))/(PP7-PP5)
      Y2CM=PP(2)+(Y2-PP6)*(PP(4)-PP(2))/(PP8-PP6)
      CALL GRDSH(1.,0.,1.)
      IF (ACHSE.gt.0) THEN
        CALL GRJMP(X1,Y1)
        CALL GRDRW(X2,Y2)
      END IF
      XDIF=X2CM-X1CM
      YDIF=Y2CM-Y1CM
      CM=SQRT(XDIF**2+YDIF**2)
      IF (CM .LE. CSIZE) RETURN
      CALL GRSCLV(PP(1),PP(2),PP(3),PP(4))
C---- RICHTUNG DER ACHSE
      PHI=ATAN2(YDIF,XDIF)
      PHIGRD=PHI*GRDFAK
C---- RICHTUNG DES TEXTES UND DER MARKIERUNGSSTRICHE
      IF (PHIGRD .LE. 0.) THEN
         ALPHA=PHIGRD+90.
         IF (LINKS) THEN
            IX=1
            XM=CSIZE*hal*COS(PHI+PI2TEL)
            YM=CSIZE*hal*SIN(PHI+PI2TEL)
         ELSE
            IX=3
            XM=CSIZE*hal*COS(PHI-PI2TEL)
            YM=CSIZE*hal*SIN(PHI-PI2TEL)
         ENDIF
      ELSE
         ALPHA=PHIGRD-90.
         IF (LINKS) THEN
            IX=3
            XM=CSIZE*hal*COS(PHI+PI2TEL)
            YM=CSIZE*hal*SIN(PHI+PI2TEL)
         ELSE
            IX=1
            XM=CSIZE*hal*COS(PHI-PI2TEL)
            YM=CSIZE*hal*SIN(PHI-PI2TEL)
         ENDIF
      ENDIF
      NETZ=ISLOGY.EQ.1234567890 .AND. YNETZ
      LETZ=ISLOGY.EQ.1234567890 .AND. YL
      IF (NETZ.OR.LETZ) THEN
         XGEGEN=PP(1)
         IF (LINKS) XGEGEN=PP(3)
      ENDIF
      IF (ISLOGY.EQ.1234567890 .AND. YINNEN) THEN
         XMPIX=-XM
         YMPIX=-YM
         IF ( LINKS ) THEN
            XM2=XM*2.
            YM2=YM*2.
         ELSE
            XM2=XM
            YM2=YM
         ENDIF
         XM=0.7*XM
         YM=0.7*YM
      ELSE
         XMPIX=XM
         YMPIX=YM
         IF ( LINKS ) THEN
            XM2=2.8*XM
            YM2=2.8*YM
         ELSE
            XM2=1.5*XM
            YM2=1.5*YM
         ENDIF
      ENDIF
      XM2PIX=2.*XMPIX
      YM2PIX=2.*YMPIX
C---- BESTIMMUNG VON DS,DSCM: ABSTAND DER ZAHLEN UND MARKIERUNGEN
      S1LOG=LOG10(S1)
      S2LOG=LOG10(S2)
      SDIF=S2LOG-S1LOG
      ISMIN=NINT(S1LOG)
      IF (10.**(ISMIN-0.001) .GT. S1) ISMIN=ISMIN-1
      ISMAX=NINT(S2LOG)
      IF (10.**(ISMAX+0.001) .LT. S2) ISMAX=ISMAX+1
C---- ABSTANDSKONTROLLE
      ABST=CM/SDIF
      X2X5= (ABST .GE. 4.*CSIZE) .AND. SUBDEC
      ISTEP=(CSIZE+CSIZE)/ABST
      IF (ISTEP .LE. 0) ISTEP=1
      ISMIN=(ISMIN-(ISTEP-1))/ISTEP*ISTEP
      ML = 2
      IF ( REAL(ISMAX-ISMIN)/ISTEP.LE.2 .AND. ABST.GE.20*CSIZE) ML = 8

      DO 100 I=ISMIN,ISMAX,ISTEP
        S=REAL(I)
        IF (X2X5) THEN
          CALL GSTXAL(IX,3)
          CALL GRCHRC(CSIZE*0.67,ALPHA,PP(16))
          DO 55 J=1,ML
            IF (INTER(S1LOG,S+FL(J),S2LOG)) THEN
              X0=X1CM+(S+FL(J)-S1LOG)/SDIF*XDIF
              Y0=Y1CM+(S+FL(J)-S1LOG)/SDIF*YDIF
              CALL GRJMP(X0,Y0)
              CALL GRDRW(X0+XMPIX,Y0+YMPIX)
              IF ( LETZ) THEN
                 CALL GRDSH(DSH1,DSH2,DSH3)
                 CALL GRJMP(X0,Y0)
                 CALL GRDRW(XGEGEN,Y0)
                 CALL GRDSH(1.,0.,1.)
              ENDIF
              IF (ACHSE.LT.2) CALL GRTXT(X0+XM*1.5,Y0+YM*1.5,1,ZL(J:J))
            ENDIF
   55     CONTINUE
        END IF
        IF (INTER(S1LOG,S,S2LOG)) THEN
          X0=X1CM+(S-S1LOG)/SDIF*XDIF
          Y0=Y1CM+(S-S1LOG)/SDIF*YDIF
          CALL GRJMP(X0,Y0)
          CALL GRDRW(X0+XM2PIX,Y0+YM2PIX)
          IF ( NETZ .AND. S.NE.S1LOG .AND. S.NE.S2LOG ) THEN
             CALL GRDSH(DSH1,DSH2,DSH3)
             CALL GRJMP(X0,Y0)
             CALL GRDRW(XGEGEN,Y0)
             CALL GRDSH(1.,0.,1.)
          ENDIF
          WRITE (ZAHL,'(I3)') I
          IST=1
          IF ( ZAHL(1:1).EQ.' ' ) IST=2
          IF ( ZAHL(2:2).EQ.' ' ) IST=3
          if (achse.lt.2) then
             IF (IX .EQ. 1) THEN
C------------- IX = 1
               CALL GRCHRC(CSIZE,ALPHA,PP(16))
               CALL GSTXAL(IX,3)
               CALL GRTXT(X0+XM2*1.25,Y0+YM2*1.25,2,'10')
               CALL GRCHRC(CSIZE*0.67,ALPHA,PP(16))
               CALL GSTXAL(IX,4)
               CALL GRTXTC(3,ZAHL(IST:3))
             ELSE
C------------- IX = 3
               CALL GRCHRC(CSIZE*0.67,ALPHA,PP(16))
               CALL GSTXAL(IX,4)
               CALL GRTXT(X0+XM2*1.25,Y0+YM2*1.25,3,ZAHL(IST:3))
               CALL GRCHRC(CSIZE,ALPHA,PP(16))
               CALL GSTXAL(IX,3)
               CALL GRTXTC(2,'10')
             END IF
          endif
        END IF
  100 CONTINUE
C---- ZURUECKSETZEN DER PLOTPARAMETER
      CALL GSTXAL(IX0,IY0)
      CALL GRCHRC(CSIZE,PP15,PP(16))
      CALL GRSCLV(PP5,PP6,PP7,PP8)
      CALL GRDSH(PP10,PP11,PP12)


      contains 
      
      function inter(X1,X,X2) result(erg)
      REAL,intent(in) :: X1,X,X2
      LOGICAL erg
      
      erg = X1-0.0001 <= X .AND. X <= X2+0.0001
      END function inter
      end subroutine GRAXLOG
