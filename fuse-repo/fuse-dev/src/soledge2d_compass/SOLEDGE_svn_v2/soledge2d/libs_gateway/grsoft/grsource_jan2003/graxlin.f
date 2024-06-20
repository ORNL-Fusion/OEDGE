C@PROCESS NOSDUMP NOGOSTMT OPT(3)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRAXLIN(VON,WO,BIS,HIER,UNT,OB,LINKS,ACHSE)
************************************************************************
*     KOORDINATENACHSE (LINEAR) IN BELIEBIGER RICHTUNG:                *
*     (VON,WO)   SIND DER ANFANG DER ACHSE                             *
*     (BIS,HIER) SIND DER ENDPUNKT DER ACHSE                           *
*     UNT,OB    SIND DIE WERTE DER ENDPUNKTE                           *
*     LINKS     GIBT INFORMATION UEBER DIE BESCHRIFTUNG MIT ZAHLEN:    *
*               .TRUE.  - BESCHRIFTUNG LINKS BZGL. DER ORIENTIERUNG    *
*               .FALSE. - BESCHRIFTUNG RECHTS BZGL. DER ORIENTIERUNG   *
*     ACHSE     1       - DIE GANZE ACHSE WIRD GEPLOTTET               *
*               0       - DIE (VORHANDENE) ACHSE WIRD NUR BESCHRIFTET  *
*               2       - DIE ACHSE WIRD gezeichnet und nur markiert   *
************************************************************************
* AUTOR: GERD GROTEN                                                   *
* COPYRIGHT: KFA JUELICH (W. GERMANY)                        1988      *
*      UPDATE 29.9.89 GROTEN (ACHSE WIRD INTEGER, SCHRIFT KANN FEHLEN)
*      UPDATE 28.6.90 BUSCH  FEHLER BEI GSTXAL KORRIGIERT
*      UPDATE 14.12.90 BUSCH  PP(9) DURCH COMMON GRPIC ERSETZT
*      update 12.11.91 Groten IMPLICIT REAL*8 weg (RS6000) , SCLC,SCLV
*      update 02.11.92 Groten MAXTEIL= N I N T (FMATEIL)
*      update 14.06.93 Groten Exponentenbeschriftung bei Enge mehr ausen
*      update 20.07.93 Groten Aenderungen zur Vermeidung von Uberlappung
************************************************************************
      REAL VON,WO,BIS,HIER,UNT,OB
      LOGICAL LINKS
      INTEGER RAHMEN
      integer achse
      INTEGER FLPIC
      REAL PP,ZEIGRO
      REAL XMAXCM,YMAXCM,XDCPIC,YDCPIC
      REAL SCLV(4),SCLC(4)
      DOUBLE PRECISION p,huepf

      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMAXCM,YMAXCM, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      SAVE /GRPP/, /GRPIC/
      EQUIVALENCE (PP(14),ZEIGRO)

      REAL XX(2),YY(2)
      LOGICAL LINKS0
      CHARACTER   FORM*6, DRUCK*9, ZEHN*5, POT*8
      DATA ZEHN /'*10  '/, POT /'     ???'/
      IF (UNT.EQ.OB) THEN
         WRITE(*,*) 'GRAXLIN: Parameter 5 = Parameter 6'
         STOP
      ENDIF
      OB0=OB
      UNT0=UNT
      VON0=VON
      WO0=WO
      BIS0=BIS
      HIER0=HIER
      LINKS0=LINKS
      IF (UNT.GT.OB) THEN
         OB0=UNT
         UNT0=OB
         VON0=BIS
         WO0=HIER
         BIS0=VON
         HIER0=WO
         LINKS0=.NOT.LINKS0
      ENDIF
      IF (ACHSE.gt.0) THEN
         XX(1)=VON0
         YY(1)=WO0
         XX(2)=BIS0
         YY(2)=HIER0
         CALL GPL(2,XX,YY)
      ENDIF

      FAX=(PP(3)-PP(1))/(PP(7)-PP(5))
      FAY=(PP(4)-PP(2))/(PP(8)-PP(6))
      HOEHE=(HIER0-WO0)*FAY
      BREITE=(BIS0-VON0)*FAX
      GLAENGE=SQRT(HOEHE**2+BREITE**2)
      FMATEIL=GLAENGE/(5.*ZEIGRO)
      IF (FMATEIL.LT.0.5) GOTO 99
      MAXTEIL=NINT(FMATEIL)
      DIST=ABS(OB0-UNT0)
      OB1=OB0+DIST*1E-5
      UNT1=UNT0-DIST*1E-5
      IDIMUN=LOG10(ABS(UNT1)*1.000001)+1000D0
      IDIMOB=LOG10(ABS(OB1) *1.000001)+1000D0
      IST=MAX(IDIMOB-1000,IDIMUN-1000)
      IF ( IST.GE.-1 .AND. IST.LE.3 ) THEN
         IEXP=0
      ELSE
         IF (IST.LT.-1) THEN
            IEXP=IST+1
         ELSE
            IEXP=IST-3
         ENDIF
      ENDIF
      FAC=10D0**(-IEXP)
      OB1=OB1*FAC
      UNT1=UNT1*FAC
      ISTEL=MAX(0,IST+1-IEXP)
      IF ( IST.LT.0 ) ISTEL=MAX(0,ISTEL-1)
      IF ( UNT0.LT.0. ) ISTEL=ISTEL+1
      VON1 =PP(1)+(VON0-PP(5))*FAX
      WO1  =PP(2)+(WO0-PP(6))*FAY
      BIS1 =PP(1)+(BIS0-PP(5))*FAX
      HIER1=PP(2)+(HIER0-PP(6))*FAY
C---------------- UEBERGANG ZU NICHTVERZERRTEN CENTIMETER-KOORDINATEN
      if (achse.lt.2) then
C------- HOLE ALTE TEXTHOEHE ZUM SICHERN
         CALL GQCHH(J,CHH)
C------- SETZE TEXTHOEHE NEU
         CALL GSCHH(ZEIGRO)
C------- HOLE ALTEN CHARACTER_EXPANSION_FACTOR
         CALL GQCHXP(J,CXP)
C------- SETZE CHARACTER_EXPANSION_FACTOR
         CALL GSCHXP(1.)
C------- HOLE ALTES TEXTALIGNMENT ZUM SICHERN
         CALL GQTXAL(J,IALH,IALL)
C------- SETZE TEXTALIGNMENT LINKS-RECHTS, HALB
C------- BEI  GSTXAL                     GKS - FEHLER (TU BERLIN)
         IF (HOEHE.LT.0D0) THEN
            IF (LINKS0) THEN
               CALL GSTXAL(1,4)
C              CALL GSTXAL(1,7)
            ELSE
               CALL GSTXAL(3,4)
C              CALL GSTXAL(3,7)
            ENDIF
            CALL GSCHUP(-BREITE,-HOEHE)
         ELSE
            IF (LINKS0) THEN
               CALL GSTXAL(3,4)
C              CALL GSTXAL(3,7)
            ELSE
               CALL GSTXAL(1,4)
C              CALL GSTXAL(1,7)
            ENDIF
            CALL GSCHUP(BREITE,HOEHE)
         ENDIF
      endif
C 10.12.90 BUSCH S. COMMON GRPIC
C---- PP(9) = FCTR = 10.5 : XMAXCM=39.5, YMAXCM=28.7
C     XMAXCM=3.7619048*PP(9)
C     YMAXCM=2.7333333*PP(9)
C---- WINDOW AUF NDC (0#0,1#0.7266) ABBILDEN
C12.11.91  altes Fenster
      CALL GQNT(1,IERR,SCLV,SCLC)
      CALL GSWN(1,0.,XMAXCM,0.,YMAXCM)
      ISTELLEN=IST-IEXP
   11    ISTELLEN=ISTELLEN-1
         HUEPF=10D0**ISTELLEN
         INTE=OB1/HUEPF+1000.
         OBEN= (INTE-1000)*HUEPF
         INTE=UNT1/HUEPF+1000.
         UNTEN=( INTE-999)*HUEPF
      IF (NINT(ABS((OBEN-UNTEN)/HUEPF)).LE.MAXTEIL ) GOTO 11
      TEIL=HUEPF*10D0
C --- 20.7.93
      IF (ABS((OBEN-UNTEN)/TEIL).GT.MAXTEIL) TEIL=TEIL*10.
      ITEIL=NINT((OBEN-UNTEN)/TEIL)
      MERKEN=0
      IF (ITEIL*5.LE.MAXTEIL) THEN
         TEIL=TEIL/5
         ITEIL=ITEIL*5
C --- 20.7.93
         IF (ITEIL.EQ.0) ITEIL=NINT((OBEN-UNTEN)/TEIL)
      ENDIF
      IF (ITEIL*2.LE.MAXTEIL) THEN
         TEIL=TEIL/2
         ITEIL=ITEIL*2
C --- 20.7.93
         IF (ITEIL.EQ.0) ITEIL=NINT((OBEN-UNTEN)/TEIL)
      ENDIF
      IF (ITEIL*5.LE.MAXTEIL) THEN
         TEIL=TEIL/5
         ITEIL=ITEIL*5
C --- 20.7.93
         IF (ITEIL.EQ.0) ITEIL=NINT((OBEN-UNTEN)/TEIL)
      ENDIF
      IF (ITEIL*2.LE.MAXTEIL) THEN
         TEIL=TEIL/2
         ITEIL=ITEIL*2
C --- 20.7.93
         IF (ITEIL.EQ.0) ITEIL=NINT((OBEN-UNTEN)/TEIL)
         MERKEN=1
      ENDIF
      HUEPF=TEIL/10D0
      INTE=OB1/HUEPF+1000D0
      OBEN= (INTE-1000)*HUEPF
      INTE=UNT1/HUEPF+1000.
      UNTEN=(INTE-999)*HUEPF
      DIST=OBEN-UNTEN
      IVORHER=ISTEL
      GRENZE=HUEPF/10D0
      HILF1=LOG10(TEIL)+1000.00000001D0
      NACHHER=MAX(0,-INT(HILF1)+1000)+MERKEN
      HILF1=(BIS1-VON1)/(OB1-UNT1)
      HILF2=(HIER1-WO1)/(OB1-UNT1)
      LDRUCK=MIN(IVORHER+NACHHER+1,9)
      WRITE(FORM,'(''(F'',I1,''.'',I1,'')'')') LDRUCK,NACHHER
      DELX1= HOEHE/GLAENGE*ZEIGRO*0.2
      DELY1=-BREITE/GLAENGE*ZEIGRO*0.2
      IF ( LINKS0 ) THEN
         DELX1=-DELX1
         DELY1=-DELY1
      ENDIF
      DELX2=2.*DELX1
      DELY2=2.*DELY1
      I=0
      P1=UNTEN
      P2=UNTEN
      NHUEPF = (OBEN+DIST*1D-5-UNTEN)/HUEPf
      DO 100 IHUEPF = 0,NHUEPF
         P = UNTEN + HUEPF*IHUEPF
         XX(1)=VON1+(P-UNT1)*HILF1
         YY(1)=WO1 +(P-UNT1)*HILF2
         XX(2)=XX(1)+DELX1
         YY(2)=YY(1)+DELY1
         CALL GPL(2,XX,YY)
C --- 20.7.93         10000 statt 1000
         INTE=(P+GRENZE)/TEIL+10000D0
         PMODUL=ABS((INTE-10000)*TEIL-P)
         IF (PMODUL.LT.GRENZE) THEN
            XX(1)=VON1+(P-UNT1)*HILF1
            YY(1)=WO1 +(P-UNT1)*HILF2
            XX(2)=XX(1)+DELX2
            YY(2)=YY(1)+DELY2
            CALL GPL(2,XX,YY)
            if (achse.lt.2) then
               WRITE(DRUCK,FORM) P
               CALL GTX(XX(2)+DELX1,YY(2)+DELY1,DRUCK(:LDRUCK))
            endif
            I=I+1
            IF (I.EQ.1) P1=P
            IF (I.EQ.2) P2=P
         ENDIF
 100  CONTINUE
      IF ( IEXP.NE.0 .and. achse.lt.2 ) THEN
         PM=(P1+P2)/2
         DELX3=2.5*DELX2
         DELY3=2.5*DELY2
C -- 14.06.93 Raus-Schieben bei Enge   (Groten)
         ENGFA=1.
         IF (GLAENGE/ITEIl.LT.1.) ENGFA=4.
         X=VON1+(PM-UNT1)*HILF1+DELX3*ENGFA
         Y=WO1 +(PM-UNT1)*HILF2+DELY3*ENGFA
         CALL GTX(X,Y,ZEHN)
         IF (LINKS0) THEN
            IF (HOEHE.GE.0.) THEN
               X=X+DELY3
               Y=Y-DELX3
            ELSE
               X=X-DELY3
               Y=Y+DELX3
            ENDIF
         ELSE
            IF (HOEHE.GE.0.) THEN
               X=X-DELY3
               Y=Y+DELX3
            ELSE
               X=X+DELY3
               Y=Y-DELX3
            ENDIF
         ENDIF
C------- SETZE TEXTHOEHE NEU
         CALL GSCHH(ZEIGRO*.625)
         WRITE(POT(6:8),'(I3)') IEXP
         CALL GTX(X,Y,POT)
      ENDIF
C---- UMRECHNUNG DES WINDOW AUF DEN GESAMTEN VIEWPORTBEREICH
C12.11.91 XLI=PP(5)-PP(1)/FAX
C12.11.91 YLI=PP(6)-PP(2)/FAY
C12.11.91 XRE=PP(7)+(XMAXCM-PP(3))/FAX
C12.11.91 YRE=PP(8)+(YMAXCM-PP(4))/FAY
C12.11.91 CALL GSWN(1,xli,xre,yli,yre)
C---- WINDOW AUF NDC (0#0,1#0.7266) ABBILDEN (USERKOORDINATEN)
C12.11.91
      CALL GSWN(1,SCLV(1),SCLV(2),SCLV(3),SCLV(4))
      if ( achse.lt.2) then
C------- SETZE ALTEN CHARACTER_EXPANSION_FACTOR ZURUECK
         CALL GSCHXP(CXP)
C------- SETZE ALTE TEXTHOEHE ZURUECK
         CALL GSCHH(CHH)
C------- SETZE ALTES TEXTALIGNMENT ZURUECK
         CALL GSTXAL(IALH,IALL)
      endif
   99 END
