C@PROCESS OPT(3) NOSDUMP NOGOSTMT FIPS(F)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C ------------------ MIT OPT(2) UEBERSETZEN ! --------------------------
C ----G R D R D M-------------------------------------------------------
C        AENDERUNG: 31.01.85 BEMASSUNG ---------------------------------
C        AENDERUNG: 12.02.85 GRDRKU UND AEQUI --------------------------
C        AENDERUNG: 13.02.85 GRDRNE ------------------------------------
C        AENDERUNG: 21.02.85 VORDERKANTE DURCHGEZOGEN-------------------
C        AENDERUNG: 26.02.85 DUENNE LINIE BEI KANTEN (IMMER) -----------
C        AENDERUNG:  7.03.85 GRDRKU ALS SYMBOL-LINIE -------------------
C        AENDERUNG: 31.07.86 GRDRAX ABGETRENNT -------------------------
C        AENDERUNG: 28.11.86 COMMON FUER GRSCAX ------------------------
C        UPDATE   :  6.11.90 SAVE  M.BUSCH
C        Update   : 12 .2 93 Busch GRJMPS durch GRPTS ersetzt
C        Update   : 18 .3 94 Groten : wegen f90 ASSIGN entfernt
C        Update   : 27.10.94 Groten : bei grfrbn : Farbe+Dicke*2000
C ----------------------------------------------------------------------
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE GRDRDM(PA,NROW,TAB,XX,YY)
      IMPLICIT NONE
C     .. Previously untyped names ..
      INTEGER IF1,IF2,IF3,IF4,IF5,NP,NROW,ll,kk
      REAL PA,PK3,TAB,XX,XYZ,YY,BREI,CAX,CAY,COP,FA,SAX,SAY,SIP,SIT,XAX,
     +     XAY,XEX,XEY,YAX,YAY,YEX,YEY,ZA1X,ZA2X,ZE2Y,ZETMA0,ZETMI0,
     +     ZETMI1
      INTEGER IRIFA
      REAL CILBER,PP,AUR,AX,AY,BA,BR,BR2,BREITE,BRH,CH,CHI,COT,DELZ,DET,
     +     DETE,DETZ,DL,DR,DSH1,DSH2,DSH3,DX,DY,DZ,EPSX,EPSY,EPSZ,ET,
     +     F00,F01,F02,F03,F10,F11,F12,PHI,TH,THETA,W3,X0,X1,X2,X3,XA,
     +     XABR,XALT,XE,XI,XI0,XIMA,XIMI,XIST,XL,XMA,XMAX,XMI,XMIN,XO,
     +     XR,XYMI,YABR,YL,YMA,YMAX,YMI,YMIN,YO,YXMI,Z0,Z1,Z1S,Z2,
     +     Z2S,Z3,ZA,ZALT,ZAS,ZE,ZES,ZET,ZET0,ZETA,ZETMI,ZM,ZMAX,ZMI,
     +     ZMIN,ZO,ZVR,X,Y,Z
      INTEGER I,I0,IALT,IB,IK,IL,INLIOL,INTLIN,INTSYM,IPEN,IQ,IQP1,
     +        IRETC,IX,J,J0,JK,JL,K,KANTEN,KANTSC,L,L1,L2,L3,L4,L5,
     +        L6,LBX,LBY,LM1,LMAX,LMP,LS,LSFA,LSFE,LSLA,LSLE,NB,NI,
     +        NJ,NK,NL,NRET,NSY,icola,icolor,lidi,ifa,ifont,jfont,ifo
      REAL ECKWUE,SICHT,WUERF,XP,YP
      INTEGER LICHT

    
C     ..
CCC   COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
CCC
CCC
      COMMON /GRPP/PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      EQUIVALENCE (PP(16),INTSYM), (PP(13),INTLIN),(PP(16),ICOLOR),
     +            (PP(9),IFONT)
CCC
CCC
      COMMON /GRCIL/CILBER(8)
CDEC$ PSECT /GRCIL/ NOSHR
      COMMON /GRDRA/SIP,SIT,COP,FA,IRIFA
CDEC$ PSECT /GRDRA/ NOSHR
      COMMON /GRSCAC/XAX,XAY,XEX,XEY,YAX,YAY,YEX,YEY,ZA1X,ZETMI0,ZETMI1,
     +       ZA2X,ZETMA0,ZE2Y,SAX,CAX,SAY,CAY,BREI
CDEC$ PSECT /GRSCAC/ NOSHR
      DIMENSION PA(16),TAB(NROW,*),XX(*),YY(*),PK3(3,NP),XYZ(3,NROW,*)
CDEB  DIMENSION PA(16),TAB(NROW,NROW),XX(NROW),YY(NROW),PK3(3,NP),
CDEB $          XYZ(3,NROW,NROW)
      DIMENSION SICHT(4,1536),LICHT(4,1536),XP(3072),YP(3072)
      DIMENSION ECKWUE(3,8),WUERF(2,8)
      INTEGER FLAWUE(3,6),KANWUE(4,12)
      LOGICAL VORWUE(6)
      EQUIVALENCE (SICHT(1,1),LICHT(1,1))
      EQUIVALENCE (SICHT(1,1),XP(1)), (SICHT(1,768),YP(1))
      INTEGER HBX,HBY,FUNKTI,RAHMEN,RICHSC
      LOGICAL SS,SX,OBEN,KURVE,AEQUI,IRICHT,JRICHT,NETZ
      CHARACTER CZ*10
      SAVE /GRCIL/,/GRDRA/,/GRSCAC/, /GRPP/
C     SAVE OBEN
      DATA ECKWUE/0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,0.,1.,
     +     1.,1.,1.,0.,1.,1./
      DATA KANWUE/1,2,1,4,2,3,4,5,4,3,3,4,1,4,4,6,5,6,1,2,6,7,2,5,8,7,2,
     +     3,5,8,2,6,2,6,1,5,3,7,3,5,4,8,3,6,1,5,1,6/
      DATA FLAWUE/1,2,5,5,6,8,8,7,4,4,3,1,2,3,6,1,5,4/
      DATA FUNKTI/1/,KANTEN/2/,KANTSC/1/,RAHMEN/3/,RICHSC/1/

C

C
C PA(1)=2.       EIN SOCKEL SOLL GEZEICHNET WERDEN
C      =1.       EIN RAHMEN SOLL GEZEICHNET WERDEN
C      =.03125   DER KANTENWUERFEL SOLL GEZEICHNET WERDEN
C      =.015625  DER KANTENWUERFEL SOLL SKALIERT   WERDEN
C      =.0078125 BEIM SKALIEREN SOLLEN DIE X-Y-ACHSEN IN 10 TEILE ZERF.
C       KOMBINIEREN DURCH ADDIEREN
C PA(5)          RETCODE
C      = 0: KEIN FEHLER FESTGESTELLT
C      = 1: ABS(PHI)>360
C      = 2: THETA<-90 ODER THETA>90
C      = 3: BREITE<=0 ODER BREITE>200
C      =10: LBX < 1
C      =11: HBX > NROW
C      =12: HBX-LBX<1
C      =13: XX NICHT MONOTON STEIGEND
C      =14: XMIN GROESSER ALS XX( LBX )
C      =15: XMAX KLEINER ALS XX( HBX )
C      =16: XMAX=XMIN
C      =17: LBX ODER HBX NICHT KONSISTENT MIT DIMENS. VON TAB
C      =20-27: WIE 10-17 ENTSPRECHEND FUER Y-KOORDINATE
C      =30: ZMAX<=ZMIN
C      =100: DER HORIZONT IST ZU LANG
C    (XMIN,XMAX)      WERTEBEREICH DER X-ACHSE
C    (YMIN,YMAX)      WERTEBEREICH DER Y-ACHSE
C    (ZMIN,ZMAX)      WERTEBEREICH DER TABELLE
C    (LBX,HBX)        KLEINSTER UND GROESSTER INDEX IN X
C    (LBY,HBY)        KLEINSTER UND GROESSTER INDEX IN Y
C    TAB(NROW,*)      TABELLE DER ZU ZEICHNENDEN FUNKTION
C    XX(*)            RASTER AUF DER X-ACHSE
C    YY(*)            RASTER AUF DER Y-ACHSE
C    (XO,YO,ZO)       MITTELPUNKT DES ZU PROJEZ. KUBUS
C    (DX,DY,DZ)       WERTE-SPANNE,WELCHE DER KUBUS UMFASST
C    SICHT(4,*)=LICHT(4,*)    SPEICHER FUER DEN HORIZONT
C    LICHT(1,*)       ZEIGER AUF DEN NAECHSTEN PUNKT
C    LICHT(2,*)       ZEIGER AUF DEN VORIGEN PUNKT
C    SICHT(3,*)       XI-WERT
C    SICHT(4,*)       ZETA-WERT
C    LSFA             ZEIGER AUF DEN ANF. DER FREIEN KETTE
C    LSFE             ZEIGER AUF DAS ENDE DER FREIEN KETTE
C    LSLA             ZEIGER AUF DEN ANF. DES HORIZONTES
C    LSLE             ZEIGER AUF DAS ENDE DES HORIZONTES
C    L3               ENDE DES RECHTEN HORIZONT-TEILES
C    L4               ANFANG DES NEUEN ZWISCHENSTUECKES
C    L5               ENDE DES NEUEN ZWISCHENSTUECKES
C    L6               ANFANG DES LINKEN HORIZONT-TEILES
C    (L1,L2,LS)       HILFSZEIGER
C
      OBEN = .TRUE.
      GO TO 70707
C ----------------------------------------------------------------------
C ----G R F R B N-------------------------------------------------------
C ----------------------------------------------------------------------
      ENTRY GRFRBN(IF1,IF2,IF3,IF4,IF5)

      IF (IF1.NE.666) FUNKTI = IF1
      IF (IF2.NE.666) KANTEN = IF2
      IF (IF3.NE.666) KANTSC = IF3
      IF (IF4.NE.666) RAHMEN = IF4
      IF (IF5.NE.666) RICHSC = IF5
      IRIFA = RICHSC
C---- RETURN
      GO TO 10101
C ----------------------------------------------------------------------
C ----G R D R K U-------------------------------------------------------
C ----------------------------------------------------------------------
      ENTRY GRDRKU(PA,NP,PK3)

      IF (NP.GT.3072) THEN
          PA(5) = 200.
          GO TO 9999
      END IF
      LBX = 1
      HBX = NP
      KURVE = .TRUE.
      AEQUI = .TRUE.
      BA = PA(1)
      NETZ = .FALSE.
      GO TO 80808
C ----------------------------------------------------------------------
C ----G R D R N E-------------------------------------------------------
C ----------------------------------------------------------------------
      ENTRY GRDRNE(PA,NROW,XYZ)

      IF (PA(16).NE.1. .AND. PA(16).NE.2.) THEN
          IRICHT = .TRUE.
          JRICHT = .TRUE.
      ELSE
          IRICHT = PA(16) .EQ. 1.
          JRICHT = PA(16) .EQ. 2.
      END IF
      AEQUI = .TRUE.
      NETZ = .TRUE.
      BA = PA(1)
      GO TO 60606
C ----------------------------------------------------------------------
C ----G R D R D U-------------------------------------------------------
C ----------------------------------------------------------------------
      ENTRY GRDRDU(PA,NROW,TAB,XX,YY)

      OBEN = .FALSE.
70707 NETZ = .FALSE.
      AEQUI = .TRUE.
      BA = PA(1)
      X = 64.*BA
      IX = X
      IF ((X-IX).LT.0.5) AEQUI = .FALSE.
60606 LBX = PA(12)
      HBX = PA(13)
      LBY = PA(14)
      HBY = PA(15)
      KURVE = .FALSE.
80808 PHI = PA(2)
      THETA = PA(3)
      BREITE = PA(4)
      XMIN = PA(6)
      XMAX = PA(7)
      YMIN = PA(8)
      YMAX = PA(9)
      ZMIN = PA(10)
      ZMAX = PA(11)
      IRIFA = RICHSC
      GO TO 3000
C
C #LINE: BEWEGEN DER ZEICHENFEDER
C
 1000 IF (ZETA.LE.BREITE) GO TO 1001
      ZETA = BREITE
      GO TO 1002
 1001 IF (ZETA.LT.0.) ZETA = 0.
 1002 IF (XI.EQ.XALT .AND. ZETA.EQ.ZALT) GO TO 2000
      IF (IPEN.EQ.3 .AND. IALT.EQ.3) GO TO 1010
      IF (IALT.NE.3) GO TO 1003
 1111 CALL GRJMP(XALT,ZALT)
      GO TO 1010
 1003 CALL GRDRW(XALT,ZALT)
 1010 XALT = XI
      ZALT = ZETA
      IALT = IPEN
C
C #MERK: VERMERKEN EINES PUNKTES IM HORIZONT
C
 2000 IF (LICHT(1,LSFA).NE.0) GO TO 2001
      PA(5) = 100.
      GO TO 9998
 2001 L1 = LS
      LS = LSFA
      LSFA = LICHT(1,LS)
      IF (L5.NE.0) GO TO 2002
      L5 = LS
      L4 = LS
      LICHT(2,LS) = 0
      GO TO 2003
 2002 LICHT(1,L5) = LS
      LICHT(2,LS) = L5
      L5 = LS
 2003 LICHT(1,LS) = 0
      LICHT(2,LSFA) = 0
      SICHT(3,LS) = XI
      SICHT(4,LS) = ZETA
      LS = L1
      GO TO (4000,4001,4002,4003,4004,4005,4012,4014,4015,4021,4022,
     +       6011,6030,6101,6120,6190,6225,6341,6500,6700,6760,
     +       6770) NRET
C
C #ANF:  PRUEFEN DER EINGABEDATEN
C
 3000 PA(5) = 0.
      CHI = PHI
      DO 3001 I = 45,315,90
          IF (ABS(CHI).EQ.I) CHI = CHI + .0166667
 3001 CONTINUE
      IF (ABS(CHI).LE.360.) GO TO 3004
      PA(5) = 1.
      GO TO 9999
 3004 IF (THETA.GE.-90. .AND. THETA.LE.90.) GO TO 3005
      PA(5) = 2.
      GO TO 9999
 3005 IF (BREITE.GT.0. .AND. BREITE.LE.200.) GO TO 3012
      PA(3) = 3.
      GO TO 9999
 3012 IF (XMAX.GT.XMIN) GO TO 3019
      PA(5) = 16.
      GO TO 9999
 3019 IF (YMAX.GT.YMIN) GO TO 3020
      PA(5) = 26.
      GO TO 9999
 3020 IF (ZMAX.GT.ZMIN) GO TO 3002
      PA(5) = 30.
      GO TO 9999
 3002 IF (KURVE) GO TO 3100
      IF (LBX.GE.1) GO TO 3008
      PA(5) = 10.
      GO TO 9999
 3008 IF (HBX.LE.NROW) GO TO 3013
      PA(5) = 11.
      GO TO 9999
 3013 IF (LBX.GE.1 .AND. HBX.LE.NROW) GO TO 3014
      PA(5) = 17.
      GO TO 9999
 3014 IF (LBY.GE.1) GO TO 3015
      PA(5) = 20.
      GO TO 9999
C RETURNCODE 21 IN FORTRAN IMPOSSIBLE TO DETECT
 3015 IF ((HBY-LBY).GE.1) GO TO 3006
      PA(5) = 22.
      GO TO 9999
 3006 IF ((HBX-LBX).GE.1) GO TO 3009
      PA(5) = 12.
      GO TO 9999
 3009 IF (NETZ) GO TO 3100
      J = LBX + 1
      DO 3010 I = J,HBX
          IF (XX(I).GE.XX(I-1)) GO TO 3010
          PA(5) = 13.
          GO TO 9999
 3010 CONTINUE
      IF (XMIN.LE.XX(LBX)) GO TO 3011
      PA(5) = 14.
      GO TO 9999
 3011 IF (XMAX.GE.XX(HBX)) GO TO 3016
      PA(5) = 15.
      GO TO 9999
 3016 J = LBY + 1
      DO 3017 I = J,HBY
          IF (YY(I).GE.YY(I-1)) GO TO 3017
          PA(5) = 23.
          GO TO 9999
 3017 CONTINUE
      IF (YMIN.LE.YY(LBY)) GO TO 3018
      PA(5) = 24.
      GO TO 9999
 3018 IF (YMAX.GE.YY(HBY)) GO TO 3100
      PA(5) = 25.
      GO TO 9999
C
C  FESTLEGEN VON KONSTANTEN UND ANFANGSWERTEN
C
 3100 FA = ATAN(1.)/45.
      icola = icolor
      INLIOL = INTLIN
      jfont = ifont
      W3 = SQRT(3.)
      XO = .5* (XMIN+XMAX)
      DX = XMAX - XMIN
      YO = .5* (YMIN+YMAX)
      DY = YMAX - YMIN
      ZO = .5* (ZMIN+ZMAX)
      DZ = ZMAX - ZMIN
      TH = FA*THETA
      CH = FA*CHI
      SIT = SIN(TH)
      SIP = SIN(CH)
      COT = COS(TH)
      COP = COS(CH)
      BREI = BREITE/W3
      F00 = (.5+ ((XO/DX*SIP-YO/DY*COP)*SIT-ZO/DZ*COT)/W3)*BREITE
      F01 = -SIP*SIT/DX*BREI
      F02 = COP*SIT/DY*BREI
      F03 = COT/DZ*BREI
      F10 = (.5- (YO/DY*SIP+XO/DX*COP)/W3)*BREITE
      F11 = COP/DX*BREI
      F12 = SIP/DY*BREI
      XABR = PP(1) + BREITE
      YABR = PP(2) + BREITE
      CALL GRSCLC(PP(1),PP(2),XABR,YABR)
      CALL GRSCLV(0.,0.,BREITE,BREITE)
      IF (KURVE .OR. NETZ) GO TO 3106
      L3 = 0
      L4 = 0
      L5 = 0
      L6 = 0
      LICHT(2,1) = 0
      LICHT(2,2) = 1
      LSFA = 1
      DO 3111 KK = 2,1535
          LICHT(1,KK-1) = KK
          LICHT(2,KK+1) = KK
 3111 CONTINUE
      LSFE = 1536
      LICHT(1,1535) = 1536
      LICHT(1,1536) = 0
      LSLA = 0
      LSLE = 0
      XALT = 0.
      ZALT = 0.
      IALT = 3
      NI = HBX - LBX + 1
      NJ = HBY - LBY + 1
      IF (SIP.LT.0.) GO TO 3101
      IQ = 1
      IF (COP.GT.0.) IQ = 0
      GO TO 3102
 3101 IQ = 3
      IF (COP.LT.0.) IQ = 2
 3102 IF (SIP.EQ.0. .AND. COP.LT.0.) IQ = 3
      IQP1 = IQ + 1
      GO TO (3107,3103,3104,3105) IQP1
 3107 I0 = LBX - 1
      IK = 1
      IL = 0
      NK = NI
      J0 = LBY - 1
      JK = 0
      JL = 1
      NL = NJ
      GO TO 3106
 3103 I0 = HBX + 1
      IK = 0
      IL = -1
      NK = NJ
      J0 = LBY - 1
      JK = 1
      JL = 0
      NL = NI
      GO TO 3106
 3104 I0 = HBX + 1
      IK = -1
      IL = 0
      NK = NI
      J0 = HBY + 1
      JK = 0
      JL = -1
      NL = NJ
      GO TO 3106
 3105 I0 = LBX - 1
      IK = 0
      IL = 1
      NK = NJ
      J0 = HBY + 1
      JK = -1
      JL = 0
      NL = NI
C
C          ZEICHNEN DES RAHMENS
C
 3106 IX = .5*BA
      IF ((BA-2*IX).LT.0.5) GO TO 3200
      lidi = RAHMEN/2000
      ifa  =  RAHMEN-lidi*2000
      CALL GRSPTS(lidi+17)
      CALL GRNWPN(ifa)
      DSH1 = PP(10)
      DSH2 = PP(11)
      DSH3 = PP(12)
      CALL GRDSH(1.,0.,1.)
      CALL GRJMP(0.,0.)
      CALL GRDRW(BREITE,0.)
      CALL GRDRW(BREITE,BREITE)
      CALL GRDRW(0.,BREITE)
      CALL GRDRW(0.,0.)
      CALL GRDSH(DSH1,DSH2,DSH3)
C
C        ZEICHNEN DES KANTENWUERFELS
C
 3200 DO 3290 I = 1,8
          X = XMIN + DX*ECKWUE(1,I)
          Y = YMIN + DY*ECKWUE(2,I)
          Z = ZMIN + DZ*ECKWUE(3,I)
          WUERF(1,I) = FXI(X,Y)
          WUERF(2,I) = FZETA(X,Y,Z)
 3290 CONTINUE
      X = 16.*BA
      IX = X
      IF ((X-IX).LT.0.5) GO TO 3300
      lidi = KANTEN/2000
      ifa  = KANTEN-lidi*2000
      CALL GRSPTS(lidi+17)
      CALL GRNWPN(ifa)
      DSH1 = PP(10)
      DSH2 = PP(11)
      DSH3 = PP(12)
      DO 3291 I = 1,6
          X1 = WUERF(1,FLAWUE(2,I)) - WUERF(1,FLAWUE(1,I))
          X2 = WUERF(1,FLAWUE(3,I)) - WUERF(1,FLAWUE(1,I))
          Z1 = WUERF(2,FLAWUE(2,I)) - WUERF(2,FLAWUE(1,I))
          Z2 = WUERF(2,FLAWUE(3,I)) - WUERF(2,FLAWUE(1,I))
          VORWUE(I) = (Z1*X2-X1*Z2) .LE. 0.
 3291 CONTINUE
      DO 3210 I = 1,12
          CALL GRJMP(WUERF(1,KANWUE(1,I)),WUERF(2,KANWUE(1,I)))
          IF (VORWUE(KANWUE(3,I)) .OR. VORWUE(KANWUE(4,I))) THEN
              CALL GRDSH(1.,0.,1.)
          ELSE
              CALL GRDSH(.5,.5,.5)
          END IF
          CALL GRDRW(WUERF(1,KANWUE(2,I)),WUERF(2,KANWUE(2,I)))
 3210 CONTINUE
      CALL GRDSH(DSH1,DSH2,DSH3)
C
C      MARKIEREN UND BESCHRIFTEN DER KANTEN (AUCH FUER GRSCAX)
C
 3300 XIMI = FXI(XMIN,YMIN)
      XIMA = XIMI
      DO 3303 I = 1,2
          Y = YMIN + DY* (I-1)
          DO 3302 J = 1,2
              X = XMIN + DX* (J-1)
              XI = FXI(X,Y)
              IF (XI.LE.XIMI) THEN
                  XIMI = XI
                  XMI = X
                  YMI = Y
              END IF
              IF (XI.GE.XIMA) THEN
                  XIMA = XI
                  XMA = X
                  YMA = Y
              END IF
 3302     CONTINUE
 3303 CONTINUE
      ZETMI0 = FZETA(XMI,YMI,ZMIN)
      ZMI = ZMIN + DZ
      ZETMI1 = FZETA(XMI,YMI,ZMI)
      ZETMA0 = FZETA(XMA,YMA,ZMIN)
      AUR = PP(15)
      BR = PP(14)
      BRH = BR*.75
      BR2 = BR*2.
      ZA1X = XIMI - BR2
      ZA2X = XIMA + BR2
      ZE2Y = ZETMA0 + (ZETMI1-ZETMI0)
      ZETMI = FZETA(XMIN,YMIN,ZMIN)
      YXMI = YMIN
      EPSX = .01* (YMAX-YMIN)
      ZETA = FZETA(XMIN,YMAX,ZMIN)
      IF (ZETA.LT.ZETMI) THEN
          YXMI = YMAX
          EPSX = -EPSX
      END IF
      XAX = FXI(XMIN,YXMI-EPSX)
      XAY = FZETA(XMIN,YXMI-EPSX,ZMIN)
      XEX = FXI(XMAX,YXMI-EPSX)
      XEY = FZETA(XMAX,YXMI-EPSX,ZMIN)
      AX = -SIP*SIT
      IF (AX.NE.0. .OR. COP.NE.0.) THEN
          AX = ATAN2(AX,COP)/FA
      END IF
      AX = AX + 90.
      IF (SIP.LT.0.) AX = AX - 180.
      SAX = SIN(AX*FA)
      CAX = COS(AX*FA)
      DETE = BR2*CAX
      DETZ = BR2*SAX
      IF (SAX*CAX.GE.0.) THEN
          DETE = -DETE
          DETZ = -DETZ
      END IF
      XAX = XAX + DETE
      XAY = XAY + DETZ
      XEX = XEX + DETE
      XEY = XEY + DETZ
      XYMI = XMIN
      EPSY = .01* (XMAX-XMIN)
      ZETA = FZETA(XMAX,YMIN,ZMIN)
      IF (ZETA.LT.ZETMI) THEN
          XYMI = XMAX
          EPSY = -EPSY
      END IF
      YAX = FXI(XYMI-EPSY,YMIN)
      YAY = FZETA(XYMI-EPSY,YMIN,ZMIN)
      YEX = FXI(XYMI-EPSY,YMAX)
      YEY = FZETA(XYMI-EPSY,YMAX,ZMIN)
      AY = COP*SIT
      IF (AY.NE.0. .OR. SIP.NE.0.) THEN
          AY = ATAN2(AY,SIP)/FA
      END IF
      AY = AY + 90.
      IF (COP.GT.0.) AY = AY - 180.
      SAY = SIN(AY*FA)
      CAY = COS(AY*FA)
      DETE = BR2*CAY
      DETZ = BR2*SAY
      IF (SAY*CAY.GE.0.) THEN
          DETE = -DETE
          DETZ = -DETZ
      END IF
      YAX = YAX + DETE
      YAY = YAY + DETZ
      YEX = YEX + DETE
      YEY = YEY + DETZ
C
C      MARKIEREN UND BESCHRIFTEN DER Z-KANTE
C
      X = 32.*BA
      IX = X
      IF ((X-IX).LT.0.5) GO TO 3999
      CALL GRCHRC(BR,0.,INTSYM)
      DELZ = (ZETMI1-ZETMI0)/10.
      EPSZ = .1
      ZETA = ZETMI0
      XI = XIMI
      lidi = KANTSC/2000
      ifa  =  KANTSC-lidi*2000
      IF ( lidi .eq. 0) THEN
        ifo = jfont
      ELSE IF( lidi .eq. 1) THEN
        ifo = -2
      ELSE
        ifo = -51
      ENDIF
      CALL GRFONT(ifo)
      CALL GRNWPN(ifa)
      LMAX = 0
      CALL GSTXAL(3,0)
      DET = -BRH
      DO 3305 I = 1,2
          ZVR = 0.
          DO 3304 J = 0,10
              CALL GRJMP(XI,ZETA)
              XIST = XI - EPSZ
              CALL GRDRW(XIST,ZETA)
              IF (ZETA.GE.ZVR) THEN
                  Z = (ZMIN* (10-J)+ZMAX*J)/10.
                  CALL GRFTOC(Z,CZ,L)
                  LMAX = MAX(LMAX,L)
                  ET = XIST + DET
                  CALL GRTXT(ET,ZETA,L,CZ)
                  ZVR = ZETA + BR2
              END IF
              ZETA = ZETA + DELZ
 3304     CONTINUE
          EPSZ = -EPSZ
          DET = BRH
          XI = XIMA
          ZETA = ZETMA0
          CALL GSTXAL(1,0)
 3305 CONTINUE
      ZA1X = ZA1X - BR*LMAX*.7 - BRH
      ZA2X = ZA2X + BR*LMAX*.7 + BRH
C
C      MARKIEREN UND BESCHRIFTEN DER X-KANTE
C
      CALL GRCHRC(BR,AX,INTSYM)
      ZVR = 0.
      IF (AEQUI) THEN
          IB = 0
          NB = 10
          X = XMIN
      ELSE
          IB = LBX
          NB = HBX
          X = XX(LBX)
      END IF
      XI0 = FXI(X,YXMI)
      ZET0 = FZETA(X,YXMI,ZMIN)
      IF (SAX*CAX*TH.GE.0.) THEN
          DETE = -BRH*CAX
          DETZ = -BRH*SAX
          CALL GSTXAL(3,0)
      ELSE
          DETE = BRH*CAX
          DETZ = BRH*SAX
          CALL GSTXAL(1,0)
      END IF

      Y = YXMI - EPSX
      LMAX = 0
      DO 3307 I = IB,NB
          IF (AEQUI) THEN
              X = (XMIN* (10-I)+XMAX*I)/10.
          ELSE
              X = XX(I)
          END IF
          XI = FXI(X,YXMI)
          ZETA = FZETA(X,YXMI,ZMIN)
          CALL GRJMP(XI,ZETA)
          XL = SQRT((XI-XI0)**2+ (ZETA-ZET0)**2)
          XI = FXI(X,Y)
          ZETA = FZETA(X,Y,ZMIN)
          CALL GRDRW(XI,ZETA)
          IF (XL.GE.ZVR) THEN
              ZVR = XL + BR2
              CALL GRFTOC(X,CZ,L)
              LMAX = MAX(LMAX,L)
              ET = XI + DETE
              ZET = ZETA + DETZ
              CALL GRTXT(ET,ZET,L,CZ)
          END IF
 3307 CONTINUE
      LMP = 2
      IF (XAX.LE.XEX) LMP = 2
      IF (SAX*CAX*TH.GE.0.) THEN
          DETE = - (LMAX+LMP)*BR*CAX*.7
          DETZ = - (LMAX+LMP)*BR*SAX*.7
      ELSE
          DETE = (LMAX+LMP)*BR*CAX*.7
          DETZ = (LMAX+LMP)*BR*SAX*.7
      END IF
      XAX = XAX + DETE
      XAY = XAY + DETZ
      XEX = XEX + DETE
      XEY = XEY + DETZ
C
C      MARKIEREN UND BESCHRIFTEN DER Y-KANTE
C
      CALL GRCHRC(BR,AY,INTSYM)
      ZVR = 0.
      IF (AEQUI) THEN
          IB = 0
          NB = 10
          Y = YMIN
      ELSE
          IB = LBY
          NB = HBY
          Y = YY(LBY)
      END IF
      XI0 = FXI(XYMI,Y)
      ZET0 = FZETA(XYMI,Y,ZMIN)
      X = XYMI - EPSY
      LMAX = 0
      IF (SAY*CAY*TH.GE.0.) THEN
          CALL GSTXAL(3,0)
          DETE = -BRH*CAY
          DETZ = -BRH*SAY
      ELSE
          CALL GSTXAL(1,0)
          DETE = BRH*CAY
          DETZ = BRH*SAY
      END IF
      DO 3309 I = IB,NB
          IF (AEQUI) THEN
              Y = (YMIN* (10-I)+YMAX*I)/10.
          ELSE
              Y = YY(I)
          END IF
          XI = FXI(XYMI,Y)
          ZETA = FZETA(XYMI,Y,ZMIN)
          CALL GRJMP(XI,ZETA)
          YL = SQRT((XI-XI0)**2+ (ZETA-ZET0)**2)
          XI = FXI(X,Y)
          ZETA = FZETA(X,Y,ZMIN)
          CALL GRDRW(XI,ZETA)
          IF (YL.GE.ZVR) THEN
              ZVR = YL + BR2
              CALL GRFTOC(Y,CZ,L)
              LMAX = MAX(LMAX,L)
              ET = XI + DETE
              ZET = ZETA + DETZ
              CALL GRTXT(ET,ZET,L,CZ)
          END IF
 3309 CONTINUE
      LMP = 0
      IF (XAX.LE.XEX) LMP = 2
      IF (SAY*CAY*TH.GE.0.) THEN
          DETE = - (LMAX+LMP)*BR*CAY*.7
          DETZ = - (LMAX+LMP)*BR*SAY*.7
      ELSE
          DETE = (LMAX+LMP)*BR*CAY*.7
          DETZ = (LMAX+LMP)*BR*SAY*.7
      END IF
      YAX = YAX + DETE
      YAY = YAY + DETZ
      YEX = YEX + DETE
      YEY = YEY + DETZ
C
C        ANFANGSFESTLEGUNG DES HORIZONTS
C
 3999 CALL GRCHRC(PP(14),AUR,INTSYM)
      lidi = FUNKTI/2000
      ifa  =  FUNKTI-lidi*2000
      CALL GRSPTS(lidi+17)
      CALL GRNWPN(ifa)
      IF (KURVE) GO TO 9000
      IF (NETZ) GO TO 9100
      XI = BREITE
      ZETA = BREITE
      NRET = 1
      GO TO 2000
 4000 I = IN(NK,NL)
      J = JN(NK,NL)
      XI = FXI(XX(I),YY(J))
      ZETA = FZETA(XX(I),YY(J),TAB(I,J))
      IPEN = 3
      NRET = 2
      GO TO 1000
 4001 IPEN = 2
      IF (BA.LT.2.) GO TO 4010
      X1 = XI
      Z1 = ZETA
      I = IN(NK,NL)
      J = JN(NK,NL)
      XI = FXI(XX(I),YY(J))
      ZETA = FZETA(XX(I),YY(J),ZMIN)
      NRET = 3
      GO TO 1000
 4002 I = IN(NK,1)
      J = JN(NK,1)
      XI = FXI(XX(I),YY(J))
      ZETA = FZETA(XX(I),YY(J),ZMIN)
      NRET = 4
      GO TO 1000
 4003 I = IN(1,1)
      J = JN(1,1)
      XI = FXI(XX(I),YY(J))
      ZETA = FZETA(XX(I),YY(J),ZMIN)
      NRET = 5
      GO TO 1000
 4004 XI = .5*XI
      ZETA = 0.
      NRET = 6
      GO TO 2000
 4005 K = NK
      L = NL - 1
      GO TO 4020
 4010 L = NL - 1
 4011 I = IN(NK,L)
      J = JN(NK,L)
      XI = FXI(XX(I),YY(J))
      ZETA = FZETA(XX(I),YY(J),TAB(I,J))
      NRET = 7
      GO TO 1000
 4012 L = L - 1
      IF (L.GE.1) GO TO 4011
      K = NK - 1
 4013 I = IN(K,1)
      J = JN(K,1)
      XI = FXI(XX(I),YY(J))
      ZETA = FZETA(XX(I),YY(J),TAB(I,J))
      NRET = 8
      GO TO 1000
 4014 K = K - 1
      IF (K.GE.1) GO TO 4013
      XI = .5*XI
      ZETA = 0.
      NRET = 9
      GO TO 2000
 4015 K = NK - 1
      L = 2
      I = IN(NK,2)
      J = JN(NK,2)
      X1 = FXI(XX(I),YY(J))
      Z1 = FZETA(XX(I),YY(J),TAB(I,J))
 4020 XI = 0.
      ZETA = 0.
      NRET = 10
      GO TO 2000
 4021 ZETA = BREITE
      NRET = 11
      GO TO 2000
 4022 LSLA = L4
      LSLE = L5
      L4 = 0
      L5 = 0
      L3 = LICHT(1,LSLA)
      L6 = LICHT(2,LSLE)
      L6 = LICHT(2,L6)
      GO TO 7000
C
C        DURCHGEHEN DES RASTERS
C
 5000 IF (K.NE.NK) GO TO 5100
      IF (L.LE.0) GO TO 5050
      I = IN(NK,L)
      J = JN(NK,L)
      X2 = FXI(XX(I),YY(J))
      Z2 = FZETA(XX(I),YY(J),TAB(I,J))
      X3 = X2
      Z3 = FZETA(XX(I),YY(J),ZMIN)
      L = L - 1
      GO TO 6000
 5050 L = 1
      K = NK - 1
 5051 I = IN(K,1)
      J = JN(K,1)
      X2 = FXI(XX(I),YY(J))
      Z2 = FZETA(XX(I),YY(J),TAB(I,J))
      X3 = X2
      Z3 = FZETA(XX(I),YY(J),ZMIN)
      K = K - 1
      GO TO 6000
 5100 IF (L.NE.1) GO TO 5150
      IF (K.GT.0) GO TO 5051
      K = NK - 1
      L = 2
 5101 I = IN(NK,L)
      J = JN(NK,L)
      X1 = FXI(XX(I),YY(J))
      Z1 = FZETA(XX(I),YY(J),TAB(I,J))
      GO TO 5151
 5150 IF (K.LE.0) GO TO 5180
 5151 I = IN(K,L)
      J = JN(K,L)
      X2 = FXI(XX(I),YY(J))
      Z2 = FZETA(XX(I),YY(J),TAB(I,J))
      LM1 = L - 1
      I = IN(K,LM1)
      J = JN(K,LM1)
      X3 = FXI(XX(I),YY(J))
      Z3 = FZETA(XX(I),YY(J),TAB(I,J))
      K = K - 1
      GO TO 6000
 5180 IF (L.GE.NL) GO TO 9998
      K = NK - 1
      L = L + 1
      GO TO 5101
C
C             DURCHGEHEN EINER MASCHE
C
 6000 LS = LSLA
 6009 LS = LICHT(1,LS)
      IF (SICHT(3,LS).GT.X1) GO TO 6009
      IF (SICHT(3,LS).EQ.X1) LS = LICHT(1,LS)
      XE = SICHT(3,LS)
      ZE = SICHT(4,LS)
      L3 = LICHT(2,LS)
      LS = L3
      XA = SICHT(3,LS)
      ZA = SICHT(4,LS)
      LS = LICHT(1,LS)
      SS = .FALSE.
      SX = .FALSE.
 6010 IF (X1.NE.X2) GO TO 6200
      IF (XA.NE.XE) GO TO 6100
      ZM = AMAX1(ZA,ZE)
      IF ((Z1.LT.ZM.AND.OBEN) .OR. Z1.EQ.ZM .OR.
     +    (Z1.GT.ZM.AND..NOT.OBEN)) THEN
          SS = .FALSE.
      ELSE
          GO TO 6020
      END IF
      IF ((Z2.GT.ZM.AND.OBEN) .OR. (Z2.LT.ZM.AND..NOT.OBEN)) THEN
          IPEN = 3
          XI = X2
          ZETA = ZM
          NRET = 12
          GO TO 1000
      ELSE
          GO TO 6190
      END IF
 6011 IPEN = 2
      ZETA = Z2
      SS = .TRUE.
      NRET = 16
      GO TO 1000
 6020 IF (SS) GO TO 6030
      IPEN = 3
      XI = X1
      ZETA = Z1
      NRET = 13
      GO TO 1000
 6030 IPEN = 2
      XI = X2
      IF ((Z2.LT.ZM.AND.OBEN) .OR. (Z2.GT.ZM.AND..NOT.OBEN)) THEN
          ZETA = ZM
          SS = .FALSE.
          NRET = 16
          GO TO 1000
      END IF
      ZETA = Z2
      SS = .TRUE.
      NRET = 16
      GO TO 1000
 6100 Z0 = (ZA* (XE-X1)+ZE* (X1-XA))/ (XE-XA)
      IF ((Z1.LT.Z0.AND.OBEN) .OR. (Z1.GT.Z0.AND..NOT.OBEN)) THEN
          IF ((Z2.GT.Z0.AND.OBEN) .OR. (Z2.LT.Z0.AND..NOT.OBEN)) THEN
              IPEN = 3
              XI = X2
              ZETA = Z0
              NRET = 14
              GO TO 1000
          ELSE
              GO TO 6109
          END IF
      ELSE
          GO TO 6110
      END IF
 6101 IPEN = 2
      ZETA = Z2
      SS = .TRUE.
      NRET = 16
      GO TO 1000
 6109 SS = .FALSE.
      GO TO 6190
 6110 IF (SS) GO TO 6120
      IPEN = 3
      XI = X1
      ZETA = Z1
      NRET = 15
      GO TO 1000
 6120 IPEN = 2
      XI = X2
      IF ((Z2.GT.Z0.AND.OBEN) .OR. (Z2.LT.Z0.AND..NOT.OBEN)) THEN
          ZETA = Z2
          SS = .TRUE.
          NRET = 16
          GO TO 1000
      END IF
      ZETA = Z0
      SS = .FALSE.
      NRET = 16
      GO TO 1000
 6190 IF (SX) THEN
          L6 = LS
          GO TO 6900
      ELSE
          GO TO 6800
      END IF
 6200 IF (XE.NE.XA) GO TO 6300
      ZM = AMAX1(ZA,ZE)
      Z0 = (Z2* (X1-XA)+Z1* (XA-X2))/ (X1-X2)
      IF ((ZM.LT.Z0.AND.OBEN) .OR. (ZM.EQ.Z0) .OR.
     +    (ZM.GT.Z0.AND..NOT.OBEN)) GO TO 6700
      ZM = AMIN1(ZA,ZE)
      IF ((ZM.GT.Z0.AND.OBEN) .OR. (ZM.EQ.Z0) .OR.
     +    (ZM.LT.Z0.AND..NOT.OBEN)) THEN
          XI = XE
          ZETA = ZE
          NRET = 20
          GO TO 2000
      END IF
      XI = XA
      ZETA = Z0
      IF ((ZA.GT.Z0.AND.OBEN) .OR. (ZA.LT.Z0.AND..NOT.OBEN)) THEN
          IPEN = 3
          SS = .TRUE.
          NRET = 20
          GO TO 1000
      END IF
      IF (.NOT.SS) GO TO 6225
      IPEN = 2
      SS = .FALSE.
      NRET = 17
      GO TO 1000
 6225 XI = XE
      ZETA = ZE
      NRET = 20
      GO TO 2000
 6300 IF (X1.GT.XA) GO TO 6301
      ZAS = (ZA* (XE-X1)+ZE* (X1-XA))/ (XE-XA)
      Z1S = Z1
      XR = X1
      GO TO 6310
 6301 Z1S = (Z1* (X2-XA)+Z2* (XA-X1))/ (X2-X1)
      ZAS = ZA
      XR = XA
 6310 IF (X2.GT.XE) GO TO 6311
      Z2S = (Z1* (X2-XE)+Z2* (XE-X1))/ (X2-X1)
      ZES = ZE
      XL = XE
      GO TO 6320
 6311 ZES = (ZA* (XE-X2)+ZE* (X2-XA))/ (XE-XA)
      Z2S = Z2
      XL = X2
 6320 DL = ZES - Z2S
      DR = ZAS - Z1S
      IF (DL*DR.GT.0.) GO TO 6400
      IF (DL.NE.0. .OR. DR.NE.0.) GO TO 6321
      X0 = XR
      Z0 = ZAS
      GO TO 6330
 6321 X0 = (XR*DL-XL*DR)/ (DL-DR)
      Z0 = (ZAS*DL-ZES*DR)/ (DL-DR)
 6330 IF (X0.GE.XL) GO TO 6331
      X0 = XL
      GO TO 6340
 6331 IF (X0.GT.XR) X0 = XR
 6340 IF ((DL.GT.DR.AND.OBEN) .OR. (DL.EQ.DR) .OR.
     +    (DL.LT.DR.AND..NOT.OBEN)) THEN
          IF (.NOT.SS) GO TO 6345
      ELSE
          GO TO 6350
      END IF
 6341 IPEN = 2
      XI = X0
      ZETA = Z0
      SS = .FALSE.
      NRET = 19
      GO TO 1000
 6345 IPEN = 3
      XI = XR
      ZETA = Z1S
      NRET = 18
      GO TO 1000
 6350 IPEN = 3
      XI = X0
      ZETA = Z0
      SS = .TRUE.
      NRET = 19
      GO TO 1000
 6400 CONTINUE
      IF ((OBEN.AND.DL.LT.0E0.AND.DR.LT.0E0) .OR.
     +    (.NOT.OBEN.AND.DL.GT.0E0.AND.DR.GT.0E0)) THEN
          IF (SS) GO TO 6500
          IPEN = 3
          XI = XR
          ZETA = Z1S
          SS = .TRUE.
          NRET = 19
          GO TO 1000
      ELSE
          SS = .FALSE.
      END IF
 6500 IF (X2.GE.XE) GO TO 6750
      IF (SS) GO TO 6700
      XI = XE
      ZETA = ZE
      NRET = 20
      GO TO 2000
 6700 XA = XE
      ZA = ZE
      LS = LICHT(1,LS)
      XE = SICHT(3,LS)
      ZE = SICHT(4,LS)
      GO TO 6200
 6750 IF (SS) THEN
          IPEN = 2
          XI = X2
          ZETA = Z2
          NRET = 21
          GO TO 1000
      END IF
 6760 IF (SX) THEN
          L6 = LS
          GO TO 6900
      END IF
 6761 IF (X2.NE.XE) GO TO 6800
      IF (SS) GO TO 6770
      XI = XE
      ZETA = ZE
      NRET = 22
      GO TO 2000
 6770 XA = XE
      ZA = ZE
      LS = LICHT(1,LS)
      XE = SICHT(3,LS)
      ZE = SICHT(4,LS)
 6800 X1 = X2
      Z1 = Z2
      X2 = X3
      Z2 = Z3
      SX = .TRUE.
      GO TO 6010
C
C        BILDEN DES NEUEN HORIZONTS
C
 6900 IF (LICHT(1,L3).EQ.L6) GO TO 6901
      LS = LICHT(2,L6)
      LICHT(1,LS) = LSFA
      LICHT(2,LSFA) = LS
      LSFA = LICHT(1,L3)
      LICHT(2,LSFA) = 0
 6901 IF (L4.NE.0) GO TO 6902
      LICHT(1,L3) = L6
      LICHT(2,L6) = L3
      GO TO 6903
 6902 LICHT(1,L3) = L4
      LICHT(2,L4) = L3
      LICHT(1,L5) = L6
      LICHT(2,L6) = L5
 6903 L4 = 0
      L5 = 0
C
C      SAEUBERUNG DES HORIZONTS
C
 7000 X0 = SICHT(3,L6)
      L2 = L3
      L1 = LICHT(1,L2)
      LS = LICHT(2,L2)
 7001 XA = SICHT(3,L1)
      ZA = SICHT(4,L1)
      XE = SICHT(3,L2)
      ZE = SICHT(4,L2)
      IF (XA.NE.XE .OR. (ZA.NE.ZE.AND.XA.NE.SICHT(3,LS))) GO TO 7050
      LICHT(2,LSFA) = L2
      LICHT(1,L2) = LSFA
      LSFA = L2
      LICHT(2,LSFA) = 0
      LICHT(1,LS) = L1
      LICHT(2,L1) = LS
      L2 = L1
      L1 = LICHT(1,L2)
      GO TO 7100
 7050 LS = L2
      L2 = L1
      L1 = LICHT(1,L2)
 7100 IF (XA.LT.X0) THEN
          GO TO 5000
      ELSE
          GO TO 7001
      END IF
C
C      ZEICHNEN EINER KURVE
C
 9000 DO 9001 I = 1,NP
          XP(I) = FXI(PK3(1,I),PK3(2,I))
          YP(I) = FZETA(PK3(1,I),PK3(2,I),PK3(3,I))
 9001 CONTINUE
      IF (PA(12).EQ.999.) THEN
          NSY = PA(13)
c        DO 9002 I=1,NP
c9002       CALL GRJMPS(XP(I),YP(I),NSY)
          CALL GRPTS(XP,YP,NP,NSY)
      ELSE IF (PA(12).EQ.1001.) THEN
          NSY = PA(13)
          CALL GRCHN(XP,YP,NP,NSY)
      ELSE
          CALL GRLN(XP,YP,NP)
      END IF
      GO TO 9990
C
C      ZEICHNEN EINES NETZES
C
 9100 IF (IRICHT) THEN
          DO 9102 J = LBY,HBY
              DO 9101 I = LBX,HBX
                  XP(I) = FXI(XYZ(1,I,J),XYZ(2,I,J))
                  YP(I) = FZETA(XYZ(1,I,J),XYZ(2,I,J),XYZ(3,I,J))
 9101         CONTINUE
              CALL GRLN(XP,YP,HBX-LBX+1)
 9102     CONTINUE
      END IF
      IF (JRICHT) THEN
          DO 9104 I = LBX,HBX
              DO 9103 J = LBY,HBY
                  XP(J) = FXI(XYZ(1,I,J),XYZ(2,I,J))
                  YP(J) = FZETA(XYZ(1,I,J),XYZ(2,I,J),XYZ(3,I,J))
 9103         CONTINUE
              CALL GRLN(XP,YP,HBY-LBY+1)
 9104     CONTINUE
      END IF
      GO TO 9990
C
C      FEHLERNACHRICHTEN
C
 9998 IF (IALT.NE.2) GO TO 9990
      CALL GRDRW(XALT,ZALT)
 9990 IF (PA(5).EQ.100.) GO TO 9999
      CILBER(1) = BREITE
      CILBER(2) = BREITE
      CILBER(8) = BREITE
      CALL GRNWPN(1)
      GO TO 777
 9999 IRETC = PA(5)
      WRITE (*,FMT='('' FEHLER IN GRDRDM,RETURNCODE ='',I4)') IRETC
  777 CALL GSTXAL(0,0)
      CALL GRSPTS(INLIOL)
      CALL GRNWPN(icola)
      CALL GRFONT(jfont)


      contains 
      
      function IN(KK,LL) result(erg)
      integer, intent(in) :: KK,LL
      integer erg
      erg = (I0+IK*KK) + IL*LL
      END function in
      
      function JN(KK,LL) result(erg)
      integer, intent(in) :: KK,LL
      integer erg
      erg =  (J0+JK*KK) + JL*LL
      END function JN

      function   FXI(X,Y) result(erg)
      REAL, intent(in) ::X,Y
      REAL  erg
      erg =   (F10+F11*X) + F12*Y
      END function FXI
    
      function FZETA(X,Y,Z)  result(erg)
      REAL, intent(in) ::X,Y,Z
      REAL  erg
      erg =   ((F00+F01*X)+F02*Y) + F03*Z
      END function FZETA
 
     

 
10101 END SUBROUTINE GRDRDM
