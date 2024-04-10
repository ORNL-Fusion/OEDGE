C@process opt(3) nosdump nogostmt fips(f)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE 27. 11.1990
C UPDATE 01. 09.1992 wegen Theta<0 .
      SUBROUTINE GRSCAX(CX,CY,CZ1,CZ2)
      CHARACTER(len=*) CX,CY,CZ1,CZ2

      CHARACTER(len=256) CHELP
      COMMON /GRSCAC/ XAX,XAY,XEX,XEY,YAX,YAY,YEX,YEY,ZA1X,ZA1Y,ZE1Y,
     $                ZA2X,ZA2Y,ZE2Y,SAX,CAX,SAY,CAY,BREI
CDEC$ PSECT /GRSCAC/ NOSHR
      SAVE /GRSCAC/

C --- XAX := X-Achse, Anfang (X-Wert)
C --- XAY := X-Achse, Anfang (Y-Wert)
C --- XEX := X-Achse, Ende   (X-Wert)
C --- XEY := X-Achse, Ende   (Y-Wert)

C --- YAX := Y-Achse, Anfang (X-Wert)
C --- YAY := Y-Achse, Anfang (Y-Wert)
C --- YEX := Y-Achse, Ende   (X-Wert)
C --- YEY := Y-Achse, Ende   (Y-Wert)

C --- ZA1X := Z-Achse 1, Anfang (X-Wert)
C --- ZA1Y := Z-Achse 1, Anfang (Y-Wert)
C --- ZE1Y := Z-Achse 1, Ende   (Y-Wert)

C --- ZA2X := Z-Achse 2, Anfang (X-Wert)
C --- ZA2Y := Z-Achse 2, Anfang (Y-Wert)
C --- ZE2Y := Z-Achse 2, Ende   (Y-Wert)

C --- SAX := Sinus der X-Achsenrichtung
C --- CAX := Cosinus der X-Achsenrichtung
C --- SAY := Sinus der Y-Achsenrichtung
C --- CAY := Cosinus der Y-Achsenrichtung

C --- BREI := Breite des quadratischen Zeichenfeldes

      COMMON /GRDRA/ SIP,SIT,COP,FA,IRIFA
CDEC$ PSECT /GRDRA/ NOSHR
      SAVE /GRDRA/

C --- SIP:= SIN(PHI)          PHI:= Winkel der Drehung der Blickrichtung
C ---                               um die Z-Achse.
C --- COP := COS(PHI)
C --- SIT := SIN(THETA)     THETA:= Winkel der nachfolgenden Drehung der
C --- FA  := PI/180                 Blickrichtung um die neue X-Achse
C --- IRIFA:= Farbe des Richtungsdreibeins.

      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      SAVE  /GRPP/

      LX=LEN(CX)
      DO 1 LX=LX,2,-1
         IF (CX(LX:LX).NE.' ') GOTO 2
   1  CONTINUE
   2  LY=LEN(CY)
      DO 3 LY=LY,2,-1
         IF (CY(LY:LY).NE.' ') GOTO 4
   3  CONTINUE
   4  LZ1=LEN(CZ1)
      DO 5 LZ1=LZ1,2,-1
         IF (CY(LZ1:LZ1).NE.' ') GOTO 6
   5  CONTINUE
   6  LZ2=LEN(CZ2)
      DO 7 LZ2=LZ2,2,-1
         IF (CY(LZ2:LZ2).NE.' ') GOTO 8
   7  CONTINUE
   8  CALL GSCHH(PP(14))
      CALL GSCHUP(-1.,0.)
      BEL=LZ1*PP(14)
      AXL=SQRT((ZA2X-ZA1X)**2+(ZE1Y-ZA1Y)**2)
      IF (BEL.GT.AXL) THEN
         EX=AXL/BEL
      ELSE
         EX=1.
      ENDIF
      CALL GSCHXP(EX)
      CALL GRTXSCN ( CZ1, LEN(CZ1), CHELP ,LHELP)
      CALL GTX(ZA1X,ZA1Y,CHELP(:LHELP))
      BEL=LZ2*PP(14)
      IF (BEL.GT.AXL) THEN
         EX=AXL/BEL
      ELSE
         EX=1.
      ENDIF
      CALL GSCHXP(EX)
      CALL GRTXSCN ( CZ2, LEN(CZ2), CHELP ,LHELP)
      CALL GTX(ZA2X+PP(14),ZA2Y,CHELP(:LHELP))
      SIC=SIP*SIT
      IF (COP.GE.0.) THEN
         XX=XAX
         YX=XAY
         CALL GSCHUP(SIC,COP)
      ELSE
         XX=XEX
         YX=XEY
         CALL GSCHUP(-SIC,-COP)
      ENDIF
      IF (SIT.GE.0) THEN
         FACX=1
         FACY=1
C 1.9.92
         IF (COP.LT.0) FACY= 2.
      ELSE
         FACX=-4.
C 1.9.92
         FACY=-4.
         IF (COP.LT.0) FACY=-5.
      ENDIF

      IF (SIC*COP.GE.0.) THEN
         XX=XX-CAX*PP(14)*FACX
         YX=YX-SAX*PP(14)*FACX
      ELSE
         XX=XX+CAX*PP(14)*FACX
         YX=YX+SAX*PP(14)*FACX
      ENDIF
      BEL=LX*PP(14)
      AXL=SQRT((XEX-XAX)**2+(XEY-XAY)**2)
      IF (BEL.GT.AXL) THEN
         EX=AXL/BEL
      ELSE
         EX=1.
      ENDIF
      CALL GSCHXP(EX)
      CALL GRTXSCN ( CX, LEN(CX), CHELP ,LHELP)
      CALL GTX(XX,YX,CHELP(:LHELP))
      SIC=COP*SIT
      IF (SIP.GE.0.) THEN
         XY=YAX
         YY=YAY
         CALL GSCHUP(-SIC,SIP)
      ELSE
         XY=YEX
         YY=YEY
         CALL GSCHUP(SIC,-SIP)
      ENDIF
      IF (SIC*SIP.LT.0) THEN
         XY=XY-CAY*PP(14)*FACY
         YY=YY-SAY*PP(14)*FACY
      ELSE
         XY=XY+CAY*PP(14)*FACY
         YY=YY+SAY*PP(14)*FACY
      ENDIF
      BEL=LY*PP(14)
      AXL=SQRT((YEX-YAX)**2+(YEY-YAY)**2)
      IF (BEL.GT.AXL) THEN
         EX=AXL/BEL
      ELSE
         EX=1.
      ENDIF
      CALL GSCHXP(EX)
      CALL GRTXSCN ( CY, LEN(CY), CHELP ,LHELP)
      CALL GTX(XY,YY,CHELP(:LHELP))
      CALL GSCHXP(1.)
      END
