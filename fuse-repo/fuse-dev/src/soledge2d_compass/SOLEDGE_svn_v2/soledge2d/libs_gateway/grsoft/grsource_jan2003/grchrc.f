C@PROCESS NOSDUMP NOGOSTMT OPT(3)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C Update : 1. 8.91 Busch
      SUBROUTINE GRCHRC(HOEHE,WINKEL,IDUMMY)
      REAL    HOEHE,WINKEL

*KFA  7.8.90 Busch, Groten
*     1) INT nur noch historisch bedingt als Argument
*        hat keine Bedeutung mehr.

      REAL      WIN,PIFAC,ANGEL
      INTEGER IDUMMY
      PARAMETER (PIFAC=6.2831853072/360)
      COMMON /SCALE/ XMAXDC,XUNITS,YUNITS
CDEC$ PSECT /SCALE/ NOSHR
C---- COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR

      SAVE /SCALE/,  /GRPP/

      PP(14)=HOEHE

C---- Der Winkel wird in das Intervall 0-360 gebracht
      ANGEL=MOD(WINKEL+360D0,360D0)

      IF (ANGEL.GT.45..AND.ANGEL.LT.135..OR.ANGEL.GT.225.AND.ANGEL
     >   .LT.315) THEN
C------- WINKEL ZWISCHEN: 45-135 UND 225-315:
         CHH=HOEHE*XUNITS
         CHXP=YUNITS/XUNITS
C------- WINKEL ZWISCHEN: 0-45 UND 135-225 UND 315-360
      ELSE
         CHH=HOEHE*YUNITS
         CHXP=XUNITS/YUNITS
      ENDIF

C---- SET CHARACTER HEIGHT
      CALL GSCHH(CHH)

C---- SET CHARACTER EXPANSION FACTOR
      CALL GSCHXP(CHXP)

      PP(15)=ANGEL


      WIN=ANGEL*PIFAC
C---- SET CHARACTER UP VECTOR
      XWIN=-SIN(WIN)
      YWIN=COS(WIN)
      CALL GSCHUP(XWIN,YWIN)

      END
