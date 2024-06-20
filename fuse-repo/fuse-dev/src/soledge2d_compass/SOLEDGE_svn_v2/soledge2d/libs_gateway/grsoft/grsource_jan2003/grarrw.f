C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM) FIPS(F)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE 8.8.1990
C UPDATE 23.10.91 Groten: Korrektur fuer jeden Massstab, nicht nur CM.
      SUBROUTINE GRARRW (XANF,YANF,XEND,YEND,ALEN,AWID,ICODE)

      REAL XANF,YANF,XEND,YEND,ALEN,AWID
      INTEGER ICODE

      REAL XA,YA,XE,YE,AL,AW,PL,EX1,EX2,EY1,EY2,XM,YM,X1,X2,Y1,Y2

      REAL PP,XMAXDC,XUNITS,YUNITS
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /SCALE/ XMAXDC,XUNITS,YUNITS
CDEC$ PSECT /SCALE/ NOSHR

      save /GRPP/,/SCALE/
C
C----- ZEICHNEN EINES PFEILES VON (XANF,YANF) NACH (XEND,YEND)
C----- XA,YA,XE,YE : CM-Masseinheiten
C
      XA = PP(1) + (XANF-PP(5))/XUNITS
      YA = PP(2) + (YANF-PP(6))/YUNITS
      XE = PP(1) + (XEND-PP(5))/XUNITS
      YE = PP(2) + (YEND-PP(6))/YUNITS
C
      PL = SQRT ( (XE-XA)**2 + (YE-YA)**2 )
C
      IF ( PL .EQ. 0. ) GOTO 99
      AL = ABS(ALEN)
      AW = ABS(AWID)
C
      CALL GRJMP (XANF,YANF)
      CALL GRDRW (XEND,YEND)
C
      EX1 = (XE-XA)/PL
      EY1 = (YE-YA)/PL
      EX2 = - EY1
      EY2 =   EX1
C
      IF ( PL .LT. AL ) THEN
         XM = XA
         YM = YA
         AW = AW * (PL/AL)
      ELSE
         XM = XA + (PL-AL)*EX1
         YM = YA + (PL-AL)*EY1
      END IF
C
      X1 = XM - .5*AW*EX2
      Y1 = YM - .5*AW*EY2
      X2 = XM + .5*AW*EX2
      Y2 = YM + .5*AW*EY2

C---- Zurueckrechnen in Benutzerkoordinaten

      X1 = PP(5)+(X1-PP(1))*XUNITS
      Y1 = PP(6)+(Y1-PP(2))*YUNITS
      X2 = PP(5)+(X2-PP(1))*XUNITS
      Y2 = PP(6)+(Y2-PP(2))*YUNITS
C
      CALL GRJMP (X1,Y1)
      CALL GRDRW (XEND,YEND)
      CALL GRDRW (X2,Y2)
C
      IF ( ICODE .NE. 0 ) CALL GRDRW (X1,Y1)
C
  99  END
