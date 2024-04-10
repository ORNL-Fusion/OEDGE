C@PROCESS OPT(3) IL(DIM) NOGOSTMT NOSDUMP
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C----------------------------------------------------------------
C UPDATE 27.9.90 M. BUSCH    Daten werden formatfrei geschrieben
C                            frueher (20A4) d.h. CRAY unvertraeglich
C UPDATE 16.12.90 G.Groten   Pruefen der Daten weggenommen
C                            einige Modernisierung am FORTRAN
C UPDATE 25. 9.91 G.GROTEN
C UPDATE  5. 5.94 G.GROTEN  Korrektur bei zwei Punkten (wird jetzt
C                           Gerade)
C----------------------------------------------------------------
      SUBROUTINE SKURF(X,Y,IST,FORM,Z)
C     AUTOR: GERD GROTEN
      REAL   X(IST),Y(IST),Z(7,IST),A(100),B(100)
      COMMON /GRCIL/ CILBER(8)
CDEC$ PSECT /GRCIL/ NOSHR
      SAVE /GRCIL/
      REAL   DIFREN
      INTEGER IO
      LOGICAL   BX,BY
      PARAMETER( DIFREN=3., IO=91)

      ISK=0
      GOTO 77
C----------------------------------------------------------------
      ENTRY SKURL(X,Y,IST,FORM,Z,KSK)
      ISK=MOD(KSK,4)

   77 BY = MOD(ISK,2).EQ.1
      BX = ISK.GE.2
      FO=ABS(FORM)
C
C     FESTLEGUNG DER INTENSITAET BZW. FARBE
C     UEBERPRUEFEN DER UEBERGABEPARAMETER
C
      ISY=FO/10.
      FO=FO-10*ISY
      IF (ISY.EQ.0) THEN
         ISY=99
         IF (FORM.LT.0.) ISY=ISY+2000
      ENDIF
      JSY=MOD(ISY,2000)
      IF( JSY.GT. 118 .AND. JSY.LT.200    .OR.
     *    JSY.GT. 318 .AND. JSY.LT.400    .OR.
     *    JSY.GT. 518 .AND. JSY.LT.600    .OR.
     *    JSY.GT. 718 .AND. JSY.LT.800    .OR.
     *    JSY.GT. 918 .AND. JSY.LT.1000   .OR.
     *    JSY.GT.1118 .AND. JSY.LT.1200   .OR.
     *    JSY.GT.1318                     .OR.
     *    IST.LT.2                        .OR.
     *    FO .LT.0.01)          GOTO 150
C
C     BESTIMUNG DER MAXIMA UND MINIMA
C
      IF (CILBER(3).EQ.0) THEN
         DX1=X(1)
         DX2=DX1
         DY1=Y(1)
         DY2=DY1
      ELSE
         DX1=CILBER(4)
         DX2=CILBER(5)
         DY1=CILBER(6)
         DY2=CILBER(7)
      ENDIF
      CILBER(3)=CILBER(3)+1

      DO 30 J=1,IST
         DX1=MIN(DX1,X(J))
         DX2=MAX(DX2,X(J))
         DY1=MIN(DY1,Y(J))
         DY2=MAX(DY2,Y(J))
   30    CONTINUE
C
C     BESTIMMUNG DER GEOMETRIE-FAKTOREN
C
      IF (DX2.EQ.DX1.AND.DY2.EQ.DY1) GOTO 150
      CILBER(4)=DX1
      CILBER(5)=DX2
      CILBER(6)=DY1
      CILBER(7)=DY2
      F1=0.
      IF (BX) THEN
         DX1=ALOG10(DX1)
         DX2=ALOG10(DX2)
      ENDIF
      IF (DX2.NE.DX1) F1=200./(DX2-DX1)
      F2=0.
      IF (BY) THEN
         DY1=ALOG10(DY1)
         DY2=ALOG10(DY2)
      ENDIF
      IF (DY2.NE.DY1) F2=200./(DY2-DY1)*FO
C
C     BERECHNUNG DES PARAMETERS KURVENLAENGE
C
      Z(1,1)=0.
      K=IST-1
      Z(2,1)=X(1)
      IF(BX) Z(2,1)=ALOG10(X(1))
      Z(3,1)=Y(1)
      IF(BY) Z(3,1)=ALOG10(Y(1))
      DO 40 I=1,K
         Z(2,I+1)=X(I+1)
         IF (BX) Z(2,I+1)=ALOG10(X(I+1))
         Z(3,I+1)=Y(I+1)
         IF (BY) Z(3,I+1)=ALOG10(Y(I+1))
         U=SQRT(((Z(2,I+1)-Z(2,I ))*F1)**2+((Z(3,I +1)-Z(3,I ))*F2)**2)
         IF(U.EQ.0) GOTO 160
         Z(1,I+1)=Z(1,I)+U
   40 CONTINUE

      IF (Z(1,IST).EQ.0.) GOTO 150
      Z(4,1)=2.
      Z(6,1)=0.
      Z(7,1)=0.
      Z(6,IST)=0.
      Z(7,IST)=0.
      Z(5,1)=0.
      DT2=Z(1,2)
      DX2=(Z(2,2)-Z(2,1))/DT2
      DY2=(Z(3,2)-Z(3,1))/DT2
      DO 50 I=2,K
         DT1=DT2
         DX1=DX2
         DY1=DY2
         DT2=Z(1,I+1)-Z(1,I)
         DX2=(Z(2,I +1)-Z(2,I ))/DT2
         DY2=(Z(3,I +1)-Z(3,I ))/DT2
         DDT=DT1+DT2
         Z(5,I)=DT2/DDT
         Z(4,I)=1.-Z(5,I)
         Z(6,I)=6.*(DX2-DX1)/DDT
         Z(7,I)=6.*(DY2-DY1)/DDT
   50 CONTINUE

      DO 60 I=2,K
         L=I-1
         FAC=Z(4,I)/Z(4,L)
         Z(6,I)=Z(6,I)-FAC*Z(6,L)
         Z(7,I)=Z(7,I)-FAC*Z(7,L)
         Z(4,I)=2.-Z(5,L)*FAC
   60 CONTINUE

      KP2=K+2
      DO 70 I=2,K
         II=KP2-I
         Z(6,II)=(Z(6,II)-Z(5,II)*Z(6,II+1))/Z(4,II)
         Z(7,II)=(Z(7,II)-Z(5,II)*Z(7,II+1))/Z(4,II)
   70 CONTINUE

      IST1=Z(1,IST)
      IF (Z(1,IST)-IST1.GT.0.) IST1=IST1+1
      IST1=(IST1+(IST+1)*DIFREN)/DIFREN
      I=1
      K=0
      L=0
      V1=0.
      U=0.
      DX2=0.
      DY2=0.
   90 IF (U.LE.Z(1,IST)) THEN
         IF (U.GE.Z(1,I) .AND. U.NE.Z(1,IST)) THEN
            I=I+1
            DX1=DX2
            DY1=DY2
            DX2=Z(6,I)/6.
            DY2=Z(7,I)/6.
            V0=V1
            V1=Z(1,I)
            U=V0
            H=V1-V0
            H2=H*H
            DX3=(Z(2,I-1)-DX1*H2)/H
            DY3=(Z(3,I-1)-DY1*H2)/H
            DX4=(Z(2,I)-DX2*H2)/H
            DY4=(Z(3,I)-DY2*H2)/H
            DX1=DX1/H
            DY1=DY1/H
            D1=DX2/H
            D2=DY2/H
         ENDIF
         K=K+1
         L=L+1
         IF (K.GT.100) THEN
            WRITE(IO)IST1,ISY,A,B
            K=1
         ENDIF
         F3=V1-U
         F4=U-V0
         F1=F3*F3
         F2=F4*F4
         A(K)=F3*(DX1*F1+DX3)+F4*(D1*F2+DX4)
         IF (BX) A(K)=10E0**A(K)
         B(K)=F3*(DY1*F1+DY3)+F4*(D2*F2+DY4)
         IF (BY) B(K)=10E0**B(K)
         U=U+DIFREN
         CILBER(4)=MIN(CILBER(4),A(K))
         CILBER(5)=MAX(CILBER(5),A(K))
         CILBER(6)=MIN(CILBER(6),B(K))
         CILBER(7)=MAX(CILBER(7),B(K))
         GOTO 90
      ENDIF

      L=L+1
      IF (L.GT.IST1) GOTO 145
      DO 140 II=L,IST1
         K=K+1
         IF (K.GT.100) THEN
            WRITE(IO)IST1,ISY,A,B
            K=1
         ENDIF
         A(K)=X(I )
         B(K)=Y(I )
  140 CONTINUE

  145 WRITE(IO)IST1,ISY,A,B
      GOTO 9999
  150 WRITE(6,210)IST,FORM
      GOTO 9999
  160 WRITE(6,220)I
      GOTO 9999
C
C     FORMATANGABEN
C
  210 FORMAT(' SKURF: FEHLER IN DER PARAMETERLISTE',I11,E16.7)
  220 FORMAT(' SKURF: DOPPELTER PUNKT, INDEX',I10)
  230 FORMAT(8A4)
 9999 END
