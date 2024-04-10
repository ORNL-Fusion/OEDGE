C@PROCESS OPT(3) IL(DIM) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C----------------------------------------------------------
C UPDATE 27.9.90 M. BUSCH    Daten werden formatfrei geschrieben
C                            frueher (20A4) d.h. CRAY unvertraeglich
C UPDATE 16.12.90 G.GROTEN   Abfrage auf Ungueltigkeit weg
C                            Modernisierung des FORTRAN
C UPDATE 25. 9.91 G.GROTEN
C     AUTOR: GERD GROTEN
C----------------------------------------------------------
      SUBROUTINE PSKURF(X,Y,IST,FORM,Z)
      REAL   X(IST),Y(IST),Z(7,IST),A(100),B(100)
      COMMON /GRCIL/ CILBER(8)
CDEC$ PSECT /GRCIL/ NOSHR
      SAVE /GRCIL/
      LOGICAL   BX,BY
      INTEGER IO
      REAL DIFREN
      PARAMETER( DIFREN = 3., IO=91)

      ISK=0
      GOTO 77
C----------------------------------------------------------
      ENTRY PSKURL(X,Y,IST,FORM,Z,KSK)
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
      IF( JSY.GT. 118 .AND. JSY.LT. 200    .OR.
     *    JSY.GT. 318 .AND. JSY.LT. 400    .OR.
     *    JSY.GT. 518 .AND. JSY.LT. 600    .OR.
     *    JSY.GT. 718 .AND. JSY.LT. 800    .OR.
     *    JSY.GT. 918 .AND. JSY.LT.1000    .OR.
     *    JSY.GT.1118 .AND. JSY.LT.1200    .OR.
     *    JSY.GT.1318                      .OR.
     *    IST.LT. 2                        .OR.
     *    FO .LT. 0.01                       ) GOTO 100
C
C     BESTIMMUNG DER MAXIMA UND MINIMA
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

      DO 11 I=1,IST
         IF (X(I).LT.DX1) DX1=X(I)
         IF (X(I).GT.DX2) DX2=X(I)
         IF (Y(I).LT.DY1) DY1=Y(I)
         IF (Y(I).GT.DY2) DY2=Y(I)
   11 CONTINUE
C
C     BESTIMMUNG DER GEOMETRIE-FAKTOREN
C
      IF (DX2.EQ.DX1.AND.DY2.EQ.DY1) GO TO 100
      CILBER(4)=DX1
      CILBER(5)=DX2
      CILBER(6)=DY1
      CILBER(7)=DY2
      F1=0.
      IF(BX) THEN
         DX1=ALOG10(DX1)
         DX2=ALOG10(DX2)
      ENDIF
      IF (DX2.NE.DX1) F1=200./(DX2-DX1)
      F2=0.
      IF(BY) THEN
         DY1=ALOG10(DY1)
         DY2=ALOG10(DY2)
      ENDIF
      IF (DY2.NE.DY1) F2=200./(DY2-DY1)*FO
      N=IST-1
      NM1=N-1
      NM2=N-2
C
C     BRECHNUNG DER SPLINES
C
      Z(1,1)=0.
      Z(2,1)=X(1)
      IF (BX) Z(2,1)=ALOG10(X(1))
      Z(3,1)=Y(1)
      IF (BY) Z(3,1)=ALOG10(Y(1))
      DO 22 I=1,N
         IP1=I+1
         Z(2,IP1)=X(IP1)
         IF (BX) Z(2,IP1)=ALOG10(X(IP1))
         Z(3,IP1)=Y(IP1)
         IF (BY) Z(3,IP1)=ALOG10(Y(IP1))
         U=SQRT( ((Z(2,IP1)-Z(2,I))*F1)**2 + ((Z(3,IP1)-Z(3,I))*F2)**2 )
         IF(U.EQ.0)GOTO 200
         Z(1,IP1)=Z(1,I)+U
   22 CONTINUE
C
C     BERECHNUNG DER MOMENTE
C
      KK=0
      KKK=1

    1 KK=KK+1
      KKK=KKK+1
      IF (N.EQ.1) THEN
         Z(KKK,1)=0.
         Z(KKK,2)=0.
         GO TO 4
      ENDIF
      H1=Z(1,2)-Z(1,1)
      DX2=H1
      IF (KK.EQ.1) Q1=(Z(2,2)-Z(2,1))/DX2
      IF (KK.EQ.2) Q1=(Z(3,2)-Z(3,1))/DX2
      DY2=Q1
      DO 33 I=1,NM1
         DX1=DX2
         DY1=DY2
         DX2=Z(1,I+2)-Z(1,I+1)
         IF (KK.EQ.1) DY2=(Z(2,I+2)-Z(2,I+1))/DX2
         IF (KK.EQ.2) DY2=(Z(3,I+2)-Z(3,I+1))/DX2
         DDX=DX1+DX2
         Z(5,I)=DX2/DDX
         Z(4,I)=1.-Z(5,I)
         Z(6,I)=6.*(DY2-DY1)/DDX
   33 CONTINUE
      IF (N.EQ.2) THEN
         Z(KKK,IST)=-Z(6,1)
         Z(KKK,N)=Z(6,1)
         Z(KKK,NM1)=-Z(6,1)
         GO TO 4
      ENDIF
      Z(7,N)=0.
      Z(7,NM1)=0.
      DDX=DX2+H1
      U=H1/DDX
      Z(6,N)=6.*(Q1-DY2)/DDX
      Z(4,N)=1.-U
      Z(7,1)=Z(4,1)
      Z(4,1)=2.
      DO 44 I=2,NM2
         N1=I-1
         FAC=Z(4,I)/Z(4,N1)
         Z(7,I)=-FAC*Z(7,N1)
         Z(6,I)=Z(6,I)-FAC*Z(6,N1)
         Z(4,I)=2.-FAC*Z(5,N1)
         FAC=U/Z(4,N1)
         U=-FAC*Z(5,N1)
         Z(7,N)=Z(7,N)-FAC*Z(7,N1)
         Z(6,N)=Z(6,N)-FAC*Z(6,N1)
   44 CONTINUE
      FAC=Z(4,NM1)/Z(4,NM2)
      Z(4,NM1)=2.-FAC*Z(5,NM2)
      Z(5,NM1)=Z(5,NM1)-FAC*Z(7,NM2)
      Z(6,NM1)=Z(6,NM1)-FAC*Z(6,NM2)
      FAC=U/Z(4,NM2)
      U=Z(4,N)-FAC*Z(5,NM2)
      Z(7,N)=2.+Z(7,N)-FAC*Z(7,NM2)
      Z(6,N)=Z(6,N)-FAC*Z(6,NM2)
      FAC=U/Z(4,NM1)
      Z(4,N)=Z(7,N)-FAC*Z(5,NM1)
      Z(6,N)=Z(6,N)-FAC*Z(6,NM1)
      FAC=Z(6,N)/Z(4,N)
      Z(KKK,IST)=FAC
      DO 55 I=1,NM1
         II=NM1-I+1
         Z(KKK,II+1)=(Z(6,II)-FAC*Z(7,II)-Z(5,II)*Z(KKK,II+2))/Z(4,II)
   55 CONTINUE
      Z(KKK,1)=Z(KKK,IST)
c    4 GO TO (1,5),KK
 4    if (kk == 1 ) then
         goto 1
      elseif (kk == 2) then 
         goto 5
      endif
C
C     SPLINE-INTERPOLATION UND SCHREIBEN AUF SYSUT1
C
    5 IST1=Z(1,IST)
      IF (Z(1,IST)-IST1.GT.0.) IST1=IST1+1
      IST1=(IST1+(IST+1)*DIFREN)/DIFREN
      I=1
      K=0
      L=0
      V1=0.
      U=0.
      DX2=Z(2,1)/6.
      DY2=Z(3,1)/6.
      DO 45 J=1,IST
         Z(4,J)=X(J)
         IF (BX) Z(4,J)=ALOG10(X(J))
         Z(5,J)=Y(J)
         IF (BY) Z(5,J)=ALOG10(Y(J))
   45 CONTINUE

C---- Loop

    6 IF (U.LT.Z(1,I)) GO TO 7
      IF (U.GT.Z(1,IST)) GO TO 9
      IF (U.EQ.Z(1,IST)) GO TO 7
      I=I+1
      DX1=DX2
      DY1=DY2
      DX2=Z(2,I)/6.
      DY2=Z(3,I)/6.
      V0=V1
      V1=Z(1,I)
      U=V0
      H=V1-V0
      H2=H*H
      DX3=(Z(4,I-1)-DX1*H2)/H
      DY3=(Z(5,I-1)-DY1*H2)/H
      DX4=(Z(4,I)-DX2*H2)/H
      DY4=(Z(5,I)-DY2*H2)/H
      DX1=DX1/H
      DY1=DY1/H
      D1=DX2/H
      D2=DY2/H

    7 K=K+1
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
      CILBER(4)=AMIN1(CILBER(4),A(K))
      CILBER(5)=AMAX1(CILBER(5),A(K))
      CILBER(6)=AMIN1(CILBER(6),B(K))
      CILBER(7)=AMAX1(CILBER(7),B(K))
      GO TO 6

    9 L=L+1
      DO 66 II=L,IST1
         K=K+1
         IF (K.GT.100) THEN
            WRITE(IO)IST1,ISY,A,B
            K=1
         ENDIF
         A(K)=X(I )
         B(K)=Y(I )
   66 CONTINUE
      WRITE(IO)IST1,ISY,A,B
      GOTO 9999
C
C     Fehlermeldungen
C
  100 WRITE(6,20)IST,FORM
   20 FORMAT(' PSKURF:  Fehler in der Parameterliste:',I11,E16.7)
      GOTO 9999
  200 WRITE(6,30)I
   30 FORMAT(' PSKURF:  Doppelter Punkt, Index',I8)

 9999 END
