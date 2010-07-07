C

      SUBROUTINE FL2O(A00,A1,A2,A3,A4,A5,A6,A7,A8,A9,INDE,X0,Y0,Z0,
     >                CX,CY,CZ,R,B0,B1,B2,B3,C0,C1,C2,C3,EPSIN,NMACH)
C
C***********************************************************************
C                                                  11. Juli    1988
C
C     Das Programm transformiert die allgemeine Gleichung
C     einer algebraischen Flaeche so, dass erkennbar ist,
C     welchen Koerper die Gleichung darstellt.
C
C     EINGABEPARAMETER :
C     ==================
C
C     A00, ... ,A9  -   KOEFFIZIENTEN DER ALGEBRAISCHEN FLAECHE
C                       A00 + A1*X + A2*Y + A3*Z + A4*X**2 + A5*Y**2
C                       + A6*Z**2 + A7*X*Y + A8*X*Z + A9*Y*Z = 0
C
C     AUSGABEPARAMETER :
C     ==================
C
C     INDE         -    INDEX DER ANGIBT, WELCHER KOERPER DURCH
C                       DIE ALGEBRAISCHE GLEICHUNG DARGESTELLT WIRD
C
C      INDE
C     -----------------------------------------
C        0
C        1
C        2
C        3
C        4  zylinder
C        5
C        6
C        7
C        8  kegel
C        9
C       10
C       11
C       12
C       13  ellipsoid
C       14  1 punkt
C
C
C     Hat INDE einen der folgenden Werte, so werden weitere
C     Parameter zurueckgegeben:
C
C     INDE = 1 :   DIE EBENENGLEICHUNG LAUTET:
C                  B0 + B1*X + B2*Y + B3*Z = 0
C
C     INDE = 2 :   DIE 1. EBENENGLEICHUNG LAUTET:
C                  B0 + B1*X + B2*Y + B3*Z = 0
C
C                  DIE 2. EBENENGLEICHUNG LAUTET:
C                  C0 + C1*X + C2*Y + C3*Z = 0
C
C     INDE = 3 :   DIE 1. EBENENGLEICHUNG LAUTET:
C                  B0 + B1*X + B2*Y + B3*Z = 0
C
C                  DIE 2. EBENENGLEICHUNG LAUTET:
C                  C0 + C1*X + C2*Y + C3*Z = 0
C
C     INDE =  4, 5, 6 ODER 7 :
C                  PARAMETERFORM DER ZYLINDERACHSE:
C                  ( X )     ( X0 )          ( CX )
C                  ( Y )  =  ( Y0 )  + MUE * ( CY )
C                  ( Z )     ( Z0 )          ( CZ )
C
C                                 T
C                  ( CX, CY, CZ )    IST DER EINHEITSRICHTUNGSVEKTOR
C                                    DER ZYLINDERACHSE
C
C                                 T
C                  ( X0, Y0, Z0 )    IST DER PUNKT DER ZYLINDERACHSE,
C                                    DER VOM URSPRUNG DEN KUERZESTEN
C                                    ABSTAND BESITZT.
C
C     INDE = 4 :   ZYLINDERRADIUS : R
C                                 T
C     INDE = 8 :   ( X0, Y0, Z0 )    IST DIE SPITZE DES KEGELS
C
C                  ( CX, CY, CZ )    IST DIE RICHTUNG DES KEGELS
C                  BEI EINEM KREISFOERMIGEN KEGEL IST DER
C                  OEFFNUNGSWINKEL: R
C     INDE = 13:   ( X0, Y0, Z0 )    IST DER MITTELPUNKT
C
C                  ( CX, CY, CZ )    SIND DIE HALBACHSEN
C                  BEI EINER KUGEL IST DER
C                  RADIUS: R
C***********************************************************************
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A00, A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >                      EPSIN
      REAL(DP), INTENT(OUT) :: X0, Y0, Z0, CX, CY, CZ, R,
     >                       B0, B1, B2, B3, C0, C1, C2, C3
      INTEGER, INTENT(IN) :: NMACH
      INTEGER, INTENT(OUT) :: INDE
      REAL(DP) :: A( 3,3 ), LAMBDA( 3 ), EV( 3,3 ), EPS, C, B(3),
     >          M(3),NUE,NORM,P(3),AMERK(3,3),RES(3),WRK(6)
     >          ,Q(5,5),DD(3),QEV(5,5)
      REAL(DP) :: AM4(4,4), LAMORI(3), A0
C
      INTEGER :: INULL, DIM, KBASIS, IW(3), I, J, ICOUNT
      LOGICAL :: NOPOS
C
C     DATA              EPS  / 5.D-10 /
      EPS=EPSIN
C
      INDE = 0
C
C     Die Koeffizienten der algebraischen Gleichung werden auf
C     die Matrix A uebertragen, so dass diese Gleichung uebergeht
C     in die Form :
C      T           T
C     X  A X + 2  B  X + C  = 0
C
C
C
      C      = A00
C
      A(1,1) = A4
      A(2,1) = A7/2.D0
      A(3,1) = A8/2.D0
      A(1,2) = A7/2.D0
      A(2,2) = A5
      A(3,2) = A9/2.D0
      A(1,3) = A8/2.D0
      A(2,3) = A9/2.D0
      A(3,3) = A6
C
      B(1) = A1/2.
      B(2) = A2/2.
      B(3) = A3/2.

      AM4(1:3,1:3) = A
      AM4(4,1:3) = B
      AM4(1:3,4) = B
      AM4(4,4) = C

      ICOUNT=1
  4   CONTINUE
C
C     PRUEFEN, OB ALLE ZEILEN DER MATRIX A GLEICH NULL SIND
C
      DIM = 0
C
      DO 2, I=1,3
         DO 2, J=1,3
            IF ( ABS( A(I,J) ) .GT. EPS ) DIM = 1
   2  CONTINUE
C
C
      IF ( DIM .EQ. 0 ) THEN
C
C     ALGEBRAISCHE FLAECHE STELLT EINE EBENE DAR
C     DIE KOEFFIZIENTEN DER EBENENGLEICHUNG WERDEN BESTIMMT
C
         B0 = 0.D0
         B1 = 0.D0
         B2 = 0.D0
         B3 = 0.D0
C
         IF ( ABS( A1 ) .GT. EPS ) THEN
            B0 = A00 / A1
            B1 = 1.D0
         ELSE IF ( ABS( A2 ) .GT. EPS ) THEN
            B0 = A00 / A2
            B2 = 1.D0
         ELSE
            B0 = A00 / A3
            B3 = 1.D0
         ENDIF
C
         INDE = 1
         GOTO 999

      ELSE
C
C        BERECHNUNG DER EIGENWERTE UND EIGENVEKTOREN VON A
C
         IF (NMACH.NE.1) CALL DEVCSF( 3,A,3,LAMBDA, EV, 3 )
         IF (NMACH.EQ.1) CALL EVCSF( 3,A,3,LAMBDA, EV, 3 )
C        WRITE (iunout,*) ' LAMBDA ',(LAMBDA(I),I=1,3)
C        WRITE (iunout,*) ' EV '
C        WRITE (iunout,*) ((EV(I,J),J=1,3),I=1,3)
         DO 4814 I=1,3
4814     IF (ABS(LAMBDA(I)).LT.EPS) LAMBDA(I)=0.
         LAMORI = LAMBDA
C
C        NORMIEREN DER EIGENVEKTOREN
C
         DO 7, I = 1,3
            NORM = SQRT( EV(1,I)**2 + EV(2,I)**2 + EV(3,I)**2 )
            DO 11, J = 1,3
               IF (NORM.GT.EPS)     EV(J,I) = EV(J,I)/NORM
  11        CONTINUE
  7      CONTINUE
C
      ENDIF
C
C     SORTIEREN DER EIGENWERTE UND DER DAZUGEHOERIGEN EIGENVEKTOREN,
C     SO DASS ZUERST DIE POSITIVEN, DANN DIE NEGATIVEN UND ZUM SCHLUSS
C     DIE EIGENWERTE, DIE NULL SIND, STEHEN.
C     WEITERHIN GILT : LAMBDA( 1 ) <= LAMBDA( 2 ) FALLS LAMBDA( 1),
C     LAMBDA( 2 ) > 0
C
C
      CALL SORT( EV, LAMBDA, INULL, NOPOS, EPS )
C     NOPOS = TRUE ==> es gibt keine positiven Eigenwerte, daher wird
C                      die gesamte Gleichung negiert
C
      IF (NOPOS.AND.ICOUNT.EQ.1) THEN
         DO 8, I = 1,3
            DO 9, J = 1,3
                  A(I,J) = -1. * A(I,J)
  9         CONTINUE
            B(I) = - B(I)
  8      CONTINUE
         C = -1. * C
         ICOUNT=2
         GOTO 4
      ELSEIF (NOPOS.AND.ICOUNT.EQ.2) THEN
C
C  FEHLER IN EVCSF ODER SORT
C
         GOTO 999
      ENDIF
C
C     Es sind INULL Eigenwerte Null
C
      DO 6, I = 1,3
         B(I) = -B(I)
 6    CONTINUE
C
C     Das Gleichungssystem Am=b ist loesbar, wenn Rang(A) = Rang(AB)
C     Das Gleichungssystem wird wird mit dem QR - Algorithmus geloest.
C     Sind die Residuen (b-A*m=RES) alle gleich Null, so gibt es
C     mindestens eine Loesung, sonst nicht.
C     IMSL - Verfahren:  DLSBRR
c
      IF (NMACH.NE.1) CALL DLSBRR(3,3,A,3,B,EPS,M,RES,KBASIS)
      IF (NMACH.EQ.1) CALL LSBRR(3,3,A,3,B,EPS,M,RES,KBASIS)
C     WRITE (iunout,*) ' IMSL  KBASIS = ',KBASIS
C     WRITE (iunout,*) ' M ',M
C
      DO 17, I = 1,3
         B(I) = -B(I)
17    CONTINUE
C
      DIM = 0
      DO 18, I = 1,3
         IF (ABS(RES(I)).GT.EPS) DIM = DIM + 1
 18   CONTINUE
      IF (DIM.EQ.0) THEN
C        Gleichungssystem ist loesbar!
C        d.h. singulaeres Gebilde ist nicht leer.
C        eine Loesung : M(I)

C                T
C        NUE := B  M + C
C
         NUE = C
         DO 20, I = 1,3
            NUE = NUE + B(I) * M(I)
   20    CONTINUE
C
C        Berechne normierten Eigenvektor EV1 = EV2 X EV3 ZU LAMBDA(1)
C
C
         EV(1,1) = EV(2,2)*EV(3,3) - EV(3,2)*EV(2,3)
         EV(2,1) = EV(3,2)*EV(1,3) - EV(3,3)*EV(1,2)
         EV(3,1) = EV(1,2)*EV(2,3) - EV(2,2)*EV(1,3)
C
         NORM = 0.D0
         DO 30, I = 1,3
            NORM = EV(I,1)* EV(I,1) + NORM
 30      CONTINUE
         NORM = SQRT(NORM)
C
         DO 40, I = 1,3
            EV(I,1) = EV(I,1) /NORM
 40      CONTINUE
C
C
C        Ermittlung, welchen Typs die algebraische Gleichung ist
C        Lambda(1) ist positiv
C
C        INDE = 0 ==> KEINE REELE FLAECHE
C
         IF (LAMBDA(2).GT. EPS) THEN
            IF (LAMBDA(3).GT. EPS) THEN
               IF (NUE.LT. 0.D0) THEN
C                 ELLIPSOID
                  INDE = 13
                  CALL KUGEL(AM4,LAMORI,X0,Y0,Z0,CX,CY,CZ,R,EPS)
               ELSEIF (ABS(NUE).LT.EPS) THEN
C                 1 PUNKT
                  INDE = 14
               ENDIF
            ELSEIF (LAMBDA(3).LT.-EPS) THEN
               IF (NUE.GT.EPS) THEN
C                 ZWEISCHALIGES HYPERBOLOID
                  INDE = 10
               ELSEIF (NUE.LT.-EPS) THEN
C                 EINSCHALIGES HYPERBOLOID
                  INDE = 9
               ELSE
C                 KEGEL
                  INDE = 8
                  CALL KEGEL(EV,LAMBDA,0._DP,M,X0,Y0,Z0,CX,CY,CZ,R,EPS)
               ENDIF
            ELSE
               IF (NUE.LT.-EPS) THEN
C                 ELLIPTISCHER ZYLINDER
                  INDE = 5
                  CALL ELLZYL(EV,LAMBDA,NUE,M,X0,Y0,Z0,CX,CY,CZ,R,INDE,
     .                        EPS)
               ELSEIF (ABS(NUE).LT.EPS) THEN
C                 1 GERADE
C                 INDE = ?????????????????????????????????
C                 X3-ACHSE
               ENDIF
            ENDIF
         ELSEIF (LAMBDA(2).LT.-EPS) THEN
            IF ((ABS(LAMBDA(3)).LT.EPS).AND.(ABS(NUE).LT.EPS)) THEN
C              2 SICH SCHNEIDENDE EBENEN
C              SCHNITTGERADE = X3-ACHSE
               INDE = 3
               CALL SCHEBE(EV,LAMBDA,0._DP,M,B0,B1,B2,B3,C0,C1,C2,C3)
            ELSEIF (ABS(LAMBDA(3)).LT.EPS) THEN
C              HYPERBOLISCHER ZYLINDER
               INDE = 6
            ENDIF
         ELSE
*           LAMBDA(2) = 0.D0
            IF (ABS(LAMBDA(3)).LT.EPS) THEN
               IF (NUE.LT.0.D0) THEN
C                 ZWEI PARALLELE EBENEN
C                 PARALLEL ZUR X2 - X3- EBENE
                  INDE = 2
                  CALL PAREBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3,
     .                        EPS)
               ELSEIF (ABS(NUE).LT.EPS) THEN
C                 DOPPELEBENE
C                 X2 - X3 -EBENE
                  INDE = 1
                  CALL DOPEBE(EV,LAMBDA,0._DP,M,B0,B1,B2,B3,EPS)
               ENDIF
            ENDIF
         ENDIF
C
C
C        Normalform :
C      lambda(1)*x1 **2 + lambda(2)*x2 **2 + lambda(3)*x3**2 + nue = 0.
C
C
C        Transformationsformel:
C        X = (ev1,ev2,ev3) * x + m
C
C        DIE ERSTEN 14 FAELLE SIND DAMIT ERSCHLAGEN
C        jetzt die Faelle mit leerem singulaerem Gebilde:
C
      ELSEIF ((ABS(LAMBDA(2)).LT.EPS).AND.(ABS(LAMBDA(3)).LT.EPS)) THEN
C        1 Fall
C        Typ: Parabolischer Zylinder
C
         INDE = 7
C        Stehen EV aufeinander senkrecht???????????????????????????
C                          T
C        Berechne P3 := EV3  * B (<>0)
C        aendere evt. Vorzeichen von ev3 so, dass P3 < 0
C
         P(3) = 0.
         DO 80, I = 1,3
            P(3) = P(3) + EV(I,3) * B(I)
  80     CONTINUE
         IF (P(3).GT.0.D0) THEN
            P(3) = -P(3)
            DO 90, I = 1,3
               EV(I,3) = -1.D0 * EV(I,3)
  90        CONTINUE
         ENDIF
C
C        Normalform:
C        Lambda(1) * X1 ** 2 + 2 * P(3) * X3 = 0.
C
C        Berechne normierten Eigenvektor EV1 = EV2 X EV3 zu LAMBDA(1)
C
         EV(1,1) = EV(2,2)*EV(3,3) - EV(3,2)*EV(2,3)
         EV(2,1) = EV(3,2)*EV(1,3) - EV(3,3)*EV(1,2)
         EV(3,1) = EV(1,2)*EV(2,3) - EV(2,2)*EV(1,3)
C
         NORM = 0.D0
         DO 95, I = 1,3
            NORM = EV(I,1)* EV(I,1) + NORM
 95      CONTINUE
         NORM = SQRT(NORM)
C
         DO 100, I = 1,3
            EV(I,1) = EV(I,1) /NORM
 100     CONTINUE
C
C        Scheitelpunkt:
C        M = (-p1/lambda(1)) * ev1 + 1./ (2.(p3) *
C            (p1**2 /lambda(1) - c) * ev3
C                     T
C        mit p1 := ev1  * b
C
         DO 110,I = 1,3
            M(I) = (-P(1) / LAMBDA(1)) * EV(I,1)
     F             + 1. / (2. * P(3)) * (P(1) ** 2 /LAMBDA(1) -C)
     F             * EV(I,3)
  110    CONTINUE
C
C        Transformationsformel:
C        X = (EV1,EV2,EV3) * X + M
C
      ELSE
C
         IF (LAMBDA(2).LT.-EPS) THEN
C           Typ: Hyperbolischer Paraboloid
            INDE = 12
         ELSE
C           Typ: Elliptischer Paraboloid
            INDE = 11
         ENDIF
C
C        Berechne P3 := EV3  * B (<>0)
C        aendere evt. Vorzeichen von ev3 so, dass P3 < 0
C
         P(3) = 0.
         DO 120, I = 1,3
            P(3) = P(3) + EV(I,3) * B(I)
  120    CONTINUE
         IF (P(3).GT.0.D0) THEN
            P(3) = -P(3)
            DO 130, I = 1,3
               EV(I,3) = -1.D0 * EV(I,3)
  130       CONTINUE
         ENDIF
C
C        Normalform:
C        lambda(1) * x1 ** 2 + lambda(2) * x2 ** 2 + 2.*p3*x3 = 0.
C
C        Berechne normierten Eigenvektor EV1 = EV2 X EV3 zu LAMBDA(1)
C
         EV(1,1) = EV(2,2)*EV(3,3) - EV(3,2)*EV(2,3)
         EV(2,1) = EV(3,2)*EV(1,3) - EV(3,3)*EV(1,2)
         EV(3,1) = EV(1,2)*EV(2,3) - EV(2,2)*EV(1,3)
C
         NORM = 0.D0
         DO 140, I = 1,3
            NORM = EV(I,1)* EV(I,1) + NORM
 140     CONTINUE
         NORM = SQRT(NORM)
C
         DO 150, I = 1,3
            EV(I,1) = EV(I,1) /NORM
 150     CONTINUE
C
C
C     Scheitelpunkt:
C     m : = (-p1/lambda(1)) * ev1 + (-p2/lambda(2)) * ev2
C         + 1./2./p3 * (p1**2 / lambda(1) + p2**2/lambda(2)-c)*ev3
C                  t                     t
C     mit p1 := ev1  * b    und p2 := ev2  *b
C

         DO 160, I = 1,2
            P(I) = 0.D0
            DO 170, J = 1,3
               P(I) = P(I) + EV(J,I) * B(J)
  170       CONTINUE
  160    CONTINUE
C
         DO 180, I = 1,3
            M(I) = (-P(1)/LAMBDA(1)) * EV(I,1) +
     F             (-P(2)/LAMBDA(2)) * EV(I,2) +
     F             1.D0/ (2. * P(3)) * ( P(1) ** 2 / LAMBDA(1)
     F                  + P(2)**2/LAMBDA(2) - C) * EV(I,3)
 180     CONTINUE
C
C        Transformationsformel:
C        X = (ev1,ev2,ev3) * x + m
      ENDIF
 999  CONTINUE
C
      END
