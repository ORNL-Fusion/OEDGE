 

      SUBROUTINE COFACT (A4,A3,I,J)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A4(4,4)
      REAL(DP), INTENT(OUT) :: A3(3,3)
      INTEGER, INTENT(IN) :: I, J
      INTEGER :: II, IS, JJ, JS

      IS = 0
      DO II=1,4
        IF (II /= I) IS = IS + 1
        JS = 0
        DO JJ=1,4
          IF (JJ /= J) JS = JS + 1
          IF ((II /= I) .AND. (JJ /= J)) A3(IS,JS) = A4(II,JJ)
        END DO
      END DO

      RETURN
      END


      SUBROUTINE DEVCSF(N,A,LDA,EVAL,EVEC,LDEVEC)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, LDA, LDEVEC
      REAL(DP), INTENT(IN) :: A(LDA,N)
      REAL(DP), INTENT(OUT) :: EVAL(N), EVEC(LDEVEC,N)
      INTEGER :: IERR, J, I, IAA
      REAL(DP), ALLOCATABLE :: EVEC1(:,:),EVAL1(:),
     .                         FV1(:),FV2(:)

      entry evcsf (n,a,lda,eval,evec,ldevec)

      ALLOCATE (EVEC1(LDA,LDA))
      ALLOCATE (EVAL1(LDA))
      ALLOCATE (FV1(LDA))
      ALLOCATE (FV2(LDA))

      CALL RS(LDA,N,A,EVAL1,1,EVEC1,FV1,FV2,IERR)
      DO 10,I=1,N
         EVAL(I) = EVAL1(N-I+1)
         EVEC(1:LDEVEC,I) = EVEC1(1:LDEVEC,N-I+1)
10    CONTINUE

      DEALLOCATE (EVEC1)
      DEALLOCATE (EVAL1)
      DEALLOCATE (FV1)
      DEALLOCATE (FV2)
      END





      SUBROUTINE DLSBRR(NRA,NCA,A,LDA,B,TOL,X,RES,KBASIS)
      USE PRECISION
      IMPLICIT NONE
      INTEGER NRA,NCA,LDA,KBASIS,I,J,IAA
      REAL(DP) A(LDA,NCA),B(NRA),X(NCA),RES(NRA),TOL
      REAL(DP), ALLOCATABLE :: Q(:,:),R(:)
      INTEGER, ALLOCATABLE :: S(:)

      ALLOCATE (Q(NRA+2,NCA+2))
      ALLOCATE (R(NRA))
      ALLOCATE (S(NRA))

      DO 10,I=1,NRA
         DO 20,J=1,NCA
            Q(I,J) = A(I,J)
20       CONTINUE
         Q(I,NCA+1) = B(I)
10    CONTINUE
      CALL MA20A(Q,RES,X,R,S,NRA+2,NRA,NCA,TOL)
      KBASIS = Q(NRA+1,NCA+2)

      DEALLOCATE (Q)
      DEALLOCATE (R)
      DEALLOCATE (S)

      END
C
C
C
      SUBROUTINE DOPEBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,EPS)
**********************************************************************
*                                                     1. JUNI 1988   *
*     Inde = 1 ===> Es liegt eine Doppelebenen zur u2,-u3 Ebene vor  *
*     Die Normalform lautet: lambda(1) * u1 ** 2       = 0           *
*     d.h.: nue=0.
*     Aus dieser Gleichung ergibt sich die Loesung:                  *
*     u1 = 0.                                                        *
*     Diese Gleichung wird ruecktransformiert und mit                *
*     Hilfe eines Punktes, der die Normalform erfuellt, kann die     *
*     Ebenengleichung     e1: bo + b1*x + b2 * y + b3 * z            *
*     aufgestellt werden.                                            *
**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: LAMBDA( 3 ), EV(3,3), EPS,
     >                       M(3), NUE
      REAL(DP), INTENT(OUT) :: B0, B1, B2, B3
      REAL(DP) :: NORM, P(3)
C
C     DATA              EPS  / 5.D-10 /
C
*     Loesung der Gleichung:
*     U(1) = 0.0
*     U(2) = 0.0
*     U(3) = 0.0
*
*     Aufstellen des Normalenvektors:
*
      B1 = M(1)
      B2 = M(2)
      B3 = M(3)
*
*     Normieren:
*
      NORM = B1**2 + B2**2 + B3**2
      NORM = SQRT(NORM)
*
      B1 = B1 / NORM
      B2 = B2 / NORM
      B3 = B3 / NORM
*
*     Ein Punkt der Ebene:   P = (0.0,1.,1.)
*     Transformation:
      P(1) = EV(1,2) + EV(1,3) + M(1)
      P(2) = EV(2,2) + EV(2,3) + M(2)
      P(3) = EV(3,2) + EV(3,3) + M(3)
*
      B0 = -(B1 * P(1) + B2 * P(2) + B3 * P(3))
*
      END

c
      function e1 (x)
      USE PRECISION
      implicit none
      real(dp), intent(in) :: x
      real(dp) :: e1, s13aaf
      integer :: ifail
      ifail = 0
      e1 = s13aaf(x, ifail)
      return
      end
C
C
C
      SUBROUTINE ELLZYL(EV,LAMBDA,NUE,M,X0,Y0,Z0,CX,CY,CZ,R,INDE,EPS)
**********************************************************************
*                                                     6. JUNI 1988   *
*     Inde = 5 ===> Es liegt ein elliptischer Zylinder vor.          *
*     Die Normalform lautet:                                         *
*     lambda(1) * u1**2 + lambda(2) * u2**2 + nue = 0.               *
*     Aus dieser Gleichung ergibt sich die Loesung:                  *
*     u1 = sqrt ( (-nue -lambda(2) *u2**2)/ lambda(1) ) Mit          *
*     Hilfe eines Punktes, der die Normalform erfuellt, kann die     *
*     Gleichung des Zylinders:                                       *
*     ( X )     ( X0 )          ( CX )                               *
*     ( Y )  =  ( Y0 )  + MUE * ( CY )                               *
*     ( Z )     ( Z0 )          ( CZ )                               *
*                                                                    *
*                   T                                                *
*     ( CX, CY, CZ )    IST DER EINHEITSRICHTUNGSVEKTOR              *
*                       DER ZYLINDERACHSE                            *
*                                                                    *
*                   T                                                *
*     ( X0, Y0, Z0 )    IST DER PUNKT DER ZYLINDERACHSE,             *
*                       DER VOM URSPRUNG DEN KUERZESTEN              *
*                       ABSTAND BESITZT.                             *
*                                                                    *
*     INDE = 4 :   ZYLINDERRADIUS : R                                *
*                                                                    *
*     aufgestellt werden.                                            *
**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: X0, Y0, Z0, CX, CY, CZ, R
      REAL(DP), INTENT(IN) :: LAMBDA(3), EV(3,3), EPS,
     >                      M(3), NUE 
      REAL(DP) :: NORM, P(3)
C
      INTEGER, INTENT(OUT) :: INDE
C
C     DATA              EPS  / 5.D-10 /
*     Zentrum:
*
      CX  = EV(1,3)
      CY  = EV(2,3)
      CZ  = EV(3,3)
*
      X0 = M(1)
      Y0 = M(2)
      Z0 = M(3)
*
      IF (ABS(LAMBDA(2)-LAMBDA(1)).LT.EPS) THEN
         INDE = 4
*        Kreisfoermiger Zylinder
         R = SQRT(- NUE/LAMBDA(1))
      ENDIF
      END
c
c
      function erf (x)
c  s15aef is a nag routine: error function
c  erf(x)=2/sqrt(pi)*int^x_0 dt exp(-(t**2))
      USE PRECISION
      implicit none
      real(dp), intent(in) :: x
      real(dp) :: erf, s15aef
      integer :: ifail
      erf=s15aef(x,ifail)
      return
      end
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
C        WRITE (6,*) ' LAMBDA ',(LAMBDA(I),I=1,3)
C        WRITE (6,*) ' EV '
C        WRITE (6,*) ((EV(I,J),J=1,3),I=1,3)
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
C     WRITE (6,*) ' IMSL  KBASIS = ',KBASIS
C     WRITE (6,*) ' M ',M
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
C
*
      FUNCTION H1RN(DUMMY)
*
*#**********************************************************************
*# RANDOM NUMBER GENERATOR AS ADVOCATED BY F. JAMES FROM PROPOSAL OF   *
*# MARSAGLIA AND ZAMAN FSU-SCRI-87-50 AND MODIFIED BY F. JAMES 1988 TO *
*# PRODUCE VECTOR OF NUMBERS.                                          *
*# ENTRIES ARE:                                                        *
*#     FUNCTION    H1RN(DUMMY)     SINGLE RANDOM NUMBER                *
*#     SUBROUTINE  H1RNV(VEC,LEN)  VECTOR OF RANDOM NUMBERS            *
*#     SUBROUTINE  H1RNIN(IJ,KL)   INITIALISE WITH SEEDS               *
*#     SUBROUTINE  H1RNIV(VEC)     INITIALISE/RESTART WITH SEED ARRAY  *
*#     SUBROUTINE  H1RNSV(VEC)     SAVE SEED ARRAY VEC(100)            *
*#                                                                     *
*# NOTE: -H1RNIN OR H1RNIV MUST BE CALLED BEFORE GENERATING ANY        *
*#        RANDOM NUMBER(S).                                            *
*#       -H1RNSV SAVES SEED ARRAY INTO VEC(100) ONLY. THE USER HAS TO  *
*#        OUTPUT IT.                                                   *
*#                                                                     *
*# CHANGED BY: G. GRINDHAMMER AT: 90/03/14                             *
*# REASON :                                                            *
*#**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: DUMMY
      REAL(DP) :: H1RN
      INTEGER :: ISEED1, ISEED2
      CHARACTER(16) :: FLAG,CHECK
      REAL(DP) :: U, C, CD, CM
      INTEGER :: I, J
      COMMON /RASET1/ U(97),C,CD,CM,I,J
      COMMON /RASET2/ FLAG
*
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST) THEN
         IF (FLAG .NE. CHECK) THEN
cpb         WRITE(6,*) ' H1RN (RANMAR): INITIALIZED WITH DEFAULT SEED'
            ISEED1      = 12345
            ISEED2      = 98765
            CALL H1RNIN(ISEED1,ISEED2)
            FIRST = .FALSE.
         ELSE
            FIRST = .FALSE.
         ENDIF
      ENDIF
*
  100 CONTINUE
      H1RN = U(I)-U(J)
      IF(H1RN .LT. 0.) H1RN = H1RN + 1.
      U(I) = H1RN
      I = I - 1
      IF( I .EQ. 0) I=97
      J = J - 1
      IF( J .EQ. 0) J=97
      C = C - CD
      IF( C .LT. 0) C = C + CM
      H1RN = H1RN-C
      IF(H1RN .LE. 0.) H1RN = H1RN + 1.
      IF(H1RN .GE. 1.) GOTO 100
      RETURN
      END
*
*
      SUBROUTINE H1RNIN(IJ,KL)
*
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN OUT) :: IJ, KL
      INTEGER :: L, II, J, K, JJ, M, I
      REAL(DP) :: S, T
      CHARACTER*16    FLAG
      REAL(DP) :: U, C, CD, CM
      INTEGER :: IP, JP
      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
      COMMON /RASET2/ FLAG
*
      IJ = IABS(IJ)
      KL = IABS(KL)
      IJ = MOD(IJ,31329)
      KL = MOD(KL,30082)
      I  = MOD(IJ/177, 177) + 2
      J  = MOD(IJ, 177)     + 2
      K  = MOD(KL/169, 178) + 1
      L  = MOD(KL, 169)
*
      DO 300 II= 1, 97
         S= 0.
         T= 0.5
         DO 250 JJ= 1,24
            M = MOD(MOD(I*J,179)*K, 179)
            I = J
            J = K
            K = M
            L = MOD(53*L+1, 169)
            IF ( MOD(L*M,64) .GE. 32) S = S + T
  250       T = 0.5*T
  300 U(II) = S
      C  =   362436./16777216.
      CD =  7654321./16777216.
      CM = 16777213./16777216.
      IP = 97
      JP = 33
*
      FLAG = 'H1RN INITIALISED'
cpb   WRITE(6,610) IJ,KL,I,J,K,L
  610 FORMAT(' H1RNIN: H1RN (RANMAR) INITIALISED WITH 2 SEEDS: ',
     >       2I7,4I4)
*
      RETURN
      END
*
*
      SUBROUTINE H1RNIV(VEC)
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: VEC(100)
      INTEGER :: IC
      CHARACTER*16    FLAG
      REAL(DP) :: U, C, CD, CM
      INTEGER :: IP, JP
      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
      COMMON /RASET2/ FLAG
*
      DO 400 IC = 1, 97
  400 U(IC) = VEC(IC)
      C  = VEC(98)
      CD =  7654321./16777216.
      CM = 16777213./16777216.
      IP = NINT(VEC(99))
      JP = NINT(VEC(100))
*
      FLAG = 'H1RN INITIALISED'
      WRITE(6,*) ' H1RNIV: H1RN (RANMAR) INITIALISED/RESTARTED WITH',
     >           ' SEED ARRAY VEC(100)'
*
      RETURN
      END
*
*
      SUBROUTINE H1RNSV(VEC)
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: VEC(100)
      INTEGER :: IC
      REAL(DP) :: U, C, CD, CM
      INTEGER :: I, J
      COMMON /RASET1/ U(97),C,CD,CM,I,J
*
      DO 10 IC = 1, 97
   10 VEC(IC) = U(IC)
      VEC(98) = C
      VEC(99) = REAL(I)
      VEC(100)= REAL(J)
      RETURN
      END
*
      SUBROUTINE H1RNV(RVEC,LEN)
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: RVEC(1)
      INTEGER, INTENT(IN) :: LEN
      REAL(DP) :: UNI
      INTEGER :: IVEC
      CHARACTER(16) :: FLAG,CHECK
      REAL(DP) :: U, C, CD, CM
      INTEGER :: I, J
      COMMON /RASET1/ U(97),C,CD,CM,I,J
      COMMON /RASET2/ FLAG
*
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST) THEN
         IF (FLAG .NE. CHECK) THEN
            WRITE(6,*) ' H1RNV (RANMAR): CALL H1RNIN OR H1RNIV BEFORE',
     >                 ' CALLING H1RN.'
            STOP
         ELSE
            FIRST = .FALSE.
         ENDIF
      ENDIF
*
      DO 200 IVEC = 1,LEN
  190    CONTINUE
         UNI = U(I)-U(J)
         IF(UNI .LT. 0.) UNI = UNI + 1.
         U(I) = UNI
         I = I - 1
         IF( I .EQ. 0) I=97
         J = J - 1
         IF( J .EQ. 0) J=97
         C = C - CD
         IF( C .LT. 0) C = C + CM
         UNI = UNI-C
         IF(UNI .LE. 0.) UNI = UNI + 1.
         IF(UNI .GE. 1.) GOTO 190
         RVEC(IVEC) = UNI
  200 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE KEGEL(EV,LAMBDA,NUE,M,X0,Y0,Z0,CX,CY,CZ,R,EPS)
**********************************************************************
*                                                    11. JULI 1988   *
*     Inde = 8 ===> Es liegt ein elliptischer Kegel vor.             *
*     Die Normalform lautet:                                         *
*     lambda(1) * u1**2 + lambda(2) * u2**2 + lambda(3) * u3**2 =0.  *
*     lambda(1) und lambda(2) sind positiv, lambda(3) ist negativ.   *
*     d.h. nue=0.
*     Aus dieser Gleichung ergibt sich :                             *
*     Der transformierte Kegel hat sein Spitze im Koordinaten-       *
*     ursprung (dieser wird ruecktransformiert), der Kegel steht     *
*     senkrecht zur Z- Achse (d.h. der Punkt (0,0,1) erfuellt        *
*     die Ebenengleichung.)                                          *
*     Falls Lambda(1) und Lambda(2) gleich sind, handelt es sich um  *
*     einen kreisfoermigen Kegel, in diesem Fall kann der Winkel     *
*     bestimmt werden. (R)                                           *
*                   T                                                *
*     ( CX, CY, CZ )    IST DER EINHEITSRICHTUNGSVEKTOR              *
*                       DER ZYLINDERACHSE                            *
*                                                                    *
*                   T                                                *
*     ( X0, Y0, Z0 )    IST DIE SPITZE DES KEGELS                    *
*                                                                    *
*     R   :   WINKEL                                                 *
*                                                                    *
**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: X0, Y0, Z0, CX, CY, CZ, R
      REAL(DP), INTENT(IN) :: LAMBDA(3), EV(3,3), EPS,
     >                      M(3), NUE 
      REAL(DP) :: NORM, P(3)
C
C     DATA              EPS  / 5.D-10 /
*
      NORM = 0.0
      R = 0.0
*
      X0  = M(1)
      Y0  = M(2)
      Z0  = M(3)
*
      CX = EV(1,3) * (-LAMBDA(3))
      CY = EV(2,3) * (-LAMBDA(3))
      CZ = EV(3,3) * (-LAMBDA(3))
      NORM = SQRT (CX ** 2 + CY ** 2 + CZ ** 2)
      CX = CX / NORM
      CY = CY / NORM
      CZ = CZ / NORM
*
      IF (ABS(LAMBDA(2)-LAMBDA(1)).LT.EPS) THEN
*        Kreisfoermiger Kegel
         R = ATAN(sqrt(-LAMBDA(3))/sqrt(LAMBDA(1)))
      ENDIF
      END
C
C
C
      SUBROUTINE KUGEL(A,LAMBDA,X0,Y0,Z0,CX,CY,CZ,R,EPS)
**********************************************************************
*                                                    17. Maerz 1995  *
*     Inde =13 ===> Es liegt ein Ellipsoid vor.                      *
*     Die Normalform lautet:                                         *
*     lambda(1)*u1**2+lambda(2)*u2**2+lambda(3)*u3**2 +nue  =0.      *
*     lambda(1),lambda(2) lambda(3)sind positiv, nue ist negativ.    *
*     Aus dieser Gleichung ergibt sich :                             *
*     Der transformierte Kegel hat sein Spitze im Koordinaten-       *
*     ursprung (dieser wird ruecktransformiert), der Kegel steht     *
*     senkrecht zur Z- Achse (d.h. der Punkt (0,0,1) erfuellt        *
*     die Ebenengleichung.)                                          *
*     Falls alle Lambda  gleich sind, handelt es sich um             *
*     eine Kegel, in diesem Fall kann der RADIUS                     *
*     bestimmt werden. (R)                                           *
*                   T                                                *
*     ( CX, CY, CZ )    SIND DIE HALBACHSEN                          *
*                                                                    *
*                   T                                                *
*     ( X0, Y0, Z0 )    IST DER MITTELPUNKT (URSPRUNG)               *
*                                                                    *
*     R   :   WINKEL                                                 *
*                                                                    *
**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: X0, Y0, Z0, CX, CY, CZ, R
      REAL(DP), INTENT(IN) :: LAMBDA(3), A(4,4), EPS
      REAL(DP) :: A3(3,3), DD, AA, D1, D2, D3, D4, ADD, SARRUS
C
C     DATA              EPS  / 5.D-10 /
cdr  noch nicht fertig
!     write (6,*) 'exit in subr. kugel '
!     call exit
*
      A3 = A(1:3,1:3)
      DD = SARRUS(A3)
*
      CALL COFACT(A,A3,1,1)
      D1 = SARRUS(A3)
*
      CALL COFACT(A,A3,2,1)
      D2 = SARRUS(A3)
*
      CALL COFACT(A,A3,3,1)
      D3 = SARRUS(A3)
*
      CALL COFACT(A,A3,4,1)
      D4 = SARRUS(A3)
*
      AA = A(1,1)*D1 - A(2,1)*D2 + A(3,1)*D3 - A(4,1)*D4

      ADD = AA / (LAMBDA(1)*LAMBDA(2)*LAMBDA(3))

      CX = SQRT(-1.D0/LAMBDA(1)*ADD)
      CY = SQRT(-1.D0/LAMBDA(2)*ADD)
      CZ = SQRT(-1.D0/LAMBDA(3)*ADD)
      R = (CX + CY + CZ) / 3.D0

      A3(1:3,1) = A(1:3,4)
      A3(1:3,2:3) = A(1:3,2:3)
      X0 = - SARRUS(A3) / DD

      A3(1:3,1) = A(1:3,1)
      A3(1:3,2) = A(1:3,4)
      A3(1:3,3) = A(1:3,3)
      Y0 = - SARRUS(A3) / DD

      A3(1:3,1) = A(1:3,1)
      A3(1:3,2) = A(1:3,2)
      A3(1:3,3) = A(1:3,4)
      Z0 = - SARRUS(A3) / DD

      WRITE (6,*) ' IN KUGEL, ELLIPSOID '
      WRITE (6,*) ' X0,Y0,Z0 ',X0,Y0,Z0
      WRITE (6,*) ' CX,CY,CZ ',CX,CY,CZ
      WRITE (6,*) ' R ',R

      END


      SUBROUTINE LSBRR(NRA,NCA,A,LDA,B,TOL,X,RES,KBASIS)
      USE PRECISION
      IMPLICIT NONE
      INTEGER NRA,NCA,LDA,KBASIS,I,J,IAA
      REAL(DP) A(LDA,NCA),B(NRA),X(NCA),RES(NRA),TOL
      write (6,*) ' lsbrr is called '
      call exit
      end
C
C
C##       MA20A          28/06/72
C NAME MA20A(R)                  CHECK
      SUBROUTINE MA20A(Q,D,A,R,S,IQ,M,N,TOLER)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IQ, M, N
      REAL(DP), INTENT(IN) :: TOLER
      REAL(DP), INTENT(OUT) :: Q(IQ,*),A(*),D(*),R(*)
      INTEGER, INTENT(OUT) :: S(*) 
      REAL(DP) :: SUM, B, PIVOT, BIG
      REAL(DP) :: MIN,MAX
      INTEGER :: OUT, KL, KR, IN, L, K, KOUNT, M2, I, J, N1,
     .           N2, M1
      LOGICAL :: STAGE, TEST
C  ***BIG MUST BE SET EQUAL TO ANY VERY LARGE REAL CONSTANT.
C  ***ITS VALUE HERE IS APPROPRIATE FOR THE IBM 370.
      DATA BIG /1.D75/
C  ***INITIALIZATION
      M2=M+2
      N2=N+2
      M1=M+1
      N1=N+1
      DO 1 J=1,N
      Q(M2,J)=J
    1 A(J)=0.
      DO 3 I=1,M
      Q(I,N2)=N+I
      D(I)=0.
      IF(Q(I,N1).GE.0) GO TO 3
      DO 2 J=1,N2
    2 Q(I,J)=-Q(I,J)
    3 CONTINUE
C  ***COMPUTE MARGINAL COSTS
      DO 5 J=1,N1
      SUM=0.
      DO 4 I=1,M
    4 SUM=SUM+Q(I,J)
    5 Q(M1,J)=SUM
C  ***STAGE I
C  ***DETERMINE VECTOR TO ENTER THE BASIS
      STAGE=.TRUE.
      KOUNT=0
      KR=1
      KL=1
    6 MAX=-1.
      DO 7 J=KR,N
      IF(ABS(Q(M2,J)).GT.N) GO TO 7
      B=ABS(Q(M1,J))
      IF(B.LE.MAX) GO TO 7
      MAX=B
      IN=J
    7 CONTINUE
      IF(Q(M1,IN).GE.0) GO TO 9
      DO 8 I=1,M2
    8 Q(I,IN)=-Q(I,IN)
C  ***DETERMINE VECTOR TO LEAVE THE BASIS
    9 K=0
      DO 10 I=KL,M
      B=Q(I,IN)
      IF(B.LE.TOLER) GO TO 10
      K=K+1
      R(K)=Q(I,N1)/B
      S(K)=I
      TEST=.TRUE.
   10 CONTINUE
   11 IF(K.GT.0) GO TO 12
      TEST=.FALSE.
      GO TO 14
   12 MIN=BIG
      DO 13 I=1,K
      IF(R(I).GE.MIN) GO TO 13
      J=I
      MIN=R(I)
      OUT=S(I)
   13 CONTINUE
      R(J)=R(K)
      S(J)=S(K)
      K=K-1
C  ***CHECK FOR LINEAR DEPENDENCE IN STAGE I
   14 IF(TEST.OR..NOT.STAGE) GO TO 16
      DO 15 I=1,M2
      B=Q(I,KR)
      Q(I,KR)=Q(I,IN)
   15 Q(I,IN)=B
      KR=KR+1
      GO TO 25
   16 IF(TEST) GO TO 17
      Q(M2,N1)=2.
      GO TO 34
   17 PIVOT=Q(OUT,IN)
      IF(Q(M1,IN)-PIVOT-PIVOT.LE.TOLER) GO TO 19
      DO 18 J=KR,N1
      B=Q(OUT,J)
      Q(M1,J)=Q(M1,J)-B-B
   18 Q(OUT,J)=-B
      Q(OUT,N2)=-Q(OUT,N2)
      GO TO 11
C  ***PIVOT ON Q(OUT,IN)
   19 DO 20 J=KR,N1
      IF(J.EQ.IN) GO TO 20
      Q(OUT,J)=Q(OUT,J)/PIVOT
   20 CONTINUE
      DO 22 I=1,M1
      IF(I.EQ.OUT) GO TO 22
      B=Q(I,IN)
      DO 21 J=KR,N1
      IF(J.EQ.IN) GO TO 21
      Q(I,J)=Q(I,J)-B*Q(OUT,J)
   21 CONTINUE
   22 CONTINUE
      DO 23 I=1,M1
      IF(I.EQ.OUT) GO TO 23
      Q(I,IN)=-Q(I,IN)/PIVOT
   23 CONTINUE
      Q(OUT,IN)=1./PIVOT
      B=Q(OUT,N2)
      Q(OUT,N2)=Q(M2,IN)
      Q(M2,IN)=B
      KOUNT=KOUNT+1
      IF(.NOT.STAGE) GO TO 26
C  ***INTERCHANGE ROWS IN STAGE I
      KL=KL+1
      DO 24 J=KR,N2
      B=Q(OUT,J)
      Q(OUT,J)=Q(KOUNT,J)
   24 Q(KOUNT,J)=B
   25 IF(KOUNT+KR.NE.N1) GO TO 6
C  ***STAGE II
      STAGE=.FALSE.
C  ***DETERMINE VECTOR TO ENTER THE BASIS
   26 MAX=-BIG
      DO 28 J=KR,N
      B=Q(M1,J)
      IF(B.GE.0) GO TO 27
      IF(B.GT.-2.) GO TO 28
      B=-B-2.
   27 IF(B.LE.MAX) GO TO 28
      MAX=B
      IN=J
   28 CONTINUE
      IF(MAX.LE.TOLER) GO TO 30
      IF(Q(M1,IN).GT.0) GO TO 9
      DO 29 I=1,M2
   29 Q(I,IN)=-Q(I,IN)
      Q(M1,IN)=Q(M1,IN)-2.
      GO TO 9
C  ***PREPARE OUTPUT
   30 L=KL-1
      DO 32 I=1,L
      IF(Q(I,N1).GE.0) GO TO 32
      DO 31 J=KR,N2
   31 Q(I,J)=-Q(I,J)
   32 CONTINUE
      Q(M2,N1)=0.
      IF(KR.NE.1) GO TO 34
      DO 33 J=1,N
      B=ABS(Q(M1,J))
      IF(B.LE.TOLER.OR.2.-B.LE.TOLER) GO TO 34
   33 CONTINUE
      Q(M2,N1)=1.
   34 DO 37 I=1,M
      K=Q(I,N2)
      B=Q(I,N1)
      IF(K.GT.0) GO TO 35
      K=-K
      B=-B
   35 IF(I.GE.KL) GO TO 36
      A(K)=B
      GO TO 37
   36 K=K-N
      D(K)=B
   37 CONTINUE
      Q(M2,N2)=KOUNT
      Q(M1,N2)=N1-KR
      SUM=0.
      DO 38 I=KL,M
   38 SUM=SUM+Q(I,N1)
      Q(M1,N1)=SUM
C     WRITE (6,*) ' A ',(A(I),I=1,M)
      RETURN
      END
C
C
C
      SUBROUTINE PAREBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3,EPS)
**********************************************************************
*                                                     1. JUNI 1988   *
*     Inde = 2 ===> Es liegen 2 parallele Ebenen zur u2,-u3 Ebene vor*
*     Die Normalform lautet: lambda(1) * u1 ** 2 + nue = 0           *
*     Aus dieser Gleichung ergeben sich die Loesungen:               *
*     u1 = sqrt ( -nue/lambda(1))   und  u1 = -sqrt(-nue/lambda(1)   *
*     Diese beiden Gleichungen werden ruecktransformiert und mit     *
*     Hilfe eines Punktes, der die Normalform erfuellt koennen die   *
*     Ebenengleichungen   e1: bo + b1*x + b2 * y + b3 * z            *
*                         e2: co + b1*x + b2 * y + b3 * z            *
*     aufgestellt werden.                                            *
**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: B0, B1, B2, B3, C0, C1, C2, C3
      REAL(DP), INTENT(IN) :: LAMBDA( 3 ), EV(3,3), EPS,
     >                      M(3), NUE
      REAL(DP) :: NORM, P(3), U(3)

C
C     DATA              EPS  / 5.D-10 /
C
*     1. Loesung der Gleichung:
*
      U(1) = SQRT(-NUE/LAMBDA(1))
      U(2) = 0.0
      U(3) = 0.0
*
*     Aufstellen des Normalenvektors:'
*
      B1 = EV(1,1) * U(1)
      B2 = EV(2,1) * U(1)
      B3 = EV(3,1) * U(1)
*
*     Normieren:
*
      NORM = B1**2 + B2**2 + B3**2
      NORM = SQRT(NORM)
*
      B1 = B1 / NORM
      B2 = B2 / NORM
      B3 = B3 / NORM
*
*     Ein Punkt der Ebene:   P = (U(1),1.,1.)
*     Transformation:
      P(1) = EV(1,1) * U(1) + EV(1,2) + EV(1,3) + M(1)
      P(2) = EV(2,1) * U(1) + EV(2,2) + EV(2,3) + M(2)
      P(3) = EV(3,1) * U(1) + EV(3,2) + EV(3,3) + M(3)
*
      B0 = -(B1 * P(1) + B2 * P(2) + B3 * P(3))
*
*
*     2. Loesung der Gleichung:
*
      U(1) = -SQRT(-NUE/LAMBDA(1))
      U(2) = 0.0
      U(3) = 0.0
*
*     Aufstellen des Normalenvektors:'
*
      C1 = EV(1,1) * U(1)
      C2 = EV(2,1) * U(1)
      C3 = EV(3,1) * U(1)
*
*     Normieren:
*
      NORM = C1**2 + C2**2 + C3**2
      NORM = SQRT(NORM)
*
      C1 = C1 / NORM
      C2 = C2 / NORM
      C3 = C3 / NORM
*
*     Ein Punkt der Ebene:   P = (U(1),1.,1.)
*     Transformation:
      P(1) = EV(1,1) * U(1) + EV(1,2) + EV(1,3) + M(1)
      P(2) = EV(2,1) * U(1) + EV(2,2) + EV(2,3) + M(2)
      P(3) = EV(3,1) * U(1) + EV(3,2) + EV(3,3) + M(3)
*
      C0 = -(C1 * P(1) + C2 * P(2) + C3 * P(3))
*
      END

      SUBROUTINE QUADEQ (A,B,C,T1,T2,IER)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A, B, C
      REAL(DP), INTENT(OUT) :: T1, T2
      INTEGER, INTENT(OUT) :: IER
      REAL(DP) :: EPS10, BH, SQ, ROOT
      DATA EPS10 /1.E-10_DP/

      IER = 0
      IF ( ABS(A) < EPS10) THEN
        WRITE (6,*) ' A = 0 IN QUADEQ '
        IER = 1
      ELSE
        BH = B/(2.D0*A)
        ROOT = BH*BH - C/A
        IF ( ROOT < 0.D0 ) THEN
          WRITE (6,*) ' PROBLEM: SQRT(-X) IN QUADEQ '
          IER = 2
        ELSE
          SQ = SQRT(ROOT)
          T1 = -BH + SQ
          T2 = -BH - SQ
!         WRITE (6,*) ' T1,T2 ',T1,T2
        END IF
      END IF

      RETURN
      END
      function ranf ()
      USE PRECISION
      implicit none
      real(dp) :: ra, dummy, h1rn, ranf

      ra=h1rn(dummy)
      ranf=ra
      return
      end
c 
c
      function ranget(iseed)
      USE PRECISION
      implicit none
      integer, intent(inout) :: iseed
      real(dp) :: ran, ranf
      integer ranget, ifirst, large
      data ifirst/0/
      save large
      if (ifirst.eq.0) then
        large=HUGE(1)
        ifirst=1
      endif
      ran=ranf()
      iseed=ran*large
      ranget=iseed
      return
      end
c
c
      function ranset(iseed)
      USE PRECISION
      implicit none
      integer, intent(in) :: iseed
      integer :: iseed1, iseed2
      real(dp) :: ranset
      iseed1=iseed
      iseed2=2000000-iseed
CPB   iseed2=iseed+1
      call h1rnin(iseed1,iseed2)
      ranset=0
      return
      end


C   NAME                -   RS
C-----------------------------------------------------------------------
C
C   PURPOSE             -   CALCULATION OF EIGENVALUES AND OPTIONALLY
C                           EIGENVECTORS OF A REAL SYMMETRIC MATRIX
C
C   COMPILER            -   VS FORTRAN
C
C   PRECISION           -   IBM DOUBLE
C
C   USAGE               -   CALL RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C
C   ARGUMENTS    NM     -   MUST BE SET TO THE ROW DIMENSION OF THE
C                           TWO-DIMENSIONAL ARRAY PARAMETERS AS DECLARED
C                           IN THE CALLING PROGRAM DIMENSION STATEMENT.
C                                                                (INPUT)
C
C                N      -   ORDER OF THE MATRIX  A.              (INPUT)
C
C                A      -   MATRIX OF DIMENSION N BY N,
C                           CONTAINS THE REAL SYMMETRIC MATRIX   (INPUT)
C
C                W      -   VECTOR OF LENGTH N,
C                           CONTAINS THE EIGENVALUES IN ASCENDING ORDER
C                                                               (OUTPUT)
C
C                MATZ   -   INTEGER VARIABLE SET EQUAL TO ZERO IF
C                           ONLY EIGENVALUES ARE DESIRED,  OTHERWISE
C                           SET TO ANY NON-ZERO INTEGER FOR BOTH
C                           EIGENVALUES AND EIGENVECTORS.        (INPUT)
C
C                Z      -   MATRIX OF DIMENSION N BY N,
C                           CONTAINS THE EIGENVECTORS IF MATZ IS
C                           NOT ZERO.                           (OUTPUT)
C
C                FV1    -   VECTOR OF LENGTH N. WORK AREA.
C
C                FV2    -   VECTOR OF LENGTH N. WORK AREA.
C
C                IERR   -   INTEGER OUTPUT VARIABLE SET EQUAL TO
C                           AN ERROR COMPLETION CODE DESCRIBED IN
C                           SECTION 2.3 OF THE REFERENCE. THE
C                           NORMAL COMPLETION CODE IS ZERO.     (OUTPUT)
C
C   REQD. ROUTINES      -   TRED1, TQLRAT AND
C                           TRED2, TQL2 ARE PROVIDED.
C
C   REFERENCES          -   B.T. SMITH ET.AL.
C                           MATRIX EIGENSYSTEM ROUTINES - EISPACK GUIDE
C                           SECOND EDITION 1976
C                           LECTURE NOTES IN COMPUTER SCIENCES 6
C
C   REMARKS             -   THIS SUBROUTINE CALLS THE RECOMMENDED
C                           SEQUENCE OF SUBROUTINES FROM THE EIGENSYSTEM
C                           SUBROUTINE PACKAGE (EISPACK) TO FIND THE
C                           EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C                           OF A REAL SYMMETRIC MATRIX.
C
C                           THE CORRESPONDING CRAY-VERSION IS A MEMBER
C                           OF THE $SCILIB.
C
C   ALGORITHM           -   THE MATRIX A IS REDUCED TO A SYMMETRIC
C                           TRIDIAGONAL MATRIX BY ORTHOGONAL TRANSFOR-
C                           MATIONS. THE EIGENVALUE PROBLEM FOR THE RE-
C                           DUCED MATRIX IS SOLVED BY THE QL ALGORITHM.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ------------------------------------------------------------------
C
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
CCCCC
CCCCC     AENDERUNG:    REAL ---> REAL(DP)
CCCCC
      USE PRECISION
      IMPLICIT NONE
C
      INTEGER :: N,NM,IERR,MATZ
      REAL(DP) :: A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC MATRIX.
C
C     ON INPUT-
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT,
C
C        N  IS THE ORDER OF THE MATRIX  A,
C
C        A  CONTAINS THE REAL SYMMETRIC MATRIX,
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT-
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     MODIFIED FOR CFT BY MIKE ESS, CRI, JULY 1980
C
C     ------------------------------------------------------------------
C
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     ********** FIND EIGENVALUES ONLY **********
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
C     ********** LAST CARD OF RS **********
      END


      FUNCTION SARRUS (A)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A(3,3)
      REAL(DP) :: SARRUS, A1, A2, A3, B1, B2, B3

      A1 = A(1,1) * A(2,2) * A(3,3)
      A2 = A(1,2) * A(2,3) * A(3,1)
      A3 = A(1,3) * A(2,1) * A(3,2)
      B1 = A(3,1) * A(2,2) * A(1,3)
      B2 = A(3,2) * A(2,3) * A(1,1)
      B3 = A(3,3) * A(2,1) * A(1,2)
      SARRUS = A1 + A2 + A3 - B1 - B2 - B3
!     WRITE (6,*) ' A1,A2,A3 ',A1,A2,A3
!     WRITE (6,*) ' B1,B2,B3 ',B1,B2,B3
!     WRITE (6,*) ' DET ',A1 + A2 + A3 - B1 - B2 - B3

      RETURN
      END
C
C
C
      SUBROUTINE SCHEBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3)
**********************************************************************
*                                                     1. JUNI 1988   *
*     Inde = 3 ===> Es liegen 2 sich schneidende Ebenen vor.         *
*     Schnittgerade: u3 - Achse                                      *
*     Die Normalform lautet:                                         *
*      lambda(1) * u1 ** 2 + lambda(2) * u2**2  = 0                  *
*     d.h. nue=0.
*     Aus dieser Gleichung ergeben sich die Loesungen:               *
*     u1 = sqrt ( -lambda(2)/lambda(1)) * u2                         *
*     u1 = - sqrt ( -lambda(2)/lambda(1)) * u2                       *
*     Fuer jede dieser Gleichung werden Normalenvektor und ein Punkt *
*     der Ebene ruecktransformiert, so dass die Ebenengleichungen    *
*                         e1: bo + b1*x + b2 * y + b3 * z            *
*                         e2: co + b1*x + b2 * y + b3 * z            *
*     aufgestellt werden.                                            *
**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: LAMBDA( 3 ), EV(3,3),
     >                      M(3), NUE
      REAL(DP), INTENT(OUT) :: B0, B1, B2, B3, C0, C1, C2, C3
      REAL(DP) :: NORM, P(3), U(3)

*
*
*     1. Loesung der Gleichung:
*     1. Normalenvektor
      U(1) = 1.0
      U(2) = -SQRT(-LAMBDA(2)/LAMBDA(1))
      U(3) = 0.0
*
      NORM = U(1)**2 + U(2)**2
      NORM = SQRT(NORM)
*
      U(1) = U(1) / NORM
      U(2) = U(2) / NORM
*
*     Ruecktransformation des Normalenvektors:'
*
      B1 = EV(1,1) * U(1) + EV(1,2) * U(2)
      B2 = EV(2,1) * U(1) + EV(2,2) * U(2)
      B3 = EV(3,1) * U(1) + EV(2,3) * U(2)
*
*     Normieren:
*
      NORM = B1**2 + B2**2 + B3**2
      NORM = SQRT(NORM)
*
      B1 = B1 / NORM
      B2 = B2 / NORM
      B3 = B3 / NORM
*
*     Ein Punkt der Ebene:   P = (sqrt(-lambda(2)/lambda(1),1.,0.)
*     Ruecktransformation:
      P(1) = EV(1,1)*SQRT(-LAMBDA(2)/LAMBDA(1)) + EV(1,2)+M(1)
      P(2) = EV(2,1)*SQRT(-LAMBDA(2)/LAMBDA(1)) + EV(2,2)+M(2)
      P(3) = EV(3,1)*SQRT(-LAMBDA(2)/LAMBDA(1)) + EV(3,2)+M(3)
*
      B0 = -(B1 * P(1) + B2 * P(2) + B3 * P(3))
*
*
*     2. Loesung der Gleichung:
*     Normalenvektor:
      U(1) = 1.0
      U(2) = SQRT(-LAMBDA(2)/LAMBDA(1))
      U(3) = 0.0

      NORM = U(1)**2 + U(2)**2
      NORM = SQRT(NORM)
*
      U(1) = U(1) / NORM
      U(2) = U(2) / NORM
*
*     Ruecktransformation des Normalenvektors:'
*
      C1 = EV(1,1) * U(1) + EV(1,2) * U(2)
      C2 = EV(2,1) * U(1) + EV(2,2) * U(2)
      C3 = EV(3,1) * U(1) + EV(3,2) * U(2)
*
*     Normieren:
*
      NORM = C1**2 + C2**2 + C3**2
      NORM = SQRT(NORM)
*
      C1 = C1 / NORM
      C2 = C2 / NORM
      C3 = C3 / NORM
*
*     Ein Punkt der Ebene:   P = (-SQRT(-LAMBDA(2)/LAMBDA(1)),1.,0.)
*     Transformation:
      P(1) = -EV(1,1)*SQRT(-LAMBDA(2)/LAMBDA(1))+EV(1,2) + M(1)
      P(2) = -EV(2,1)*SQRT(-LAMBDA(2)/LAMBDA(1))+EV(2,2) + M(2)
      P(3) = -EV(3,1)*SQRT(-LAMBDA(2)/LAMBDA(1))+EV(3,2) + M(3)
*
      C0 = -(C1 * P(1) + C2 * P(2) + C3 * P(3))
      END
C
C
C
      SUBROUTINE SORT ( EV, LAMBDA, INULL, NOPOS, EPS )
C
C     SORTIEREN DREIER EIGENWERTE, SO DASS ZUERST DIE POSITIVEN,
C     DANN DIE NEGATIVEN UND DANN DIE EIGENWERTE, DIE NULL SIND
C     KOMMEN.
C
C     UEBERGABEPARAMETER :
C     EV         :    MATRIX, IN DENEN SPALTENWEISE DIE EIGENVEKTOREN
C                     STEHEN
C     LAMBDA     :    VEKTOR MIT DEN DREI EIGENWERTEN
C     INULL      :    ANZAHL DER EIGENWERTE, DIE NULL SIND
C     NOPOS      :    = TRUE: KEIN POSITIVER EIGENWERT
C     EPS        :    GENAUIGKEITSSCHRANKE
C***********************************************************************
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: EV( 3,3 ), LAMBDA( 3 )
      REAL(DP), INTENT(IN) :: EPS
      INTEGER, INTENT(OUT) :: INULL
      LOGICAL, INTENT(OUT) :: NOPOS

      REAL(DP) :: HILFL( 3 ), HILFEV( 3,3 ), MERK
      INTEGER :: J, L, I, K
C
C
C     AUF HILF WERDEN ZUNAECHST DIE EIGENWERTE, DIE POSITIV SIND,
C     ZUGEWIESEN.
C
C
      NOPOS = .TRUE.
      J = 0
      DO 10, I=1,3
         IF ( LAMBDA( I ) .GT. EPS ) THEN
            NOPOS = .FALSE.
            HILFL( J+1 ) = LAMBDA( I )
            DO 15, K = 1,3
               HILFEV( K,J+1  ) = EV( K,I )
   15       CONTINUE
            J = J + 1
         ENDIF
   10 CONTINUE
C
C     DIE POSITIVEN EIGENWERTE UND DIE DAZUGEHOERIGEN EIGENVEKTOREN
C     WERDEN SO SORTIERT, DASS LAMBDA( 1 ) <= LAMBDA( 2 ) IST
C
      IF ( J .GE. 2 ) THEN
C
         DO 17, K = 1,J-1
            DO 17, L = K+1,J
C
               IF ( HILFL( K ) .GT. HILFL( L ) ) THEN
                  DO 16, I =1,3
                     MERK         = HILFEV( I,K )
                     HILFEV( I,K )= HILFEV( I,L )
                     HILFEV( I,L )= MERK
   16             CONTINUE
                  MERK       = HILFL( K )
                  HILFL( K ) = HILFL( L )
                  HILFL( L ) = MERK
               ENDIF
   17    CONTINUE
C
      ENDIF
C
C
C     NUN DIE NEGATIVEN EIGENWERTE
      DO 20, I=1,3
C
C
         IF ( ( ABS( LAMBDA( I ) ) .GT. EPS ) .AND.
     >        (      LAMBDA( I )   .LT. EPS ) ) THEN
C
            HILFL( J+1 ) = LAMBDA( I )
            DO 25, K = 1,3
               HILFEV( K,J+1 ) = EV( K,I )
   25       CONTINUE
            J = J + 1
         ENDIF
   20 CONTINUE
C
C     Die Eigenwerte, die Null sind ans Ende
      DO 30, I = 1,3
         IF (ABS(LAMBDA(I)).LT.EPS) THEN
            HILFL(J+1) = LAMBDA(I)
            DO 35, K = 1,3
               HILFEV( K,J+1 ) = EV(K,I)
   35       CONTINUE
            J = J + 1
         ENDIF
   30 CONTINUE
C
C     HILF AUF LAMBDA UEBERSCHREIBEN
      DO 40, I=1,3
         LAMBDA( I ) = HILFL( I )
         DO 45, K = 1,3
            EV( K,I ) = HILFEV( K,I )
   45    CONTINUE
   40 CONTINUE
C
      INULL = 3 - J
      END
*DK EIRMAT
c
c  SUBROUTINE SPLINE(X,Y,N,A,B,C,D)
c  SUBROUTINE SORT ( EV, LAMBDA, INULL, NOPOS, EPS )
c  SUBROUTINE PAREBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3,EPS)
c  SUBROUTINE DOPEBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,EPS)
C  SUBROUTINE SCHEBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3)
C  SUBROUTINE ELLZYL(EV,LAMBDA,NUE,M,X0,Y0,Z0,CX,CY,CZ,R,INDE,EPS)
C  SUBROUTINE KEGEL (EV,LAMBDA,NUE,M,X0,Y0,Z0,CX,CY,CZ,R,EPS)
C  SUBROUTINE KUGEL (EV,LAMBDA,NUE,M,X0,Y0,Z0,CX,CY,CZ,R,EPS)
C  SUBROUTINE FL2O(A00,A1,A2,A3,A4,A5,A6,A7,A8,A9,INDE,X0,Y0,Z0,
C .                CX,CY,CZ,R,B0,B1,B2,B3,C0,C1,C2,C3,EPSIN,NMACH)
C
C  SUBROUTINE MA20A(Q,D,A,R,S,IQ,M,N,TOLER)
C  SUBROUTINE EA03A(A,B,N,ND,E)
C
c  24.4.95: imsl routinen fuer spline (iqhscu,....) raus, ersetzt durch : SPLINE
c  24.4.95: imsl routine mmdei (ellipt. integral) raus, redundant
c  24.4.95: inter raus, (interpolation in rechtecknetzen???), redundant
C
C*DK SPLINE
      SUBROUTINE SPLINE(X,Y,N,A,B,C,D)
C
C  P(X)=A(I)+B(I)*(X-X(I))+C(I)*(X-X(I))**2+D(I)*(X-X(I))**3
C  FUER X(I) <= X < X(I+1)       (D.H. A(I)=Y(I))
C
C  P(X)=A(1)+B(1)*(X-X(1))
C  FUER  X <= X(1)
C
C  P(X)=A(N)+B(N)*(X-X(N))
C  FUER X(N) <= X
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X(*), Y(*)
      REAL(DP), INTENT(OUT) :: A(*), B(*), C(*), D(*)
      INTEGER, INTENT(IN) :: N
      REAL(DP) :: FAC, DX1, DX2, DY1, DY2, DDX, H
      INTEGER :: I, J, K, NMH, NM1

      NM1=N-1
C  SETZE C-ARRAY (LOESE TRIDIAGONALE MATRIX)
      A(1)=2.E0
      C(1)=0.E0
      C(N)=0.E0
      D(1)=0.E0
      B(1)=0.E0
      DX2=X(2)-X(1)
      DY2=(Y(2)-Y(1))/DX2
      DO 1 I=2,NM1
      DX1=DX2
      DY1=DY2
      DX2=X(I+1)-X(I)
      DY2=(Y(I+1)-Y(I))/DX2
      DDX=DX1+DX2
      B(I)=DX2/DDX
      A(I)=1.E0-B(I)
      D(I)=6.E0*(DY2-DY1)/DDX
    1 CONTINUE
      DO 2 I=2,NM1
      J=I-1
      FAC=A(I)/A(J)
      D(I)=D(I)-FAC*D(J)
      A(I)=2.E0-B(J)*FAC
    2 CONTINUE
      K=N
      NMH=N
      DO 3 I=2,NM1
      K=K-1
      J=NMH-1
      C(J)=(D(K)-B(K)*C(NMH))/A(K)
      NMH=J
    3 CONTINUE
C  SETZE A,B,C UND D-ARRAY UND STEIGUNG BEI X=X(N)
      DO 4 J=2,NM1
4     C(J)=C(J)/2.
      DO 5 J=1,NM1
      H=X(J+1)-X(J)
      A(J)=Y(J)
      B(J)=(Y(J+1)-Y(J))/H
      B(J)=B(J)-H/3.*(C(J+1)+2.*C(J))
      D(J)=(C(J+1)-C(J))/(3.*H)
5     CONTINUE
      A(N)=Y(N)
      B(N)=B(NM1)+2.*C(NM1)*(X(N)-X(NM1))+3.*D(NM1)*(X(N)-X(NM1))**2
      RETURN
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      USE PRECISION
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      REAL(DP) D(N),E(N),Z(NM,N)
      REAL(DP) B,C,F,G,H,P,R,S,MACHEP,SMACH
      REAL(DP) DSQRT,DABS,DSIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT:
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1;
C
C        E HAS BEEN DESTROYED;
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
      MACHEP = SMACH ( 1 )
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))
         IF (B .LT. H) B = H
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::
         DO 110 M = L, N
            IF (ABS(E(M)) .LE. B) GO TO 120
C     :::::::::: E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ::::::::::
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     :::::::::: FORM SHIFT ::::::::::
         P = (D(L+1) - D(L)) / (2.0D0 * E(L))
         R = DSQRT(P*P+1.0D0)
         H = D(L) - E(L) / (P + SIGN(R,P))
C
         DO 140 I = L, N
  140    D(I) = D(I) - H
C
         F = F + H
C     :::::::::: QL TRANSFORMATION ::::::::::
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         MML = M - L
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = DSQRT(C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 160
  150       C = P / E(I)
            R = DSQRT(C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     :::::::::: FORM VECTOR ::::::::::
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         E(L) = S * P
         D(L) = C * P
         IF (ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = L
 1001 RETURN
C     :::::::::: LAST CARD OF TQL2 ::::::::::
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TQLRAT(N,D,E2,IERR)
CCCCC
CCCCC     AENDERUNG:    REAL ---> REAL(DP)
CCCCC
      USE PRECISION
      IMPLICIT NONE
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL(DP) D(N),E2(N)
      REAL(DP) B,C,F,G,H,P,R,S,MACHEP,SMACH
C     REAL SQRT,ABS,SIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     MODIFIED FOR CFT BY MIKE ESS, CRI, JULY 1980
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
C
C     MACHEP IS THE SMALLEST MACHINE REPRESENTABLE REAL NUMBER SUCH
C      THAT 1.+MACHEP .NE. 1.
C
      MACHEP = SMACH( 1 )
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0
      B = 0.0
      E2(N) = 0.0
C
      DO 290 L = 1, N
         J = 0
C        H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
         H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
         IF (B .GT. H) GO TO 105
         B = H
         C = B * B
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
C        S = SQRT(E2(L))
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0 * S)
C        R = SQRT(P*P+1.0)
         R = SQRT(P*P+1.0_DP)
C        D(L) = S / (P + SIGN(R,P))
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     ********** RATIONAL QL TRANSFORMATION **********
         G = D(M)
         IF (G .EQ. 0.0) G = B
         H = G
         S = 0.0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
         IF (H .EQ. 0.0) GO TO 210
C        IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0) GO TO 130
  210    P = D(L) + F
C     ********** ORDER EIGENVALUES **********
         IF (L .EQ. 1) GO TO 250
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQLRAT **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
      USE PRECISION
      IMPLICIT NONE
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL(DP) A(NM,N),D(N),E(N),E2(N)
      REAL(DP) F,G,H,SCALE
      REAL(DP) DSQRT,DABS,DSIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT:
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
  100 D(I) = A(I,I)
C     :::::::::: FOR I=N STEP -1 UNTIL 1 DO -- ::::::::::
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K))
C
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290
C
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            H = H + A(I,K) * A(I,K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = A(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 270
         F = 0.0D0
C
         DO 240 J = 1, L
            G = 0.0D0
C     :::::::::: FORM ELEMENT OF A*U ::::::::::
            DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
C     :::::::::: FORM ELEMENT OF P ::::::::::
  220       E(J) = G / H
            F = F + E(J) * A(I,J)
  240    CONTINUE
C
         H = F / (H + H)
C     :::::::::: FORM REDUCED A ::::::::::
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - H * F
            E(J) = G
C
            DO 260 K = 1, J
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE
C
  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)
C
  290    H = D(I)
         D(I) = A(I,I)
         A(I,I) = H
  300 CONTINUE
C
      RETURN
C     :::::::::: LAST CARD OF TRED1 ::::::::::
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C
      USE PRECISION
      IMPLICIT NONE
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL(DP) A(NM,N),D(N),E(N),Z(NM,N)
      REAL(DP) F,G,H,HH,SCALE
      REAL(DP) DSQRT,DABS,DSIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT:
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;
C
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
C          PRODUCED IN THE REDUCTION;
C
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 320
C     :::::::::: FOR I=N STEP -1 UNTIL 2 DO -- ::::::::::
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
C
         IF (SCALE .NE. 0.0_DP) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
C
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0D0
C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / (SCALE * H)
            G = 0.0D0
C     :::::::::: FORM ELEMENT OF A*U ::::::::::
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     :::::::::: FORM ELEMENT OF P ::::::::::
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     :::::::::: FORM REDUCED A ::::::::::
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
C
         DO 280 K = 1, L
  280    Z(I,K) = SCALE * Z(I,K)
C
  290    D(I) = H
  300 CONTINUE
C
  320 D(1) = 0.0D0
      E(1) = 0.0D0
C     :::::::::: ACCUMULATION OF TRANSFORMATION MATRICES ::::::::::
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0D0) GO TO 380
C
         DO 360 J = 1, L
            G = 0.0D0
C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0D0
         IF (L .LT. 1) GO TO 500
C
         DO 400 J = 1, L
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  400    CONTINUE
C
  500 CONTINUE
C
      RETURN
C     :::::::::: LAST CARD OF TRED2 ::::::::::
      END
