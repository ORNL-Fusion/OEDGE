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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(*),Y(*),A(*),B(*),C(*),D(*)
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
C
C
      SUBROUTINE SORT ( EV, LAMBDA, INULL, NOPOS, EPS )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLEPRECISION   EV( 3,3 ), LAMBDA( 3 ), EPS,
     >                  HILFL( 3 ), HILFEV( 3,3 ), MERK
C
      INTEGER           INULL
      LOGICAL           NOPOS
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
C
C
C
      SUBROUTINE PAREBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLE PRECISION              B0,B1,B2,B3,C0,C1,C2,C3
      DOUBLEPRECISION   LAMBDA( 3 ), EV(3,3),EPS,
     >                  M(3),NUE,NORM,P(3),U(3)
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
C
C
C
      SUBROUTINE DOPEBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLEPRECISION   LAMBDA( 3 ), EV(3,3),EPS,
     >                  M(3),NUE,NORM,P(3)
      DOUBLE PRECISION              B0,B1,B2,B3
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
C
C
C
      SUBROUTINE SCHEBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLEPRECISION   LAMBDA( 3 ), EV(3,3),
     >                  M(3),NUE,NORM,P(3),U(3)
      DOUBLE PRECISION              B0,B1,B2,B3,C0,C1,C2,C3

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
      SUBROUTINE ELLZYL(EV,LAMBDA,NUE,M,X0,Y0,Z0,CX,CY,CZ,R,INDE,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLE PRECISION              X0, Y0, Z0, CX, CY, CZ, R
      DOUBLEPRECISION   LAMBDA(3), EV(3,3), EPS,
     >                  M(3), NUE, NORM, P(3)
C
      INTEGER           INDE
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
C
C
C
      SUBROUTINE KEGEL(EV,LAMBDA,NUE,M,X0,Y0,Z0,CX,CY,CZ,R,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLE PRECISION              X0, Y0, Z0, CX, CY, CZ, R
      DOUBLEPRECISION   LAMBDA(3), EV(3,3), EPS,
     >                  M(3), NUE, NORM
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLE PRECISION              X0, Y0, Z0, CX, CY, CZ, R
      DOUBLEPRECISION   LAMBDA(3), A(4,4), EPS
      DOUBLE PRECISION A3(3,3), DD, AA, D1, D2, D3, D4, ADD
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

      SUBROUTINE COFACT (A4,A3,I,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A4(4,4),A3(3,3)

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

      FUNCTION SARRUS (A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(3,3)

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
      SUBROUTINE FL2O(A00,A1,A2,A3,A4,A5,A6,A7,A8,A9,INDE,X0,Y0,Z0,
     >                CX,CY,CZ,R,B0,B1,B2,B3,C0,C1,C2,C3,EPSIN,NMACH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLEPRECISION   A( 3,3 ), LAMBDA( 3 ), EV( 3,3 ), EPS, C, B(3),
     >                  M(3),NUE,NORM,P(3),AMERK(3,3),RES(3),WRK(6)
     >                 ,Q(5,5),DD(3),QEV(5,5)
      DOUBLE PRECISION     A00,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,X0,Y0,Z0,
     >                  CX,CY,CZ,R,B0,B1,B2,B3,C0,C1,C2,C3
      DOUBLE PRECISION AM4(4,4),LAMORI(3)
C
      INTEGER           INDE, INULL, DIM, KBASIS, IW(3)
      LOGICAL           NOPOS
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
                  CALL KEGEL(EV,LAMBDA,0.D0,M,X0,Y0,Z0,CX,CY,CZ,R,EPS)
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
               CALL SCHEBE(EV,LAMBDA,0.D0,M,B0,B1,B2,B3,C0,C1,C2,C3)
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
                  CALL DOPEBE(EV,LAMBDA,0.D0,M,B0,B1,B2,B3,EPS)
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
C
C##       MA20A          28/06/72
C NAME MA20A(R)                  CHECK
      SUBROUTINE MA20A(Q,D,A,R,S,IQ,M,N,TOLER)
c slmod begin -f90
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION SUM
      DOUBLE PRECISION MIN,MAX
      INTEGER OUT,S(*)
      LOGICAL STAGE,TEST
      DOUBLE PRECISION Q(IQ,*),A(*),D(*),R(*)

      DATA BIG /1.D75/
c
c      IMPLICIT REAL*8 (A-H,O-Z)
c      REAL*8 SUM
c      REAL*8 MIN,MAX
c      INTEGER OUT,S(*)
c      LOGICAL STAGE,TEST
c      REAL*8 Q(IQ,*),A(*),D(*),R(*)
C  ***BIG MUST BE SET EQUAL TO ANY VERY LARGE REAL CONSTANT.
C  ***ITS VALUE HERE IS APPROPRIATE FOR THE IBM 370.
c      DATA BIG /1.E75/
c slmod end
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
C##       EA03A          20/03/79
C NAME EA03A(R)                  CHECK
      SUBROUTINE EA03A(A,B,N,ND,E)
      IMPLICIT REAL*8 (A-H,O-Z)
C STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
      DIMENSION A(1),B(1),ARMAX(200),JRMAX(200)
C PURPOSE FINDS ALL THE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C         BY JACOBI'S METHOD (CLASSICAL METHOD)
C
C ARGUMENTS
C         A(A TWO DIMENSIONAL ARRAY TO THE USER) MUST CONTAIN THE MATRIX
C          TREATED IN THE FIRST N ROWS AND COLUMNS OF THE ARRAY A.
C              ********* THE ORIGINAL MATRIX WILL BE DESTROYED
C              *********** THE EIGENVALUES WILL BE FOUND IN A(I,I) I=1,2
C         N THE DIMENSION OF THE MATRIX
C         B THE EIGENVECTORS WILL BE FOUND IN THE COLUMNS OF B
C         ND THE FIRST DIMENSION OF THE ARRAYS A AND B IN THE CALLING PR
C         MUST BE DIMENSION A(ND,  ),B(ND,  )
C         E ACCURACY CONTROL<<<< THE PROCESS CHOOSES THE OFF DIAGONAL EL
C         WITH THE LARGEST MODULUS AND PERFORMS AN ORTHOGONAL TRANSFORMA
C         REDUCE THE ELEMENT TO ZERO.THE PROCESS IS REPEATED UNTIL THE L
C         MAGNITUDE AMONG OFF DIAGONAL ELEMENT IS LESS THAN E.
C *******'A IS ASSUMED TO BE SYMMETRIC' ******
C***********************************************************************
C
C FINDS THE OFF DIAGONAL ELEMENT WITH THE MAXIMUM MODULUS
C
C***********************************************************************
      NDN = ND*N
      DO 1 K=1,NDN
      B(K) = 0.0
    1 CONTINUE
      DO 2 K=1,N
      KK = K*(ND+1)-ND
      ARMAX(K) = 0.0
      B(KK) = 1.0
      DO 3 L=K,N
      IF(L-K)4,3,4
    4 KL = K+ND*(L-1)
      Y = ABS(A(KL))
      IF(ARMAX(K)-Y)5,3,3
    5 ARMAX(K) = Y
      JRMAX(K) = L
    3 CONTINUE
    2 CONTINUE
   11 AMAX = 0.0
      DO 6 K=1,N
      Y = ABS(ARMAX(K))
      IF(AMAX-Y)7,6,6
    7 AMAX = Y
      I = K
    6 CONTINUE
      IF(E-AMAX)8,9,9
    8 NDI = ND*(I-1)
      J = JRMAX(I)
C***********************************************************************
C
C   (K)  (K-1)         (K-1)         (K)            MAX MODULUS IN P,Q
C  A   =     COS X + A     SIN X = A               POSITION PERFORMS
C   I,P  I,P           I,Q           P,I            AN ORTHOGONAL
C                                         I#P,Q    TRANSFORMATION TO
C   (K)   (K-1)        (K-1)         (K)            REDUCE ELEMENT TO 0
C  A   =-A     SIN X+ A     COS X = A               X=ANGLE OF ROTATION
C   I,Q   I,P          I,Q           Q,I
C
C   (K)  (K-1)   2      (K-1)              (K-1)   2        (K)    (K)
C  A   =A     COS X +2 A     COS X SIN X +A     SIN X      A   =A
C   P,P  P,P            P,Q                Q,Q              P,Q  Q,P
C
C   (K)  (K-1)   2      (K-1)              (K-1)   2
C  A   =A     SIN X -2 A     COS X SIN X +A     COS X
C   Q,Q  P,P            P,Q                Q,Q
C
      NDJ = ND*(J-1)
      II = I+NDI
      JJ = J+NDJ
      IJ = I+NDJ
      JI = J+NDI
      AII = A(II)
      AJJ = A(JJ)
      AIJ = A(IJ)
      Y = 2.0*AIJ
      X = AII-AJJ
      T = SIGN(1.0D0,X)*Y/(ABS(X)+SQRT(X**2+Y**2))
      TSQ = T**2
      C = 1.0/SQRT(1.0+TSQ)
      TY = T*Y
      S = T*C
      CSQ = C**2
      A(II) = CSQ*(AII+TY+AJJ*TSQ)
      A(JJ) = CSQ*(AJJ-TY+AII*TSQ)
      A(IJ) = 0.0
      A(JI) = 0.0
      DO 10 K=1,N
      JTES = (K-I)*(K-J)
      NDK = ND*(K-1)
      KI = K+NDI
      KJ = K+NDJ
      IF(JTES)13,12,13
   13 JK = J+NDK
      IK = I+NDK
      A(KI) = C*A(IK)+S*A(JK)
      A(KJ) =-S*A(IK)+C*A(JK)
      A(JK) = A(KJ)
      A(IK) = A(KI)
   12 X = B(KI)
      B(KI) = C*X+S*B(KJ)
      B(KJ) =-S*X+C*B(KJ)
   10 CONTINUE
C***********************************************************************
C
C FINDS THE MAXIMUM MODULUS OFF DIAGONAL ELEMENT BY MODIFIC ATION OF
C PREVIOUS INFORMATION OF OFF DIAGONAL ELEMENTS
C
C***********************************************************************
      ARMAX(I) = 0.0
      DO 14 K=1,N
      IF(K-I)15,14,15
   15 IK = I+ND*(K-1)
      Y = ABS(A(IK))
      IF(ARMAX(I)-Y)16,14,14
   16 ARMAX(I) = Y
      JRMAX(I) = K
   14 CONTINUE
      ARMAX(J) = 0.0
      DO 17 K=1,N
      IF(K-J)18,17,18
   18 JK = J+ND*(K-1)
      Y = ABS(A(JK))
      IF(ARMAX(J)-Y)19,17,17
   19 ARMAX(J) = Y
      JRMAX(J) = K
   17 CONTINUE
      DO 20 K=1,N
      ITES = (K-I)*(K-J)
      KI = K+NDI
      KJ = K+NDJ
      IF(ITES)21,20,21
   21 X = ABS(A(KI))
      Y = ABS(A(KJ))
      JR = J
      IF(X-Y)22,22,23
   23 Y = X
      JR = I
   22 IF(ARMAX(K)-Y)24,20,20
   24 ARMAX(K) = Y
      JRMAX(K) = JR
   20 CONTINUE
      GO TO 11
    9 RETURN
      END

      SUBROUTINE QUADEQ (A,B,C,T1,T2,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA EPS10 /1.E-10/

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
