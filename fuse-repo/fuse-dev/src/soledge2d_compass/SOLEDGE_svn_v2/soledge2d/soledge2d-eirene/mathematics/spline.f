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
      SUBROUTINE EIRENE_SPLINE(X,Y,N,A,B,C,D)
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
      USE EIRMOD_PRECISION
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
