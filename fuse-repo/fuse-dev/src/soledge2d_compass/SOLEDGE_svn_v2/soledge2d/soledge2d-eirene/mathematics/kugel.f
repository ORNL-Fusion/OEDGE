C
C
C
      SUBROUTINE EIRENE_KUGEL(A,LAMBDA,X0,Y0,Z0,CX,CY,CZ,R,EPS)
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
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: X0, Y0, Z0, CX, CY, CZ, R
      REAL(DP), INTENT(IN) :: LAMBDA(3), A(4,4), EPS
      REAL(DP) :: A3(3,3), DD, AA, D1, D2, D3, D4, ADD, EIRENE_SARRUS
C
C     DATA              EPS  / 5.D-10 /
cdr  noch nicht fertig
!     write (iunout,*) 'exit in subr. kugel '
!     call exit_own(1)
*
      A3 = A(1:3,1:3)
      DD = EIRENE_SARRUS(A3)
*
      CALL EIRENE_COFACT(A,A3,1,1)
      D1 = EIRENE_SARRUS(A3)
*
      CALL EIRENE_COFACT(A,A3,2,1)
      D2 = EIRENE_SARRUS(A3)
*
      CALL EIRENE_COFACT(A,A3,3,1)
      D3 = EIRENE_SARRUS(A3)
*
      CALL EIRENE_COFACT(A,A3,4,1)
      D4 = EIRENE_SARRUS(A3)
*
      AA = A(1,1)*D1 - A(2,1)*D2 + A(3,1)*D3 - A(4,1)*D4
 
      ADD = AA / (LAMBDA(1)*LAMBDA(2)*LAMBDA(3))
 
      CX = SQRT(-1.D0/LAMBDA(1)*ADD)
      CY = SQRT(-1.D0/LAMBDA(2)*ADD)
      CZ = SQRT(-1.D0/LAMBDA(3)*ADD)
      R = (CX + CY + CZ) / 3.D0
 
      A3(1:3,1) = A(1:3,4)
      A3(1:3,2:3) = A(1:3,2:3)
      X0 = - EIRENE_SARRUS(A3) / DD
 
      A3(1:3,1) = A(1:3,1)
      A3(1:3,2) = A(1:3,4)
      A3(1:3,3) = A(1:3,3)
      Y0 = - EIRENE_SARRUS(A3) / DD
 
      A3(1:3,1) = A(1:3,1)
      A3(1:3,2) = A(1:3,2)
      A3(1:3,3) = A(1:3,4)
      Z0 = - EIRENE_SARRUS(A3) / DD
 
      WRITE (iunout,*) ' IN KUGEL, ELLIPSOID '
      WRITE (iunout,*) ' X0,Y0,Z0 ',X0,Y0,Z0
      WRITE (iunout,*) ' CX,CY,CZ ',CX,CY,CZ
      WRITE (iunout,*) ' R ',R
 
      END
