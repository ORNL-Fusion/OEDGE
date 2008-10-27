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
