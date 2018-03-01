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
