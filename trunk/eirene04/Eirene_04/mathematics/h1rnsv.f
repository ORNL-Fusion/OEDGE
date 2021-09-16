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
