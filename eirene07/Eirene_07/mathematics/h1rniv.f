*
*
      SUBROUTINE H1RNIV(VEC)
*
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
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
      WRITE(iunout,*) 
     >           ' H1RNIV: H1RN (RANMAR) INITIALISED/RESTARTED WITH',
     >           ' SEED ARRAY VEC(100)'
*
      RETURN
      END
