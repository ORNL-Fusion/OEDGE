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
