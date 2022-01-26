C
C
      SUBROUTINE INVERT(A,AI)

      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A(3,3)
      REAL(DP), INTENT(OUT) :: AI(3,3)
      REAL(DP) :: DET
C
      DET = A(1,1) * A(2,2) * A(3,3) + A(1,2) * A(2,3) * A(3,1)
     .    + A(1,3) * A(2,1) * A(3,2) - A(1,3) * A(2,2) * A(3,1)
     .    - A(1,2) * A(2,1) * A(3,3) - A(1,1) * A(2,3) * A(3,2)
C
C     WRITE (*,*) DET
C
      AI(1,1) = (A(2,2) * A(3,3) - A(2,3) * A(3,2)) / DET
      AI(1,2) = (A(1,3) * A(3,2) - A(1,2) * A(3,3)) / DET
      AI(1,3) = (A(1,2) * A(2,3) - A(2,2) * A(1,3)) / DET
      AI(2,1) = (A(2,3) * A(3,1) - A(2,1) * A(3,3)) / DET
      AI(2,2) = (A(1,1) * A(3,3) - A(1,3) * A(3,1)) / DET
      AI(2,3) = (A(1,3) * A(2,1) - A(1,1) * A(2,3)) / DET
      AI(3,1) = (A(2,1) * A(3,2) - A(2,2) * A(3,1)) / DET
      AI(3,2) = (A(1,2) * A(3,1) - A(1,1) * A(3,2)) / DET
      AI(3,3) = (A(1,1) * A(2,2) - A(1,2) * A(2,1)) / DET
C
      RETURN
      END
