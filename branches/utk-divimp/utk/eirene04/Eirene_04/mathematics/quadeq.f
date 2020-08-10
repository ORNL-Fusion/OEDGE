
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
