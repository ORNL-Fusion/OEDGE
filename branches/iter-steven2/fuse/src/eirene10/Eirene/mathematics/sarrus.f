

      FUNCTION SARRUS (A)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A(3,3)
      REAL(DP) :: SARRUS, A1, A2, A3, B1, B2, B3

      A1 = A(1,1) * A(2,2) * A(3,3)
      A2 = A(1,2) * A(2,3) * A(3,1)
      A3 = A(1,3) * A(2,1) * A(3,2)
      B1 = A(3,1) * A(2,2) * A(1,3)
      B2 = A(3,2) * A(2,3) * A(1,1)
      B3 = A(3,3) * A(2,1) * A(1,2)
      SARRUS = A1 + A2 + A3 - B1 - B2 - B3
!     WRITE (iunout,*) ' A1,A2,A3 ',A1,A2,A3
!     WRITE (iunout,*) ' B1,B2,B3 ',B1,B2,B3
!     WRITE (iunout,*) ' DET ',A1 + A2 + A3 - B1 - B2 - B3

      RETURN
      END
