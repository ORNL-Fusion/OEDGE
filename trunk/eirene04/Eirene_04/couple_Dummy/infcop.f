*DK INFCOP
      SUBROUTINE INFCOP
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1, I2, I3
      ENTRY IF0COP
      ENTRY IF1COP
      ENTRY IF2COP(I1)
      ENTRY IF3COP(I1,I2,I3)
      ENTRY IF4COP
      END