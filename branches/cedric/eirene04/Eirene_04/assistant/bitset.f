C
C
      SUBROUTINE BITSET(IBIT,N1LOW,N1UP,IROW,JCOL,ISET,NBITS)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N1LOW, N1UP, IROW, JCOL, ISET, NBITS
      INTEGER, INTENT(INOUT) :: IBIT(N1LOW:N1UP,*)
      INTEGER :: JELEM, JB

      JELEM=JCOL/NBITS
      JB=MOD(JCOL,NBITS)
      IF (JB == 0) THEN
        JB=NBITS-1
      ELSE
        JELEM=JELEM+1
        JB=JB-1
      END IF

      IF (ISET == 0) IBIT(IROW,JELEM) = IBCLR(IBIT(IROW,JELEM),JB)
      IF (ISET == 1) IBIT(IROW,JELEM) = IBSET(IBIT(IROW,JELEM),JB)

      RETURN
      END