C
C
      SUBROUTINE SWRDRM (F,NROW,IX,IY,ISW)

      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NROW, IX, IY, ISW
      REAL(SP), INTENT(INOUT) :: F(NROW,*)
      REAL(SP) :: H
      INTEGER :: I, IS, J

      IF (MOD(ISW,2).EQ.1) THEN
        IS=IX/2
        DO 1 I=1,IS
        DO 1 J=1,IY
          H=F(I,J)
          F(I,J)=F(IX-I+1,J)
1         F(IX-I+1,J)=H
      ENDIF
      IF (ISW.GE.2) THEN
        IS=IY/2
        DO 2 J=1,IS
        DO 2 I=1,IX
          H=F(I,J)
          F(I,J)=F(I,IY-J+1)
2         F(I,IY-J+1)=H
      ENDIF
C
      RETURN
      END
