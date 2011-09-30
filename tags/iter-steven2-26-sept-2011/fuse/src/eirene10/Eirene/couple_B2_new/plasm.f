C
*//PLASM//
C=======================================================================
C          S U B R O U T I N E   P L A S M
C=======================================================================
      SUBROUTINE PLASM(KARD,NDIMX,NDIMY,NDIMF,N,M,NF,DUMMY)

      USE PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KARD, NDIMX, NDIMY, NDIMF, N, M, NF
      REAL(DP), INTENT(INOUT) :: DUMMY(0:N+1,0:M+1,NF)
      INTEGER :: ND1, LIM, IF, III, IX, IY

      ND1 = NDIMX + 2
      LIM = (ND1/5)*5 - 4
      DUMMY(0:N+1,0:M+1,NF)=0._DP
      DO    110  IF = 1,NDIMF
      DO    110  IY = 0,NDIMY+1
      DO    100  IX = 1,LIM,5
100     READ(KARD,910,END=500) (DUMMY(-1+IX-1+III,IY,IF),III = 1,5)
        IF( (LIM+4).EQ.ND1 )     GOTO 110
        READ(KARD,910,END=500) (DUMMY(-1+IX,IY,IF),IX = LIM+5,ND1)
110   CONTINUE
500   RETURN
910   FORMAT(5(E16.8))
*//END PLASM//
      END
