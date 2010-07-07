C
C
*//NEUTR//
C=======================================================================
C          S U B R O U T I N E   N E U T R
C=======================================================================
      SUBROUTINE NEUTR(KARD,NDIMX,NDIMY,NDIMF,DUMMY,LDMX,LDMY,LDMF,
     .                 LDNS,IS)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KARD,NDIMX,NDIMY,NDIMF,LDMX,LDMY,LDMF,
     .                       LDNS,IS
      INTEGER :: ND1,LIM,IX,IY,III,IF
      REAL(DP), INTENT(IN) :: DUMMY(0:LDMX+1,0:LDMY+1,LDMF,LDNS)
C
      ND1 = NDIMX
      LIM = (ND1/5)*5 - 4
      DO  500  IF = 1,NDIMF
        DO  110  IY = 1,NDIMY
          DO  100  IX = 1,LIM,5
  100     WRITE(KARD,910) (DUMMY(IX-1+III,IY,IF,IS),III = 1,5)
          IF( (LIM+4).EQ.ND1 )   GOTO 110
          WRITE(KARD,910) (DUMMY(IX,IY,IF,IS),IX = LIM+5,ND1)
  110   CONTINUE
  500 CONTINUE
      RETURN
  910 FORMAT(5(E16.8))
*//END NEUTR//
      END
