C
C
C----------------------------------------------------------------------*
C FUNCTION IEXP10
C----------------------------------------------------------------------*
      FUNCTION EIRENE_IEXP10(ZAHL)
C
C  FUNKTION, WELCHE DEN 10-ER-EXPONENTEN EINER ZAHL
C  BERECHNET. (AUSNAHME: IEXP10(0)=0)
C
C  EINGABE: ZAHL             REAL
C  AUSGABE: IEXP10(ZAHL)     INTEGER
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: ZAHL
      REAL(DP) :: ZHL2
      INTEGER :: I, EIRENE_IEXP10
 
      I=0
      IF( ZAHL.EQ.0) THEN
         EIRENE_IEXP10=I
         RETURN
      ENDIF
      ZHL2=ABS(ZAHL)
  10  IF (ZHL2.LT.1) THEN
         I=I-1
         ZHL2=ZHL2*10
         GOTO 10
      ENDIF
  20  IF (ZHL2.GE.10) THEN
         I=I+1
         ZHL2=ZHL2/10
         GOTO 20
      ENDIF
      EIRENE_IEXP10=I
      RETURN
      END
