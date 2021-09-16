C
C
C----------------------------------------------------------------------*
C SUBROUTINE ANPSGL                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE ANPSGL(MIN,MAX,MINL,MAXL,IERR)
C
C  DIESES PROGRAMM ERHAELT ZWEI INTERVALLGRENZEN (MIN,MAX) UND
C  BERECHNET WELCHE ZEHNERPOTENZ ALS LINKE INTERVALLGRENZE UND
C  WELCHE ALS RECHT IN FRAGE KOMMT (MINL,MAXL).
C  DAS INTERVALL (MINL,MAXL) EIGNET SICH DANN ZUR LOGARTHMISCHEN
C  SKALIERUNG EINER KOORDINATENACHSE.
C
C  EINGABE :  MIN   DEC DBLE(6)
C             MAX   DEC DBLE(6)
C  AUSGABE :  MINL  BIN FIXED(15)
C             MAXL  BIN FIXED(15)                                  *
C
      USE PRECISION

      IMPLICIT NONE

      REAL(SP), INTENT(IN) ::  MIN, MAX
      INTEGER, INTENT(INOUT) :: IERR
      INTEGER, INTENT(OUT) :: MINL, MAXL
      INTEGER :: IEXP10
C
      IF (MIN.LE.0.OR.MAX.LE.0) THEN
         WRITE(6,*)  '*** PARAMETERFEHLER IN ANPSGL ***'
         WRITE(6,*)  '***    MIN<=0  ODER  MAX<=0    ***'
         WRITE(6,*)  'MIN =',MIN,', MAX =',MAX
         IERR=IERR+1
         RETURN
      ENDIF
      IF (MIN.GT.MAX) THEN
         WRITE(6,*)  '*** PARAMETERFEHLER IN ANPSGL ***'
         WRITE(6,*)  '***        MIN  >  MAX         ***'
         WRITE(6,*)  'MIN= ',MIN,' MAX= ',MAX
         IERR=IERR+1
         RETURN
      ENDIF
      IF (MAX-MIN.GE.0) THEN
         MINL=IEXP10(REAL(MIN,KIND=DP))
         MAXL=IEXP10(REAL(MAX,KIND=DP))+1
      ENDIF
C
      RETURN
      END
