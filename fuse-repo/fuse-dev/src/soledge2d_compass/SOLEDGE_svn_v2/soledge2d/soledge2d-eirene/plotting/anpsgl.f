C
C
C----------------------------------------------------------------------*
C SUBROUTINE ANPSGL                                                    *
C----------------------------------------------------------------------*
      SUBROUTINE EIRENE_ANPSGL(MIN,MAX,MINL,MAXL,IERR)
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
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      REAL(SP), INTENT(IN) ::  MIN, MAX
      INTEGER, INTENT(INOUT) :: IERR
      INTEGER, INTENT(OUT) :: MINL, MAXL
      INTEGER :: EIRENE_IEXP10
C
      IF (MIN.LE.0.OR.MAX.LE.0) THEN
         WRITE(iunout,*)  '*** PARAMETERFEHLER IN ANPSGL ***'
         WRITE(iunout,*)  '***    MIN<=0  ODER  MAX<=0    ***'
         WRITE(iunout,*)  'MIN =',MIN,', MAX =',MAX
         IERR=IERR+1
         RETURN
      ENDIF
      IF (MIN.GT.MAX) THEN
         WRITE(iunout,*)  '*** PARAMETERFEHLER IN ANPSGL ***'
         WRITE(iunout,*)  '***        MIN  >  MAX         ***'
         WRITE(iunout,*)  'MIN= ',MIN,' MAX= ',MAX
         IERR=IERR+1
         RETURN
      ENDIF
      IF (MAX-MIN.GE.0) THEN
         MINL=EIRENE_IEXP10(REAL(MIN,KIND=DP))
         MAXL=EIRENE_IEXP10(REAL(MAX,KIND=DP))+1
      ENDIF
C
      RETURN
      END
