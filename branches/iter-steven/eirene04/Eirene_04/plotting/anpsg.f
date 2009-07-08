C
C
C----------------------------------------------------------------------*
C SUBROUTINE ANPSG                                                     *
C----------------------------------------------------------------------*
      SUBROUTINE ANPSG(MIN,MAX,INTNR,STPSZ,IERR)
C
C  DIESES PROGRAMM ERHAELT ZWEI INTERVALLGRENZEN (MIN,MAX) UND
C  FORMT DIESE MEIST "KRUMMEN" ZAHLEN UM IN SOLCHE, DIE SICH
C  ZUR ACHSENBESCHRIFTUNG EIGNEN.
C  INTNR GIBT GIBT ANZAHL DER TEILINTERVALLE ZWISCHEN MIN UND MAX
C  AN, DIE SICH AUS DER UMFORMUNG ERGEBEN.
C  STPSZ IST DIE LAENGE EINES SOLCHEN TEILINTERVALLES.
C
C  EINGABE :  MIN   DEC DBLE(6)
C             MAX   DEC DBLE(6)
C  AUSGABE :  MIN   DEC DBLE(6)
C             MAX   DEC DBLE(6)
C             INTNR BIN FIXED(15)
C             STPSZ DEC DBLE(6)
C
      USE PRECISION

      IMPLICIT NONE
      REAL(SP), INTENT(INOUT) :: MAX, MIN
      REAL(SP), INTENT(OUT) :: STPSZ
      INTEGER, INTENT(INOUT) :: IERR
      INTEGER, INTENT(OUT) :: INTNR
      REAL(DP) :: LKS, RTS, MDDL, STPS2, DIFF
      INTEGER :: IPT, IEXP10
C
      DIFF=MAX-MIN
      IF (DIFF.EQ.0) THEN
         IPT=IEXP10(REAL(MAX,KIND=DP))
         MIN=MIN-10.**(IPT-1)
         MAX=MAX+10.**(IPT-1)
         DIFF=MAX-MIN
      ENDIF
      IF (DIFF.GT.0) THEN
         STPS2=DIFF/10
         IPT=IEXP10(STPS2)
         STPS2=STPS2/10.**IPT
         STPS2=AINT(STPS2+0.5)*10.**IPT
         IPT=IEXP10(DIFF)
         MDDL=(MAX+MIN)*0.5/10.**(IPT-1)
         MDDL=AINT(MDDL)*10.**(IPT-1)
         LKS=MDDL
         RTS=MDDL
         INTNR=0
    5    IF (LKS.GT.MIN) THEN
            INTNR=INTNR+1
            LKS=LKS-STPS2
            GOTO 5
         ENDIF
   10    IF (RTS.LT.MAX) THEN
            INTNR=INTNR+1
            RTS=RTS+STPS2
            GOTO 10
         ENDIF
         IPT=IEXP10(LKS)
         MIN=LKS+SIGN(1._DP,LKS)*10.D0**(-14+IPT)
         IPT=IEXP10(RTS)
         MAX=RTS+SIGN(1._DP,RTS)*10.D0**(-14+IPT)
         IPT=IEXP10(STPS2)
         STPSZ=STPS2+10.**(-14+IPT)
      ELSE
         WRITE(6,*)  '------------------------------------'
         WRITE(6,*)  'PARAMETERFEHLER IN ANPSG: MIN > MAX'
         WRITE(6,*)  'MIN =',MIN,', MAX =',MAX
         IERR=IERR+1
         RETURN
      ENDIF
C
      RETURN
      END
