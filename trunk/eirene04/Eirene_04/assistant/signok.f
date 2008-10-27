


C-----------------------------------------------------------------------
            SUBROUTINE SIGNOK(AUSDRU,AKTLEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     DAS UNTERPROGRAMM UEBERPRUEFT AUSDRU AUF ZULAESSIGE ZEICHEN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(33), PARAMETER :: 
     P                  GULTIG='ABCDEFGHIJKLMNOPQRSTUVWXYZ+-*/^()'

C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(IN) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN

C
C     HILFSVARIABLE :
C
         INTEGER :: I


      I=1
C
C     WHILE : SOLANGE KEIN FEHLER AUFGETRETEN UND STRING NOCH NICHT
C             VOLLSTAENDIG DURCHLAUFEN
C
10    IF (ERROR .EQ. 0  .AND.  I .LE. AKTLEN) THEN
         IF (INDEX(GULTIG,AUSDRU(I:I)) .GT. 0) THEN
            I=I+1
         ELSE
C
C           UNGUELTIGES ZEICHEN IN AUSDRU GEFUNDEN
C
            ERROR=1
         ENDIF
         GOTO 10
      ENDIF
C
C     ENDWHILE
C
C     ENDE VON SIGNOK
C
      END
