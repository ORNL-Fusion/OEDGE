 
 
 
C-----------------------------------------------------------------------
           SUBROUTINE EIRENE_OPRAND(AUSDRU,AKTLEN,OANDEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     UEBERPRUEFEN AUF ZULAESSIGE OPERANDEN
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
 
C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(26), PARAMETER ::
     P                  BUCHST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU
 
         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD
 
C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN
 
C
C     AUSGABEPARAMETER :
C
         INTEGER, INTENT(OUT) :: OANDEN
C           : ANZAHL DER OPERANDEN IN AUSDRU
 
C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS
 
 
      OANDEN=0
      AUSDRU(AKTLEN+1:AKTLEN+2)='  '
 
      I=1
C
C     WHILE-1 : SOLANGE AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C               UND KEIN FEHLER AUFGETRETEN
C
20    IF (I .LE. AKTLEN  .AND.  ERROR .EQ. 0 ) THEN
         POS= INDEX (BUCHST, AUSDRU(I:I) )
C
C        WHILE-2 : SOLANGE ANFANG VON EINEM OPERANDEN GEFUNDEN
C                  UND AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C                  UND KEIN FEHLER AUFGETRETEN
C
10       IF (POS .GT. 0  .AND. I .LE. AKTLEN .AND. ERROR .EQ. 0 ) THEN
            IF (INDEX(BUCHST, AUSDRU(I+1:I+1)) .GT. 0  .AND.
     >          INDEX(BUCHST, AUSDRU(I+2:I+2)) .GT. 0 ) THEN
C
C              OPERAND BESTEHT AUS MEHR ALS 2 BUCHSTABEN
C               => ZU LANG
C
               ERROR=6
            ELSE
C
C              GEHE WEITER IM STRING
C
               I=I+2
               POS=INDEX( BUCHST, AUSDRU(I:I) )
               OANDEN=OANDEN+1
            ENDIF
            GOTO 10
         ENDIF
C
C        ENDWHILE-2
C
C        KEIN ANFANG VON EINEM OPERANDEN GEFUNDEN, GEHE WEITER IM STRING
C
         I=I+1
         GOTO 20
      ENDIF
C
C     ENDWHILE-1
C
C     ENDE VON OPRAND
C
      END
