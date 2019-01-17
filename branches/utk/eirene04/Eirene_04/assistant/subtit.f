


C-----------------------------------------------------------------------
               SUBROUTINE SUBTIT(AUSDRU, AKTLEN, HILFE)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     ELIMINATION VON BLANKS UND SUBTITUTION VON '**' DURCH '^'
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     EIN/AUSGABEPARAMETER :
C
         INTEGER, INTENT(INOUT) :: AKTLEN
C           : AKTUELLE LAENGE VON AUSDRU

         CHARACTER(*), INTENT(INOUT) :: AUSDRU
C           : AUSDRUCK, DER IM UNTERPROGRAMM ZERLEGT WIRD
 
         CHARACTER(*), INTENT(INOUT) :: HILFE
C           : HILFSSTRING, DER FUER ZUWEISUNGEN BENOETIGT WIRD

C
C     HILFSVARIABLEN :
C
         INTEGER :: I, POS

C
C     ELIMINATION VON BLANKS
C
      I=1
C
C     WHILE-1 : SOLANGE AUSDRU NOCH NICHT VOLLSTAENDIG DURCHLAUFEN
C
30    IF (I .LT. AKTLEN) THEN
C
C        ENTFERNEN VON BLANKS
C
         IF (AUSDRU(I:I) .EQ. ' ') THEN
            IF ( I .EQ. 1) THEN
               HILFE=AUSDRU(I+1:AKTLEN) // ' '
            ELSE
               HILFE=AUSDRU(1:I-1) // AUSDRU(I+1:AKTLEN) // ' '
            ENDIF
            AUSDRU=HILFE
            AKTLEN=AKTLEN-1
         ELSE
            I=I+1
         ENDIF
         GOTO 30
      ENDIF
C
C     ENDWHILE-1
C
C
C     SUBTITUTION VON '**' DURCH '^'
C
      POS=INDEX(AUSDRU,'**')
C
C     WHILE-2 : SOLANGE AUSDRU NOCH '**' ENTHAELT
C
10    IF ( POS .GT. 1 .AND. POS .LT. AKTLEN-1) THEN
C
C        ERSETZEN VON '**' DURCH '^'
C
         HILFE=AUSDRU(1:POS-1) // '^' // AUSDRU(POS+2:AKTLEN)
         AUSDRU=HILFE
         AKTLEN=AKTLEN-1
         POS=INDEX(AUSDRU,'**')
         GOTO 10
      ENDIF
C
C     ENDWHILE-2
C
C     ENDE VON SUBTIT
C
      END
