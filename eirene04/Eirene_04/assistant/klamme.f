


C-----------------------------------------------------------------------
             SUBROUTINE KLAMME(AUSDRU,AKTLEN,ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     UEBERPRUEFUNG AUF KORREKTE KLAMMERUNG
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

C
C     KONSTANTENDEKLARATION :
C
         INTEGER, PARAMETER :: KLAMAX=3
C           : MAXIMAL ZULAESSIGE ANZAHL VON KLAMMERSCHACHTELUNGEN

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
C     LOKALE VARIABLEN :
C
         INTEGER :: ANZAHL
C           : VARIABLE, DIE JE NACH IHREM WERT ANGIBT
C             = 0  ANZAHL DER OEFFNENDEN KLAMMERN IST GLEICH DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN
C             > 0  ANZAHL DER OEFFNENDEN KLAMMERN IST GROESSER DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN
C             < 0  ANZAHL DER OEFFNENDEN KLAMMERN IST KLEINER DER
C                  ANZAHL DER SCHLIESSENDEN KLAMMERN

         INTEGER :: MINKLA
C           : MINIMUM VON ANZAHL

         INTEGER :: MAXKLA
C           : MAXIMUM VON ANZAHL

         INTEGER :: KLAMON
C           : POSITION DER ZULETZ GEFUNDENEN OEFFNENDEN KLAMMER

         INTEGER :: KLAMOF
C           : POSITION DER ZULETZ GEFUNDENEN SCHLIESSENDEN KLAMMER

C
C     HILFSVARIABLE :
C
         INTEGER :: I


      ANZAHL=0
      MAXKLA=0
      MINKLA=0
      KLAMON=0
      KLAMOF=-1

      DO 10, I=1,AKTLEN
         IF (AUSDRU(I:I) .EQ. '(' .AND. ERROR .EQ. 0 ) THEN
            ANZAHL=ANZAHL+1
            MAXKLA=MAX(MAXKLA,ANZAHL)
            KLAMON=I

            IF (KLAMON .EQ. KLAMOF+1) THEN
C
C              ZWISCHEN ZWEI KLAMMERAUSDRUECKEN STEHT KEIN OPERATOR
C
               ERROR=2
            ENDIF
         ELSEIF (AUSDRU(I:I) .EQ. ')' .AND.  ERROR .EQ. 0 ) THEN
            ANZAHL=ANZAHL-1
            MINKLA=MIN(MINKLA,ANZAHL)
            KLAMOF=I

            IF (KLAMOF .EQ. KLAMON+1) THEN
C
C              DER KLAMMERAUSDRUCK IST LEER
C
               ERROR=3
            ENDIF
         ENDIF
10       CONTINUE

      IF ((ANZAHL .NE. 0  .OR. MINKLA .LT. 0) .AND. ERROR .EQ. 0) THEN
C
C        FALSCHE KLAMMERSETZUNG
C
         ERROR=4
      ELSEIF (MAXKLA .GT. KLAMAX .AND. ERROR .EQ. 0) THEN
C
C        ES WURDEN MEHR ALS KLAMAX KLAMMERN GESCHACHTELT
C
         ERROR=5
      ENDIF
C
C     ENDE VON KLAMME
C
      END
