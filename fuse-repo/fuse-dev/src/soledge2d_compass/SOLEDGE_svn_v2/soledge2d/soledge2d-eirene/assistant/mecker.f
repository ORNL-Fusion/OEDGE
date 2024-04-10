 
 
 
C-----------------------------------------------------------------------
                SUBROUTINE EIRENE_MECKER(ERROR)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     AUSGABE DER REGELVERLETZUNGEN
C
C-----------------------------------------------------------------------
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
 
C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: ERROR
C           : FEHLERVARIABLE: > 0, FALLS EIN FEHLER AUFGETRETEN
 
 
      IF     (ERROR .EQ.  1) THEN
        WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT EIN UNGUELTIGES ZEICHEN.'
      ELSEIF (ERROR .EQ.  2) THEN
        WRITE(iunout,*)
     >  'ZWISCHEN ZWEI KLAMMERAUSDRUECKEN BEFINDET SICH',
     >  ' KEIN OPERATOR.'
      ELSEIF (ERROR .EQ.  3) THEN
         WRITE(iunout,*) 'EIN KLAMMERAUSDRUCK IST LEER.'
      ELSEIF (ERROR .EQ.  4) THEN
         WRITE(iunout,*) 'DER AUSDRUCK IST FALSCH GEKLAMMERT.'
      ELSEIF (ERROR .EQ.  5) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT MEHR ALS 3 INEINANDER',
     >              'GESCHACHTELTE KLAMMERN.'
      ELSEIF (ERROR .EQ.  6) THEN
        WRITE(iunout,*)
     >  'EIN OPERAND BESTEHT AUS MEHR ALS ZWEI BUCHSTABEN.'
      ELSEIF (ERROR .EQ.  7) THEN
         WRITE(iunout,*) 'ZWEI OPERATOREN STEHEN NEBENEINANDER.'
      ELSEIF (ERROR .EQ.  8) THEN
         WRITE(iunout,*) 'NACH EINEM OPERATOR FOLGT EINE SCHLIESSENDE',
     >              ' KLAMMER.'
      ELSEIF (ERROR .EQ.  9) THEN
        WRITE(iunout,*)
     >  'DAS ERSTES ZEICHEN DES AUSDRUCKS IST OPERATOR, ',
     >  'UND KEIN PRAEFIX.'
      ELSEIF (ERROR .EQ. 10) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENDET MIT EINEM OPERATOR.'
      ELSEIF (ERROR .EQ. 11) THEN
         WRITE(iunout,*) 'DAS ERSTE ZEICHEN IN EINEM KLAMMERAUSDRUCK',
     >              ' IST EIN OPERATOR.'
      ELSEIF (ERROR .EQ. 12) THEN
         WRITE(iunout,*) 'DER AUSDRUCK ENTHAELT MEHR ALS 15 OPERATOREN.'
      ELSEIF (ERROR .EQ. 13) THEN
         WRITE(iunout,*) 'DIE (ANZAHL DER OPERANDEN)-1 IST UNGLEICH ',
     >              ' DER (ANZAHL DER OPERATOREN).'
      ENDIF
C
C     ENDE VON MECKER
C
      END
