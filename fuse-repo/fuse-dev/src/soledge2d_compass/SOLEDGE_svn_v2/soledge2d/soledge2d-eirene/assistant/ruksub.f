 
 
 
C-----------------------------------------------------------------------
        SUBROUTINE EIRENE_RUKSUB(TEIL, IPART, PART, IARITH, ARITH,
     .  HILFE)
C-----------------------------------------------------------------------
C
C     FUNKTION:
C
C     ERSETZEN VON '^' DURCH '**' UND EINFUEGEN VON 'Z' FUER
C     ZERLEGUNG
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
 
C
C     KONSTANTENDEKLARATION :
C
         CHARACTER(10), PARAMETER :: ZIFFER='0123456789'
 
C
C     EINGABEPARAMETER :
C
         INTEGER, INTENT(IN) :: TEIL
C           : AKTUELLE ANZAHL DER ZERLEGUNGEN
 
C
C     EIN/AUSGABEPARAMETER :
C
         CHARACTER(*), INTENT(INOUT) :: PART(*)
C           : FELD VON STRINGS, AUF DENEN DIE EINZELNEN
C             ELEMENTARZERLEGUNGEN FESTGEHALTEN WERDEN
 
         INTEGER, INTENT(INOUT) :: IPART(*)
C           : AKTUELLE LAENGEN VON PART(ZMAX)
 
         CHARACTER(*), INTENT(INOUT) :: ARITH(*)
C           : FELD VON STRINGS, AUF DENEN DIE TEIL-TE GENERATION
C             VON AUSDRU FESTGEHALTEN WIRD
 
         INTEGER, INTENT(INOUT) :: IARITH(*)
C           : AKTUELLE LAENGEN VON ARITH(ZMAX)
 
         CHARACTER(*), INTENT(INOUT) :: HILFE
C           : HILFSSTRING, DER FUER ZUWEISUNGEN BENOETIGT WIRD
 
C
C     HILFSVARIABLEN :
C
         INTEGER :: I, J
 
      DO 10, J=1,TEIL
C
C        ERSETZEN VON '^' DURCH '**' IN PART(J)
C
         I=2
20       IF (I .LT. IPART(J)) THEN
            IF (PART(J)(I:I) .EQ. '^') THEN
CPB            HILFE=PART(J)(1:I-1)//'**'//PART(J)(I+1:IPART(J))//' '
CPB            I=I+1
CPB            IPART(J)=IPART(J)+1
CPB            PART(J)=HILFE
            ENDIF
            I=I+1
            GOTO 20
         ENDIF
C
C        ERSETZEN VON '^' DURCH '**' IN ARITH(J)
C
         I=2
30       IF (I .LT. IARITH(J)) THEN
            IF (ARITH(J)(I:I) .EQ. '^') THEN
CPB            HILFE=ARITH(J)(1:I-1)//'**'//ARITH(J)(I+1:IARITH(J))//' '
CPB            I=I+1
CPB            IARITH(J)=IARITH(J)+1
CPB            ARITH(J)=HILFE
            ENDIF
            I=I+1
            GOTO 30
         ENDIF
C
C        EINFUEGEN  VON 'Z' FUER ZERLEGUNG IN PART(J)
C
         I=1
40       IF (I .LT. IPART(J)) THEN
            IF (INDEX(ZIFFER,PART(J)(I:I)) .GT. 0) THEN
               IF (I .EQ. 1) THEN
                  HILFE='Z'//PART(J)(I:IPART(J))//' '
               ELSE
                  HILFE=PART(J)(1:I-1)//'Z'//
     >                   PART(J)(I:IPART(J))//' '
               ENDIF
 
               I=I+2
               IPART(J)=IPART(J)+1
               PART(J)=HILFE
            ENDIF
            I=I+1
            GOTO 40
         ENDIF
C
C        EINFUEGEN VON 'Z' FUER ZERLEGUNG IN ARITH(J)
C
         I=1
50       IF (I .LT. IARITH(J)) THEN
            IF (INDEX(ZIFFER,ARITH(J)(I:I)) .GT. 0) THEN
               IF (I .EQ. 1) THEN
                  HILFE='Z'//ARITH(J)(I:IARITH(J))//' '
               ELSE
                  HILFE=ARITH(J)(1:I-1)//'Z'//
     >                   ARITH(J)(I:IARITH(J))//' '
               ENDIF
 
               I=I+2
               IARITH(J)=IARITH(J)+1
               ARITH(J)=HILFE
            ENDIF
            I=I+1
            GOTO 50
         ENDIF
 
10       CONTINUE
C
C     ENDE VON RUKSUB
C
      END
