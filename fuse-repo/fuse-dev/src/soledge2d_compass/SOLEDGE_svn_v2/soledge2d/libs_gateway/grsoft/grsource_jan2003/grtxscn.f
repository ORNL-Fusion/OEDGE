C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C   GLIGKS: mit \\a,\\o,\\u,\\s,@A,@O,@U
C           Alle Umlaute mit @ moeglich (16.9.96)
C   Update: Busch 24.394
C           mit ^a,^u,^o,^s Kleinbuchstaben
C           die doppelten Backslashes sollen nur uebergangsweise
C           erlaubt sein und werden in der Doku nicht mehr erwaehnt
C   Update: 16.9.96 ^ aus Abfrage streichen, da das ^ selbst sonst nicht
C           ausgegeben werden kann
C   XGKS:   keine Umlaute
C
C***********************************************************************
      SUBROUTINE GRTXSCN(STR1, L1, STR2, L2)
      CHARACTER STR1*(*), STR2*(*)
c
      CHARACTER BS
      CHARACTER ASCII(7), GKSTYP*7
      INTEGER GERMAN(7)
      LOGICAL ESCAPE
      COMMON /GRGKS/ GKSTYP
CDEC$ PSECT /GRGKS/ NOSHR
      SAVE /GRGKS/
      DATA ASCII /'A','O','U','a','o','u','s'/
      DATA GERMAN /196,214,220,228,246,252,223/
c
      IF (GKSTYP .EQ. 'GLIGKS') THEN
        I = 0
        L2 = 0
        BS = CHAR(92)
        ESCAPE = .FALSE.
        DO 100, I = 1, L1
          IF (ESCAPE) THEN
            DO 200, J = 1, 7
              IF (STR1(I:I) .EQ. ASCII(J)) THEN
                L2 = L2 + 1
                STR2(L2:L2) = CHAR(GERMAN(J))
              END IF

  200       CONTINUE
            ESCAPE = .FALSE.

C Busch 16.9.96
C         ELSE IF (STR1(I:I).EQ.BS .OR. STR1(I:I).EQ.'@'
c    $            .OR. STR1(I:I).EQ.'^' ) THEN
          ELSE IF (STR1(I:I).EQ.BS .OR. STR1(I:I).EQ.'@') THEN
            ESCAPE = .TRUE.

          ELSE
            L2 = L2 + 1
            STR2(L2:L2) = STR1(I:I)
          END IF

  100   CONTINUE
        STR2(L2+1:L2+1) = CHAR(0)
CJHeinen
      ELSE
        L2 = L1
        STR2(:L2) = STR1(:L1)
      END IF

      END
