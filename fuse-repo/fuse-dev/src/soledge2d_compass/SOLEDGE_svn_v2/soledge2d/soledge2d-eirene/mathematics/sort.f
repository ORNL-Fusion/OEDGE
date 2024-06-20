C
C
C
      SUBROUTINE EIRENE_SORT ( EV, LAMBDA, INULL, NOPOS, EPS )
C
C     SORTIEREN DREIER EIGENWERTE, SO DASS ZUERST DIE POSITIVEN,
C     DANN DIE NEGATIVEN UND DANN DIE EIGENWERTE, DIE NULL SIND
C     KOMMEN.
C
C     UEBERGABEPARAMETER :
C     EV         :    MATRIX, IN DENEN SPALTENWEISE DIE EIGENVEKTOREN
C                     STEHEN
C     LAMBDA     :    VEKTOR MIT DEN DREI EIGENWERTEN
C     INULL      :    ANZAHL DER EIGENWERTE, DIE NULL SIND
C     NOPOS      :    = TRUE: KEIN POSITIVER EIGENWERT
C     EPS        :    GENAUIGKEITSSCHRANKE
C***********************************************************************
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: EV( 3,3 ), LAMBDA( 3 )
      REAL(DP), INTENT(IN) :: EPS
      INTEGER, INTENT(OUT) :: INULL
      LOGICAL, INTENT(OUT) :: NOPOS
 
      REAL(DP) :: HILFL( 3 ), HILFEV( 3,3 ), MERK
      INTEGER :: J, L, I, K
C
C
C     AUF HILF WERDEN ZUNAECHST DIE EIGENWERTE, DIE POSITIV SIND,
C     ZUGEWIESEN.
C
C
      NOPOS = .TRUE.
      J = 0
      DO 10, I=1,3
         IF ( LAMBDA( I ) .GT. EPS ) THEN
            NOPOS = .FALSE.
            HILFL( J+1 ) = LAMBDA( I )
            DO 15, K = 1,3
               HILFEV( K,J+1  ) = EV( K,I )
   15       CONTINUE
            J = J + 1
         ENDIF
   10 CONTINUE
C
C     DIE POSITIVEN EIGENWERTE UND DIE DAZUGEHOERIGEN EIGENVEKTOREN
C     WERDEN SO SORTIERT, DASS LAMBDA( 1 ) <= LAMBDA( 2 ) IST
C
      IF ( J .GE. 2 ) THEN
C
         DO 17, K = 1,J-1
            DO 17, L = K+1,J
C
               IF ( HILFL( K ) .GT. HILFL( L ) ) THEN
                  DO 16, I =1,3
                     MERK         = HILFEV( I,K )
                     HILFEV( I,K )= HILFEV( I,L )
                     HILFEV( I,L )= MERK
   16             CONTINUE
                  MERK       = HILFL( K )
                  HILFL( K ) = HILFL( L )
                  HILFL( L ) = MERK
               ENDIF
   17    CONTINUE
C
      ENDIF
C
C
C     NUN DIE NEGATIVEN EIGENWERTE
      DO 20, I=1,3
C
C
         IF ( ( ABS( LAMBDA( I ) ) .GT. EPS ) .AND.
     >        (      LAMBDA( I )   .LT. EPS ) ) THEN
C
            HILFL( J+1 ) = LAMBDA( I )
            DO 25, K = 1,3
               HILFEV( K,J+1 ) = EV( K,I )
   25       CONTINUE
            J = J + 1
         ENDIF
   20 CONTINUE
C
C     Die Eigenwerte, die Null sind ans Ende
      DO 30, I = 1,3
         IF (ABS(LAMBDA(I)).LT.EPS) THEN
            HILFL(J+1) = LAMBDA(I)
            DO 35, K = 1,3
               HILFEV( K,J+1 ) = EV(K,I)
   35       CONTINUE
            J = J + 1
         ENDIF
   30 CONTINUE
C
C     HILF AUF LAMBDA UEBERSCHREIBEN
      DO 40, I=1,3
         LAMBDA( I ) = HILFL( I )
         DO 45, K = 1,3
            EV( K,I ) = HILFEV( K,I )
   45    CONTINUE
   40 CONTINUE
C
      INULL = 3 - J
      END
