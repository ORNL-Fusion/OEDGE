C
C
      FUNCTION EIRENE_IDEZ(I,J,IZIF)
C
C   I IS AN INTEGER WITH  N DEZIMALS, THE FUNCTION RETURNS
C   IDEZ, THE J.TH DEZIMAL
C   J MUST BE .LE. IZIF
C      IF N .LT. IZIF AND J .GT. N  , IDEZ=0
C      IF N .GT. IZIF AND J .EQ. IZIF  , IDEZ=PART OF I ON THE LEFT OF
C      THE DEZIMAL J, INCLUDING THE J TH DEZIMAL
C      E.G.: I=1234,IZIF=3, J=1: IDEZ=4
C                           J=2: IDEZ=3
C                           J=3: IDEZ=12
C                           J=4: ERROR ,J GT IZIF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, J, IZIF
      INTEGER :: EIRENE_IDEZ, K, IZ, IQ, IT, IDIF
 
      IDIF=IZIF-J+1
      IF (IDIF.LE.0) THEN
        CALL EIRENE_MASAGE
     .  ('ERROR IN FUNCTION EIRENE_IDEZ                         ')
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      IZ=I
      DO 1 K=1,IDIF
        IT=10**(IZIF-K)
        IQ=IZ/IT
        IZ=IZ-IQ*IT
1     CONTINUE
      EIRENE_IDEZ=IQ
      RETURN
      END
