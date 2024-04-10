C
      FUNCTION EIRENE_SMACH(I)
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      REAL(DP)  EIRENE_SMACH
      IF(I.EQ.1) THEN
        EIRENE_smach=epsilon(1.d0)
      ELSEIF(I.EQ.2) THEN
        EIRENE_smach=tiny(1.d0)
      ELSEIF(I.EQ.3) THEN
        EIRENE_smach=huge(1.d0)
      ELSE
        WRITE(iunout,'(A)')
     .  ' SMACH CALLED EIRENE_WITH PARAMETER JOB NOT EQUAL TO 1,2 OR 3.'
      ENDIF
      RETURN
      END
