C
      FUNCTION SMACH(I)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      REAL(DP)  SMACH
      IF(I.EQ.1) THEN
        smach=epsilon(1.d0)
      ELSEIF(I.EQ.2) THEN
        smach=tiny(1.d0)
      ELSEIF(I.EQ.3) THEN
        smach=huge(1.d0)
      ELSE
        WRITE(6,'(A)')
     @   ' SMACH CALLED WITH PARAMETER JOB NOT EQUAL TO 1,2 OR 3.'
      ENDIF
      RETURN
      END
