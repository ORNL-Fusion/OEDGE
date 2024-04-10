C
C
      SUBROUTINE EIRENE_EXIT_OWN (ICC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ICC
c      CALL GREND
c      STOP

c yannick, try to implement clean crash, mpi-wise ...

      include 'mpif.h'

c this may not be the master process, or another one might have crashed before ... 
      stop



      END SUBROUTINE EIRENE_EXIT_OWN
