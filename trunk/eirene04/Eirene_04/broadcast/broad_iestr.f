

      SUBROUTINE BROAD_IESTR(IESTR)
      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INTEGER, INTENT(IN) :: IESTR
      INTEGER :: IER

      CALL MPI_BARRIER (MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (IESTR,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BARRIER (MPI_COMM_WORLD,IER)

      RETURN
      END