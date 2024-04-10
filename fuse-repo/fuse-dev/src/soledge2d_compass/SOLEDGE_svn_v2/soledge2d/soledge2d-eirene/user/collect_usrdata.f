

      SUBROUTINE EIRENE_COLLECT_USRDATA
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comsou
      use styx2eirene
      use eirmod_cpes
      IMPLICIT NONE
      integer :: ier
      INCLUDE 'mpif.h'

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      call mpi_gather(timpara,1,MPI_REAL8,time_para,
     .    1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      RETURN
      END
