C EIRENE06 COMPILATION
C ===== SOURCE: broadcast.f
      SUBROUTINE BROADCAST
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: broad_iestr.f

!      SUBROUTINE BROAD_USR
!      IMPLICIT NONE
!      RETURN
!      END

      SUBROUTINE BROAD_IESTR(IESTR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IESTR
      RETURN
      END
C ===== SOURCE: broadref.f

      SUBROUTINE BROADREF
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: broadsput.f

      SUBROUTINE BROADSPUT(ES,M2M1,ETF,ETH,Q,N1,N2)
      USE PRECISION
      USE PARMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN):: N1, N2
      REAL(DP), INTENT(IN OUT) :: ES(N1), M2M1(N1,N2), ETF(N1,N2),
     .                            ETH(N1,N2), Q(N1,N2)
      RETURN
      END
C ===== SOURCE: calstr.f

      SUBROUTINE CALSTR
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: colsum.f

      SUBROUTINE COLSUM
      IMPLICIT NONE
      RETURN
      END
C ===== SOURCE: mpi_barrier.f

      SUBROUTINE MPI_BARRIER (MPI_COM,IER)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MPI_COM
      INTEGER, INTENT(OUT) :: IER
      IER=0
      RETURN
      END
C ===== SOURCE: mpi_comm_rank.f

      SUBROUTINE MPI_COMM_RANK (MPI_COM,MY_PE,IER)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MPI_COM
      INTEGER, INTENT(OUT) :: MY_PE, IER
      MY_PE=0
      IER=0
      RETURN
      END
C ===== SOURCE: mpi_comm_size.f

      SUBROUTINE MPI_COMM_SIZE (MPI_COM,NPRS,IER)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MPI_COM
      INTEGER, INTENT(OUT) :: NPRS, IER
      NPRS=1
      IER=0
      RETURN
      END
C ===== SOURCE: mpi_finalize.f

      SUBROUTINE MPI_FINALIZE (IER)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: IER
      IER=0
      RETURN
      END
C ===== SOURCE: mpi_init.f

      SUBROUTINE MPI_INIT (IER)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: IER
      IER=0
      RETURN
      END
