subroutine exit_eirene_coupling()
  use eirmod_cpes
  implicit none
  integer*4 :: ier

  include 'mpif.h'

  if (my_pe==0) call grend
  call mpi_barrier(mpi_comm_world,ier)
  call mpi_finalize(ier)

  
end subroutine exit_eirene_coupling
