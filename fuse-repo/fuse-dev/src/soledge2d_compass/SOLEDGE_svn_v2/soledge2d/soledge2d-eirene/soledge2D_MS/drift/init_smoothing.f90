subroutine init_smoothing()
  use all_variables, only : global_parameters, zones, global_variables
  use Msmoothing_vars
  implicit none
  CSC_smoothing%pastix_comm=0 ! no mpi
  call find_delta_r()
  call init_implicit_pastix(CSC_smoothing)
end subroutine init_smoothing
