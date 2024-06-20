subroutine init_vorticity()
  use all_variables, only : global_parameters, zones, global_variables
  use Mvorticity_vars
  use Mdefinitions
  implicit none
  call allocate_vorticity()
  call init_vorticity_fields()
  CSC_vort%pastix_comm=0 ! no mpi
!  call init_implicit_pastix(CSC_vort)
  global_variables%dt_vort=0.d0 ! reset (incremented in step_in_soledge)
end subroutine init_vorticity
