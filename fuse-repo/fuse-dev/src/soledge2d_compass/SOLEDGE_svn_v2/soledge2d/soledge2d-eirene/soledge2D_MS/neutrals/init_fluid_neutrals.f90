subroutine init_fluid_neutrals()
  use all_variables, only : global_parameters, zones, global_variables
  use Mneutral_vars
  use Mdefinitions
  implicit none
  call allocate_neutrals()
  call init_fluid_neutrals_fields()
  CSC_neutral%pastix_comm=0 ! no mpi
  call init_implicit_pastix(CSC_neutral)
end subroutine init_fluid_neutrals
