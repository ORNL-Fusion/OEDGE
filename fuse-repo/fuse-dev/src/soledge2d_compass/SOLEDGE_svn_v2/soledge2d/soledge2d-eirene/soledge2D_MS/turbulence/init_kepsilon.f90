subroutine init_kepsilon()
  use all_variables, only : global_parameters, zones, kepsilon_param
  use Mturbulence
  implicit none
  call allocate_kepsilon()
  call read_kepsilon_param()
  call init_kepsilon_fields()
end subroutine init_kepsilon
