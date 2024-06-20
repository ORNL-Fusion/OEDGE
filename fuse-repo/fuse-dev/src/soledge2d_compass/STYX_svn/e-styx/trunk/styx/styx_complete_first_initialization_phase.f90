subroutine styx_complete_first_initialization_phase
  use all_variables, only : global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  implicit none

  ! electron temperature
  allocate(eiv_e%T(Neir_cells))
  eiv_e%T=0._dp
  allocate(vpar_tri(Neir_cells,global_parameters%N_ions))
  allocate(pflux_in(Neir_cells,3,Ntor_cells))
  
  ! to be handle properly in the multifluid case ...
  !nplsi=1
  !!!!!!!!!!! now use the actual value from the input file !!!!

! initialisation of plot library

  call GRSTRT(35,8)

end subroutine styx_complete_first_initialization_phase
