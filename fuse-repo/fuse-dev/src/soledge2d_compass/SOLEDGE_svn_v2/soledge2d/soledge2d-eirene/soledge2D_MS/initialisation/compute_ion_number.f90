subroutine compute_ion_number()
  use all_variables, only : global_parameters
  use Mlog_message
  integer*4 n
  global_parameters%N_ions=0
  do n=1,global_parameters%N_species
     global_parameters%N_ions=global_parameters%N_ions&
          +global_parameters%element_list(n)%Z
  end do
end subroutine compute_ion_number
