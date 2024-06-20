subroutine allocate_polynoms(i)
  use all_variables, only : global_parameters
  implicit none
  integer*4,intent(in) :: i
  integer*4 :: Z
  integer*4 :: k
  integer*4 :: degree
  integer*4 :: read_polynom_degree
  Z=global_parameters%element_list(i)%Z
  degree=read_polynom_degree(global_parameters%element_list(i)%Symbol)
  global_parameters%element_list(i)%amdata_polynom_degree=degree
  do k=1,Z
     allocate(global_parameters%element_list(i)%amdatas(k)%recombination_rate_polynom%coefficients(1:degree*(degree+1)/2))
     allocate(global_parameters%element_list(i)%amdatas(k)%ionization_rate_polynom%coefficients(1:degree*(degree+1)/2))
     allocate(global_parameters%element_list(i)%amdatas(k)%line_excitation_polynom%coefficients(1:degree*(degree+1)/2))
     allocate(global_parameters%element_list(i)%amdatas(k)%line_recombination_polynom%coefficients(1:degree*(degree+1)/2))
  end do
end subroutine allocate_polynoms
