subroutine read_polynoms_coefficients(element_index)
  use all_variables, only : global_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: element_index
  character(50) :: filename
  integer*4 :: Z,degree
  integer*4 :: k
  Type(TAm_polynom),allocatable :: polynom_buffer(:)
  Z=global_parameters%element_list(element_index)%Z
  degree=global_parameters%element_list(element_index)%amdata_polynom_degree
  allocate(polynom_buffer(1:Z))
  do k=1,Z
     allocate(polynom_buffer(k)%coefficients(degree*(degree+1)/2))
  end do
  !ionization
  filename='scd_'//trim(global_parameters%element_list(element_index)%Symbol)//'.dat'
  call read_amdata_file(filename,polynom_buffer,Z,degree)
  do k=1,Z
     global_parameters%element_list(element_index)%amdatas(k)%ionization_rate_polynom=polynom_buffer(k)
  end do
  !recombination
  filename='acd_'//trim(global_parameters%element_list(element_index)%Symbol)//'.dat'
  call read_amdata_file(filename,polynom_buffer,Z,degree)
  do k=1,Z
     global_parameters%element_list(element_index)%amdatas(k)%recombination_rate_polynom=polynom_buffer(k)
  end do
  !line excitation
  filename='plt_'//trim(global_parameters%element_list(element_index)%Symbol)//'.dat'
  call read_amdata_file(filename,polynom_buffer,Z,degree)
  do k=1,Z
     global_parameters%element_list(element_index)%amdatas(k)%line_excitation_polynom=polynom_buffer(k)
  end do
  !line recombination and Bremsstrahlung
  filename='prb_'//trim(global_parameters%element_list(element_index)%Symbol)//'.dat'
  call read_amdata_file(filename,polynom_buffer,Z,degree)
  do k=1,Z
     global_parameters%element_list(element_index)%amdatas(k)%line_recombination_polynom=polynom_buffer(k)
  end do
  deallocate(polynom_buffer)
end subroutine read_polynoms_coefficients
