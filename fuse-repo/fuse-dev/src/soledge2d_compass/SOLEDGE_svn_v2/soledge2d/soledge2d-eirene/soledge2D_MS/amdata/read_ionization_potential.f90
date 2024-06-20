subroutine read_ionization_potential(element_index)
  use all_variables, only : global_parameters
  use Mreaders
  implicit none
  integer*4,intent(in) :: element_index
  character(50) :: filename
  integer*4 :: k
  filename='ADAS_ionization_potentials_'//trim(global_parameters%element_list(element_index)%Symbol)
  open(unit=20,file=trim(filename),status='unknown')
  call skip_line(20,7)
  do k=1,global_parameters%element_list(element_index)%Z
     read(20,'(E14.7)') global_parameters%element_list(element_index)%amdatas(k)%ionization_potential
  end do
end subroutine read_ionization_potential
