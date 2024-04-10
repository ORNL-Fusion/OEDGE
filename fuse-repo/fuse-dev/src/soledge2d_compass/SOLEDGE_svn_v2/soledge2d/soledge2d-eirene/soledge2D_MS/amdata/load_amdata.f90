subroutine load_amdata()
  use all_variables
  use Mlog_message
  implicit none
  integer*4 :: i
  integer*4 :: Z
  do i=1,global_parameters%N_species
     select case(global_parameters%element_list(i)%symbol)
     case('H ','D ','T ')
	Z=global_parameters%element_list(i)%Z
        allocate(global_parameters%element_list(i)%amdatas(1:Z))     
	call read_ionization_potential(i)
     case('He','Li','Be','B ','C ','N ','O ','Ne')
        Z=global_parameters%element_list(i)%Z
        allocate(global_parameters%element_list(i)%amdatas(1:Z))     
        call allocate_polynoms(i)
        call read_polynoms_coefficients(i)
        call read_ionization_potential(i)
     case default
        call write_log_message(0,'Unrecognised species : '//&
             global_parameters%element_list(i)%symbol)
        stop
     end select
  end do
end subroutine load_amdata
