subroutine check_element_list()
  use all_variables, only : global_parameters
  use Mlog_message
  use Mphysics
  implicit none
  integer*4 i
  global_parameters%element_list(0)%name='electron'
  global_parameters%element_list(0)%mass=0.
  global_parameters%element_list(0)%mass2=m_e/m_u
  global_parameters%element_list(0)%Z=-1
  do i=1,global_parameters%N_species
     select case(global_parameters%element_list(i)%symbol)
     case('H ')
        global_parameters%element_list(i)%name='Hydrogen'
        global_parameters%element_list(i)%mass=1.008
        global_parameters%element_list(i)%Z=1
     case('D ')
        global_parameters%element_list(i)%name='Deuterium'
        global_parameters%element_list(i)%mass=2.014
        global_parameters%element_list(i)%Z=1
     case('T ')
        global_parameters%element_list(i)%name='Tritium'
        global_parameters%element_list(i)%mass=3.016
        global_parameters%element_list(i)%Z=1
     case('He')
        global_parameters%element_list(i)%name='Helium'
        global_parameters%element_list(i)%mass=4.003
        global_parameters%element_list(i)%Z=2
     case('Li')
        global_parameters%element_list(i)%name='Lithium'
        global_parameters%element_list(i)%mass=6.941
        global_parameters%element_list(i)%Z=3
     case('Be')
        global_parameters%element_list(i)%name='Beryllium'
        global_parameters%element_list(i)%mass=9.012
        global_parameters%element_list(i)%Z=4
     case('B ')
        global_parameters%element_list(i)%name='Boron'
        global_parameters%element_list(i)%mass=10.811
        global_parameters%element_list(i)%Z=5
     case('C ')
        global_parameters%element_list(i)%name='Carbon'
        global_parameters%element_list(i)%mass=12.011
        global_parameters%element_list(i)%Z=6
     case('N ')
        global_parameters%element_list(i)%name='Nitrogen'
        global_parameters%element_list(i)%mass=14.001
        global_parameters%element_list(i)%Z=7
     case('O ')
        global_parameters%element_list(i)%name='Oxygen'
        global_parameters%element_list(i)%mass=15.999
        global_parameters%element_list(i)%Z=8
     case('Ne')
        global_parameters%element_list(i)%name='Neon'
        global_parameters%element_list(i)%mass=20.180
        global_parameters%element_list(i)%Z=10
     case default
        call write_log_message(0,'Unrecognised species : '//&
             global_parameters%element_list(i)%symbol)
        stop
     end select
     global_parameters%element_list(i)%mass2=global_parameters%element_list(i)%mass
  end do
end subroutine check_element_list
