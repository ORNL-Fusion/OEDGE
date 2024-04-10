module Mlog_message
  
  integer*4,parameter ::  msg_soledge_started = 1
  integer*4,parameter ::  msg_soledge_ended = 2
  integer*4,parameter ::  msg_element_list = 3
  integer*4,parameter ::  msg_species_list = 4
  integer*4,parameter ::  msg_mesh_loaded = 5
  integer*4,parameter ::  msg_simulation_crashed = 6
  integer*4,parameter ::  msg_loop_started = 7

contains

  subroutine write_log_message(msg_number,msg)
    implicit none
    integer*4,intent(in) :: msg_number
    character(len=*),intent(in),optional :: msg
    select case(msg_number) 
    case(0)
       write(20,*) msg
       write(*,*) msg
    case(msg_soledge_started)
       call write_start_message()
    case(msg_soledge_ended)
       call write_end_message()
    case(msg_element_list)
       call write_element_list()
    case(msg_species_list)
       call write_species_list()
    case(msg_mesh_loaded)
       call write_mesh_loaded()
    case(msg_simulation_crashed)
       call write_crash()
    case(msg_loop_started)
       call write_loop_begin_message()
    end select
  end subroutine write_log_message


end module Mlog_message
