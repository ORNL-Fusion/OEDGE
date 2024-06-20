module Mreaders

  
  contains

  subroutine skip_line(file,line_nbr)
    implicit none
    integer,intent(in) :: file
    integer*4,intent(in) :: line_nbr
    integer*4 :: i
    do i=1,line_nbr
       read(file,*)
    end do
  end subroutine skip_line

  subroutine  parse_line_elements(String,N_elements)
    use all_variables, only : global_parameters
    use Mphysics
    use Mlog_message
    implicit none
    character(128),intent(in) :: String
    integer*4,intent(in) :: N_elements
    integer :: pos
    character(128) :: buffer
    integer*4 i
    character(2) :: dummy_symbol
    buffer=String
    do i=1,N_elements-1
       pos=index(buffer,",")
       if(pos.ge.2) then
          dummy_symbol=adjustl(trim(buffer(1:pos-1)))
          global_parameters%element_list(i)%symbol=dummy_symbol
          buffer=adjustl(buffer(pos+1:))
       else
          call write_log_message(0,'Unrecognised species list')
          stop
       end if
    end do
    !last species
    dummy_symbol=adjustl(trim(buffer))
    global_parameters%element_list(N_elements)%symbol=dummy_symbol
  end subroutine parse_line_elements

  subroutine  parse_line_elements_puff(String,N_puffs,nstrat0)
    use all_variables, only : global_parameters
    use Mphysics
    use Mlog_message
    use styx2eirene, only : puffs
    implicit none
    character(128),intent(in) :: String
    integer*4,intent(in) :: N_puffs,nstrat0
    character(128) :: buffer
    integer :: pos
    integer*4 i
    character(2) :: dummy_symbol
    buffer=String
    do i=1,N_puffs-1
       pos=index(buffer,",")
       if(pos.ge.2) then
          dummy_symbol=adjustl(trim(buffer(1:pos-1)))
          puffs(nstrat0+i)%species=dummy_symbol
          buffer=adjustl(buffer(pos+1:))
       else
          call write_log_message(0,'Unrecognised species list for gas puffs')
          stop
       end if
    end do
    !last species
    dummy_symbol=adjustl(trim(buffer))
    puffs(nstrat0+N_puffs)%species=dummy_symbol
  end subroutine parse_line_elements_puff

  subroutine parse_line_integer(String,N_integer,int_list)
    implicit none
    character(128),intent(in) :: String
    integer*4,intent(in) :: N_integer
    integer*4,intent(out) :: int_list(1:N_integer)
    integer :: pos
    integer*4 :: i
    character(128) :: buffer
    pos=index(String,'=')
    buffer=String(pos+1:)
    read(buffer,*) int_list
  end subroutine parse_line_integer

  subroutine parse_line_float(String,N_float,float_list)
    implicit none
    character(128),intent(in) :: String
    integer*4,intent(in) :: N_float
    real*8,intent(out) :: float_list(1:N_float)
    integer :: pos
    integer*4 :: i
    character(128) :: buffer
    pos=index(String,'=')
    buffer=String(pos+1:)
    read(buffer,*) float_list
  end subroutine parse_line_float

  subroutine parse_line_logical(String,N_logical,logical_list)
    implicit none
    character(128),intent(in) :: String
    integer*4,intent(in) :: N_logical
    logical,intent(out) :: logical_list(1:N_logical)
    integer :: pos
    integer*4 :: i
    character(128) :: buffer
    pos=index(String,'=')
    buffer=String(pos+1:)
    read(buffer,*) logical_list
  end subroutine parse_line_logical


end module Mreaders
