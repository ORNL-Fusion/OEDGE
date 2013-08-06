module read_data
  use error_handling

  implicit none

  integer,parameter :: line_length=512
  character,parameter :: line_form='(a512)'

  interface read_data

     module procedure read_str_data, read_int_data, read_r_data, read_r8_data


  end interface



contains

subroutine read_str_data(unit, str)
  implicit none
  
  integer :: unit
  character*(*) str
  character*line_length :: line

  call get_next_line(unit,line)





end subroutine read_str_data


subroutine get_next_line(unit,line,ios)
implicit none

integer :: unit,ios
character*(*) :: line

10 read(unit,line_form,iostat=ios) line

if (ios.eq.0) then 

    if (line(1:1) .eq. '#'.or. line(1:1).eq. '$') goto 10























end module read_data
