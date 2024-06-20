integer*4 function read_polynom_degree(Symbol)
  use Mreaders
  implicit none
  character(2),intent(in) :: Symbol
  character(50) :: filename
  filename='scd_'//trim(Symbol)//'.dat'
  open(unit=20,file=trim(filename),status='unknown')
  call skip_line(20,11)
  read(20,1) read_polynom_degree
1 format(20X,I3)
  close(20)
end function read_polynom_degree
