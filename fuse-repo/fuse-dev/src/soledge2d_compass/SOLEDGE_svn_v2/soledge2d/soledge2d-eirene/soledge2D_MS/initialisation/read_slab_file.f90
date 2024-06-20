subroutine read_slab_file()
  use all_variables, only : reference_parameters
  use Mreaders
  implicit none
  integer :: finput
  finput=111
  open(unit=finput,file='SLAB.txt',status='unknown')
  call skip_line(finput,4)
  read(finput,1)  reference_parameters%geometry%qref
  call skip_line(finput,2)
  read(finput,2)  reference_parameters%geometry%R0
  call skip_line(finput,2)
  read(finput,2)  reference_parameters%geometry%Rm0
  call skip_line(finput,2)
  read(finput,3)  reference_parameters%geometry%rs0
  call skip_line(finput,2)
  read(finput,1)  reference_parameters%geometry%Bpol0
  call skip_line(finput,2)
  read(finput,1)  reference_parameters%geometry%Btor0
1 format(7X,F5.2)
2 format(4X,F5.2)
3 format(10X,F5.2)
  close(finput)
end subroutine read_slab_file
