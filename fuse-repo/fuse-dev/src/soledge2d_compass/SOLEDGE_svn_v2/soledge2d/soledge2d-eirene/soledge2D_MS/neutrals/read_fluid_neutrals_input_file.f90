subroutine read_fluid_neutrals_input_file()
  use Mneutral_vars
  use Mreaders
  integer :: finput
  finput=106
  open(unit=finput,file='fluid_neutrals',status='unknown')
  call skip_line(finput,1)
  read(finput,100) FN_diffusivity
  call skip_line(finput,2)
  read(finput,101) FN_recycling_coefficient
  close(finput)
100 format(5X,ES9.2E2)
101 format(5X,F5.3)
end subroutine read_fluid_neutrals_input_file
