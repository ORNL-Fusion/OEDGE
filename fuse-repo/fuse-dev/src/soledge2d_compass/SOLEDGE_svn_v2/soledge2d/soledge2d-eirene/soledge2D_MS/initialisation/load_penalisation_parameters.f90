subroutine load_penalisation_parameters()
  use all_variables, only : penalisation_parameters, global_variables, flags
  implicit none
  open(unit=100,file='penalisation_parameters',status='unknown')
  read(100,1) penalisation_parameters%eta
  read(100,3) penalisation_parameters%dump_pen
  read(100,4) penalisation_parameters%keep
  read(100,3) penalisation_parameters%gain
  read(100,2) global_variables%min_density
  read(100,2) global_variables%min_temperature
  read(100,5) flags%use_weno
  read(100,6) global_variables%epsilon_weno
  penalisation_parameters%eta2=penalisation_parameters%eta*1.d-2
1 format(6X,F7.4)
2 format(8X,F7.4)
3 format(11X,F7.4)
4 format(7X,I6)
5 format(11X,L1)
6 format(15X,F7.4)
end subroutine load_penalisation_parameters

