subroutine compute_coulomb_logarithm(zone,STEP)
  use all_variables, only : global_parameters, transport_parameters
  use Mzone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n
  do n=0,global_parameters%N_ions
     zone%species(n)%var(STEP)%log_Lambda=transport_parameters%Coulomb_log
  end do
end subroutine compute_coulomb_logarithm
