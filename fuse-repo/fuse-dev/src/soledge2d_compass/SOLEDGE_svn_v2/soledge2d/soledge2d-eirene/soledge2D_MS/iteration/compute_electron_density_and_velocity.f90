subroutine compute_electron_density_and_velocity(zone,step)
  use all_variables, only : global_parameters
  use Mzone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: step
  integer*4 :: n
  real*8,parameter :: eps = 1.d-10
  zone%species(0)%var(step)%density=0.
  zone%species(0)%var(step)%velocity=0.
  zone%species(0)%var(step)%Gamma=0.
  do n=1,global_parameters%N_ions
     zone%species(0)%var(step)%density=zone%species(0)%var(step)%density&
          +zone%species(n)%var(step)%density*zone%species(n)%charge
     zone%species(0)%var(step)%Gamma=zone%species(0)%var(step)%Gamma&
          +zone%species(n)%var(step)%Gamma*zone%species(n)%charge
  end do
  zone%species(0)%var(step)%velocity=zone%species(0)%var(step)%Gamma&
       /(zone%species(0)%var(step)%density+eps)
end subroutine compute_electron_density_and_velocity

