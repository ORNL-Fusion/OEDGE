subroutine update_radial_diffusivity_feedback(zone)
  use all_variables, only : reference_parameters, global_parameters
  use MradialFeedback
  use Mzone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: i,j,k,n
  integer*4 :: Nx,Nz
  real*8 :: x,D,keep,D0
  real*8 :: interpolate_feedback
  k=zone%number
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  keep=radialFeedbackData%keep
  D0=(reference_parameters%geometry%rs0**2.)/reference_parameters%fields%tau0
  do i=1,Nx
     x=zone%mesh%x(i,1)
     D=interpolate_feedback(x,1) ! D (and nu)
     do n=1,global_parameters%N_ions
        zone%species(n)%transport_perp%D_p(i,:)=D/D0!(keep*zone%species(n)%transport_perp%D_p(i,:)+D)/(keep+1.D0)
        zone%species(n)%transport_perp%nu_p(i,:)=D/D0!zone%species(n)%transport_perp%D_p(i,:)
     end do
     D=interpolate_feedback(x,3) ! chii
     do n=1,global_parameters%N_ions
        zone%species(n)%transport_perp%chi_p(i,:)=D/D0!(keep*zone%species(n)%transport_perp%chi_p(i,:)+D)/(keep+1.D0)
     end do
     D=interpolate_feedback(x,2) ! chie
     zone%species(0)%transport_perp%chi_p(i,:)=D/D0!(keep*zone%species(0)%transport_perp%chi_p(i,:)+D)/(keep+1.D0)
  end do
end subroutine update_radial_diffusivity_feedback



