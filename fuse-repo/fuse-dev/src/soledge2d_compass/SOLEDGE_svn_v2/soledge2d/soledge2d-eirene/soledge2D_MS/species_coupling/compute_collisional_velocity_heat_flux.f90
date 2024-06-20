subroutine compute_collisional_velocity_heat_flux(zone)
  use all_variables, only : global_parameters, reference_parameters
  use Mphysics
  use MZone
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4 :: Nx,Nz,i,j,n,m
  real*8 :: density_mass
  real*8,parameter :: coef=0.71D0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=0,Nz+1
        do n=0,global_parameters%N_ions
           do m=0,global_parameters%N_ions
              density_mass=(zone%species(n)%var(1)%density(i,j)*zone%species(m)%var(1)%density(i,j))&
                   /(zone%species(n)%var(1)%density(i,j)/zone%species(n)%element%mass2&
                   +zone%species(m)%var(1)%density(i,j)/zone%species(m)%element%mass2)
              zone%species(n)%coupling_terms%qu(i,j,m)=coef*zone%species(n)%var(1)%temperature(i,j)&
                   *zone%species(n)%charge**2*density_mass/zone%species(n)%element%mass2&
                   *(zone%species(n)%var(1)%velocity(i,j)-zone%species(m)%var(1)%velocity(i,j))
           end do
        end do
     end do
  end do
end subroutine compute_collisional_velocity_heat_flux
