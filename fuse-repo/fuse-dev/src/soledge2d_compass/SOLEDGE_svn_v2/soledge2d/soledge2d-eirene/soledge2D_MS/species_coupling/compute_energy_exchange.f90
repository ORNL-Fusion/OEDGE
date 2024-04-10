subroutine compute_energy_exchange(zone)
  use all_variables, only : global_parameters, reference_parameters
  use Mphysics
  use MZone
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4 :: Nx,Nz,i,j,n,m
  real*8 :: density_mass,mass
  real*8,parameter :: coef=0.71D0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        do n=0,global_parameters%N_ions
           do m=0,global_parameters%N_ions
              density_mass=(zone%species(n)%var(1)%density(i,j)*zone%species(m)%var(1)%density(i,j))&
                   /(zone%species(n)%var(1)%density(i,j)/zone%species(n)%element%mass2&
                   +zone%species(m)%var(1)%density(i,j)/zone%species(m)%element%mass2)
              mass=(zone%species(n)%element%mass2*zone%species(m)%element%mass2)&
                   /(zone%species(n)%element%mass2+zone%species(m)%element%mass2)
              zone%species(n)%coupling_terms%Q(i,j,m)=-1.5D0*density_mass/mass/&
                   zone%species(n)%coupling_terms%tau(i,j,m)*reference_parameters%fields%tau0&
                   *(zone%species(n)%var(1)%temperature(i,j)-zone%species(m)%var(1)%temperature(i,j))
              zone%species(n)%coupling_terms%Q(i,j,m)=zone%species(n)%coupling_terms%Q(i,j,m)&
                   -zone%species(n)%coupling_terms%R(i,j,m)&
                   *(zone%species(n)%var(1)%velocity(i,j)-zone%species(m)%var(1)%velocity(i,j))&
                   *zone%species(m)%element%mass2/(zone%species(n)%element%mass2+zone%species(m)%element%mass2)
           end do
        end do
     end do
  end do  
end subroutine compute_energy_exchange
