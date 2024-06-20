subroutine compute_collision_times(zone)
  use all_variables, only : global_parameters, reference_parameters, global_variables
  use Mphysics
  use MZone
  implicit none
  type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz,i,j,n,m
  real*8 :: Adim,density
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  Adim=3.D0*sqrt(2.D0)*pi**1.5D0*epsilon_0**2/(eV**4)
  Adim=Adim*(eV*reference_parameters%fields%T0eV)**1.5D0/reference_parameters%fields%n0&
       *sqrt(m_u)
  Adim=Adim/12.D0 ! assume log Lambda is 12
  do i=1,Nx
     do j=1,Nz
        do n=0,global_parameters%N_ions
           do m=0,global_parameters%N_ions
!!$              density=2.D0*(zone%species(n)%var(1)%density(i,j)*zone%species(m)%var(1)%density(i,j))&
!!$                   /(zone%species(n)%var(1)%density(i,j)+zone%species(m)%var(1)%density(i,j))
              density=(zone%species(n)%var(1)%density(i,j)+zone%species(m)%var(1)%density(i,j))/2.D0
              zone%species(n)%coupling_terms%tau(i,j,m)=Adim&
                   *zone%species(n)%element%mass2*zone%species(m)%element%mass2&
                   /(zone%species(n)%charge**2*zone%species(m)%charge**2)&
                   /density&
                   *(zone%species(n)%var(1)%temperature(i,j)/zone%species(n)%element%mass2&
                   +zone%species(m)%var(1)%temperature(i,j)/zone%species(m)%element%mass2)**1.5D0
              zone%species(n)%coupling_terms%tau(i,j,m)=max(zone%species(n)%coupling_terms%tau(i,j,m),&
                   global_variables%dt*reference_parameters%fields%tau0)
           end do
        end do
     end do
  end do
end subroutine compute_collision_times
