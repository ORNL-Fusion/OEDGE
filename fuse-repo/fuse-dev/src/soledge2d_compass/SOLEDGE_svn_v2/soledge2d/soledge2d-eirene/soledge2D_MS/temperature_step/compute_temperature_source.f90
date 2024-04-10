subroutine compute_temperature_source(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4 :: i,j,n
  integer*4 :: Nx,Nz
  real*8 :: dt
  real*8 :: m_i
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  dt=global_variables%dt
  do i=1,Nx
     do j=1,Nz
        do n=0,global_parameters%N_ions
           m_i=zone%species(n)%element%mass
           zone%species(n)%tridiag%S(i,j)=1.5d0*zone%species(n)%var(1)%density(i,j)*zone%species(n)%var(1)%temperature(i,j)/dt&
                -0.5d0*m_i*(zone%species(n)%var(2)%Gamma(i,j)**2.D0/zone%species(n)%var(2)%density(i,j)&
                -zone%species(n)%var(1)%Gamma(i,j)**2.D0/zone%species(n)%var(1)%density(i,j))/dt
        end do
     end do
  end do
end subroutine compute_temperature_source
