subroutine compute_temperature_source_implicit_correction(zone)
  use all_variables, only : global_parameters
  use MZone
  use Moperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: i,j,n
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  !ions
  do n=1,global_parameters%N_ions
     Fluxes=0.D0
     !east
     do i=1,Nx
        do j=1,Nz
           Fluxes(i,j,3)=zone%species(n)%implicit_coefs%east_density(i,j)&
                *(zone%species(n)%var(1)%temperature(i,j)+zone%species(n)%var(1)%temperature(i,j+1))*0.5D0&
                *(zone%species(n)%var(2)%density(i,j+1)-zone%species(n)%var(2)%density(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))*2.5d0
           Fluxes(i,j,4)=zone%species(n)%implicit_coefs%west_density(i,j)&
                *(zone%species(n)%var(1)%temperature(i,j)+zone%species(n)%var(1)%temperature(i,j-1))*0.5D0&
                *(zone%species(n)%var(2)%density(i,j)-zone%species(n)%var(2)%density(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))*2.5d0

           Fluxes(i,j,3)=Fluxes(i,j,3)+zone%species(n)%implicit_coefs%east_density(i,j)&
                *(zone%species(n)%var(1)%velocity(i,j)**2.d0+zone%species(n)%var(1)%velocity(i,j+1)**2.d0)*0.5D0&
                *(zone%species(n)%var(2)%density(i,j+1)-zone%species(n)%var(2)%density(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))*0.5d0*zone%species(n)%element%mass
           Fluxes(i,j,4)=Fluxes(i,j,4)+zone%species(n)%implicit_coefs%west_density(i,j)&
                *(zone%species(n)%var(1)%velocity(i,j)**2.d0+zone%species(n)%var(1)%velocity(i,j-1)**2.d0)*0.5D0&
                *(zone%species(n)%var(2)%density(i,j)-zone%species(n)%var(2)%density(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))*0.5d0*zone%species(n)%element%mass

           Fluxes(i,j,3)=Fluxes(i,j,3)+zone%species(n)%implicit_coefs%east_velocity(i,j)&
                *(zone%species(n)%var(1)%Gamma(i,j)+zone%species(n)%var(1)%Gamma(i,j+1))*0.5D0&
                *(zone%species(n)%var(2)%velocity(i,j+1)-zone%species(n)%var(2)%velocity(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))*zone%species(n)%element%mass
           zone%species(n)%fluxes%fluxE(i,j,3)=zone%species(n)%fluxes%fluxE(i,j,3)-Fluxes(i,j,3)
           Fluxes(i,j,4)=Fluxes(i,j,4)+zone%species(n)%implicit_coefs%west_velocity(i,j)&
                *(zone%species(n)%var(1)%Gamma(i,j)+zone%species(n)%var(1)%Gamma(i,j-1))*0.5D0&
                *(zone%species(n)%var(2)%velocity(i,j)-zone%species(n)%var(2)%velocity(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))*zone%species(n)%element%mass
           zone%species(n)%fluxes%fluxE(i,j,4)=zone%species(n)%fluxes%fluxE(i,j,4)-Fluxes(i,j,4)

        end do
     end do
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%sources%SE=zone%species(n)%sources%SE+Source
  end do
  !electrons
  Fluxes=0.D0
  do n=1,global_parameters%N_ions
     !east
     do i=1,Nx
        do j=1,Nz
           Fluxes(i,j,3)=Fluxes(i,j,3)+zone%species(n)%implicit_coefs%east_density(i,j)&
                *(zone%species(0)%var(1)%temperature(i,j)+zone%species(0)%var(1)%temperature(i,j+1))*0.5D0&
                *(zone%species(n)%var(2)%density(i,j+1)-zone%species(n)%var(2)%density(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))*2.5d0*zone%species(n)%charge
           zone%species(0)%fluxes%fluxE(i,j,3)=zone%species(0)%fluxes%fluxE(i,j,3)-Fluxes(i,j,3)
           Fluxes(i,j,4)=Fluxes(i,j,4)+zone%species(n)%implicit_coefs%west_density(i,j)&
                *(zone%species(0)%var(1)%temperature(i,j)+zone%species(0)%var(1)%temperature(i,j-1))*0.5D0&
                *(zone%species(n)%var(2)%density(i,j)-zone%species(n)%var(2)%density(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))*2.5d0*zone%species(n)%charge
           zone%species(0)%fluxes%fluxE(i,j,4)=zone%species(0)%fluxes%fluxE(i,j,4)-Fluxes(i,j,4)
        end do
     end do
  end do
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%species(0)%sources%SE=zone%species(0)%sources%SE+Source
  deallocate(Fluxes,Source)
end subroutine compute_temperature_source_implicit_correction
