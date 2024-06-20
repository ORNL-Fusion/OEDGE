subroutine compute_implicit_fluxes(zone)
  use all_variables, only : global_parameters, global_variables, reference_parameters
  use MZone
  use Mphysics
  implicit none
  Type(TZone),intent(inout) :: zone
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: kappa(:,:)
  integer*4 :: i,j,Nx,Nz,n
  real*8 :: rs0,dt,R0,Teps
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(kappa(0:Nx+1,0:Nz+1))
  Teps=global_variables%Teps
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  do n=1,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           !density
           Fluxes(i,j,3)=zone%species(n)%implicit_coefs%east_density(i,j)&
                *(zone%species(n)%var(2)%density(i,j+1)-zone%species(n)%var(2)%density(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
           zone%species(n)%fluxes%fluxn(i,j,3)=zone%species(n)%fluxes%fluxn(i,j,3)-Fluxes(i,j,3)
           Fluxes(i,j,4)=zone%species(n)%implicit_coefs%west_density(i,j)&
                *(zone%species(n)%var(2)%density(i,j)-zone%species(n)%var(2)%density(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
           zone%species(n)%fluxes%fluxn(i,j,4)=zone%species(n)%fluxes%fluxn(i,j,4)-Fluxes(i,j,4)
           !velocity 
           Fluxes(i,j,3)=zone%species(n)%implicit_coefs%east_velocity(i,j)&
                *(zone%species(n)%var(1)%density(i,j+1)+zone%species(n)%var(1)%density(i,j))*0.5d0&
                *(zone%species(n)%var(2)%velocity(i,j+1)-zone%species(n)%var(2)%velocity(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
           zone%species(n)%fluxes%fluxG(i,j,3)=zone%species(n)%fluxes%fluxG(i,j,3)-Fluxes(i,j,3)
           Fluxes(i,j,4)=zone%species(n)%implicit_coefs%west_velocity(i,j)&
                *(zone%species(n)%var(1)%density(i,j)+zone%species(n)%var(1)%density(i,j-1))*0.5d0&
                *(zone%species(n)%var(2)%velocity(i,j)-zone%species(n)%var(2)%velocity(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
           zone%species(n)%fluxes%fluxG(i,j,4)=zone%species(n)%fluxes%fluxG(i,j,4)-Fluxes(i,j,4)
        end do
     end do
  end do
  do n=0,global_parameters%N_ions
     kappa=zone%species(n)%transport_para%kappa/&
          zone%species(n)%var(1)%log_Lambda
     do i=1,Nx
        do j=1,Nz
           !temperature perp diffusion 
           Fluxes(i,j,3)=zone%species(n)%implicit_coefs%east_temperature(i,j)&
                *(zone%species(n)%var(1)%density(i,j+1)+zone%species(n)%var(1)%density(i,j))*0.5d0&
                *(zone%species(n)%var(2)%temperature(i,j+1)-zone%species(n)%var(2)%temperature(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
           zone%species(n)%fluxes%fluxE(i,j,3)=zone%species(n)%fluxes%fluxE(i,j,3)-Fluxes(i,j,3)
           Fluxes(i,j,4)=zone%species(n)%implicit_coefs%west_temperature(i,j)&
                *(zone%species(n)%var(1)%density(i,j)+zone%species(n)%var(1)%density(i,j-1))*0.5d0&
                *(zone%species(n)%var(2)%temperature(i,j)-zone%species(n)%var(2)%temperature(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
           zone%species(n)%fluxes%fluxE(i,j,4)=zone%species(n)%fluxes%fluxE(i,j,4)-Fluxes(i,j,4)
           !parallel diffusion
           Fluxes(i,j,3)=(kappa(i,j)*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0&
                +kappa(i,j+1)*(max(zone%species(n)%var(1)%temperature(i,j+1),Teps))**2.5d0)*0.5d0&
                *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5d0&
                *(zone%species(n)%var(2)%temperature(i,j+1)-zone%species(n)%var(2)%temperature(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
                *zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)
           zone%species(n)%fluxes%fluxE(i,j,3)=zone%species(n)%fluxes%fluxE(i,j,3)-Fluxes(i,j,3)
           Fluxes(i,j,4)=(kappa(i,j)*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0&
                +kappa(i,j-1)*(max(zone%species(n)%var(1)%temperature(i,j-1),Teps))**2.5d0)*0.5d0&
                *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5d0&
                *(zone%species(n)%var(2)%temperature(i,j)-zone%species(n)%var(2)%temperature(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
                *zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)
           zone%species(n)%fluxes%fluxE(i,j,4)=zone%species(n)%fluxes%fluxE(i,j,4)-Fluxes(i,j,4)
        end do
     end do
  end do
end subroutine compute_implicit_fluxes
