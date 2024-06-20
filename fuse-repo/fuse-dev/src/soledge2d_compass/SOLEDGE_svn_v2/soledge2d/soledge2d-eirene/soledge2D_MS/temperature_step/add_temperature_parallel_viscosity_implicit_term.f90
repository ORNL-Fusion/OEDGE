subroutine add_temperature_parallel_viscosity_implicit_term(zone)
  use all_variables, only : global_parameters, global_variables
  use MZone
  use Moperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: i,j,n
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  real*8,allocatable :: nu(:,:)
  real*8 :: Teps
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  allocate(nu(0:Nx+1,0:Nz+1))
  Teps=global_variables%Teps
  !ions
  do n=1,global_parameters%N_ions
     nu=zone%species(n)%transport_para%nu/&
          zone%species(n)%var(1)%log_Lambda
     Fluxes=0.D0
     !east
     do i=1,Nx
        do j=1,Nz
           Fluxes(i,j,3)=Fluxes(i,j,3)+0.5D0*&
                (nu(i,j+1)*(max(zone%species(n)%var(1)%temperature(i,j+1),Teps))**2.5d0&
                +nu(i,j)*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0)*0.5d0&
                *(zone%metric_coefficients%G(i,j+1)+zone%metric_coefficients%G(i,j))*0.5D0&
                *(zone%species(n)%var(2)%velocity(i,j+1)**2-zone%species(n)%var(2)%velocity(i,j)**2)&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))*(1.d0-zone%masks%chi4(i,j+1))*(1.d0-zone%masks%chi4(i,j))
           zone%species(n)%fluxes%fluxE(i,j,3)=zone%species(n)%fluxes%fluxE(i,j,3)-Fluxes(i,j,3)
           Fluxes(i,j,4)=Fluxes(i,j,4)+0.5D0*&
                (nu(i,j-1)*(max(zone%species(n)%var(1)%temperature(i,j-1),Teps))**2.5d0&
                +nu(i,j)*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0)*0.5d0&
                *(zone%metric_coefficients%G(i,j-1)+zone%metric_coefficients%G(i,j))*0.5D0&
                *(zone%species(n)%var(2)%velocity(i,j)**2-zone%species(n)%var(2)%velocity(i,j-1)**2)&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))*(1.d0-zone%masks%chi4(i,j-1))*(1.d0-zone%masks%chi4(i,j))
           zone%species(n)%fluxes%fluxE(i,j,4)=zone%species(n)%fluxes%fluxE(i,j,4)-Fluxes(i,j,4)
        end do
     end do
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%sources%SE=zone%species(n)%sources%SE+Source
  end do
  deallocate(Fluxes,Source,nu)
end subroutine add_temperature_parallel_viscosity_implicit_term
