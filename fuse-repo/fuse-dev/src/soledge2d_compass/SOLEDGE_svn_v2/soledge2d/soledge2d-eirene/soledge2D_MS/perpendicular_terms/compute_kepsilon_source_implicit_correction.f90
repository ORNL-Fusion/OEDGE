subroutine compute_kepsilon_source_implicit_correction(zone)
  use all_variables, only : global_parameters
  use MZone
  use Moperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: i,j
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  Fluxes=0.D0
  !k
  do i=1,Nx
     do j=1,Nz
        Fluxes(i,j,3)=zone%species(1)%implicit_coefs%east_density(i,j)&
             *(zone%kepsilon(1)%k(i,j)+zone%kepsilon(1)%k(i,j+1))*0.5D0&
             *(zone%species(1)%var(2)%density(i,j+1)-zone%species(1)%var(2)%density(i,j))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
        Fluxes(i,j,4)=zone%species(1)%implicit_coefs%west_density(i,j)&
             *(zone%kepsilon(1)%k(i,j)+zone%kepsilon(1)%k(i,j-1))*0.5D0&
             *(zone%species(1)%var(2)%density(i,j)-zone%species(1)%var(2)%density(i,j-1))&
             /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
     end do
  end do
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%kepsilon(1)%Sk=zone%kepsilon(1)%Sk+Source
  Fluxes=0.D0
  !epsilon
  do i=1,Nx
     do j=1,Nz
        Fluxes(i,j,3)=zone%species(1)%implicit_coefs%east_density(i,j)&
             *(zone%kepsilon(1)%epsilon(i,j)+zone%kepsilon(1)%epsilon(i,j+1))*0.5D0&
             *(zone%species(1)%var(2)%density(i,j+1)-zone%species(1)%var(2)%density(i,j))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
        Fluxes(i,j,4)=zone%species(1)%implicit_coefs%west_density(i,j)&
             *(zone%kepsilon(1)%epsilon(i,j)+zone%kepsilon(1)%epsilon(i,j-1))*0.5D0&
             *(zone%species(1)%var(2)%density(i,j)-zone%species(1)%var(2)%density(i,j-1))&
             /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
     end do
  end do
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%kepsilon(1)%Sepsilon=zone%kepsilon(1)%Sepsilon+Source
  deallocate(Fluxes,Source)
end subroutine compute_kepsilon_source_implicit_correction
