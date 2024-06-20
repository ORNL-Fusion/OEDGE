subroutine compute_velocity_source_implicit_correction(zone)
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
  do n=1,global_parameters%N_ions
     Fluxes=0.D0
     !east
     do i=1,Nx
        do j=1,Nz
           Fluxes(i,j,3)=zone%species(n)%implicit_coefs%east_density(i,j)&
                *(zone%species(n)%var(1)%velocity(i,j)+zone%species(n)%var(1)%velocity(i,j+1))*0.5D0&
                *(zone%species(n)%var(2)%density(i,j+1)-zone%species(n)%var(2)%density(i,j))&
                /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
           zone%species(n)%fluxes%fluxG(i,j,3)=zone%species(n)%fluxes%fluxG(i,j,3)-Fluxes(i,j,3)
           Fluxes(i,j,4)=zone%species(n)%implicit_coefs%west_density(i,j)&
                *(zone%species(n)%var(1)%velocity(i,j)+zone%species(n)%var(1)%velocity(i,j-1))*0.5D0&
                *(zone%species(n)%var(2)%density(i,j)-zone%species(n)%var(2)%density(i,j-1))&
                /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
           zone%species(n)%fluxes%fluxG(i,j,4)=zone%species(n)%fluxes%fluxG(i,j,4)-Fluxes(i,j,4)
        end do
     end do
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%sources%SG=zone%species(n)%sources%SG+Source
  end do
  deallocate(Fluxes,Source)
end subroutine compute_velocity_source_implicit_correction
