subroutine compute_collisional_velocity_heat_flux_divergence(zone)
  use all_variables, only : global_parameters, flags, reference_parameters
  use Mzone
  use Moperator
  use Mphysics
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: n,m
  integer*4 :: Nx,Nz
  real*8,allocatable :: Fluxes_E(:,:,:)
  real*8,allocatable :: Source(:,:)
  integer*4 :: i,j
  real*8 :: R0, rs0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes_E(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  do n=0,global_parameters%n_ions
     Fluxes_E=0.D0
     do m=0,global_parameters%n_ions
        do i=1,Nx
           do j=1,Nz
              Fluxes_E(i,j,3)=Fluxes_E(i,j,3)+(zone%species(n)%coupling_terms%qu(i,j,m)+&
                   zone%species(n)%coupling_terms%qu(i,j+1,m))*0.5D0&
                   *zone%metric_coefficients%sinepitch_east(i,j)*(2.d0*pi*R0/rs0)
              Fluxes_E(i,j,4)=Fluxes_E(i,j,4)+(zone%species(n)%coupling_terms%qu(i,j,m)+&
                   zone%species(n)%coupling_terms%qu(i,j-1,m))*0.5D0&
                   *zone%metric_coefficients%sinepitch_west(i,j)*(2.d0*pi*R0/rs0)
           end do
        end do
     end do
     zone%species(n)%fluxes%fluxE=zone%species(n)%fluxes%fluxE+Fluxes_E
     Source=divergence(zone,Fluxes_E,Nx,Nz)
     zone%species(n)%sources%SE=zone%species(n)%sources%SE-Source
  end do
  deallocate(Fluxes_E,Source)
end subroutine compute_collisional_velocity_heat_flux_divergence
