subroutine compute_explicit_pseudo_convective_perp_source_termsE(zone)
  use all_variables, only : zones, global_parameters
  use MZone
  use MOperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: Nx,Nz
  real*8,allocatable :: Field(:,:)
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  integer*4 :: i,j
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  allocate(Field(0:Nx+1,0:Nz+1))
  !compute source term for ion temperature equation
  Field=5.D0/2.D0*zone%species(0)%var(1)%temperature
  Fluxes=0.D0
  do i=1,Nx
     do j=1,Nz
           !north
           Fluxes(i,j,1)=Fluxes(i,j,1)+0.5D0*(Field(i+1,j)+Field(i,j))&
                *zone%species(0)%fluxes%fluxn(i,j,1)
           !south
           Fluxes(i,j,2)=Fluxes(i,j,2)+0.5D0*(Field(i-1,j)+Field(i,j))&
                *zone%species(0)%fluxes%fluxn(i,j,2)
           !east
           Fluxes(i,j,3)=Fluxes(i,j,3)+0.5D0*(Field(i,j+1)+Field(i,j))&
                *zone%species(0)%fluxes%fluxn(i,j,3)
           !west
           Fluxes(i,j,4)=Fluxes(i,j,4)+0.5D0*(Field(i,j-1)+Field(i,j))&
                *zone%species(0)%fluxes%fluxn(i,j,4)
     end do
  end do
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%species(0)%fluxes%fluxE=zone%species(0)%fluxes%fluxE+Fluxes
  zone%species(0)%sources%SE=zone%species(0)%sources%SE-Source
  deallocate(Fluxes,Source,Field)
end subroutine compute_explicit_pseudo_convective_perp_source_termsE
