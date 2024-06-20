subroutine compute_explicit_pseudo_convective_perp_source_termsI(zone)
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
  do n=1,global_parameters%N_ions
     !compute source term for ion temperature equation
     Field=5.D0/2.D0*zone%species(n)%var(1)%temperature+&
          1.D0/2.D0*zone%species(n)%var(1)%velocity*&
          zone%species(n)%var(1)%velocity*zone%species(n)%element%mass
     do i=1,Nx
        do j=1,Nz
           !north
           Fluxes(i,j,1)=0.5D0*(Field(i+1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,1)
           !south
           Fluxes(i,j,2)=0.5D0*(Field(i-1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,2)
           !east
           Fluxes(i,j,3)=0.5D0*(Field(i,j+1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,3)
           !west
           Fluxes(i,j,4)=0.5D0*(Field(i,j-1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,4)
        end do
     end do
     Field=zone%species(n)%var(1)%velocity*zone%species(n)%element%mass !friction contribution
     do i=1,Nx
        do j=1,Nz
           !north
           Fluxes(i,j,1)=Fluxes(i,j,1)+0.5D0*(Field(i+1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxG(i,j,1)
           !south
           Fluxes(i,j,2)=Fluxes(i,j,2)+0.5D0*(Field(i-1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxG(i,j,2)
           !east
           Fluxes(i,j,3)=Fluxes(i,j,3)+0.5D0*(Field(i,j+1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxG(i,j,3)
           !west
           Fluxes(i,j,4)=Fluxes(i,j,4)+0.5D0*(Field(i,j-1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxG(i,j,4)
        end do
     end do
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%fluxes%fluxE=zone%species(n)%fluxes%fluxE+Fluxes
     zone%species(n)%sources%SE=zone%species(n)%sources%SE-Source
     !compute source term for velocity equation
     Field=zone%species(n)%var(1)%velocity
     do i=1,Nx
        do j=1,Nz
           !north
           Fluxes(i,j,1)=0.5D0*(Field(i+1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,1)
           !south
           Fluxes(i,j,2)=0.5D0*(Field(i-1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,2)
           !east
           Fluxes(i,j,3)=0.5D0*(Field(i,j+1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,3)
           !west
           Fluxes(i,j,4)=0.5D0*(Field(i,j-1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,4)
        end do
     end do
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%fluxes%fluxG=zone%species(n)%fluxes%fluxG+Fluxes
     zone%species(n)%sources%SG=zone%species(n)%sources%SG-Source
  end do
  deallocate(Fluxes,Source,Field)
end subroutine compute_explicit_pseudo_convective_perp_source_termsI
