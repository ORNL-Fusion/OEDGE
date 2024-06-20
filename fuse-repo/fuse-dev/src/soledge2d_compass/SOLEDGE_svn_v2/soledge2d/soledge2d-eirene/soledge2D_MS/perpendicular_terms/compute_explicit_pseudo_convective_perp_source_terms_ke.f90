subroutine compute_explicit_pseudo_convective_perp_source_terms_ke(zone)
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
  !compute source term for k equation
  Field=zone%kepsilon(1)%k
  Fluxes=0.D0
  do i=1,Nx
     do j=1,Nz
        do n=1,global_parameters%N_ions
           !north
           Fluxes(i,j,1)=Fluxes(i,j,1)+0.5D0*(Field(i+1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,1)*zone%species(n)%charge
           !south
           Fluxes(i,j,2)=Fluxes(i,j,2)+0.5D0*(Field(i-1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,2)*zone%species(n)%charge
           !east
           Fluxes(i,j,3)=Fluxes(i,j,3)+0.5D0*(Field(i,j+1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,3)*zone%species(n)%charge
           !west
           Fluxes(i,j,4)=Fluxes(i,j,4)+0.5D0*(Field(i,j-1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,4)*zone%species(n)%charge
        end do
     end do
  end do
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%kepsilon(1)%Sk=zone%kepsilon(1)%Sk-Source
  !compute source term for epsilon equation
  Field=zone%kepsilon(1)%epsilon
  Fluxes=0.D0
  do i=1,Nx
     do j=1,Nz
        do n=1,global_parameters%N_ions
           !north
           Fluxes(i,j,1)=Fluxes(i,j,1)+0.5D0*(Field(i+1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,1)*zone%species(n)%charge
           !south
           Fluxes(i,j,2)=Fluxes(i,j,2)+0.5D0*(Field(i-1,j)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,2)*zone%species(n)%charge
           !east
           Fluxes(i,j,3)=Fluxes(i,j,3)+0.5D0*(Field(i,j+1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,3)*zone%species(n)%charge
           !west
           Fluxes(i,j,4)=Fluxes(i,j,4)+0.5D0*(Field(i,j-1)+Field(i,j))&
                *zone%species(n)%fluxes%fluxn(i,j,4)*zone%species(n)%charge
        end do
     end do
  end do
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%kepsilon(1)%Sepsilon=zone%kepsilon(1)%Sepsilon-Source
  deallocate(Fluxes,Source,Field)
end subroutine compute_explicit_pseudo_convective_perp_source_terms_ke
