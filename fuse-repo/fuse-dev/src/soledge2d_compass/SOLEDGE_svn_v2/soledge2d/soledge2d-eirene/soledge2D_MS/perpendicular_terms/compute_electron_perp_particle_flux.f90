subroutine compute_electron_perp_particle_flux(zone)
  use all_variables, only : zones, global_parameters
  use MZone
  use MOperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: Nx,Nz
  integer*4 :: i,j
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  !compute source term for ion temperature equation
  zone%species(0)%fluxes%fluxn=0.d0
  do i=1,Nx
     do j=1,Nz
        do n=1,global_parameters%N_ions
           !north
           zone%species(0)%fluxes%fluxn(i,j,1)=zone%species(0)%fluxes%fluxn(i,j,1)+&
                zone%species(n)%fluxes%fluxn(i,j,1)*zone%species(n)%charge
           !south
           zone%species(0)%fluxes%fluxn(i,j,2)=zone%species(0)%fluxes%fluxn(i,j,2)+&
                zone%species(n)%fluxes%fluxn(i,j,2)*zone%species(n)%charge
           !east
           zone%species(0)%fluxes%fluxn(i,j,3)=zone%species(0)%fluxes%fluxn(i,j,3)+&
                zone%species(n)%fluxes%fluxn(i,j,3)*zone%species(n)%charge
           !west
           zone%species(0)%fluxes%fluxn(i,j,4)=zone%species(0)%fluxes%fluxn(i,j,4)+&
                zone%species(n)%fluxes%fluxn(i,j,4)*zone%species(n)%charge
        end do
     end do
  end do
end subroutine compute_electron_perp_particle_flux
