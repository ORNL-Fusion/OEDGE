subroutine allocate_fluxes_and_sources(zone)
  use all_variables, only : global_parameters
  use MZone 
  implicit none
  type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: n
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     allocate(zone%species(n)%sources%Sn(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%SG(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%SE(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%Sn_am(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%SG_am(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%SE_am(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%rad(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%Sn_n(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%Sn_G(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%Sn_E(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%Volumic_sources_n(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%Volumic_sources_G(1:Nx,1:Nz))
     allocate(zone%species(n)%sources%Volumic_sources_E(1:Nx,1:Nz))
     !third dimension for fluxes NSEW=1234
     allocate(zone%species(n)%fluxes%fluxn(1:Nx,1:Nz,1:4))
     allocate(zone%species(n)%fluxes%fluxG(1:Nx,1:Nz,1:4))
     allocate(zone%species(n)%fluxes%fluxE(1:Nx,1:Nz,1:4))
     zone%species(n)%sources%Sn_n=0.d0
     zone%species(n)%sources%Sn_G=0.d0
     zone%species(n)%sources%Sn_E=0.d0
  end do
end subroutine allocate_fluxes_and_sources
