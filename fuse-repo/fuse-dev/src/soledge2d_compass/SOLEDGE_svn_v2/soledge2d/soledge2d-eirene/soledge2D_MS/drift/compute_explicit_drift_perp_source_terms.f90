subroutine compute_explicit_drift_perp_source_terms(zone)
  use all_variables, only : zones, global_parameters, reference_parameters
  use MZone
  use MZone
  use MOperator
  use MDrift_perp
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: n,i,j
  integer*4 :: Nx,Nz
  real*8,allocatable :: v_psi(:,:),v_theta(:,:)
  real*8,allocatable :: Field(:,:)
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  real*8 :: R0,c0,rs0
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  c0=reference_parameters%fields%c0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  allocate(v_psi(0:Nx+1,0:Nz+1))
  allocate(v_theta(0:Nx+1,0:Nz+1))
  allocate(Field(0:Nx+1,0:Nz+1))
  do n=1,global_parameters%N_ions

     !! Curvature drift
     v_psi=2.d0*pi*R0/c0*(zone%species(n)%drifts%uBp)
     v_theta=2.d0*pi*R0/c0*(zone%species(n)%drifts%uBt)/(2.d0*pi)
     !! compute source term for density equation
     Field=zone%species(n)%var(1)%density
     Fluxes=explicit_perp_drift(Field,v_psi,v_theta,zone,Nx,Nz)
     call set_perp_drift_on_core(zone,Fluxes,Nx,Nz)
     call set_perp_drift_on_the_wall(zone,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%fluxes%fluxn=zone%species(n)%fluxes%fluxn-Fluxes
     zone%species(n)%sources%Sn=zone%species(n)%sources%Sn+Source
     !! compute source term for velocity equation
     Field=zone%species(n)%var(1)%Gamma
     Fluxes=explicit_perp_drift(Field,v_psi,v_theta,zone,Nx,Nz)
     call set_perp_drift_on_core(zone,Fluxes,Nx,Nz)
     call set_perp_drift_on_the_wall(zone,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%sources%SG=zone%species(n)%sources%SG+Source
	
!!$     !! ExB drift
!!$     v_psi=2.d0*pi*R0/c0*(zone%species(n)%drifts%uEp)
!!$     v_theta=2.d0*pi*R0/c0*(zone%species(n)%drifts%uEt)/(2.d0*pi)
!!$     !! compute source term for density equation
!!$     Field=zone%species(n)%var(1)%density
!!$     Fluxes=explicit_perp_drift(Field,v_psi,v_theta,zone,Nx,Nz)
!!$     call set_perp_drift_on_core(zone,Fluxes,Nx,Nz)
!!$     call set_perp_drift_on_the_wall(zone,Fluxes,Nx,Nz)
!!$     Source=divergence(zone,Fluxes,Nx,Nz)
!!$     zone%species(n)%fluxes%fluxn=zone%species(n)%fluxes%fluxn-Fluxes
!!$     zone%species(n)%sources%Sn=zone%species(n)%sources%Sn+Source


  end do

!!$  do n=0,global_parameters%N_ions
!!$     !! compute source term for temperature equation
!!$     !! Curvature drift
!!$     v_psi=2.d0*pi*R0/c0*(zone%species(n)%drifts%uBp)
!!$     v_theta=2.d0*pi*R0/c0*(zone%species(n)%drifts%uBt)/(2.d0*pi)
!!$     Field=2.5d0*zone%species(n)%var(1)%temperature*zone%species(n)%var(1)%density+0.5d0*global_parameters%element_list(n)%mass*&
!!$		zone%species(n)%var(1)%velocity*zone%species(n)%var(1)%velocity*zone%species(n)%var(1)%density
!!$     Fluxes=explicit_perp_drift(Field,v_psi,v_theta,zone,Nx,Nz)
!!$     call set_perp_drift_on_core(zone,Fluxes,Nx,Nz)
!!$     call set_perp_drift_on_the_wall(zone,Fluxes,Nx,Nz)
!!$     Source=divergence(zone,Fluxes,Nx,Nz)
!!$     zone%species(n)%sources%SE=zone%species(n)%sources%SE+Source
!!$
!!$  end do

  deallocate(Fluxes,Source,v_psi,v_theta,Field)
end subroutine compute_explicit_drift_perp_source_terms
