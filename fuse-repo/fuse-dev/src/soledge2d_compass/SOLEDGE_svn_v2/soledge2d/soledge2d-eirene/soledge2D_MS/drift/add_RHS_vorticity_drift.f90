subroutine add_RHS_vorticity_drift(zone)
  use all_variables,only : global_variables, global_parameters, reference_parameters, drift_flags
  use MZone
  use Moperator
  use MDrift_perp
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz,i,j
  real*8, allocatable :: Fluxes(:,:,:)
  real*8, allocatable :: Source(:,:) , v_psi(:,:) , v_theta(:,:), Field(:,:)
  real*8 :: R0,c0,rs0
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  c0=reference_parameters%fields%c0
  Nx = zone%mesh%Nx
  Nz = zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  allocate(v_psi(0:Nx+1,0:Nz+1))
  allocate(v_theta(0:Nx+1,0:Nz+1))
  allocate(Field(0:Nx+1,0:Nz+1))
!!$  !! compute source term for jExB = - W * ExB
!!$  v_psi=2.d0*pi*(zone%species(1)%drifts%uEp)
!!$  v_theta=0.d0!2.d0*pi*R0/c0*(zone%species(1)%drifts%uEt)/(2.d0*pi)
!!$  Field=zone%electric_fields(1)%vorticity
!!$  Fluxes=explicit_perp_drift(Field,v_psi,v_theta,zone,Nx,Nz)
!!$!  call set_perp_drift_on_core(zone,Fluxes,Nx,Nz)
!!$!  call set_perp_drift_on_the_wall(zone,Fluxes,Nx,Nz)
!!$  call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
!!$!  call set_perp_flux_to_zero_on_core(zone,Fluxes,Nx,Nz)
!!$
!!$  if(drift_flags%jExBW) then
!!$     zone%species(1)%drifts%jExB = Fluxes
!!$  else
!!$     zone%species(1)%drifts%jExB = 0.D0
!!$  end if
!!$
!!$  !! compute source term for jBxDB = - W * BxDB
!!$  v_psi=2.d0*pi*(zone%species(1)%drifts%uBp)
!!$  v_theta=0.d0!2.d0*pi*R0/c0*(zone%species(1)%drifts%uBt)/(2.d0*pi)
!!$  Field=zone%electric_fields(1)%vorticity
!!$  Fluxes=explicit_perp_drift(Field,v_psi,v_theta,zone,Nx,Nz)
!!$!  call set_perp_drift_on_core(zone,Fluxes,Nx,Nz)
!!$!  call set_perp_drift_on_the_wall(zone,Fluxes,Nx,Nz)
!!$  call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
!!$!  call set_perp_flux_to_zero_on_core(zone,Fluxes,Nx,Nz)
!!$
!!$  if(drift_flags%jgradBW) then
!!$     zone%species(1)%drifts%jBxDB = Fluxes
!!$  else
!!$     zone%species(1)%drifts%jBxDB = 0.D0
!!$  end if

  !! compute source term j* = n e (BxgradB_i - BxgradB_e)
  v_psi=-reference_parameters%fields%tau0/(global_parameters%element_list(1)%mass*reference_parameters%fields%W0)*reference_parameters%fields%n0*ev*&
	(zone%species(1)%drifts%uBp-zone%species(0)%drifts%uBp)
  v_theta=-reference_parameters%fields%tau0/(global_parameters%element_list(1)%mass*reference_parameters%fields%W0)*reference_parameters%fields%n0*ev/(2.d0*pi)*&
	(zone%species(1)%drifts%uBt-zone%species(0)%drifts%uBt)
  Field=zone%species(1)%var(1)%density
  Fluxes=explicit_perp_drift(Field,v_psi,v_theta,zone,Nx,Nz)
!  call set_perp_drift_on_core(zone,Fluxes,Nx,Nz)
!  call set_perp_drift_on_the_wall(zone,Fluxes,Nx,Nz)
!  call set_jdiam_on_the_wall(zone,Fluxes,Nx,Nz)
!  call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
!  call set_perp_flux_to_zero_on_core(zone,Fluxes,Nx,Nz)

  if(drift_flags%jdiam) then
     zone%species(1)%drifts%jdiam = Fluxes
  else
     zone%species(1)%drifts%jdiam = 0.D0
  end if

  Fluxes=zone%species(1)%drifts%jExB(1:Nx,1:Nz,1:4)+zone%species(1)%drifts%jBxDB(1:Nx,1:Nz,1:4)
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%electric_fields(1)%RHS=zone%electric_fields(1)%RHS+Source*global_variables%dt

  Fluxes=zone%species(1)%drifts%jdiam(1:Nx,1:Nz,1:4)
  Source=divergence(zone,Fluxes,Nx,Nz)
!  call correct_jdiam_source(zone,Source,Nx,Nz)
  zone%electric_fields(1)%RHS=zone%electric_fields(1)%RHS+Source*global_variables%dt

!!$  do i=1,Nx
!!$     write(2000+zone%number,100) (Source(i,j),j=1,Nz)
!!$  end do
!!$100 format(512es15.7)

  deallocate(Fluxes,Source,v_psi,v_theta,Field)
end subroutine add_RHS_vorticity_drift
