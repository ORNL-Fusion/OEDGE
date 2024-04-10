subroutine add_RHS_pe_Ohm_term(zone)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables
  use Moperator
  use MZone
  use MDefinitions
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: i,j,k
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  real*8,allocatable :: D(:,:)
  real*8 :: eta_para0,delta,beta
  real*8 :: rs0, R0
  real*8 :: Teps
  Teps=global_variables%Teps
  eta_para0=8.e-4*reference_parameters%fields%T0eV**(-1.5d0)
  delta=(reference_parameters%fields%tau0*reference_parameters%geometry%rs0)/&
       (reference_parameters%fields%n0*m_u*eta_para0*(2.d0*pi*reference_parameters%geometry%R0))

!  delta=(reference_parameters%fields%tau0*reference_parameters%geometry%rs0)/&
!       (reference_parameters%fields%n0*m_u*eta_para0*(reference_parameters%geometry%rs0))
  
!  delta=0.D0
  
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0

  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  Fluxes=0.D0
  allocate(D(0:Nx+1,0:Nz+1))
  do i=0,Nx+1
     do j=0,Nz+1
        D(i,j)=max(zone%species(0)%var(STEP_NEW)%temperature(i,j),Teps)**(1.5d0) 
     end do
  end do
  do i=1,Nx
     do j=1,Nz
        !east
        Fluxes(i,j,3)=(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5D0&
             *(zone%species(0)%var(1)%density(i,j+1)*zone%species(0)%var(1)%temperature(i,j+1)&
             -zone%species(0)%var(1)%density(i,j)*zone%species(0)%var(1)%temperature(i,j))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
        Fluxes(i,j,3)=Fluxes(i,j,3)/(0.5D0*(zone%species(0)%var(1)%density(i,j+1)+zone%species(0)%var(1)%density(i,j)))
        !adding friction force
        Fluxes(i,j,3)=Fluxes(i,j,3)+0.71d0*(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5D0&
             *(zone%species(0)%var(1)%temperature(i,j+1)-zone%species(0)%var(1)%temperature(i,j))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
        Fluxes(i,j,3)=Fluxes(i,j,3)*delta*(D(i,j)+D(i,j+1))*0.5D0
        !west
        Fluxes(i,j,4)=(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5D0&
             *(zone%species(0)%var(1)%density(i,j)*zone%species(0)%var(1)%temperature(i,j)&
             -zone%species(0)%var(1)%density(i,j-1)*zone%species(0)%var(1)%temperature(i,j-1))&
             /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
        Fluxes(i,j,4)=Fluxes(i,j,4)/(0.5D0*(zone%species(0)%var(1)%density(i,j)+zone%species(0)%var(1)%density(i,j-1)))
        !adding friction force
        Fluxes(i,j,4)=Fluxes(i,j,4)+0.71d0*(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5D0&
             *(zone%species(0)%var(1)%temperature(i,j)-zone%species(0)%var(1)%temperature(i,j-1))&
             /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
        Fluxes(i,j,4)=Fluxes(i,j,4)*delta*(D(i,j)+D(i,j-1))*0.5D0

        Fluxes(i,j,3)=Fluxes(i,j,3)*zone%metric_coefficients%sinepitch_east(i,j)
        Fluxes(i,j,4)=Fluxes(i,j,4)*zone%metric_coefficients%sinepitch_west(i,j)
     end do
  end do
  zone%electric_fields(1)%j_parallel=Fluxes ! first contribution to Ohm current 
  !  (the phi one is implicit -- see after inversion -- add_phi_Ohm_term.f90)
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%electric_fields(1)%RHS=zone%electric_fields(1)%RHS+Source*global_variables%dt_vort
  deallocate(D,Fluxes,Source)
end subroutine add_RHS_pe_Ohm_term
