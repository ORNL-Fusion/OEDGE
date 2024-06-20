subroutine set_core_boundary_conditions_north(zone,STEP)
#include "compile_opt.inc"
  use all_variables, only : boundary_conditions, global_parameters&
       ,reference_parameters, element_variables, ballooning_parameters
  use Meirene_vars
  use Mfeedback_control
  use MZone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n,n_element
  integer*4 :: Nx,Nz
  real*8 :: total_flux
  real*8,allocatable :: phi_bc(:),temp1(:),temp2(:),temp3(:)
  real*8,allocatable :: D(:),cpp(:),v_pinch(:)
  real*8,allocatable :: grad_density(:),density(:),chi(:)
  real*8,allocatable :: total_density(:)
  real*8 :: rs0,c0,n0,T0eV,tau0
  integer*4 :: m
  real*8 :: Puff
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  n0=reference_parameters%fields%n0
  tau0=reference_parameters%fields%tau0
  T0eV=reference_parameters%fields%T0eV
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(total_density(1:Nz))
  allocate(phi_bc(1:Nz),temp1(1:Nz),temp2(1:Nz),temp3(1:Nz))
  allocate(D(1:Nz),cpp(1:Nz),v_pinch(1:Nz),chi(1:Nz))
  allocate(density(1:Nz),grad_density(1:Nz))
  !compute sum of ions density / used to compute concentration latter
  total_density=0.d0
  do n=1,global_parameters%N_ions
     total_density=total_density+zone%species(n)%var(STEP)%density(Nx,1:Nz)
  end do
  do n=1,global_parameters%N_ions
     n_element=global_parameters%ions_list(n,1)
     ! set density BC to the completely ionized ion for the element
     ! for the other ions, copy density (grad approx 0)
     if(zone%species(n)%charge.lt.global_parameters%element_list(n_element)%Z) then
        zone%species(n)%var(STEP)%density(Nx+1,1:Nz) = zone%species(n)%var(STEP)%density(Nx,1:Nz)
        zone%species(n)%var(STEP)%Gamma(Nx+1,1:Nz) = zone%species(n)%var(STEP)%Gamma(Nx,1:Nz)
        select case(boundary_conditions%BCT_model(n_element))
        case(0,1)
           zone%species(n)%var(STEP)%temperature(Nx+1,1:Nz) = zone%species(n)%var(STEP)%temperature(Nx,1:Nz)
        case(2) !flux proportional to concentration
           total_flux=boundary_conditions%BCTi(1) ! flux to share on species 1 BC
           phi_bc=total_flux*zone%species(1)%transport_perp%D_p(Nx,1:Nz)/&
                ballooning_parameters%core_weighted_integral
           phi_bc=phi_bc*zone%species(n)%var(STEP)%density(Nx,1:Nz)/total_density
           D=(zone%species(n)%transport_perp%D_p(Nx,1:Nz)&
                +zone%species(n)%transport_perp%D_p(Nx+1,1:Nz))*0.5
           chi=(zone%species(n)%transport_perp%chi_p(Nx,1:Nz)&
                +zone%species(n)%transport_perp%chi_p(Nx+1,1:Nz))*0.5
           cpp=(zone%metric_coefficients%cpp(Nx,1:Nz)&
                +zone%metric_coefficients%cpp(Nx+1,1:Nz))*0.5
           v_pinch=(zone%species(n)%transport_perp%v_pinch(Nx,1:Nz)&
                +zone%species(n)%transport_perp%v_pinch(Nx+1,1:Nz))*0.5
           grad_density=(zone%species(n)%var(STEP)%density(Nx+1,1:Nz)-zone%species(n)%var(STEP)%density(Nx,1:Nz))&
                /(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))
           density=(zone%species(n)%var(STEP)%density(Nx+1,1:Nz)+zone%species(n)%var(STEP)%density(Nx,1:Nz))*0.5
           temp1=-D*sqrt(cpp)*grad_density*rs0/tau0*rs0+v_pinch*rs0/n0*density
           temp2=2.5*temp1/2.-chi*density*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0
           temp3=2.5*temp1/2.+chi*density*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0
           zone%species(n)%var(STEP)%temperature(Nx+1,1:Nz)=-zone%species(n)%var(STEP)%temperature(Nx,1:Nz)*temp3/temp2&
                -phi_bc/temp2*rs0/(n0*eV*T0eV)
        end select
     else
        !############## boundary conditions for density ########################################
        select case(boundary_conditions%BCn_model(n_element))
        case(0) !dirichlet
           zone%species(n)%var(STEP)%density(Nx+1,1:Nz) = boundary_conditions%BCn(n_element)/n0
        case(1,2) !flux
           if(boundary_conditions%BCn_model(n_element).eq.1) then
              total_flux=boundary_conditions%BCn(n_element)*reference_parameters%geometry%Score
           else
              if(eirene_vars%feedback.eq.2) then
                 call puff_feedback(Control_data,Puff,error_data)
                 total_flux=element_variables(n_element)%total_flux_core_nonionized+Puff
              else
                 total_flux=element_variables(n_element)%total_flux_core_nonionized
              end if
           end if
           phi_bc=total_flux*zone%species(1)%transport_perp%D_p(Nx,1:Nz)/&
                ballooning_parameters%core_weighted_integral
           D=(zone%species(n)%transport_perp%D_p(Nx,1:Nz)&
                +zone%species(n)%transport_perp%D_p(Nx+1,1:Nz))*0.5
           cpp=(zone%metric_coefficients%cpp(Nx,1:Nz)&
                +zone%metric_coefficients%cpp(Nx+1,1:Nz))*0.5
           v_pinch=(zone%species(n)%transport_perp%v_pinch(Nx,1:Nz)&
                +zone%species(n)%transport_perp%v_pinch(Nx+1,1:Nz))*0.5
           temp1=-D*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0&
                +v_pinch*rs0/(2.*n0)
           temp2=D*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0&
                +v_pinch*rs0/(2.*n0)
           zone%species(n)%var(STEP)%density(Nx+1,1:Nz)=-zone%species(n)%var(STEP)%density(Nx,1:Nz)*temp2/temp1&
                -phi_bc*rs0/n0/temp1
        case(3) !dirichlet percentage of main ion
           zone%species(n)%var(STEP)%density(Nx+1,1:Nz) = boundary_conditions%BCn(n_element)*&
                zone%species(1)%var(STEP)%density(Nx+1,1:Nz)
        end select
        !#######################################################################################
        !############## boundary conditions for Gamma #######################################
        zone%species(n)%var(STEP)%Gamma(Nx+1,1:Nz)=zone%species(n)%var(STEP)%Gamma(Nx,1:Nz)
        !#######################################################################################
        !############## boundary conditions for temperature ####################################
        select case(boundary_conditions%BCT_model(n_element))
        case(0) !dirichlet
           zone%species(n)%var(STEP)%temperature(Nx+1,1:Nz)=boundary_conditions%BCTi(n_element)/T0eV
        case(1) !flux
           total_flux=boundary_conditions%BCTi(n_element)
           phi_bc=total_flux*zone%species(1)%transport_perp%D_p(Nx,1:Nz)/&
                ballooning_parameters%core_weighted_integral
           D=(zone%species(n)%transport_perp%D_p(Nx,1:Nz)&
                +zone%species(n)%transport_perp%D_p(Nx+1,1:Nz))*0.5
           chi=(zone%species(n)%transport_perp%chi_p(Nx,1:Nz)&
                +zone%species(n)%transport_perp%chi_p(Nx+1,1:Nz))*0.5
           cpp=(zone%metric_coefficients%cpp(Nx,1:Nz)&
                +zone%metric_coefficients%cpp(Nx+1,1:Nz))*0.5
           v_pinch=(zone%species(n)%transport_perp%v_pinch(Nx,1:Nz)&
                +zone%species(n)%transport_perp%v_pinch(Nx+1,1:Nz))*0.5
           grad_density=(zone%species(n)%var(STEP)%density(Nx+1,1:Nz)-zone%species(n)%var(STEP)%density(Nx,1:Nz))&
                /(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))
           density=(zone%species(n)%var(STEP)%density(Nx+1,1:Nz)+zone%species(n)%var(STEP)%density(Nx,1:Nz))*0.5
           temp1=-D*sqrt(cpp)*grad_density*rs0/tau0*rs0+v_pinch*rs0/n0*density
           temp2=2.5*temp1/2.-chi*density*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0
           temp3=2.5*temp1/2.+chi*density*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0
           zone%species(n)%var(STEP)%temperature(Nx+1,1:Nz)=-zone%species(n)%var(STEP)%temperature(Nx,1:Nz)*temp3/temp2&
                -phi_bc/temp2*rs0/(n0*eV*T0eV)
        case(2) !flux proportional to concentration
           total_flux=boundary_conditions%BCTi(1) ! flux to share on species 1 BC
           phi_bc=total_flux*zone%species(1)%transport_perp%D_p(Nx,1:Nz)/&
                ballooning_parameters%core_weighted_integral
           phi_bc=phi_bc*zone%species(n)%var(STEP)%density(Nx,1:Nz)/total_density
           D=(zone%species(n)%transport_perp%D_p(Nx,1:Nz)&
                +zone%species(n)%transport_perp%D_p(Nx+1,1:Nz))*0.5
           chi=(zone%species(n)%transport_perp%chi_p(Nx,1:Nz)&
                +zone%species(n)%transport_perp%chi_p(Nx+1,1:Nz))*0.5
           cpp=(zone%metric_coefficients%cpp(Nx,1:Nz)&
                +zone%metric_coefficients%cpp(Nx+1,1:Nz))*0.5
           v_pinch=(zone%species(n)%transport_perp%v_pinch(Nx,1:Nz)&
                +zone%species(n)%transport_perp%v_pinch(Nx+1,1:Nz))*0.5
           grad_density=(zone%species(n)%var(STEP)%density(Nx+1,1:Nz)-zone%species(n)%var(STEP)%density(Nx,1:Nz))&
                /(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))
           density=(zone%species(n)%var(STEP)%density(Nx+1,1:Nz)+zone%species(n)%var(STEP)%density(Nx,1:Nz))*0.5
           temp1=-D*sqrt(cpp)*grad_density*rs0/tau0*rs0+v_pinch*rs0/n0*density
           temp2=2.5*temp1/2.-chi*density*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0
           temp3=2.5*temp1/2.+chi*density*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0
           zone%species(n)%var(STEP)%temperature(Nx+1,1:Nz)=-zone%species(n)%var(STEP)%temperature(Nx,1:Nz)*temp3/temp2&
                -phi_bc/temp2*rs0/(n0*eV*T0eV)              
        end select
        !#######################################################################################
     end if
  end do
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   for electrons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !############## boundary conditions for density     ##########################################
  zone%species(0)%var(STEP)%density(Nx+1,1:Nz)=0
  do n=1,global_parameters%N_ions
     zone%species(0)%var(STEP)%density(Nx+1,1:Nz)=zone%species(0)%var(STEP)%density(Nx+1,1:Nz)&
          +zone%species(n)%var(STEP)%density(Nx+1,1:Nz)*zone%species(n)%charge
  end do
  !#############################################################################################
  !############## boundary conditions for temperature ##########################################
  select case(boundary_conditions%BCT_model(1)) !same model for electrons and for main ions
  case(0) !dirichlet
     zone%species(0)%var(STEP)%temperature(Nx+1,1:Nz)=boundary_conditions%BCTe/T0eV
  case(1,2) !flux
     total_flux=boundary_conditions%BCTe
     phi_bc=total_flux*zone%species(1)%transport_perp%D_p(Nx,1:Nz)/&    !one uses main ion ballooning
          ballooning_parameters%core_weighted_integral
     chi=(zone%species(0)%transport_perp%chi_p(Nx,1:Nz)&
          +zone%species(0)%transport_perp%chi_p(Nx+1,1:Nz))*0.5
     cpp=(zone%metric_coefficients%cpp(Nx,1:Nz)&
          +zone%metric_coefficients%cpp(Nx+1,1:Nz))*0.5
     temp1=0.D0
     do m=1,global_parameters%N_ions
        D=(zone%species(m)%transport_perp%D_p(Nx,1:Nz)&
             +zone%species(m)%transport_perp%D_p(Nx+1,1:Nz))*0.5
        v_pinch=(zone%species(m)%transport_perp%v_pinch(Nx,1:Nz)&
             +zone%species(m)%transport_perp%v_pinch(Nx+1,1:Nz))*0.5
        grad_density=(zone%species(m)%var(STEP)%density(Nx+1,1:Nz)-zone%species(m)%var(STEP)%density(Nx,1:Nz))&
             /(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))
        density=(zone%species(m)%var(STEP)%density(Nx+1,1:Nz)+zone%species(m)%var(STEP)%density(Nx,1:Nz))*0.5
        temp1=temp1+(-D*sqrt(cpp)*grad_density*rs0/tau0*rs0+v_pinch*rs0/n0*density)*zone%species(m)%charge
     end do
     density=(zone%species(0)%var(STEP)%density(Nx+1,1:Nz)+zone%species(0)%var(STEP)%density(Nx,1:Nz))*0.5
     temp2=2.5*temp1/2.-chi*density*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0
     temp3=2.5*temp1/2.+chi*density*sqrt(cpp)/(zone%mesh%x(Nx+1,1:Nz)-zone%mesh%x(Nx,1:Nz))*rs0/tau0*rs0
     zone%species(0)%var(STEP)%temperature(Nx+1,1:Nz)=-zone%species(0)%var(STEP)%temperature(Nx,1:Nz)*temp3/temp2&
          -phi_bc/temp2*rs0/(n0*eV*T0eV)
  end select
  !#############################################################################################
  deallocate(phi_bc,temp1,temp2,D,cpp,v_pinch,chi,density,grad_density,temp3)
end subroutine set_core_boundary_conditions_north
