subroutine compute_transport_coefficients()
  implicit none
  call compute_perp_DD_diffusivities()
  call compute_para_DD_diffusivities()
  call MD_broadcast_transport_coefficients()
end subroutine compute_transport_coefficients


subroutine compute_perp_DD_diffusivities()
  use all_variables, only : zones, global_parameters, ballooning_parameters,&
       reference_parameters, transport_parameters
  use Mphysics
  implicit none
  integer*4 :: k
  integer*4 :: n
  integer*4 :: n_element
  reference_parameters%geometry%A=2.*pi*reference_parameters%geometry%R0&
       /reference_parameters%geometry%rs0
  !compute reference uniform diffusivities
  do k=1,global_parameters%N_Zones
     !for the electrons
     zones(k)%species(0)%transport_perp%chi_p=transport_parameters%chie_p/&
          (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
     zones(k)%species(0)%transport_perp%chi_t=transport_parameters%chie_t/&
          (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
     zones(k)%species(0)%transport_perp%zeta_p=transport_parameters%zeta_p/&
          (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
     zones(k)%species(0)%transport_perp%zeta_t=transport_parameters%zeta_t/&
          (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
     !for the ions
     do n=1,global_parameters%N_ions
        n_element=global_parameters%ions_list(n,1)
        zones(k)%species(n)%transport_perp%D_p=transport_parameters%Dn_p(n_element)/&
             (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
        zones(k)%species(n)%transport_perp%D_t=transport_parameters%Dn_t(n_element)/&
             (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
        zones(k)%species(n)%transport_perp%nu_p=transport_parameters%nu_p(n_element)/&
             (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
        zones(k)%species(n)%transport_perp%nu_t=transport_parameters%nu_t(n_element)/&
             (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
        zones(k)%species(n)%transport_perp%chi_p=transport_parameters%chii_p(n_element)/&
             (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
        zones(k)%species(n)%transport_perp%chi_t=transport_parameters%chii_t(n_element)/&
             (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
     end do
  end do
  select case(ballooning_parameters%ballooning_model) 
  case(0)
     call set_no_ballooning()
  case(1)
     call set_gaussian_ballooning()
  case(2)
     call load_ballooning_from_file()
  case(4)
     call load_transport_coefficients()
  end select
  call compute_ballooning_weighted_Score()
end subroutine compute_perp_DD_diffusivities


subroutine compute_para_DD_diffusivities()
  use all_variables, only : zones,global_parameters,reference_parameters,transport_parameters
  use Mphysics
  implicit none
  integer*4 :: k,n
  real*8 :: n0,T0eV,R0,tau0
  real*8 :: m_i
  n0=reference_parameters%fields%n0
  T0eV=reference_parameters%fields%T0eV
  tau0=reference_parameters%fields%tau0
  R0=reference_parameters%geometry%R0
  do k=1,global_parameters%N_zones
     !kappa0 for electrons
     zones(k)%species(0)%transport_para%flux_limiter=transport_parameters%flux_limiter(0)
     !kappa = kappa0*(Te^5/2)/log_Lambda/Zeff^2
     zones(k)%species(0)%transport_para%kappa0=3.16*1.09d16*(eV*eV)/(m_e*1000.**(1.5))
     !dedimensionalised
     zones(k)%species(0)%transport_para%kappa0=zones(k)%species(0)%transport_para%kappa0&
          *(T0eV)**(2.5)*tau0/((2.*pi*R0)**2.*eV*n0)
     zones(k)%species(0)%transport_para%gamma=4.5d0
     !ions
     do n=1,global_parameters%N_ions
        zones(k)%species(n)%transport_para%flux_limiter=transport_parameters%flux_limiter(global_parameters%ions_list(n,1))
        zones(k)%species(n)%transport_para%flux_limiter_nu=transport_parameters%flux_limiter_nu(global_parameters%ions_list(n,1))
        !kappa0 for ions
        !kappa = kappa0*(Ti^5/2)/log_Lambda/Zeff^2
        m_i=zones(k)%species(n)%element%mass ! in atomic mass unit
        zones(k)%species(n)%transport_para%kappa0=3.9*6.6d17*sqrt(m_i)&
             *(eV*eV)/(m_i*m_u*1000.**(1.5)*zones(k)%species(n)%charge**2.)
        !dedimensionalised
        zones(k)%species(n)%transport_para%kappa0=zones(k)%species(n)%transport_para%kappa0&
             *(T0eV)**(2.5)*tau0/((2.*pi*R0)**2.*eV*n0)
        !nu0 for ions
        !nu = nu0*(Ti^5/2)/log_Lambda/Zeff^2
        zones(k)%species(n)%transport_para%nu0=0.96*6.6d17*sqrt(m_i)&
             *eV/(1000.**(1.5)*zones(k)%species(n)%charge**2.)
        !dedimensionalised
        zones(k)%species(n)%transport_para%nu0=zones(k)%species(n)%transport_para%nu0&
             *(T0eV)**(2.5)*tau0/((2.*pi*R0)**2.*n0*m_i*m_u)
        zones(k)%species(n)%transport_para%gamma=transport_parameters%gamma_i
     end do
  end do
end subroutine compute_para_DD_diffusivities
