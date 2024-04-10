subroutine add_temperature_penalisation_terms(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4 :: i,j,n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        do n=0,global_parameters%N_ions
           call add_temperature_penalisation_terms_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_temperature_penalisation_terms


subroutine add_temperature_penalisation_terms_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters, flags, transport_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: dt
  real*8 :: m_i
  REAL*8 :: grad1,grad2,etaT,Teps
  real*8 :: nu,nu_m1,nu_p1,Gem, Gep,R0,c0,rs0
  real*8 :: gamma, delta_e
  dt=global_variables%dt
  etaT=penalisation_parameters%eta**2
  Teps=global_variables%Teps
  if(n.ge.1) then
     nu=zone%species(n)%transport_para%nu(i,j)/&
          zone%species(n)%var(1)%log_Lambda(i,j)
     nu_p1=zone%species(n)%transport_para%nu(i,j+1)/&
          zone%species(n)%var(1)%log_Lambda(i,j+1)
     nu_m1=zone%species(n)%transport_para%nu(i,j-1)/&
          zone%species(n)%var(1)%log_Lambda(i,j-1)
  else
     nu=0.D0
     nu_m1=0.D0
     nu_p1=0.D0
  end if
  R0 = reference_parameters%geometry%R0
  rs0 = reference_parameters%geometry%rs0
  c0 = reference_parameters%fields%c0

  Gem = zone%species(n)%var(2)%density(i,j-1)*&
       (zone%species(n)%drifts%uEt(i,j-1)+zone%species(n)%drifts%uBt(i,j-1))&
       *sqrt(zone%metric_coefficients%ctt(i,j-1))/zone%metric_coefficients%G(i,j-1)&
       *(2.d0*pi*R0/rs0)
  Gep = zone%species(n)%var(2)%density(i,j+1)*&
       (zone%species(n)%drifts%uEt(i,j+1)+zone%species(n)%drifts%uBt(i,j+1))&
       *sqrt(zone%metric_coefficients%ctt(i,j+1))/zone%metric_coefficients%G(i,j+1)&
       *(2.d0*pi*R0/rs0)


  m_i=zone%species(n)%element%mass
  if(n.eq.0) then
     delta_e=transport_parameters%delta_e
     gamma=2.d0/(1.D0-delta_e)-0.5d0*log(&
          (2.D0*pi*m_e/(m_u*global_parameters%element_list(1)%mass))*&
          (1+zone%species(1)%var(1)%temperature(i,j-1)/zone%species(0)%var(1)%temperature(i,j-1))&
          /(1.D0-delta_e)**2.D0)
     gamma=min(max(gamma,2.5d0),40.d0)
  else
     gamma=zone%species(n)%transport_para%gamma
  end if
  grad1=(gamma-2.5D0)*(zone%species(n)%var(2)%Gamma(i,j-1)+3.d0/5.d0*Gem)&
       +0.5D0*(nu*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0&
       +nu_m1*(max(zone%species(n)%var(1)%temperature(i,j-1),Teps))**2.5d0)*0.5D0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5D0&
       *(zone%species(n)%var(2)%velocity(i,j)**2-zone%species(n)%var(2)%velocity(i,j-1)**2)&
       /(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
  grad1=grad1/(-(zone%metric_coefficients%G(i,j-1)+zone%metric_coefficients%G(i,j))*0.5d0*&
       (zone%species(n)%transport_para%kappa(i,j-1)/zone%species(n)%var(1)%log_lambda(i,j-1)&
       *(max(zone%species(n)%var(1)%temperature(i,j-1),Teps))**1.5d0&
       +zone%species(n)%transport_para%kappa(i,j)/zone%species(n)%var(1)%log_lambda(i,j)&
       *(max(zone%species(n)%var(1)%temperature(i,j),Teps))**1.5d0)*0.5D0)&
       *(1.D0-zone%masks%chi2(i,j-1))
  if(n.eq.0) then
     delta_e=transport_parameters%delta_e
     gamma=2.d0/(1.D0-delta_e)-0.5d0*log(&
          (2.D0*pi*m_e/(m_u*global_parameters%element_list(1)%mass))*&
          (1+zone%species(1)%var(1)%temperature(i,j+1)/zone%species(0)%var(1)%temperature(i,j+1))&
          /(1.D0-delta_e)**2.D0)
     gamma=min(max(gamma,2.5d0),40.d0)
  else
     gamma=zone%species(n)%transport_para%gamma
  end if
  grad2=(gamma-2.5D0)*(zone%species(n)%var(2)%Gamma(i,j+1)+3.d0/5.d0*Gep)&
       +0.5D0*(nu*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0&
       +nu_p1*(max(zone%species(n)%var(1)%temperature(i,j+1),Teps))**2.5d0)*0.5D0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5D0&
       *(zone%species(n)%var(2)%velocity(i,j+1)**2-zone%species(n)%var(2)%velocity(i,j)**2)&
       /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
  grad2=grad2/(-(zone%metric_coefficients%G(i,j+1)+zone%metric_coefficients%G(i,j))*0.5d0*&
       (zone%species(n)%transport_para%kappa(i,j+1)/zone%species(n)%var(1)%log_lambda(i,j+1)&
       *(max(zone%species(n)%var(1)%temperature(i,j+1),Teps))**1.5d0&
       +zone%species(n)%transport_para%kappa(i,j)/zone%species(n)%var(1)%log_lambda(i,j)&
       *(max(zone%species(n)%var(1)%temperature(i,j),Teps))**1.5d0)*0.5d0)&
       *(1.D0-zone%masks%chi2(i,j+1))
  zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)&
       +zone%masks%chi4(i,j)/penalisation_parameters%eta2*global_variables%min_temperature&
       +zone%masks%chi1(i,j)/etaT*grad1&
       -zone%masks%chi3(i,j)/etaT*grad2
  zone%species(n)%tridiag%a(i,j)=zone%species(n)%tridiag%a(i,j)&
       -zone%masks%chi1(i,j)/etaT/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)&
       +(zone%masks%chi1(i,j)/etaT/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1)))&
       +(zone%masks%chi3(i,j)/etaT/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j)))&
       +1.5d0*zone%masks%chi4(i,j)/penalisation_parameters%eta2
  zone%species(n)%tridiag%c(i,j)=zone%species(n)%tridiag%c(i,j)&
       -zone%masks%chi3(i,j)/etaT/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
end subroutine add_temperature_penalisation_terms_point
