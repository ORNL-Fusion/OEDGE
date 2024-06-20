subroutine add_uEt_penalisation_terms(zone)
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
           call add_uEt_penalisation_terms_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_uEt_penalisation_terms


subroutine add_uEt_penalisation_terms_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters, flags
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: dt
  real*8 :: m_i
  REAL*8 :: grad1,grad2,etaT,Teps
  real*8 :: nu,nu_m1,nu_p1,Gem, Gep,R0,c0
  dt=global_variables%dt
  etaT=penalisation_parameters%eta**2
  Teps=global_variables%Teps
  R0 = reference_parameters%geometry%R0
  c0 = reference_parameters%fields%c0
  grad1=0.d0
  grad2=0.d0
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
end subroutine add_uEt_penalisation_terms_point
