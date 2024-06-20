subroutine compute_implicit_tridiag_coefficients_FIELD_ke(n_ion,FIELD)
  use all_variables, only : global_parameters, zones, reference_parameters, kepsilon_param
  use Mdefinitions
  implicit none
  integer*4,intent(in) :: n_ion
  integer*4,intent(in) :: FIELD
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  real*8,allocatable :: ctt(:,:),cpt(:,:),cpp(:,:),G(:,:)
  real*8,allocatable :: D_p(:,:),D_t(:,:)
  real*8 :: east,west
  real*8 :: A
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(ctt(0:Nx+1,0:Nz+1),cpt(0:Nx+1,0:Nz+1))
     allocate(cpp(0:Nx+1,0:Nz+1),G(0:Nx+1,0:Nz+1))
     allocate(D_p(0:Nx+1,0:Nz+1),D_t(0:Nx+1,0:Nz+1))
     A=reference_parameters%geometry%A
     ctt=zones(k)%metric_coefficients%ctt
     cpt=zones(k)%metric_coefficients%cpt
     cpp=zones(k)%metric_coefficients%cpp
     G=zones(k)%metric_coefficients%G
     select case(FIELD)
     case(DENSITY_FIELD)
        D_p=zones(k)%species(n_ion)%transport_perp%D_p
        D_t=zones(k)%species(n_ion)%transport_perp%D_t
     case(VELOCITY_FIELD)
        D_p=zones(k)%species(n_ion)%transport_perp%nu_p
        D_t=zones(k)%species(n_ion)%transport_perp%nu_t
     case(TEMPERATURE_FIELD)
        D_p=zones(k)%species(n_ion)%transport_perp%chi_p
        D_t=zones(k)%species(n_ion)%transport_perp%chi_t
     case(K_FIELD)
        D_p=zones(k)%kepsilon(1)%mu_t/kepsilon_param%sigma_k
        D_t=0.d0!zones(k)%kepsilon(1)%mu_t/kepsilon_param%sigma_k
     case(EPSILON_FIELD)
        D_p=zones(k)%kepsilon(1)%mu_t/kepsilon_param%sigma_epsilon
        D_t=0.d0!zones(k)%kepsilon(1)%mu_t/kepsilon_param%sigma_epsilon
     end select
     do i=1,Nx
        do j=1,Nz
           east = &
                ((D_t(i,j)*(ctt(i,j)-1.d0/A**2.d0*G(i,j)**2.d0)&
                +(D_p(i,j)-D_t(i,j))*cpt(i,j)**2.d0/cpp(i,j))&
                +(D_t(i,j+1)*(ctt(i,j+1)-1.d0/A**2.d0*G(i,j+1)**2.d0)&
                +(D_p(i,j+1)-D_t(i,j+1))*cpt(i,j+1)**2.d0/cpp(i,j+1)))*0.5d0&
                /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j+1))))
           west = &
                ((D_t(i,j)*(ctt(i,j)-1.d0/A**2.d0*G(i,j)**2.d0)&
                +(D_p(i,j)-D_t(i,j))*cpt(i,j)**2.d0/cpp(i,j))&
                +(D_t(i,j-1)*(ctt(i,j-1)-1.d0/A**2.d0*G(i,j-1)**2.d0)&
                +(D_p(i,j-1)-D_t(i,j-1))*cpt(i,j-1)**2.d0/cpp(i,j-1)))*0.5d0&
                /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j-1))))

           !set flux to zero on wall
           east = east * (1.D0-Zones(k)%masks%chi2(i,j+1)*(1.d0-Zones(k)%masks%chi2(i,j)))
           west = west * (1.D0-(1.d0-Zones(k)%masks%chi2(i,j+1))*Zones(k)%masks%chi2(i,j))
           west = west * (1.D0-Zones(k)%masks%chi2(i,j-1)*(1.d0-Zones(k)%masks%chi2(i,j)))
           east = east * (1.D0-(1.d0-Zones(k)%masks%chi2(i,j-1))*Zones(k)%masks%chi2(i,j))

           select case(FIELD)
           case(DENSITY_FIELD)
              zones(k)%species(n_ion)%implicit_coefs%east_density(i,j)=east
              zones(k)%species(n_ion)%implicit_coefs%west_density(i,j)=west
           case(VELOCITY_FIELD)
              zones(k)%species(n_ion)%implicit_coefs%east_velocity(i,j)=east
              zones(k)%species(n_ion)%implicit_coefs%west_velocity(i,j)=west
           case(TEMPERATURE_FIELD)
              zones(k)%species(n_ion)%implicit_coefs%east_temperature(i,j)=east
              zones(k)%species(n_ion)%implicit_coefs%west_temperature(i,j)=west
           case(K_FIELD)
              zones(k)%kepsilon(1)%implicit_coefs%east_k(i,j)=east
              zones(k)%kepsilon(1)%implicit_coefs%west_k(i,j)=west
           case(EPSILON_FIELD)
              zones(k)%kepsilon(1)%implicit_coefs%east_epsilon(i,j)=east
              zones(k)%kepsilon(1)%implicit_coefs%west_epsilon(i,j)=west
           end select
        end do
     end do
     deallocate(ctt,cpt,cpp,G,D_p,D_t)
  end do
end subroutine compute_implicit_tridiag_coefficients_FIELD_ke
