subroutine recompute_diffusion_coefficients(zone)
  use all_variables, only : global_parameters, reference_parameters, kepsilon_param
  use Mzone
  use MPhysics
  implicit none
  Type(TZone), intent(inout) :: zone
  integer*4 :: n  
  zone%species(0)%transport_perp%chi_p=zone%kepsilon(1)%mu_t/kepsilon_param%sigma_T
  do n=1,global_parameters%N_ions
     zone%species(n)%transport_perp%D_p=zone%kepsilon(1)%mu_t/kepsilon_param%sigma_n
     zone%species(n)%transport_perp%nu_p=zone%kepsilon(1)%mu_t/kepsilon_param%sigma_v
     zone%species(n)%transport_perp%chi_p=zone%kepsilon(1)%mu_t/kepsilon_param%sigma_T
  end do
end subroutine recompute_diffusion_coefficients
