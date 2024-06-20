subroutine compute_implicit_tridiag_coefficients()
  use all_variables, only : global_parameters
  use Mdefinitions
  implicit none
  integer*4 :: n
  do n=1,global_parameters%n_ions
     call compute_implicit_tridiag_coefficients_FIELD(n,DENSITY_FIELD)
     call compute_implicit_tridiag_coefficients_FIELD(n,VELOCITY_FIELD)
  end do
  do n=0,global_parameters%n_ions
     call compute_implicit_tridiag_coefficients_FIELD(n,TEMPERATURE_FIELD)
  end do
  call compute_implicit_tridiag_coefficients_FIELD(n,K_FIELD)
  call compute_implicit_tridiag_coefficients_FIELD(n,EPSILON_FIELD)
end subroutine compute_implicit_tridiag_coefficients
