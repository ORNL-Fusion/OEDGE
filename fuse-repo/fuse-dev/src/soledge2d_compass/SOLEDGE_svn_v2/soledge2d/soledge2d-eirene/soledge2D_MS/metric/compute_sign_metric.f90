subroutine compute_sign_metric()
  use all_variables, only : zones, global_variables, flags
  implicit none
  global_variables%sign_metric=sign(1.d0,sum(zones(1)%metric_coefficients%G(1:zones(1)%mesh%Nx,1:zones(1)%mesh%Nz)))
end subroutine compute_sign_metric
