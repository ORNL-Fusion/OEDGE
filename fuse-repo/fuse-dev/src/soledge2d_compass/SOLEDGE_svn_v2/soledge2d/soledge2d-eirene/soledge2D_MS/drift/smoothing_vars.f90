module Msmoothing_vars

  use Mlist
  use Mpastix_solve

  implicit none

  integer*4 :: mat_smoothing_nnz    ! number of non-zero terms in neutral matrix
  type(cell), pointer :: smoothing_mat => null()
  Type(cell), pointer :: sm_ptr
  Type(cell), pointer :: ptr_prev
  Type(CSC) :: CSC_smoothing

  real*8 :: delta_r !m (smallest mesh size)
  real*8 :: smoothing_diffusivity ! M2/s

end module Msmoothing_vars
