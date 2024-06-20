module Mneutral_vars

  use Mlist
  use Mpastix_solve

  implicit none

  integer*4 :: mat_neutral_nnz    ! number of non-zero terms in neutral matrix
  type(cell), pointer :: neutral_mat => null()
  Type(cell), pointer :: nm_ptr
  Type(cell), pointer :: ptr_prev
  Type(CSC) :: CSC_neutral

  real*8 :: FN_diffusivity ! M2/s
  real*8 :: FN_recycling_coefficient 

end module Mneutral_vars
