module Mvorticity_vars

  use Mlist
  use Mpastix_solve

  implicit none

  integer*4 :: mat_vort_nnz    ! number of non-zero terms in vorticity matrix
  type(cell), pointer :: vorticity_mat => null()
  Type(cell), pointer :: vm_ptr
  Type(cell), pointer :: ptr_prev
  Type(CSC) :: CSC_vort

  integer*4,parameter :: ZERO_FLUX = 1
  integer*4,parameter :: DIRICHLET_TEST = 2
  integer*4,parameter :: LAMBDA_TE = 1
  integer*4,parameter :: BOHM_CURRENT = 2

  
end module Mvorticity_vars
