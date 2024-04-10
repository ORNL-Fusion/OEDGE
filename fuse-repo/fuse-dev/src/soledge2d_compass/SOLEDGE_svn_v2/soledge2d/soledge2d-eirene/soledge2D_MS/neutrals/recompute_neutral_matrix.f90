subroutine recompute_neutral_matrix()
  use all_variables, only : global_parameters, zones
  use Mneutral_vars
  implicit none

  call compute_laplacian()
  call add_neutral_matrix_BC()
  call fill_implicit_pastix(CSC_neutral,mat_neutral_nnz,neutral_mat)
  call free_all(neutral_mat)
  call analyze_implicit_pastix(CSC_neutral)
end subroutine recompute_neutral_matrix
