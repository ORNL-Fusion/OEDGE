subroutine prepare_penalisation_masks()
  use all_variables, only : global_parameters, zones
  implicit none
  call extra_chis()
  call MD_broadcast_masks()
  call compute_perp_masks()
  call save_masks()
  call find_pts_around_penwall()
end subroutine prepare_penalisation_masks
