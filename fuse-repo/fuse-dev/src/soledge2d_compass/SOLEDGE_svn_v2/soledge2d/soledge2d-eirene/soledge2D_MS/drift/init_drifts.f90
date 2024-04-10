subroutine init_drifts()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: k,i,j,n
!!$  do k=1,global_parameters%N_zones
!!$     call find_ExB_OK_points(zones(k))
!!$  end do
!!$  call broadcast_OK_points()
!!$  do k=1,global_parameters%N_zones
!!$     call find_ref_point_for_NOK_points(zones(k))
!!$     call find_ref_point_for_NOK_points_step2(zones(k))
!!$     call detect_problem_with_ref_points(zones(k))
!!$  end do
  call init_smoothing()
  call compute_smooth_mat()
end subroutine init_drifts
