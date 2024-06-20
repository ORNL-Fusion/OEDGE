subroutine update_all_fields()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: k
  do k=1,global_parameters%N_zones
     call update_fields(zones(k))
  end do
end subroutine update_all_fields
