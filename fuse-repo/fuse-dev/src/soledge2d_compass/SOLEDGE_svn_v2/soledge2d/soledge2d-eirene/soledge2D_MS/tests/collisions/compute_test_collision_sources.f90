subroutine compute_test_collision_sources(part)
  use test_var
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4,intent(in) :: part
  integer*4 :: k
  do k=1,global_parameters%N_zones
     call compute_test_collision_sourceN(k)
     call compute_test_collision_sourceG(k)
     select case(part)
     case(1)
        call compute_test_collision_sourceT1(k)
     case(2)
        call compute_test_collision_sourceT2(k)
     case(3)
        call compute_test_collision_sourceT3(k)
     end select
  end do
end subroutine compute_test_collision_sources
