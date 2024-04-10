subroutine MD_broadcast_kepsilon(zone,STEP)
  use all_variables, only : zones
  use MZone
  use MDefinitions
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: Nx,Nz,Nx_N,Nz_N
  integer*4 :: North,South,East,West
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  North=zone%Neighbors(1)
  select case(North)
  case(-1)
     !     call set_core_boundary_conditions_north_ke(zone,STEP)
     call set_continuity_boundary_conditions_north_ke(zone,STEP)
  case(-2)
     !     call set_wall_boundary_conditions_north_ke(zone,STEP)
     call set_continuity_boundary_conditions_north_ke(zone,STEP)
  case(-4)
     call set_continuity_boundary_conditions_north_ke(zone,STEP)
  case default
     call gather_neighboring_values_ke(zone,N_NORTH,STEP)
  end select
  South=zone%Neighbors(2)
  select case(South)
  case(-1)
     !     call set_core_boundary_conditions_south_ke(zone,STEP)
     call set_continuity_boundary_conditions_south_ke(zone,STEP)
  case(-2)
     !     call set_wall_boundary_conditions_south_ke(zone,STEP)
     call set_continuity_boundary_conditions_south_ke(zone,STEP)
  case(-4)
     call set_continuity_boundary_conditions_south_ke(zone,STEP)
  case default
     call gather_neighboring_values_ke(zone,N_SOUTH,STEP)
  end select
  East=zone%Neighbors(3)
  select case(East)
  case(-3)
!     call set_wall_boundary_conditions_east_ke(zone,STEP)
  case default
     call gather_neighboring_values_ke(zone,N_EAST,STEP)
  end select
  West=zone%Neighbors(4)
  select case(West)
  case(-3)
!     call set_wall_boundary_conditions_west_ke(zone,STEP)
  case default
     call gather_neighboring_values_ke(zone,N_WEST,STEP)
  end select
end subroutine MD_broadcast_kepsilon


