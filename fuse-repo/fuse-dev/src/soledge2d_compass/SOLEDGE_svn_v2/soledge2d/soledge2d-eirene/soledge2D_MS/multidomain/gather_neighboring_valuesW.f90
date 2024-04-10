subroutine gather_neighboring_valuesW(zone,SIDE,STEP)
  use all_variables, only : zones, global_parameters
  use Mzone
  implicit none
  integer*4 :: n
  Type(TZone) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: SIDE
  integer*4 :: Neighbor
  integer*4 :: Nx,Nz,Nx_N,Nz_N
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  select case(SIDE)
  case(N_NORTH)
     Neighbor=zone%Neighbors(1)
     zone%electric_fields(STEP)%vorticity(Nx+1,1:Nz)=zones(Neighbor)%electric_fields(STEP)%vorticity(1,1:Nz)
  case(N_SOUTH)
     Neighbor=zone%Neighbors(2)
     Nx_N=zones(Neighbor)%mesh%Nx
     zone%electric_fields(STEP)%vorticity(0,1:Nz)=zones(Neighbor)%electric_fields(STEP)%vorticity(Nx_N,1:Nz)
  case(N_EAST)
     Neighbor=zone%Neighbors(3)
     zone%electric_fields(STEP)%vorticity(1:Nx,Nz+1)=zones(Neighbor)%electric_fields(STEP)%vorticity(1:Nx,1)
  case(N_WEST)
     Neighbor=zone%Neighbors(4)
     Nz_N=zones(Neighbor)%mesh%Nz
     zone%electric_fields(STEP)%vorticity(1:Nx,0)=zones(Neighbor)%electric_fields(STEP)%vorticity(1:Nx,Nz_N)
  end select
end subroutine gather_neighboring_valuesW


