subroutine MD_broadcast_corners_ke(zone,STEP)
#include "compile_opt.inc"
  use all_variables, only : zones,global_parameters
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n
  integer*4 :: North,South,East,West
  integer*4 :: Nx_N,Nz_N
  North=zone%Neighbors(1)
  if(North.gt.0) then
     Nz_N=zones(North)%mesh%Nz
     zone%kepsilon(STEP)%cornersK(1,1,1)=zones(North)%kepsilon(STEP)%k(1,0)
     zone%kepsilon(STEP)%cornersK(1,2,1)=zones(North)%kepsilon(STEP)%k(1,Nz_N+1)
     zone%kepsilon(STEP)%cornersEpsilon(1,1,1)=zones(North)%kepsilon(STEP)%epsilon(1,0)
     zone%kepsilon(STEP)%cornersEpsilon(1,2,1)=zones(North)%kepsilon(STEP)%epsilon(1,Nz_N+1)
  end if
  South=zone%Neighbors(2)
  if(South.gt.0) then
     Nz_N=zones(South)%mesh%Nz
     Nx_N=zones(South)%mesh%Nx
     zone%kepsilon(STEP)%cornersK(2,1,1)=zones(South)%kepsilon(STEP)%k(Nx_N,0)
     zone%kepsilon(STEP)%cornersK(2,2,1)=zones(South)%kepsilon(STEP)%k(Nx_N,Nz_N+1)
     zone%kepsilon(STEP)%cornersEpsilon(2,1,1)=zones(South)%kepsilon(STEP)%epsilon(Nx_N,0)
     zone%kepsilon(STEP)%cornersEpsilon(2,2,1)=zones(South)%kepsilon(STEP)%epsilon(Nx_N,Nz_N+1)

  end if
  East=zone%Neighbors(3)
  if(East.gt.0) then
     Nx_N=zones(East)%mesh%Nx
     Nz_N=zones(East)%mesh%Nz
     zone%kepsilon(STEP)%cornersK(1,2,2)=zones(East)%kepsilon(STEP)%k(Nx_N+1,1)
     zone%kepsilon(STEP)%cornersK(2,2,2)=zones(East)%kepsilon(STEP)%k(0,1)
     zone%kepsilon(STEP)%cornersEpsilon(1,2,2)=zones(East)%kepsilon(STEP)%epsilon(Nx_N+1,1)
     zone%kepsilon(STEP)%cornersEpsilon(2,2,2)=zones(East)%kepsilon(STEP)%epsilon(0,1)
  end if
  West=zone%Neighbors(4)
  if(West.gt.0) then
     Nx_N=zones(West)%mesh%Nx
     Nz_N=zones(West)%mesh%Nz
     zone%kepsilon(STEP)%cornersK(1,1,2)=zones(West)%kepsilon(STEP)%k(Nx_N+1,Nz_N)
     zone%kepsilon(STEP)%cornersK(2,1,2)=zones(West)%kepsilon(STEP)%k(0,Nz_N)
     zone%kepsilon(STEP)%cornersEpsilon(1,1,2)=zones(West)%kepsilon(STEP)%epsilon(Nx_N+1,Nz_N)
     zone%kepsilon(STEP)%cornersEpsilon(2,1,2)=zones(West)%kepsilon(STEP)%epsilon(0,Nz_N)
  end if
  if(North.lt.0) then
     zone%kepsilon(STEP)%cornersK(1,1,1)=zone%kepsilon(STEP)%cornersK(1,1,2)
     zone%kepsilon(STEP)%cornersK(1,2,1)=zone%kepsilon(STEP)%cornersK(1,2,2)
     zone%kepsilon(STEP)%cornersEpsilon(1,1,1)=zone%kepsilon(STEP)%cornersEpsilon(1,1,2)
     zone%kepsilon(STEP)%cornersEpsilon(1,2,1)=zone%kepsilon(STEP)%cornersEpsilon(1,2,2)
  end if
  if(South.lt.0) then
     zone%kepsilon(STEP)%cornersK(2,1,1)=zone%kepsilon(STEP)%cornersK(2,1,2)
     zone%kepsilon(STEP)%cornersK(2,2,1)=zone%kepsilon(STEP)%cornersK(2,2,2)
     zone%kepsilon(STEP)%cornersEpsilon(2,1,1)=zone%kepsilon(STEP)%cornersEpsilon(2,1,2)
     zone%kepsilon(STEP)%cornersEpsilon(2,2,1)=zone%kepsilon(STEP)%cornersEpsilon(2,2,2)
  end if
  if(East.lt.0) then
     zone%kepsilon(STEP)%cornersK(1,2,2)=zone%kepsilon(STEP)%cornersK(1,2,1)
     zone%kepsilon(STEP)%cornersK(2,2,2)=zone%kepsilon(STEP)%cornersK(2,2,1)
     zone%kepsilon(STEP)%cornersEpsilon(1,2,2)=zone%kepsilon(STEP)%cornersEpsilon(1,2,1)
     zone%kepsilon(STEP)%cornersEpsilon(2,2,2)=zone%kepsilon(STEP)%cornersEpsilon(2,2,1)
  end if
  if(West.lt.0) then
     zone%kepsilon(STEP)%cornersK(1,1,2)=zone%kepsilon(STEP)%cornersK(1,1,1)
     zone%kepsilon(STEP)%cornersK(2,1,2)=zone%kepsilon(STEP)%cornersK(2,1,1)
     zone%kepsilon(STEP)%cornersEpsilon(1,1,2)=zone%kepsilon(STEP)%cornersEpsilon(1,1,1)
     zone%kepsilon(STEP)%cornersEpsilon(2,1,2)=zone%kepsilon(STEP)%cornersEpsilon(2,1,1)
  end if
end subroutine MD_broadcast_corners_ke
