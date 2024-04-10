subroutine broadcast_OK_points()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: k
  integer*4 :: Nx,Nz,Nx_N,Nz_N
  integer*4 :: North,South,East,West
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     !North boundary
     North=zones(k)%Neighbors(1)
     if(North.lt.0) then
        zones(k)%driftsExtra%OK_points(Nx+1,1:Nz)=zones(k)%driftsExtra%OK_points(Nx,1:Nz)
     else
        zones(k)%DriftsExtra%OK_points(Nx+1,1:Nz)=zones(North)%DriftsExtra%OK_points(1,1:Nz)
     end if
     !South boundary
     South=zones(k)%Neighbors(2)
     if(South.lt.0) then
        zones(k)%DriftsExtra%OK_points(0,1:Nz)=zones(k)%DriftsExtra%OK_points(1,1:Nz)
     else
        Nx_N=zones(South)%mesh%Nx
        zones(k)%DriftsExtra%OK_points(0,1:Nz)=zones(South)%DriftsExtra%OK_points(Nx_N,1:Nz)
     end if
     !East boundary
     East=zones(k)%Neighbors(3)
     if(East.lt.0) then
        zones(k)%DriftsExtra%OK_points(1:Nx,Nz+1)=zones(k)%DriftsExtra%OK_points(1:Nx,Nz)
     else
        zones(k)%DriftsExtra%OK_points(1:Nx,Nz+1)=zones(East)%DriftsExtra%OK_points(1:Nx,1)
     end if
     !West boundary
     West=zones(k)%Neighbors(4)
     if(West.lt.0) then
        zones(k)%DriftsExtra%OK_points(1:Nx,0)=zones(k)%DriftsExtra%OK_points(1:Nx,1)
     else
        Nz_N=zones(West)%mesh%Nz
        zones(k)%DriftsExtra%OK_points(1:Nx,0)=zones(West)%DriftsExtra%OK_points(1:Nx,Nz_N)
     end if
  end do
end subroutine broadcast_OK_points
