subroutine MD_broadcast_mesh()
  use all_variables, only : zones, global_parameters
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
        zones(k)%mesh%x(Nx+1,:)=zones(k)%mesh%xmax
        zones(k)%mesh%z(Nx+1,1:Nz)=zones(k)%mesh%z(Nx,1:Nz)
        zones(k)%mesh%z_minus_1half(Nx+1,1:Nz)=zones(k)%mesh%z_minus_1half(Nx,1:Nz)
        zones(k)%mesh%z_plus_1half(Nx+1,1:Nz)=zones(k)%mesh%z_plus_1half(Nx,1:Nz)
     else
        zones(k)%mesh%x(Nx+1,1:Nz)=zones(k)%mesh%xmax+zones(North)%mesh%x(1,1:Nz)-zones(North)%mesh%xmin
        zones(k)%mesh%z(Nx+1,1:Nz)=zones(North)%mesh%z(1,1:Nz)
        zones(k)%mesh%z_minus_1half(Nx+1,1:Nz)=zones(North)%mesh%z_minus_1half(1,1:Nz)
        zones(k)%mesh%z_plus_1half(Nx+1,1:Nz)=zones(North)%mesh%z_plus_1half(1,1:Nz)
        !corners
        zones(k)%mesh%x(Nx+1,0)=zones(k)%mesh%xmax+zones(North)%mesh%x(1,1)-zones(North)%mesh%xmin
        zones(k)%mesh%x(Nx+1,Nz+1)=zones(k)%mesh%xmax+zones(North)%mesh%x(1,Nz)-zones(North)%mesh%xmin
     end if
     !South boundary
     South=zones(k)%Neighbors(2)
     if(South.lt.0) then
        zones(k)%mesh%x(0,:)=zones(k)%mesh%xmin
        zones(k)%mesh%z(0,1:Nz)=zones(k)%mesh%z(1,1:Nz)
        zones(k)%mesh%z_minus_1half(0,1:Nz)=zones(k)%mesh%z_minus_1half(1,1:Nz)
        zones(k)%mesh%z_plus_1half(0,1:Nz)=zones(k)%mesh%z_plus_1half(1,1:Nz)
     else
        Nx_N=zones(South)%mesh%Nx
        zones(k)%mesh%x(0,1:Nz)=zones(k)%mesh%xmin+zones(South)%mesh%x(Nx_N,1:Nz)-zones(South)%mesh%xmax
        zones(k)%mesh%z(0,1:Nz)=zones(South)%mesh%z(Nx_N,1:Nz)
        zones(k)%mesh%z_minus_1half(0,1:Nz)=zones(South)%mesh%z_minus_1half(Nx_N,1:Nz)
        zones(k)%mesh%z_plus_1half(0,1:Nz)=zones(South)%mesh%z_plus_1half(Nx_N,1:Nz)
        !corners
        zones(k)%mesh%x(0,0)=zones(k)%mesh%xmin+zones(South)%mesh%x(Nx_N,1)-zones(South)%mesh%xmax
        zones(k)%mesh%x(0,Nz+1)=zones(k)%mesh%xmin+zones(South)%mesh%x(Nx_N,Nz)-zones(South)%mesh%xmax
     end if
     !East boundary
     East=zones(k)%Neighbors(3)
     if(East.lt.0) then
        zones(k)%mesh%z_plus_1half(1:Nx,Nz+1)=2.*zones(k)%mesh%z_plus_1half(1:Nx,Nz)&
             -zones(k)%mesh%z_plus_1half(1:Nx,Nz-1)
        zones(k)%mesh%z(1:Nx,Nz+1)=(3.*zones(k)%mesh%z_plus_1half(1:Nx,Nz)-zones(k)%mesh%z_plus_1half(1:Nx,Nz-1))/2.
        zones(k)%mesh%z_minus_1half(1:Nx,Nz+1)=zones(k)%mesh%z_plus_1half(1:Nx,Nz)
        zones(k)%mesh%x(1:Nx,Nz+1)=zones(k)%mesh%x(1:Nx,Nz)
        zones(k)%mesh%x_minus_1half(1:Nx,Nz+1)=zones(k)%mesh%x_minus_1half(1:Nx,Nz)
        zones(k)%mesh%x_plus_1half(1:Nx,Nz+1)=zones(k)%mesh%x_plus_1half(1:Nx,Nz)
     else
        zones(k)%mesh%z_plus_1half(1:Nx,Nz+1)=zones(k)%mesh%zmax+zones(East)%mesh%z_plus_1half(1:Nx,1)-zones(East)%mesh%zmin
        zones(k)%mesh%z(1:Nx,Nz+1)=zones(k)%mesh%zmax+zones(East)%mesh%z(1:Nx,1)-zones(East)%mesh%zmin
        zones(k)%mesh%z_minus_1half(1:Nx,Nz+1)=zones(k)%mesh%zmax+zones(East)%mesh%z_minus_1half(1:Nx,1)-zones(East)%mesh%zmin
        zones(k)%mesh%x(1:Nx,Nz+1)=zones(East)%mesh%x(1:Nx,1)
        zones(k)%mesh%x_minus_1half(1:Nx,Nz+1)=zones(East)%mesh%x_minus_1half(1:Nx,1)
        zones(k)%mesh%x_plus_1half(1:Nx,Nz+1)=zones(East)%mesh%x_plus_1half(1:Nx,1)
     end if
     !West boundary
     West=zones(k)%Neighbors(4)
     if(West.lt.0) then
        zones(k)%mesh%z_plus_1half(1:Nx,0)=zones(k)%mesh%z_minus_1half(1:Nx,1)
        zones(k)%mesh%z(1:Nx,0)=(3.*zones(k)%mesh%z_minus_1half(1:Nx,1)-zones(k)%mesh%z_minus_1half(1:Nx,2))/2.
        zones(k)%mesh%z_minus_1half(1:Nx,0)=2.*zones(k)%mesh%z_minus_1half(1:Nx,1)-zones(k)%mesh%z_minus_1half(1:Nx,2)
        zones(k)%mesh%x(1:Nx,0)=zones(k)%mesh%x(1:Nx,1)
        zones(k)%mesh%x_minus_1half(1:Nx,0)=zones(k)%mesh%x_minus_1half(1:Nx,1)
        zones(k)%mesh%x_plus_1half(1:Nx,0)=zones(k)%mesh%x_plus_1half(1:Nx,1)
     else
        Nz_N=zones(West)%mesh%Nz
        zones(k)%mesh%z_plus_1half(1:Nx,0)=zones(k)%mesh%zmin+zones(West)%mesh%z_plus_1half(1:Nx,Nz_N)-zones(West)%mesh%zmax
        zones(k)%mesh%z(1:Nx,0)=zones(k)%mesh%zmin+zones(West)%mesh%z(1:Nx,Nz_N)-zones(West)%mesh%zmax
        zones(k)%mesh%z_minus_1half(1:Nx,0)=zones(k)%mesh%zmin+zones(West)%mesh%z_minus_1half(1:Nx,Nz_N)-zones(West)%mesh%zmax
        zones(k)%mesh%x(1:Nx,0)=zones(West)%mesh%x(1:Nx,Nz_N)
        zones(k)%mesh%x_minus_1half(1:Nx,0)=zones(West)%mesh%x_minus_1half(1:Nx,Nz_N)
        zones(k)%mesh%x_plus_1half(1:Nx,0)=zones(West)%mesh%x_plus_1half(1:Nx,Nz_N)
     end if
  end do
end subroutine MD_broadcast_mesh
