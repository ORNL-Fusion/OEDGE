subroutine MD_broadcast_vpinch()
  use all_variables, only : zones, global_parameters
  implicit none
  integer*4 :: k,n
  integer*4 :: Nx,Nz,Nx_N,Nz_N
  integer*4 :: North,South,East,West
  do k=1,global_parameters%N_zones
     do n=0,global_parameters%N_ions
        Nx=zones(k)%mesh%Nx
        Nz=zones(k)%mesh%Nz
        !North boundary
        North=zones(k)%Neighbors(1)
        if(North.lt.0) then
           zones(k)%species(n)%transport_perp%v_pinch(Nx+1,:)=zones(k)%species(n)%transport_perp%v_pinch(Nx,:)
        else
           zones(k)%species(n)%transport_perp%v_pinch(Nx+1,:)=zones(North)%species(n)%transport_perp%v_pinch(1,:)
        end if
        !South boundary
        South=zones(k)%Neighbors(2)
        if(South.lt.0) then
           zones(k)%species(n)%transport_perp%v_pinch(0,:)=zones(k)%species(n)%transport_perp%v_pinch(1,:)
        else
           Nx_N=zones(South)%mesh%Nx
           zones(k)%species(n)%transport_perp%v_pinch(0,:)=zones(South)%species(n)%transport_perp%v_pinch(Nx_N,:)
        end if
        !East boundary
        East=zones(k)%Neighbors(3)
        if(East.lt.0) then
           zones(k)%species(n)%transport_perp%v_pinch(:,Nz+1)=zones(k)%species(n)%transport_perp%v_pinch(:,Nz)
        else
           zones(k)%species(n)%transport_perp%v_pinch(:,Nz+1)=zones(East)%species(n)%transport_perp%v_pinch(:,1)
        end if
        !West boundary
        West=zones(k)%Neighbors(4)
        if(West.lt.0) then
           zones(k)%species(n)%transport_perp%v_pinch(:,0)=zones(k)%species(n)%transport_perp%v_pinch(:,1)
        else
           Nz_N=zones(West)%mesh%Nz
           zones(k)%species(n)%transport_perp%v_pinch(:,0)=zones(West)%species(n)%transport_perp%v_pinch(:,Nz_N)
        end if
     end do
  end do
end subroutine MD_broadcast_vpinch
