subroutine MD_broadcast_corners_W(zone,STEP)
  use all_variables, only : zones,global_parameters
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n
  integer*4 :: North,South,East,West
  integer*4 :: Nx_N,Nz_N
  do n=0,global_parameters%N_ions
     North=zone%Neighbors(1)
     if(North.gt.0) then
        Nz_N=zones(North)%mesh%Nz
        zone%electric_fields(STEP)%cornersW(1,1,1)=zones(North)%electric_fields(STEP)%vorticity(1,0)
        zone%electric_fields(STEP)%cornersW(1,2,1)=zones(North)%electric_fields(STEP)%vorticity(1,Nz_N+1)
     end if
     South=zone%Neighbors(2)
     if(South.gt.0) then
        Nz_N=zones(South)%mesh%Nz
        Nx_N=zones(South)%mesh%Nx
        zone%electric_fields(STEP)%cornersW(2,1,1)=zones(South)%electric_fields(STEP)%vorticity(Nx_N,0)
        zone%electric_fields(STEP)%cornersW(2,2,1)=zones(South)%electric_fields(STEP)%vorticity(Nx_N,Nz_N+1)
     end if
     East=zone%Neighbors(3)
     if(East.gt.0) then
        Nx_N=zones(East)%mesh%Nx
        Nz_N=zones(East)%mesh%Nz
        zone%electric_fields(STEP)%cornersW(1,2,2)=zones(East)%electric_fields(STEP)%vorticity(Nx_N+1,1)
        zone%electric_fields(STEP)%cornersW(2,2,2)=zones(East)%electric_fields(STEP)%vorticity(0,1)
     end if
     West=zone%Neighbors(4)
     if(West.gt.0) then
        Nx_N=zones(West)%mesh%Nx
        Nz_N=zones(West)%mesh%Nz
        zone%electric_fields(STEP)%cornersW(1,1,2)=zones(West)%electric_fields(STEP)%vorticity(Nx_N+1,Nz_N)
        zone%electric_fields(STEP)%cornersW(2,1,2)=zones(West)%electric_fields(STEP)%vorticity(0,Nz_N)
     end if     
     if(North.lt.0) then
        zone%electric_fields(STEP)%cornersW(1,1,1)=zone%electric_fields(STEP)%cornersW(1,1,2)
        zone%electric_fields(STEP)%cornersW(1,2,1)=zone%electric_fields(STEP)%cornersW(1,2,2)
     end if
     if(South.lt.0) then
        zone%electric_fields(STEP)%cornersW(2,1,1)=zone%electric_fields(STEP)%cornersW(2,1,2)
        zone%electric_fields(STEP)%cornersW(2,2,1)=zone%electric_fields(STEP)%cornersW(2,2,2)
     end if
     if(East.lt.0) then
        zone%electric_fields(STEP)%cornersW(1,2,2)=zone%electric_fields(STEP)%cornersW(1,2,1)
        zone%electric_fields(STEP)%cornersW(2,2,2)=zone%electric_fields(STEP)%cornersW(2,2,1)
     end if
     if(West.lt.0) then
        zone%electric_fields(STEP)%cornersW(1,1,2)=zone%electric_fields(STEP)%cornersW(1,1,1)
        zone%electric_fields(STEP)%cornersW(2,1,2)=zone%electric_fields(STEP)%cornersW(2,1,1)
     end if     
 end do
end subroutine MD_broadcast_corners_W
