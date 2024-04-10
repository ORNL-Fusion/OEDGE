subroutine MD_broadcast_corners_phi(zone,STEP)
  use all_variables, only : zones,global_parameters
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n
  integer*4 :: North,South,East,West
  integer*4 :: Nx_N,Nz_N
  integer*4 :: i,j,k,Nx,Nz
  k=zone%number

  North=zone%Neighbors(1)
  South=zone%Neighbors(2)
  East=zone%Neighbors(3)
  West=zone%Neighbors(4)

  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz

  !treating corner NE
  if(North.gt.0) then
     if(East.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(1,2,1)=0.5D0*(&
             zones(North)%electric_fields(STEP)%phi(1,Nz+1)&
             +zones(East)%electric_fields(STEP)%phi(Nx+1,1))
     else
        zone%electric_fields(STEP)%cornersPhi(1,2,1)=(&
             zones(North)%electric_fields(STEP)%phi(1,Nz+1))
     end if
  else
     if(East.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(1,2,1)=(&
             zones(East)%electric_fields(STEP)%phi(Nx+1,1))
     else !real corner
        zone%electric_fields(STEP)%cornersPhi(1,2,1)=0.5D0*(&
             zone%electric_fields(STEP)%phi(Nx+1,Nz)&
             +zone%electric_fields(STEP)%phi(Nx,Nz+1))
     end if
  end if

  !treating corner NW
  if(North.gt.0) then
     if(West.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(1,1,1)=0.5D0*(&
             zones(North)%electric_fields(STEP)%phi(1,0)&
             +zones(West)%electric_fields(STEP)%phi(Nx+1,Zones(West)%mesh%Nz))
     else
        zone%electric_fields(STEP)%cornersPhi(1,1,1)=(&
             zones(North)%electric_fields(STEP)%phi(1,0))
     end if
  else
     if(West.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(1,1,1)=(&
             zones(West)%electric_fields(STEP)%phi(Nx+1,Zones(West)%mesh%Nz))
     else !real corner
        zone%electric_fields(STEP)%cornersPhi(1,1,1)=0.5D0*&
             (zone%electric_fields(STEP)%phi(Nx+1,1)&
             +zone%electric_fields(STEP)%phi(Nx,0))
     end if
  end if

  !treating corner SE
  if(South.gt.0) then
     if(East.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(2,2,1)=0.5D0*(&
             zones(South)%electric_fields(STEP)%phi(Zones(South)%mesh%Nx,Nz+1)&
             +zones(East)%electric_fields(STEP)%phi(0,1))
     else
        zone%electric_fields(STEP)%cornersPhi(2,2,1)=(&
             zones(South)%electric_fields(STEP)%phi(Zones(South)%mesh%Nx,Nz+1))
     end if
  else
     if(East.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(2,2,1)=(&
             zones(East)%electric_fields(STEP)%phi(0,1))
     else !real corner
        zone%electric_fields(STEP)%cornersPhi(2,2,1)=0.5D0*&
             (zone%electric_fields(STEP)%phi(0,Nz)&
             +zone%electric_fields(STEP)%phi(1,Nz+1))
     end if
  end if

  !treating corner SW
  if(South.gt.0) then
     if(West.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(2,1,1)=0.5D0*(&
             zones(South)%electric_fields(STEP)%phi(Zones(South)%mesh%Nx,0)&
             +zones(West)%electric_fields(STEP)%phi(0,Zones(West)%mesh%Nz))
     else
        zone%electric_fields(STEP)%cornersPhi(2,1,1)=(&
             zones(South)%electric_fields(STEP)%phi(Zones(South)%mesh%Nx,0))
     end if
  else
     if(West.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(2,1,1)=(&
             zones(West)%electric_fields(STEP)%phi(0,Zones(West)%mesh%Nz))
     else !real corner
        zone%electric_fields(STEP)%cornersPhi(2,1,1)=0.5d0*(&
             zone%electric_fields(STEP)%phi(1,0)&
             +zone%electric_fields(STEP)%phi(0,1))
     end if
  end if

  zone%electric_fields(STEP)%cornersPhi(2,1,2)=zone%electric_fields(STEP)%cornersPhi(2,1,1)
  zone%electric_fields(STEP)%cornersPhi(1,1,2)=zone%electric_fields(STEP)%cornersPhi(1,1,1)
  zone%electric_fields(STEP)%cornersPhi(1,2,2)=zone%electric_fields(STEP)%cornersPhi(1,2,1)
  zone%electric_fields(STEP)%cornersPhi(2,2,2)=zone%electric_fields(STEP)%cornersPhi(2,2,1)
  
!!$
!!$  
!!$  do n=0,global_parameters%N_ions
!!$     North=zone%Neighbors(1)
!!$     if(North.gt.0) then
!!$        Nz_N=zones(North)%mesh%Nz
!!$        zone%electric_fields(STEP)%cornersPhi(1,1,1)=zones(North)%electric_fields(STEP)%phi(1,0)
!!$        zone%electric_fields(STEP)%cornersPhi(1,2,1)=zones(North)%electric_fields(STEP)%phi(1,Nz_N+1)
!!$     end if
!!$     South=zone%Neighbors(2)
!!$     if(South.gt.0) then
!!$        Nz_N=zones(South)%mesh%Nz
!!$        Nx_N=zones(South)%mesh%Nx
!!$        zone%electric_fields(STEP)%cornersPhi(2,1,1)=zones(South)%electric_fields(STEP)%phi(Nx_N,0)
!!$        zone%electric_fields(STEP)%cornersPhi(2,2,1)=zones(South)%electric_fields(STEP)%phi(Nx_N,Nz_N+1)
!!$     end if
!!$     East=zone%Neighbors(3)
!!$     if(East.gt.0) then
!!$        Nx_N=zones(East)%mesh%Nx
!!$        Nz_N=zones(East)%mesh%Nz
!!$        zone%electric_fields(STEP)%cornersPhi(1,2,2)=zones(East)%electric_fields(STEP)%phi(Nx_N+1,1)
!!$        zone%electric_fields(STEP)%cornersPhi(2,2,2)=zones(East)%electric_fields(STEP)%phi(0,1)
!!$     end if
!!$     West=zone%Neighbors(4)
!!$     if(West.gt.0) then
!!$        Nx_N=zones(West)%mesh%Nx
!!$        Nz_N=zones(West)%mesh%Nz
!!$        zone%electric_fields(STEP)%cornersPhi(1,1,2)=zones(West)%electric_fields(STEP)%phi(Nx_N+1,Nz_N)
!!$        zone%electric_fields(STEP)%cornersPhi(2,1,2)=zones(West)%electric_fields(STEP)%phi(0,Nz_N)
!!$     end if     
!!$     zone%electric_fields(STEP)%cornersPhi(1,1,2)=zone%electric_fields(STEP)%cornersPhi(1,1,1)
!!$     zone%electric_fields(STEP)%cornersPhi(1,2,2)=zone%electric_fields(STEP)%cornersPhi(1,2,1)
!!$     zone%electric_fields(STEP)%cornersPhi(2,1,2)=zone%electric_fields(STEP)%cornersPhi(2,1,1)
!!$     zone%electric_fields(STEP)%cornersPhi(2,2,2)=zone%electric_fields(STEP)%cornersPhi(2,2,1)
!!$     if(North.lt.0) then
!!$        zone%electric_fields(STEP)%cornersPhi(1,1,1)=zone%electric_fields(STEP)%cornersPhi(1,1,2)
!!$        zone%electric_fields(STEP)%cornersPhi(1,2,1)=zone%electric_fields(STEP)%cornersPhi(1,2,2)
!!$     end if
!!$     if(South.lt.0) then
!!$        zone%electric_fields(STEP)%cornersPhi(2,1,1)=zone%electric_fields(STEP)%cornersPhi(2,1,2)
!!$        zone%electric_fields(STEP)%cornersPhi(2,2,1)=zone%electric_fields(STEP)%cornersPhi(2,2,2)
!!$     end if
!!$     if(East.lt.0) then
!!$        zone%electric_fields(STEP)%cornersPhi(1,2,2)=zone%electric_fields(STEP)%cornersPhi(1,2,1)
!!$        zone%electric_fields(STEP)%cornersPhi(2,2,2)=zone%electric_fields(STEP)%cornersPhi(2,2,1)
!!$     end if
!!$     if(West.lt.0) then
!!$        zone%electric_fields(STEP)%cornersPhi(1,1,2)=zone%electric_fields(STEP)%cornersPhi(1,1,1)
!!$        zone%electric_fields(STEP)%cornersPhi(2,1,2)=zone%electric_fields(STEP)%cornersPhi(2,1,1)
!!$     end if     
!!$ end do

end subroutine MD_broadcast_corners_phi
