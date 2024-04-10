subroutine compute_corner_W(zone,STEP)
  use all_variables, only : global_parameters,zones
  use MZones
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: i,j,k,Nx,Nz
  k=zone%number
  integer*4 :: North,South,East,West
  zone%electric_fields(STEP)%cornersPhi

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
             zones(North)%electric_fiels(STEP)%phi(1,Nz+1)&
             +zones(East)%electric_fiels(STEP)%phi(Nx+1,1))
     else
        zone%electric_fields(STEP)%cornersPhi(1,2,1)=(&
             zones(North)%electric_fiels(STEP)%phi(1,Nz+1))
     end if
  else
     if(East.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(1,2,1)=(&
             zones(East)%electric_fiels(STEP)%phi(Nx+1,1))
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
             zones(North)%electric_fiels(STEP)%phi(1,0)&
             +zones(West)%electric_fiels(STEP)%phi(Nx+1,Zones(West)%mesh%Nz))
     else
        zone%electric_fields(STEP)%cornersPhi(1,1,1)=(&
             zones(North)%electric_fiels(STEP)%phi(1,0))
     end if
  else
     if(West.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(1,1,1)=(&
             zones(West)%electric_fiels(STEP)%phi(Nx+1,Zones(West)%mesh%Nz))
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
             zones(South)%electric_fiels(STEP)%phi(Zones(South)%mesh%Nx,Nz+1)&
             +zones(East)%electric_fiels(STEP)%phi(0,1))
     else
        zone%electric_fields(STEP)%cornersPhi(2,2,1)=(&
             zones(South)%electric_fiels(STEP)%phi(Zones(South)%mesh%Nx,Nz+1))
     end if
  else
     if(East.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(2,2,1)=(&
             zones(East)%electric_fiels(STEP)%phi(0,1))
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
             zones(South)%electric_fiels(STEP)%phi(Zones(South)%mesh%Nx,0)&
             +zones(West)%electric_fiels(STEP)%phi(0,Zones(West)%mesh%Nz))
     else
        zone%electric_fields(STEP)%cornersPhi(2,1,1)=(&
             zones(South)%electric_fiels(STEP)%phi(Zones(South)%mesh%Nx,0))
     end if
  else
     if(West.gt.0) then
        zone%electric_fields(STEP)%cornersPhi(2,1,1)=(&
             zones(West)%electric_fiels(STEP)%phi(0,Zones(West)%mesh%Nz))
     else !real corner
        zone%electric_fields(STEP)%cornersPhi(2,1,1)=0.5d0*(&
             zone%electric_fields(STEP)%phi(1,0)&
             +zone%electric_fields(STEP)%phi(0,1))
     end if
  end if

end subroutine compute_corner_W
