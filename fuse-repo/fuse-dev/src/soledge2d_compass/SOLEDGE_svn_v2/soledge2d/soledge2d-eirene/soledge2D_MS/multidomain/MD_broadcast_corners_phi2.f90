subroutine MD_broadcast_corners_phi2(zone,STEP)
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

  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  
  !NE
  zone%electric_fields(STEP)%cornersPhi(1,2,1)=zone%electric_fields(STEP)%phi(Nx+1,Nz+1)
  zone%electric_fields(STEP)%cornersPhi(1,2,2)=zone%electric_fields(STEP)%phi(Nx+1,Nz+1)

  !NW
  zone%electric_fields(STEP)%cornersPhi(1,1,1)=zone%electric_fields(STEP)%phi(Nx+1,0)
  zone%electric_fields(STEP)%cornersPhi(1,1,2)=zone%electric_fields(STEP)%phi(Nx+1,0)

  !SE
  zone%electric_fields(STEP)%cornersPhi(2,2,1)=zone%electric_fields(STEP)%phi(0,Nz+1)
  zone%electric_fields(STEP)%cornersPhi(2,2,2)=zone%electric_fields(STEP)%phi(0,Nz+1)

  !SW
  zone%electric_fields(STEP)%cornersPhi(2,1,1)=zone%electric_fields(STEP)%phi(0,0)
  zone%electric_fields(STEP)%cornersPhi(2,1,2)=zone%electric_fields(STEP)%phi(0,0)

  
end subroutine MD_broadcast_corners_phi2
