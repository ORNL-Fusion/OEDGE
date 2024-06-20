subroutine vorticity_continuity_penalisation(zone,STEP)
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        if(zone%masks%chi1(i,j).eq.1) then
           zone%electric_fields(STEP)%vorticity(i,j)=zone%electric_fields(STEP)%vorticity(i,j-1)
        end if
        if(zone%masks%chi3(i,j).eq.1) then
           zone%electric_fields(STEP)%vorticity(i,j)=zone%electric_fields(STEP)%vorticity(i,j+1)
        end if
        if(zone%masks%chi5(i,j).eq.1) then
           zone%electric_fields(STEP)%vorticity(i,j)=zone%electric_fields(STEP)%vorticity(i+1,j)
        end if
        if(zone%masks%chi6(i,j).eq.1) then
           zone%electric_fields(STEP)%vorticity(i,j)=zone%electric_fields(STEP)%vorticity(i-1,j)
        end if
     end do
  end do
end subroutine vorticity_continuity_penalisation
