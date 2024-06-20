subroutine set_perp_drift_on_core(zone,Fluxes,Nx,Nz)
  use MZone
  implicit none
  Type(Tzone),intent(inout) :: zone
  real*8 :: Fluxes(1:Nx,1:Nz,1:4)
  integer*4 :: i,j
  integer*4,intent(in) :: Nx,Nz
  !north
  if(zone%Neighbors(1).eq.-1) then
     do j=1,Nz
        Fluxes(Nx,j,1)=Fluxes(Nx,j,2)
     end do
  end if
  !south
  if(zone%Neighbors(2).eq.-1) then
     do j=1,Nz
        Fluxes(1,j,2)=Fluxes(1,j,1)
     end do
  end if
end subroutine set_perp_drift_on_core
