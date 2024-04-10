subroutine set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
  use MZone
  implicit none
  Type(Tzone),intent(inout) :: zone
  real*8 :: Fluxes(1:Nx,1:Nz,1:4)
  integer*4 :: i,j
  integer*4,intent(in) :: Nx,Nz
  do i=1,Nx
     do j=1,Nz
        !zero radial flux at wall interface
!!$        Fluxes(i,j,1)=Fluxes(i,j,1)*(1.D0-zone%masks%chi2(i+1,j)*(1.d0-zone%masks%chi2(i,j)))
!!$        Fluxes(i,j,2)=Fluxes(i,j,2)*(1.D0-(1.d0-zone%masks%chi2(i-1,j))*zone%masks%chi2(i,j))
!!$        Fluxes(i,j,2)=Fluxes(i,j,2)*(1.D0-zone%masks%chi2(i-1,j)*(1.d0-zone%masks%chi2(i,j)))
!!$        Fluxes(i,j,1)=Fluxes(i,j,1)*(1.D0-(1.d0-zone%masks%chi2(i-1,j))*zone%masks%chi2(i,j))
!!$        Fluxes(i,j,3)=Fluxes(i,j,3)*(1.D0-zone%masks%chi2(i,j+1)*(1.d0-zone%masks%chi2(i,j)))
!!$        Fluxes(i,j,4)=Fluxes(i,j,4)*(1.D0-(1.d0-zone%masks%chi2(i,j+1))*zone%masks%chi2(i,j))
!!$        Fluxes(i,j,4)=Fluxes(i,j,4)*(1.D0-zone%masks%chi2(i,j-1)*(1.D0-zone%masks%chi2(i,j)))
!!$        Fluxes(i,j,3)=Fluxes(i,j,3)*(1.D0-(1.d0-zone%masks%chi2(i,j-1))*zone%masks%chi2(i,j))
        Fluxes(i,j,1)=Fluxes(i,j,1)*(1.D0-zone%masks%chi2(i+1,j))
        Fluxes(i,j,2)=Fluxes(i,j,2)*(1.D0-zone%masks%chi2(i-1,j))
        Fluxes(i,j,3)=Fluxes(i,j,3)*(1.D0-zone%masks%chi2(i,j+1))
        Fluxes(i,j,4)=Fluxes(i,j,4)*(1.D0-zone%masks%chi2(i,j-1))
     end do
  end do
  !north
  if(zone%Neighbors(1).eq.-2) then
     do j=1,Nz
        Fluxes(Nx,j,1)=0.D0
     end do
  end if
  !south
  if(zone%Neighbors(2).eq.-2) then
     do j=1,Nz
        Fluxes(1,j,2)=0.D0
     end do
  end if
  !east
  if(zone%Neighbors(3).eq.-3) then
     do i=1,Nx
        Fluxes(i,Nz,3)=0.D0
     end do
  end if
  !west
  if(zone%Neighbors(4).eq.-3) then
     do i=1,Nx
        Fluxes(i,1,4)=0.D0
     end do
  end if
end subroutine set_perp_flux_to_zero_on_the_wall
