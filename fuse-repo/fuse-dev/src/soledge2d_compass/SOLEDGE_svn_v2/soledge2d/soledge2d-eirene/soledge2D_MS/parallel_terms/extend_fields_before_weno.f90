subroutine extend_fields_before_weno(k,n,Nx,Nz,density,Gamma,T,uEt,uBt)
  use all_variables, only : zones
  implicit none
  integer*4,intent(in) :: k,Nx,Nz,n
  real*8,intent(out) :: density(1:Nx,-1:Nz+2) 
  real*8,intent(out) :: Gamma(1:Nx,-1:Nz+2)
  real*8,intent(out) :: T(1:Nx,-1:Nz+2)
  real*8,intent(out) :: uEt(1:Nx,-1:Nz+2)
  real*8,intent(out) :: uBt(1:Nx,-1:Nz+2)
  integer*4 :: East,West
  ! fill in the center part
  density(1:Nx,0:Nz+1)=zones(k)%species(n)%var(1)%density(1:Nx,0:Nz+1)
  Gamma(1:Nx,0:Nz+1)=zones(k)%species(n)%var(1)%Gamma(1:Nx,0:Nz+1)
  T(1:Nx,0:Nz+1)=zones(k)%species(n)%var(1)%temperature(1:Nx,0:Nz+1)
  uEt(1:Nx,0:Nz+1)=zones(k)%species(n)%drifts%uEt(1:Nx,0:Nz+1)
  uBt(1:Nx,0:Nz+1)=zones(k)%species(n)%drifts%uBt(1:Nx,0:Nz+1)
  ! east
  East=zones(k)%Neighbors(3)
  if(East.gt.0) then
     density(1:Nx,Nz+2)=zones(East)%species(n)%var(1)%density(1:Nx,2)
     Gamma(1:Nx,Nz+2)=zones(East)%species(n)%var(1)%Gamma(1:Nx,2)
     T(1:Nx,Nz+2)=zones(East)%species(n)%var(1)%temperature(1:Nx,2)
     uEt(1:Nx,Nz+2)=zones(East)%species(n)%drifts%uEt(1:Nx,2)
     uBt(1:Nx,Nz+2)=zones(East)%species(n)%drifts%uBt(1:Nx,2)
  else 
     density(1:Nx,Nz+2)=density(1:Nx,Nz+1)
     Gamma(1:Nx,Nz+2)=Gamma(1:Nx,Nz+1)
     T(1:Nx,Nz+2)=T(1:Nx,Nz+1)
     uEt(1:Nx,Nz+2)=uEt(1:Nx,Nz+1)
     uBt(1:Nx,Nz+2)=uBt(1:Nx,Nz+1)
  end if
  ! west
  West=zones(k)%Neighbors(4)
  if(West.gt.0) then
     density(1:Nx,-1)=zones(West)%species(n)%var(1)%density(1:Nx,zones(West)%mesh%Nz-1)
     Gamma(1:Nx,-1)=zones(West)%species(n)%var(1)%Gamma(1:Nx,zones(West)%mesh%Nz-1)
     T(1:Nx,-1)=zones(West)%species(n)%var(1)%temperature(1:Nx,zones(West)%mesh%Nz-1)
     uEt(1:Nx,-1)=zones(West)%species(n)%drifts%uEt(1:Nx,zones(West)%mesh%Nz-1)
     uBt(1:Nx,-1)=zones(West)%species(n)%drifts%uBt(1:Nx,zones(West)%mesh%Nz-1)
  else
     density(1:Nx,-1)=density(1:Nx,0)
     Gamma(1:Nx,-1)=Gamma(1:Nx,0)
     T(1:Nx,-1)=T(1:Nx,0)
     uEt(1:Nx,-1)=uEt(1:Nx,0)
     uBt(1:Nx,-1)=uBt(1:Nx,0)
  end if
end subroutine extend_fields_before_weno
