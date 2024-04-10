subroutine correct_BC_sources(zone,Fluxes_n,Fluxes_G,Fluxes_E,Nx,Nz,n_ion,fluxm,fluxp)
  use all_variables, only : reference_parameters, global_variables, flags
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: Nx,Nz,n_ion
  real*8,intent(inout) :: Fluxes_n(1:Nx,1:Nz,1:4)
  real*8,intent(inout) :: Fluxes_G(1:Nx,1:Nz,1:4)
  real*8,intent(inout) :: Fluxes_E(1:Nx,1:Nz,1:4)
  real*8,intent(in) :: fluxm(1:Nx,0:Nz,1:3)
  real*8,intent(in) :: fluxp(1:Nx,0:Nz,1:3)
  integer*4 :: i
  !east
  if(zone%Neighbors(3).lt.0) then
     do i=1,Nx
        if(Fluxes_n(i,Nz,3).lt.0.D0) then
           Fluxes_n(i,Nz,3)=fluxm(i,Nz,1)/(0.5d0*(zone%metric_coefficients%jacobian(i,Nz)+zone%metric_coefficients%jacobian(i,Nz+1)))&
                /(0.5d0*(sqrt(zone%metric_coefficients%ctt(i,Nz))+sqrt(zone%metric_coefficients%ctt(i,Nz+1))))
           Fluxes_G(i,Nz,3)=fluxm(i,Nz,2)/(0.5d0*(zone%metric_coefficients%jacobian(i,Nz)+zone%metric_coefficients%jacobian(i,Nz+1)))&
                /(0.5d0*(sqrt(zone%metric_coefficients%ctt(i,Nz))+sqrt(zone%metric_coefficients%ctt(i,Nz+1))))
           Fluxes_E(i,Nz,3)=fluxm(i,Nz,3)/(0.5d0*(zone%metric_coefficients%jacobian(i,Nz)+zone%metric_coefficients%jacobian(i,Nz+1)))&
                /(0.5d0*(sqrt(zone%metric_coefficients%ctt(i,Nz))+sqrt(zone%metric_coefficients%ctt(i,Nz+1))))
           Fluxes_n(i,Nz,3)=0.d0
           Fluxes_G(i,Nz,3)=0.d0
           Fluxes_E(i,Nz,3)=0.d0
        end if
     end do
  end if
  !west
  if(zone%Neighbors(4).lt.0) then
     do i=1,Nx
        if(Fluxes_n(i,1,4).gt.0.D0) then
           Fluxes_n(i,1,4)=fluxp(i,0,1)/(0.5d0*(zone%metric_coefficients%jacobian(i,1)+zone%metric_coefficients%jacobian(i,0)))&
                /(0.5d0*(sqrt(zone%metric_coefficients%ctt(i,1))+sqrt(zone%metric_coefficients%ctt(i,0))))
           Fluxes_G(i,1,4)=fluxp(i,0,2)/(0.5d0*(zone%metric_coefficients%jacobian(i,1)+zone%metric_coefficients%jacobian(i,0)))&
                /(0.5d0*(sqrt(zone%metric_coefficients%ctt(i,1))+sqrt(zone%metric_coefficients%ctt(i,0))))
           Fluxes_E(i,1,4)=fluxp(i,0,3)/(0.5d0*(zone%metric_coefficients%jacobian(i,1)+zone%metric_coefficients%jacobian(i,0)))&
                /(0.5d0*(sqrt(zone%metric_coefficients%ctt(i,1))+sqrt(zone%metric_coefficients%ctt(i,0))))
           Fluxes_n(i,1,4)=0.d0
           Fluxes_G(i,1,4)=0.d0
           Fluxes_E(i,1,4)=0.d0
        end if
     end do
  end if
end subroutine correct_BC_sources
