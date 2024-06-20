subroutine compute_diffareas()
  use all_variables, only : zones, global_parameters, reference_parameters
  use Mphysics
  implicit none
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  real*8 :: rs0,R0
  real*8 :: v1R,v1Z,Rmid
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           !north
           zones(k)%metric_coefficients%dS_north_PU(i,j)=0.5d0*(zones(k)%metric_coefficients%jacobian(i,j)&
                +zones(k)%metric_coefficients%jacobian(i+1,j))&
                *0.5D0*(sqrt(zones(k)%metric_coefficients%cpp(i,j))+sqrt(zones(k)%metric_coefficients%cpp(i+1,j)))&
                *2.D0*pi*(zones(k)%mesh%z_plus_1half(i,j)-zones(k)%mesh%z_minus_1half(i,j))*2.D0*pi/rs0
           zones(k)%metric_coefficients%dS_north_DD(i,j)=zones(k)%metric_coefficients%dS_north_PU(i,j)/(2.d0*pi*rs0)
           !south
           zones(k)%metric_coefficients%dS_south_PU(i,j)=0.5d0*(zones(k)%metric_coefficients%jacobian(i,j)&
                +zones(k)%metric_coefficients%jacobian(i-1,j))&
                *0.5D0*(sqrt(zones(k)%metric_coefficients%cpp(i,j))+sqrt(zones(k)%metric_coefficients%cpp(i-1,j)))&
                *2.D0*pi*(zones(k)%mesh%z_plus_1half(i,j)-zones(k)%mesh%z_minus_1half(i,j))*2.D0*pi/rs0
           zones(k)%metric_coefficients%dS_south_DD(i,j)=zones(k)%metric_coefficients%dS_south_PU(i,j)/(2.d0*pi*rs0)
           !east
           zones(k)%metric_coefficients%dS_east_PU(i,j)=0.5d0*(zones(k)%metric_coefficients%jacobian(i,j)&
                +zones(k)%metric_coefficients%jacobian(i,j+1))&
                *0.5D0*(sqrt(zones(k)%metric_coefficients%ctt(i,j))+sqrt(zones(k)%metric_coefficients%ctt(i,j+1)))&
                *2.D0*pi*(zones(k)%mesh%x_plus_1half(i,j)-zones(k)%mesh%x_minus_1half(i,j))*2.D0*pi/rs0
           zones(k)%metric_coefficients%dS_east_DD(i,j)=zones(k)%metric_coefficients%dS_east_PU(i,j)/(2.d0*pi*rs0)
           !west
           zones(k)%metric_coefficients%dS_west_PU(i,j)=0.5d0*(zones(k)%metric_coefficients%jacobian(i,j)&
                +zones(k)%metric_coefficients%jacobian(i,j-1))&
                *0.5D0*(sqrt(zones(k)%metric_coefficients%ctt(i,j))+sqrt(zones(k)%metric_coefficients%ctt(i,j-1)))&
                *2.D0*pi*(zones(k)%mesh%x_plus_1half(i,j)-zones(k)%mesh%x_minus_1half(i,j))*2.D0*pi/rs0
           zones(k)%metric_coefficients%dS_west_DD(i,j)=zones(k)%metric_coefficients%dS_west_PU(i,j)/(2.d0*pi*rs0)
        end do
     end do
  end do
end subroutine compute_diffareas
