subroutine compute_dvol()
  use all_variables, only : zones, global_parameters, reference_parameters
  use Mphysics
  implicit none
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  real*8 :: rs0,R0
  real*8 :: v1R,v1Z,Rmid,v2R,v2Z
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           zones(k)%metric_coefficients%dvol_PU(i,j)=zones(k)%metric_coefficients%Jacobian(i,j)*&
                (zones(k)%mesh%z_plus_1half(i,j)-zones(k)%mesh%z_minus_1half(i,j))*2.D0*pi*&
                (zones(k)%mesh%x_plus_1half(i,j)-zones(k)%mesh%x_minus_1half(i,j))*2.D0*pi
           Zones(k)%metric_coefficients%dvol_DD(i,j)=zones(k)%metric_coefficients%dvol_PU(i,j)/(2.d0*pi*reference_parameters%geometry%rs0**2)
        end do
     end do
  end do
end subroutine compute_dvol
