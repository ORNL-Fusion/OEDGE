subroutine compute_diffareas_slab()
  use all_variables, only : zones, global_parameters, reference_parameters
  use Mphysics
  implicit none
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  real*8 :: rs0,R0,Rm0
  real*8 :: v1R,v1Z,Rmid
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  Rm0=reference_parameters%geometry%Rm0
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     zones(k)%metric_coefficients%dS_north_PU=2.D0*pi*R0*&
          (Zones(k)%mesh%z_plus_1half(1:Nx,1:Nz)-&
          Zones(k)%mesh%z_minus_1half(1:Nx,1:Nz))*2.D0*pi*Rm0
     zones(k)%metric_coefficients%dS_north_DD=zones(k)%metric_coefficients%dS_north_PU/(2.d0*pi*rs0)

     zones(k)%metric_coefficients%dS_south_PU=2.D0*pi*R0*&
          (Zones(k)%mesh%z_plus_1half(1:Nx,1:Nz)-&
          Zones(k)%mesh%z_minus_1half(1:Nx,1:Nz))*2.D0*pi*Rm0
     zones(k)%metric_coefficients%dS_south_DD=zones(k)%metric_coefficients%dS_south_PU/(2.d0*pi*rs0)

     zones(k)%metric_coefficients%dS_east_PU=2.D0*pi*R0*&
          (Zones(k)%mesh%x_plus_1half(1:Nx,1:Nz)-&
          Zones(k)%mesh%x_minus_1half(1:Nx,1:Nz))*rs0
     zones(k)%metric_coefficients%dS_east_DD=zones(k)%metric_coefficients%dS_east_PU/(2.d0*pi*rs0)

     zones(k)%metric_coefficients%dS_west_PU=2.D0*pi*R0*&
          (Zones(k)%mesh%x_plus_1half(1:Nx,1:Nz)-&
          Zones(k)%mesh%x_minus_1half(1:Nx,1:Nz))*rs0
     zones(k)%metric_coefficients%dS_west_DD=zones(k)%metric_coefficients%dS_west_PU/(2.d0*pi*rs0)
  end do
end subroutine compute_diffareas_slab
