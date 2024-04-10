subroutine compute_dvol_slab()
  use all_variables, only : zones, global_parameters, reference_parameters
  use Mphysics
  implicit none
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  real*8 :: rs0,R0,Rm0
  real*8 :: v1R,v1Z,Rmid,v2R,v2Z
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  Rm0=reference_parameters%geometry%Rm0
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     zones(k)%metric_coefficients%dvol_PU=abs(&
          2.*pi*R0&
          *(Zones(k)%mesh%x_plus_1half(1:Nx,1:Nz)-&
          Zones(k)%mesh%x_minus_1half(1:Nx,1:Nz))*rs0&
          *2.D0*pi*(Zones(k)%mesh%z_plus_1half(1:Nx,1:Nz)-&
          Zones(k)%mesh%z_minus_1half(1:Nx,1:Nz))*Rm0)
     zones(k)%metric_coefficients%dvol_DD=zones(k)%metric_coefficients%dvol_PU/(2.d0*pi*reference_parameters%geometry%rs0**2)
  end do
end subroutine compute_dvol_slab
