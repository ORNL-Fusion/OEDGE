subroutine compute_sinepitch_slab()
  use all_variables, only : zones, global_parameters, reference_parameters
  use Mphysics
  implicit none
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  real*8 :: rs0,R0
  real*8 :: v1R,v1Z,Rmid,v2R,v2Z
  real*8 :: dSR,dSZ
  real*8 :: Br,Bz,Bphi,norm
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     zones(k)%metric_coefficients%sinepitch_east=&
          0.5*(zones(k)%metric_coefficients%G(1:Nx,1:Nz)/sqrt(zones(k)%metric_coefficients%ctt(1:Nx,1:Nz))*(rs0/(2.D0*pi*R0))+&
          zones(k)%metric_coefficients%G(1:Nx,2:Nz+1)/sqrt(zones(k)%metric_coefficients%ctt(1:Nx,2:Nz+1))*(rs0/(2.D0*pi*R0)))
     zones(k)%metric_coefficients%sinepitch_west=&
          0.5*(zones(k)%metric_coefficients%G(1:Nx,1:Nz)/sqrt(zones(k)%metric_coefficients%ctt(1:Nx,1:Nz))*(rs0/(2.D0*pi*R0))+&
          zones(k)%metric_coefficients%G(1:Nx,0:Nz-1)/sqrt(zones(k)%metric_coefficients%ctt(1:Nx,0:Nz-1))*(rs0/(2.D0*pi*R0)))
  end do
end subroutine compute_sinepitch_slab
