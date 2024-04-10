subroutine compute_metric_coefficients_slab()
  use all_variables, only : zones, global_parameters, reference_parameters
  use Mphysics
  implicit none
  integer*4 :: i,j,k
  real*8 :: DVOL
  real*8 :: R0,rs0,Rm0,qref
  integer*4 :: Nx,Nz
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  Rm0=reference_parameters%geometry%Rm0
  qref=reference_parameters%geometry%qref
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     zones(k)%metric_coefficients%Jacobian=R0*rs0*Rm0
     zones(k)%metric_coefficients%ctt=(1.d0/Rm0)**(2.d0)*&
          (rs0/(2.d0*pi))**(2.d0)
     zones(k)%metric_coefficients%cpp=1.d0
     zones(k)%metric_coefficients%cpt=0.d0
     zones(k)%metric_coefficients%c_tt=1.d0/zones(k)%metric_coefficients%ctt
     zones(k)%metric_coefficients%c_pt=0.d0
     zones(k)%metric_coefficients%c_pp=1.d0/zones(k)%metric_coefficients%cpp
     zones(k)%metric_coefficients%G=1.d0/sqrt((qref*R0)**(2.d0)+Rm0**2.d0)*R0
     zones(k)%metric_coefficients%divb=0.D0       
     zones(k)%mesh%Rgeom=R0
     zones(k)%mesh%Bphi=qref*R0*zones(k)%mesh%B/sqrt(Rm0**2.d0+qref**2.d0*R0**2.d0)
  end do
  call compute_diffareas_slab()
  call compute_dvol_slab()
  call compute_sinepitch_slab()
  call MD_broadcast_metric()
  call compute_Score()
  call save_metric()
end subroutine compute_metric_coefficients_slab
