subroutine allocate_metric_arrays(zone)
  use MZone 
  implicit none
  type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(zone%metric_coefficients%cpp(0:Nx+1,0:Nz+1))
  allocate(zone%metric_coefficients%cpt(0:Nx+1,0:Nz+1))
  allocate(zone%metric_coefficients%ctt(0:Nx+1,0:Nz+1))
  allocate(zone%metric_coefficients%c_pp(0:Nx+1,0:Nz+1))
  allocate(zone%metric_coefficients%c_pt(0:Nx+1,0:Nz+1))
  allocate(zone%metric_coefficients%c_tt(0:Nx+1,0:Nz+1))
  allocate(zone%metric_coefficients%Jacobian(0:Nx+1,0:Nz+1))
  allocate(zone%metric_coefficients%G(0:Nx+1,0:Nz+1))
  allocate(zone%metric_coefficients%dS_north_PU(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dS_south_PU(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dS_east_PU(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dS_west_PU(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dS_north_DD(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dS_south_DD(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dS_east_DD(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dS_west_DD(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dvol_PU(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%dvol_DD(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%sinepitch_east(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%sinepitch_west(1:Nx,1:Nz))
  allocate(zone%metric_coefficients%divb(1:Nx,1:Nz))
  allocate(zone%jacobian%dPdR(0:Nx+1,0:Nz+1))
  allocate(zone%jacobian%dPdZ(0:Nx+1,0:Nz+1))
  allocate(zone%jacobian%dTdR(0:Nx+1,0:Nz+1))
  allocate(zone%jacobian%dTdZ(0:Nx+1,0:Nz+1))
  allocate(zone%jacobian%dRdP(0:Nx+1,0:Nz+1))
  allocate(zone%jacobian%dRdT(0:Nx+1,0:Nz+1))
  allocate(zone%jacobian%dZdP(0:Nx+1,0:Nz+1))
  allocate(zone%jacobian%dZdT(0:Nx+1,0:Nz+1))
  !set to non zero  value
  zone%metric_coefficients%cpp(0:Nx+1,0:Nz+1)=1.d-15
  zone%metric_coefficients%cpt(0:Nx+1,0:Nz+1)=1.d-15
  zone%metric_coefficients%ctt(0:Nx+1,0:Nz+1)=1.d-15
  zone%metric_coefficients%c_pp(0:Nx+1,0:Nz+1)=1.d-15
  zone%metric_coefficients%c_pt(0:Nx+1,0:Nz+1)=1.d-15
  zone%metric_coefficients%c_tt(0:Nx+1,0:Nz+1)=1.d-15
  zone%metric_coefficients%Jacobian(0:Nx+1,0:Nz+1)=1.d-15
  zone%metric_coefficients%G(0:Nx+1,0:Nz+1)=1.d-15
  zone%jacobian%dPdR(0:Nx+1,0:Nz+1)=1.d-15
  zone%jacobian%dPdZ(0:Nx+1,0:Nz+1)=1.d-15
  zone%jacobian%dTdR(0:Nx+1,0:Nz+1)=1.d-15
  zone%jacobian%dTdZ(0:Nx+1,0:Nz+1)=1.d-15
  zone%jacobian%dRdP(0:Nx+1,0:Nz+1)=1.d-15
  zone%jacobian%dZdP(0:Nx+1,0:Nz+1)=1.d-15
  zone%jacobian%dRdT(0:Nx+1,0:Nz+1)=1.d-15
  zone%jacobian%dZdT(0:Nx+1,0:Nz+1)=1.d-15
end subroutine allocate_metric_arrays
