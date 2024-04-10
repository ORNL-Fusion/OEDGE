subroutine compute_power_logs(zone,pow_logNe,pow_logTe,max_degree,Nx,Nz,Te_min,Te_max,Ne_min,Ne_max)
  use all_variables, only : reference_parameters, global_variables
  use Mzone
  type(Tzone),intent(in) :: zone
  real*8,intent(inout) :: pow_logNe(0:max_degree-1,1:Nx,1:Nz)
  real*8,intent(inout) :: pow_logTe(0:max_degree-1,1:Nx,1:Nz)
  integer*4,intent(in) :: max_degree
  integer*4,intent(in) :: Nx,Nz
  real*8,intent(in) :: Te_min,Te_max,Ne_min,Ne_max
  integer*4 :: k
  real*8 :: logTe(1:Nx,1:Nz)
  real*8 :: logNe(1:Nx,1:Nz)
  logNe=log10(min(max(zone%species(0)%var(1)%density(1:Nx,1:Nz)*reference_parameters%fields%n0,Ne_min),Ne_max))
  logTe=log10(min(max(zone%species(0)%var(1)%temperature(1:Nx,1:Nz)*reference_parameters%fields%T0eV,Te_min),Te_max))
  pow_logNe(0,:,:)=1.D0
  pow_logTe(0,:,:)=1.D0
  do k=1,max_degree-1
     pow_logNe(k,:,:)=pow_logNe(k-1,:,:)*logNe
     pow_logTe(k,:,:)=pow_logTe(k-1,:,:)*logTe
  end do
end subroutine compute_power_logs
