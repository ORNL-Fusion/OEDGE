subroutine apply_epsilon_boundary_conditions(flux_surface,n_ion)
  use all_variables, only : global_variables, zones
  use Mflux_surface
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz,offset,i_psi
  integer*4 :: n,izone
  Nz=flux_surface%Nz
  i_psi=flux_surface%properties%i_psi
  if(.not.flux_surface%properties%is_periodic) then
     !west BC : n(0)=n(1)
     flux_surface%tridiag%b(1)=flux_surface%tridiag%b(1)+&
          flux_surface%tridiag%a(1)
     flux_surface%tridiag%a(1)=0.d0
     !east BC : n(Nz+1)=n(Nz)
     flux_surface%tridiag%b(Nz)=flux_surface%tridiag%b(Nz)+&
          flux_surface%tridiag%c(Nz)
     flux_surface%tridiag%c(Nz)=0.d0
  end if
end subroutine apply_epsilon_boundary_conditions
