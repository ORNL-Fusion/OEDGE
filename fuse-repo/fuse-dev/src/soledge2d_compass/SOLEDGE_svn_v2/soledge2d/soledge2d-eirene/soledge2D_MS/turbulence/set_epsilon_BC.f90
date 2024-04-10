subroutine set_epsilon_BC(flux_surface)
  use Mflux_surface
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4 :: Nz
  Nz=flux_surface%Nz
  if(flux_surface%properties%is_periodic) then
     flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(1)
     flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(Nz)
  else
     !BC
     flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)
     flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(1)
  end if
end subroutine set_epsilon_BC
