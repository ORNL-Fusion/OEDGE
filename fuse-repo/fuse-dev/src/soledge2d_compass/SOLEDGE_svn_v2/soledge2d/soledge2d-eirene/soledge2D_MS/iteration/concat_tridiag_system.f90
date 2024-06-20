subroutine concat_tridiag_system(flux_surface,n_ion)
  use all_variables, only : global_variables, zones
  use Mflux_surface
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz,offset,i_psi
  integer*4 :: n,izone
  Nz=flux_surface%Nz
  offset=1
  i_psi=flux_surface%properties%i_psi
  do n=1,flux_surface%properties%n_zones
     izone=flux_surface%properties%zones(n)
     flux_surface%tridiag%a(offset:offset-1+Zones(izone)%mesh%Nz)=&
          zones(izone)%species(n_ion)%tridiag%a(i_psi,1:zones(izone)%mesh%Nz)
     flux_surface%tridiag%b(offset:offset-1+Zones(izone)%mesh%Nz)=&
          zones(izone)%species(n_ion)%tridiag%b(i_psi,1:zones(izone)%mesh%Nz)
     flux_surface%tridiag%c(offset:offset-1+Zones(izone)%mesh%Nz)=&
          zones(izone)%species(n_ion)%tridiag%c(i_psi,1:zones(izone)%mesh%Nz)
     flux_surface%tridiag%S(offset:offset-1+Zones(izone)%mesh%Nz)=&
          zones(izone)%species(n_ion)%tridiag%S(i_psi,1:zones(izone)%mesh%Nz)
     offset=offset+Zones(izone)%mesh%Nz
  end do
end subroutine concat_tridiag_system
