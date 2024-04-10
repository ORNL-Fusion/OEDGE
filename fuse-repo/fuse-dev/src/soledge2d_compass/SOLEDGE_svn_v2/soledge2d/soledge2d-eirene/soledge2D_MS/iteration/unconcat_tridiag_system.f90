subroutine unconcat_tridiag_system(flux_surface,n_ion,FIELD)
  use all_variables, only : global_variables, zones
  use Mflux_surface
  use Mdefinitions
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz,offset,i_psi
  integer*4 :: n,izone
  integer*4 :: FIELD
  integer*4 :: i,j,k
  Nz=flux_surface%Nz
  i_psi=flux_surface%properties%i_psi
  select case(FIELD)
  case(DENSITY_FIELD)
     offset=0
     do k=1,flux_surface%properties%n_zones
        izone=flux_surface%properties%zones(k)
        do j=0,zones(izone)%mesh%Nz+1
           zones(izone)%species(n_ion)%var(2)%density(i_psi,j)=&
                flux_surface%tridiag%buffer(j+offset)
        end do
        offset=offset+zones(izone)%mesh%Nz
     end do
  case(VELOCITY_FIELD)
     offset=0
     do k=1,flux_surface%properties%n_zones
        izone=flux_surface%properties%zones(k)
        do j=0,zones(izone)%mesh%Nz+1
           zones(izone)%species(n_ion)%var(2)%Gamma(i_psi,j)=&
                flux_surface%tridiag%buffer(j+offset)
        end do
        offset=offset+zones(izone)%mesh%Nz
     end do
  case(TEMPERATURE_FIELD)
     offset=0
     do k=1,flux_surface%properties%n_zones
        izone=flux_surface%properties%zones(k)
        do j=0,zones(izone)%mesh%Nz+1
           zones(izone)%species(n_ion)%var(2)%temperature(i_psi,j)=&
                flux_surface%tridiag%buffer(j+offset)
        end do
        offset=offset+zones(izone)%mesh%Nz
     end do
  case(K_FIELD)
     offset=0
     do k=1,flux_surface%properties%n_zones
        izone=flux_surface%properties%zones(k)
        do j=0,zones(izone)%mesh%Nz+1
           zones(izone)%kepsilon(2)%k(i_psi,j)=&
                flux_surface%tridiag%buffer(j+offset)
        end do
        offset=offset+zones(izone)%mesh%Nz
     end do
  case(EPSILON_FIELD)
     offset=0
     do k=1,flux_surface%properties%n_zones
        izone=flux_surface%properties%zones(k)
        do j=0,zones(izone)%mesh%Nz+1
           zones(izone)%kepsilon(2)%epsilon(i_psi,j)=&
                flux_surface%tridiag%buffer(j+offset)
        end do
        offset=offset+zones(izone)%mesh%Nz
     end do
  case(UET_FIELD)
     offset=0
     do k=1,flux_surface%properties%n_zones
        izone=flux_surface%properties%zones(k)
        do j=0,zones(izone)%mesh%Nz+1
           zones(izone)%species(n_ion)%drifts%uEt(i_psi,j)=&
                flux_surface%tridiag%buffer(j+offset)*zones(izone)%metric_coefficients%G(i_psi,j)/sqrt(zones(izone)%metric_coefficients%ctt(i_psi,j))
        end do
        offset=offset+zones(izone)%mesh%Nz
     end do
  case(UBT_FIELD)
     offset=0
     do k=1,flux_surface%properties%n_zones
        izone=flux_surface%properties%zones(k)
        do j=0,zones(izone)%mesh%Nz+1
           zones(izone)%species(n_ion)%drifts%uBt(i_psi,j)=&
                flux_surface%tridiag%buffer(j+offset)*zones(izone)%metric_coefficients%G(i_psi,j)/sqrt(zones(izone)%metric_coefficients%ctt(i_psi,j))
        end do
        offset=offset+zones(izone)%mesh%Nz
     end do
  end select
end subroutine unconcat_tridiag_system
