subroutine compute_electric_field(zone)
  use all_variables, only : global_parameters
  use MZone
  implicit none
  type(TZone),intent(inout) :: zone
  integer*4 :: i,j,m
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        zone%electric_fields(1)%E(i,j) = -1.D0/zone%species(0)%var(1)%density(i,j)&
             *zone%metric_coefficients%G(i,j)&
             *(zone%species(0)%var(1)%temperature(i,j+1)*zone%species(0)%var(1)%density(i,j+1)&
             -zone%species(0)%var(1)%temperature(i,j-1)*zone%species(0)%var(1)%density(i,j-1))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))
        do m=0,global_parameters%N_ions
           zone%electric_fields(1)%E(i,j) = zone%electric_fields(1)%E(i,j) &
                + zone%species(0)%coupling_terms%R(i,j,m)/zone%species(0)%var(1)%density(i,j)
        end do

        zone%electric_fields(1)%Etheta(i,j)=-(zone%electric_fields(1)%phi(i,j+1)-zone%electric_fields(1)%phi(i,j-1))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))*sqrt(zone%metric_coefficients%ctt(i,j))

        zone%electric_fields(1)%Epsi(i,j)=-(zone%electric_fields(1)%phi(i+1,j)-zone%electric_fields(1)%phi(i-1,j))&
             /(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))*sqrt(zone%metric_coefficients%cpp(i,j))

     end do
  end do
end subroutine compute_electric_field
