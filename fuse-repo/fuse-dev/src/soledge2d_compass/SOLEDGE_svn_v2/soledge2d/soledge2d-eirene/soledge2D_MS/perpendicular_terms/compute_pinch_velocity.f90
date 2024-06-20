subroutine compute_pinch_velocity(zone)
  use all_variables, only : global_parameters, transport_parameters, reference_parameters
  use MZone
  implicit none
  Type(TZone), intent(inout) :: zone
  integer*4 :: i,j,Nx,Nz
  real*8 :: lambda_T_inv, lambda_B_inv
  integer*4 :: nspe, n
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        lambda_T_inv=((zone%species(0)%var(1)%temperature(i+1,j)-zone%species(0)%var(1)%temperature(i-1,j))&
             /(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))*&
             zone%metric_coefficients%cpp(i,j)&
             +(zone%species(0)%var(1)%temperature(i,j+1)-zone%species(0)%var(1)%temperature(i,j-1))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))*&
             zone%metric_coefficients%cpt(i,j))/zone%species(0)%var(1)%temperature(i,j)
        lambda_B_inv=((zone%mesh%B(i+1,j)-zone%mesh%B(i-1,j))&
             /(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))*&
             zone%metric_coefficients%cpp(i,j)&
             +(zone%mesh%B(i,j+1)-zone%mesh%B(i,j-1))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))*&
             zone%metric_coefficients%cpt(i,j))/zone%mesh%B(i,j)
        do n=1,global_parameters%N_ions
           nspe=global_parameters%ions_list(n,1)
           zone%species(n)%transport_perp%v_pinch(i,j)=transport_parameters%v_pinch(nspe)*&
                (2.d0*lambda_B_inv+lambda_T_inv)*reference_parameters%geometry%R0
           zone%species(n)%transport_perp%v_pinch(i,j)=zone%species(n)%transport_perp%v_pinch(i,j)*(1.d0-zone%masks%chi2(i,j))
           !DD
           zone%species(n)%transport_perp%v_pinch(i,j)=zone%species(n)%transport_perp%v_pinch(i,j)&
                /(reference_parameters%geometry%rs0/reference_parameters%fields%tau0)
        end do
        !quick fix
        zone%species(0)%transport_perp%v_pinch(i,j)=zone%species(1)%transport_perp%v_pinch(i,j)
     end do
  end do
end subroutine compute_pinch_velocity
