subroutine compute_tridiag_system_temperature(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters, drift_flags
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4 :: n,m
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           call compute_tridiag_matrix_implicit_temperature_point(zone,i,j,n)
           call add_temperature_matrix_perp_diffusion_point(zone,i,j,n)
           call add_temperature_matrix_para_diffusion_point(zone,i,j,n)
           call add_temperature_penalisation_terms_point(zone,i,j,n)
           call add_temperature_explicit_source_point(zone,i,j,n)
           call add_temperature_am_source_point(zone,i,j,n)
           call add_temperature_Eparallel_term_point(zone,i,j,n)
           call add_temperature_coupling_terms_point(zone,i,j,n)
           call add_temperature_neutrals_terms_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine compute_tridiag_system_temperature
