subroutine compute_tridiag_system_velocity(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4 :: n,m
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: rs0,dt,R0,Teps
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  Teps=global_variables%Teps
  dt=global_variables%dt
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=1,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           call compute_tridiag_matrix_implicit_velocity_point(zone,i,j,n)
           call add_velocity_matrix_perp_diffusion_point(zone,i,j,n)
           call add_velocity_matrix_para_diffusion_point(zone,i,j,n)
           call add_velocity_explicit_source_point(zone,i,j,n)
           call add_velocity_am_source_point(zone,i,j,n)
           call add_velocity_divb_term_point(zone,i,j,n)
           call add_velocity_Eparallel_term_point(zone,i,j,n)
           call add_velocity_penalisation_terms_point(zone,i,j,n)
           call add_velocity_coupling_terms_point(zone,i,j,n)
           call add_velocity_neutrals_terms_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine compute_tridiag_system_velocity
