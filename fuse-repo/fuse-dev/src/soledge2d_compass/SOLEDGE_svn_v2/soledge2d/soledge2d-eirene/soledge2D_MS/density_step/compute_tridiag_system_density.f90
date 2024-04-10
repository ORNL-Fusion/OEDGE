subroutine compute_tridiag_system_density(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4 :: n
  real*8 :: dvol,ds_east,ds_west
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: rs0,dt
  rs0=reference_parameters%geometry%rs0
  dt=global_variables%dt
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=1,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           call compute_tridiag_matrix_implicit_density_point(zone,i,j,n)
           call add_density_matrix_perp_diffusion_point(zone,i,j,n)
           call add_density_explicit_source_point(zone,i,j,n)
           call add_density_am_source_point(zone,i,j,n)
           call add_density_penalisation_terms_point(zone,i,j,n)
           call add_density_neutrals_terms_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine compute_tridiag_system_density
