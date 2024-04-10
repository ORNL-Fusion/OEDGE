subroutine read_geometry_slab()
  use all_variables, only : global_parameters,zones,reference_parameters
  use Mzone
  implicit none
  integer*4 :: k
  integer*4 :: Nx,Nz
  do k=1,global_parameters%N_zones
     zones(k)%mesh%Rgeom=reference_parameters%geometry%R0
     zones(k)%mesh%B=reference_parameters%geometry%Btor0
  end do
end subroutine read_geometry_slab

