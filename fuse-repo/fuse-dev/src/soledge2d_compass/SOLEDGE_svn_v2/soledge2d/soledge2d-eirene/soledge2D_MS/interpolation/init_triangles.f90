subroutine init_triangles()
  use all_variables, only : interp_data2
  implicit none
  call load_triangle_mesh()
  call compute_interpolation_coefficients()
  call find_closest_point_ultimate_interp()
end subroutine init_triangles
