subroutine allocate_interpolated_data2()
  use all_variables, only : global_parameters,Interp_data2
  implicit none
  allocate(Interp_data2%knots_R(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_Z(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_B(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_Br(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_Bz(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_Bphi(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_interp_points(1:Interp_data2%N_knots))
  allocate(Interp_data2%tri_knots(1:Interp_data2%N_triangles,1:3))
  allocate(Interp_data2%type_face(1:Interp_data2%N_triangles,1:3))
  ! input fields
  allocate(Interp_data2%knots_density(1:Interp_data2%N_knots,0:global_parameters%N_ions,1:1))
  allocate(Interp_data2%knots_velocity(1:Interp_data2%N_knots,0:global_parameters%N_ions,1:1))
  allocate(Interp_data2%knots_temperature(1:Interp_data2%N_knots,0:global_parameters%N_ions,1:1))
  allocate(Interp_data2%tri_fluxn(1:Interp_data2%N_triangles,1:3,0:global_parameters%N_ions,1:1))
  allocate(Interp_data2%tri_fluxG(1:Interp_data2%N_triangles,1:3,0:global_parameters%N_ions,1:1))
  allocate(Interp_data2%tri_fluxE(1:Interp_data2%N_triangles,1:3,0:global_parameters%N_ions,1:1))
  ! output fields
  allocate(Interp_data2%tri_Sn(1:Interp_data2%N_triangles,0:global_parameters%N_species,1:1))
  allocate(Interp_data2%tri_SG(1:Interp_data2%N_triangles,0:global_parameters%N_species,1:1))
  allocate(Interp_data2%tri_SE(1:Interp_data2%N_triangles,0:global_parameters%N_species,1:1))
  allocate(Interp_data2%tri_Nn(1:Interp_data2%N_triangles,0:global_parameters%N_species,1:1))
  allocate(Interp_data2%tri_Nm(1:Interp_data2%N_triangles,0:global_parameters%N_species,1:1))
  allocate(Interp_data2%tri_Tn(1:Interp_data2%N_triangles,0:global_parameters%N_species,1:1))
  allocate(Interp_data2%tri_Tm(1:Interp_data2%N_triangles,0:global_parameters%N_species,1:1))
  allocate(Interp_data2%tri_Srad(1:Interp_data2%N_triangles,0:global_parameters%N_species,1:1))
  allocate(Interp_data2%Puff(1:global_parameters%N_species))
  allocate(Interp_data2%neutral_outflux(1:global_parameters%N_species))
  !miscellaneous
  allocate(Interp_data2%knots_radiation(1:Interp_data2%N_knots,0:global_parameters%N_ions))
  allocate(Interp_data2%knots_Zeff(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_phi(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_vorticity(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_Epara(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_k(1:Interp_data2%N_knots))
  allocate(Interp_data2%knots_epsilon(1:Interp_data2%N_knots))
end subroutine allocate_interpolated_data2

subroutine allocate_wall_data()
  use all_variables, only : Interp_data2
  implicit none
  allocate(Interp_data2%wall_data%wall_triangle_face(1:Interp_data2%wall_data%N_triangles_on_wall))
  allocate(Interp_data2%wall_data%wall_triangle_surf(1:Interp_data2%wall_data%N_triangles_on_wall))
  allocate(Interp_data2%wall_data%wall_triangle_sabs(1:Interp_data2%wall_data%N_triangles_on_wall))
  allocate(Interp_data2%wall_data%back_interp_on_wall(1:Interp_data2%wall_data%N_triangles_on_wall,1:4))
  allocate(Interp_data2%wall_data%s2d_to_use(1:Interp_data2%wall_data%N_triangles_on_wall,1:4,1:3))
  allocate(Interp_data2%wall_data%weight_s2d_to_use(1:Interp_data2%wall_data%N_triangles_on_wall,1:4))
end subroutine allocate_wall_data
