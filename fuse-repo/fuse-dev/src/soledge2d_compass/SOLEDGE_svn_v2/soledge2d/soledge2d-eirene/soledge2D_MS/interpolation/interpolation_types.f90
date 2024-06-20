module MInterpolation_types

  type zone_temp
     real*8,allocatable:: val(:,:,:)
  end type zone_temp

  type :: Interpolation_temp
     real*8,allocatable :: knots_val(:,:)
     type(zone_temp),allocatable :: zones(:)
  end type Interpolation_temp
  
  Type :: TInterp_points
     integer*4 :: pass
     integer*4 :: n_soledge
     integer*4 :: n_eirene
     integer*4 :: sol(8,3) ! soledge points for interpolation (n_soledge*3 [i,j,k]) - 8 points max
     integer*4 :: eir(3) ! eirene points for interpolation (n_eirene) ! 3 points max
     real*8 :: interp_matrix(3,3) ! z=a*R+b*Z+c
     integer*4 :: type_ultimate_point ! 1: Soledge ; 2: Eirene
     integer*4 :: ultimate_coords(3)
  end type TInterp_points

  Type :: Twall_data
     integer*4 :: N_triangles_on_wall
     integer*4,allocatable :: wall_triangle_face(:)
     integer*4,allocatable :: wall_triangle_surf(:)
     integer*4,allocatable :: wall_triangle_sabs(:)
     integer*4,allocatable :: back_interp_on_wall(:,:) ! ntri_on_wall * 4 (4 ntri,i,j,k)
     integer*4,allocatable :: s2d_to_use(:,:,:) !ntri_on_wall * 4 * 3 (4NSEW, 3ijk)
     real*8,allocatable :: weight_s2d_to_use(:,:) !ntri_on_wall * 4 (4NSEW)
  end type Twall_data
  
  Type :: Tzone_data_tri
     integer*4,allocatable :: triangles(:,:,:)
     integer*4,allocatable :: num_triangles(:,:)
  end type Tzone_data_tri
  
  Type :: TInterpolated_data2
     integer*4 :: N_triangles
     integer*4 :: N_knots
     ! next 2 variables : ensure styx compatibility
     integer*4 :: n_tor
     real*8 :: ang_max
     real*8,allocatable :: knots_R(:), knots_Z(:)
     real*8,allocatable :: knots_Br(:), knots_Bz(:), knots_Bphi(:), knots_B(:)
     type(Tinterp_points),allocatable :: knots_interp_points(:)
     type(Twall_data) :: wall_data
     integer*4,allocatable :: tri_knots(:,:) ! ntriangle * 3
     integer*4,allocatable :: type_face(:,:) ! ntriangle * 3
     real*8,allocatable :: knots_density(:,:,:)
     real*8,allocatable :: knots_velocity(:,:,:)
     real*8,allocatable :: knots_temperature(:,:,:)
     real*8,allocatable :: knots_radiation(:,:)
     real*8,allocatable :: knots_Zeff(:)
     real*8,allocatable :: knots_phi(:)
     real*8,allocatable :: knots_vorticity(:)
     real*8,allocatable :: knots_Epara(:)
     real*8,allocatable :: knots_k(:)
     real*8,allocatable :: knots_epsilon(:)
     real*8,allocatable :: tri_fluxn(:,:,:,:)
     real*8,allocatable :: tri_fluxG(:,:,:,:)
     real*8,allocatable :: tri_fluxE(:,:,:,:)
     real*8,allocatable :: tri_Sn(:,:,:)
     real*8,allocatable :: tri_SG(:,:,:)
     real*8,allocatable :: tri_SE(:,:,:)
     real*8,allocatable :: tri_Nn(:,:,:)
     real*8,allocatable :: tri_Nm(:,:,:)
     real*8,allocatable :: tri_Tn(:,:,:)
     real*8,allocatable :: tri_Tm(:,:,:)
     real*8,allocatable :: tri_Srad(:,:,:)
     Type(Tzone_data_tri),allocatable :: zones_data(:)
     real*8,allocatable :: Puff(:)
     real*8,allocatable :: neutral_outflux(:)
  end type TInterpolated_data2

end module MInterpolation_types
