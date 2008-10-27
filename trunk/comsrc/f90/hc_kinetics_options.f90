module hc_kinetics_options

  integer :: hc_kinetics_opt  ! overall option for hc kinetics implementation - 0 original code - 1 updated kinetics
  integer :: hc_vperp_opt     ! vperp treatment option - 0 vperp =const, 1 vperpequal vpara, 2 vperp -> tperp 3 vperp diffuses similarly to vpara
  integer :: hc_mtc_3D_opt    ! 3D momentum transfer collision option - particle path changes by 90 degrees into plane perpendicular to current path - 0 off, 1 on
  !integer :: hc_tperp_opt     ! option for adjusting the particle temperature

  ! reaction_type
  integer,parameter :: nn_reaction = 1, ni_reaction =2, in_reaction = 3, ii_reaction=4
  ! reaction_kind
  integer,parameter :: e_reaction = 1, p_reaction=2
  ! options for mapping velocities between ion and neutral
  integer,parameter :: map_to_3D=1, map_to_bfield=2
  ! options for mapping to and from temperature and velocity values
  integer,parameter :: conv_to_t=1, conv_to_v=2


  ! kinetics debug option
  !logical,parameter :: debug_kinetics = .true.
  logical,parameter :: debug_kinetics = .false.


end module hc_kinetics_options
