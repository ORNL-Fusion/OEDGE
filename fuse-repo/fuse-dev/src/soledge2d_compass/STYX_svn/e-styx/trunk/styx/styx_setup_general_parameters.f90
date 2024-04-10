subroutine styx_setup_general_parameters
  use all_variables, only : global_parameters,Interp_data2
  use eirmod_precision
  use styx2eirene
  implicit none

  numiter=global_parameters%N_iterations

! use all tallies (a bit slower, but otherwise internal balance in EIRENE cannot be checked)
  all_tal=.true.

! EIRENE uses cylindrical geometry
! use of this options implies consistency with metric coefficients in styx2D, not implemented
!  levgeo_styx=1

! EIRENE uses toroidal geometry 
  levgeo_styx=2

  NTRI_styx=Interp_data2%n_triangles

! Option for use with 3D fluid code
#ifdef TK3X
  is_3D=.true.
  fluid_code='tokam3X'
#endif
#ifdef S2D
  is_3D=.false.
  fluid_code='soledge2D'
#endif

  if (is_3D) then
    ! testing : full torus with 10 cells
    Ntor_cells=Interp_data2%n_tor
    Neir_cells=(NTRI_styx+1)*(Ntor_cells+1)
    ang_max=Interp_Data2%ang_max
  else
    Ntor_cells=1
    Neir_cells=NTRI_styx+1 ! +1 = test
    ! number of toroidal segments in EIRENE (torus approximation) 
    NTTRA_styx=30
    ang_max=360.d0
  endif

! major radius in cm (zero if the mesh coordinates include R, as here)
  ROA_styx=0.d0 

! number of recycling/recombination strata
  Nrecyc = global_parameters%n_species
  Nrecomb= global_parameters%n_species
  NSTRATA=Nrecyc+Nrecomb+NPUFFS

  n_call_eir=0

   

end subroutine styx_setup_general_parameters
