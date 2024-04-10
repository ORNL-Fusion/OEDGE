subroutine set_test_soledge_drift_parameters(element,temperature)
  use all_variables
  implicit none
  character(len=2),intent(in) :: element
  real*8,intent(in) :: temperature
  !####################################################
  !soledge parameters
  flags%is_SLAB=.false.
  flags%is_to_the_centre=.false.
  flags%is_Pen=.false.
  flags%solve_temperature=.true.
  flags%turbulence_model=0
  flags%solve_phi=.true.
  flags%solve_drift=.true.
  flags%neutral_model=0
  flags%restart=.false.
  global_parameters%N_iterations=1
  global_parameters%CFL=0.2d0
  !####################################################
  !chemistry block
  global_parameters%N_species=1
  call allocate_species_properties(global_parameters%N_species)
  global_parameters%element_list(1)%symbol=element
  call check_element_list()
  call compute_ion_number()
  call affect_element_and_charge_to_species()
  !###################################################
  !BC
  boundary_conditions%BCn_model(1)=0
  boundary_conditions%BCn(1)=1.d19
  boundary_conditions%BCT_model(1)=0
  boundary_conditions%BCTe=temperature
  boundary_conditions%BCTi(1)=temperature
  !###################################################
  !Transport parameters
  transport_parameters%Dn_p(1)=1.D0
  transport_parameters%Dn_t(1)=0.D0
  transport_parameters%nu_p(1)=1.D0
  transport_parameters%nu_t(1)=0.D0
  transport_parameters%chie_p=1.D0
  transport_parameters%chie_t=0.D0
  transport_parameters%chii_p(1)=1.D0
  transport_parameters%chii_t(1)=0.D0
  transport_parameters%v_pinch(1)=0.D0
  ballooning_parameters%ballooning_model=0
  ballooning_parameters%zbal=0.D0
  ballooning_parameters%minmaxbal=1.D0
  ballooning_parameters%sigmabal=1.D0
  transport_parameters%Coulomb_log=12.D0
  transport_parameters%Flux_limiter=0.D0
  !###################################################
  !reference parameters
  reference_parameters%fields%n0=0.D0
  reference_parameters%fields%T0eV=0.D0
end subroutine set_test_soledge_drift_parameters
