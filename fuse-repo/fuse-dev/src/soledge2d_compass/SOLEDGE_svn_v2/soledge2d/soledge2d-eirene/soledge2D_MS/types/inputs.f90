module Minputs

  use Mphysics
  implicit none

  Type :: Tflags
     logical :: is_SLAB
     logical :: is_Pen
     logical :: restart
     logical :: is_to_the_centre
     logical :: solve_temperature
     logical :: use_triangles
     logical :: non_zero_parallel_viscosity
     integer*4 :: neutral_model
     integer*4 :: turbulence_model
     logical :: radialFeedback
     logical :: use_weno
  end type Tflags

  type :: Tdrift_flags
     logical :: vorticity_restart
     logical :: solve_phi
     logical :: solve_drift
     integer*4 :: N_solve_phi
     integer*4 :: BC_type_phi
     logical :: jdiam
     logical :: jadvW
     logical :: jdiffW
     logical :: jExBW
     logical :: jgradBW
     logical :: use_gradB_radial
     logical :: use_gradB_poloidal
     logical :: use_ExB_radial
     logical :: use_ExB_poloidal
     real*8 :: eta_perp0
     real*8 :: Mach_lim
     real*8 :: Mach_lim_rad
     logical :: reverse_B
  end type Tdrift_flags

  Type :: Tglobals
     integer*4 :: N_Zones
     integer*4 :: N_Megazones
     integer*4 :: N_species
     integer*4 :: N_ions
     integer*4 :: N_iterations
     integer*4 :: N_save
     real*8 :: CFL
     Type(Telement),allocatable :: element_list(:)
     integer*4,allocatable :: ions_list(:,:) ! first column:  element number
                                                  ! second column: charge number
     integer*4,allocatable :: ind_ion_0(:) !ion index of single charge ions - size=n_species
     integer*4 :: N_psi_surfaces
     integer*4 :: N_phi_iterations
  end type Tglobals

  Type :: Tboundary_Conditions
     integer*4,allocatable :: BCn_model(:)
     integer*4,allocatable :: BCT_model(:)
     real*8,allocatable :: BCn(:)       ! Core Boundary condition for density (either value or flux)
     real*8 :: BCTe      ! Core Boundary condition for electron energy (either value or flux)
     real*8,allocatable :: BCTi(:)      ! Core Boundary condition for ions energy (either value or flux)
  end type Tboundary_Conditions

  Type :: TTransport_Parameters
     real*8,allocatable :: Dn_t(:),Dn_p(:) ! Reference crossfield mass diffusivities
     real*8,allocatable :: nu_t(:),nu_p(:) ! Reference crossfield mass viscosities
     real*8 :: chie_t,chie_p ! Reference crossfield electron heat conductivities
     real*8,allocatable :: chii_t(:),chii_p(:) ! Reference crossfield ions heat conductivities
     real*8 :: zeta_t,zeta_p ! Reference crossfield vorticity diffusivities
     real*8,allocatable :: v_pinch(:) ! Reference v_pinch
     real*8 :: Coulomb_log
     real*8 :: delta_e      ! secondary electron emission
     real*8 :: gamma_i      ! ion sheath transmission factor
     real*8,allocatable :: Flux_limiter(:)
     real*8,allocatable :: Flux_limiter_nu(:)
  end type TTransport_Parameters

  Type :: Tballooning
     integer*4 :: ballooning_model
     real*8 :: zbal      ! gaussian peak position
     real*8 :: minmaxbal ! ratio between gaussian minimum and maximum
     real*8 :: sigmabal  ! gaussian width
     real*8 :: core_weighted_integral ! int_core ball*dS
  end type Tballooning

end module Minputs
