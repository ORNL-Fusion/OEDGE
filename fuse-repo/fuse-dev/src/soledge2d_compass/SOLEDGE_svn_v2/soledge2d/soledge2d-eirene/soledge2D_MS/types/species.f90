module Mspecies

  use Mphysics
  use Mglobal_variables
  use Mdrifts				!### Leybros Robin ###
  implicit none

  ! This module contains the variable for each specie
  ! Specie number 0 is reserved for electrons

  type :: Tcoupling_terms
     real*8,allocatable :: tau(:,:,:)
     real*8,allocatable :: R(:,:,:)
     real*8,allocatable :: Q(:,:,:)
     real*8,allocatable :: qu(:,:,:)
  end type Tcoupling_terms

  type :: Ttransport_perp
     real*8,allocatable :: D_p(:,:)
     real*8,allocatable :: nu_p(:,:)
     real*8,allocatable :: chi_p(:,:)
     real*8,allocatable :: zeta_p(:,:)
     real*8,allocatable :: D_t(:,:)
     real*8,allocatable :: nu_t(:,:)
     real*8,allocatable :: chi_t(:,:)
     real*8,allocatable :: zeta_t(:,:)
     real*8,allocatable :: v_pinch(:,:)
  end type Ttransport_perp
  
  type :: Ttransport_para
     real*8,allocatable :: kappa0(:,:)
     real*8,allocatable :: kappa(:,:)
     real*8,allocatable :: nu0(:,:)
     real*8,allocatable :: nu(:,:)
     real*8 :: gamma
     real*8 :: delta_e 
     real*8 :: flux_limiter
     real*8 :: flux_limiter_nu
  end type Ttransport_para

  type :: TVar_species
     real*8,allocatable :: density(:,:)
     real*8,allocatable :: velocity(:,:)
     real*8,allocatable :: temperature(:,:)
     real*8,allocatable :: Mach(:,:)
     real*8,allocatable :: Gamma(:,:)
     real*8,allocatable :: log_Lambda(:,:)
  end type TVar_species

  type :: TCorners
!    NW(1,1,1)  NE(1,2,1)       WN(1,1,2)  EN(1,2,2)
!    SW(2,1,1)  SE(2,2,1)       WS(2,1,2)  ES(2,2,2)
     real*8 :: density(2,2,2) 
     real*8 :: Gamma(2,2,2) 
     real*8 :: velocity(2,2,2) 
     real*8 :: temperature(2,2,2)
  end type TCorners

  type :: TAm_vars
     real*8,allocatable :: ionization_rate_coefficient(:,:)
     real*8,allocatable :: recombination_rate_coefficient(:,:)
     real*8,allocatable :: radiation_function_excitation(:,:)
     real*8,allocatable :: radiation_function_recombination_bremsstrahlung(:,:)
  end type TAm_vars

  type :: TSources
     real*8,allocatable :: Sn(:,:)
     real*8,allocatable :: SG(:,:)
     real*8,allocatable :: SE(:,:)
     real*8,allocatable :: Sn_am(:,:)
     real*8,allocatable :: SG_am(:,:)
     real*8,allocatable :: SE_am(:,:)
     real*8,allocatable :: rad(:,:)
     real*8,allocatable :: Sn_n(:,:)
     real*8,allocatable :: Sn_G(:,:)
     real*8,allocatable :: Sn_E(:,:)
     real*8,allocatable :: Volumic_sources_n(:,:)
     real*8,allocatable :: Volumic_sources_G(:,:)
     real*8,allocatable :: Volumic_sources_E(:,:)
  end type TSources

  type :: Tfluxes
     real*8,allocatable :: Fluxn(:,:,:)
     real*8,allocatable :: FluxG(:,:,:)
     real*8,allocatable :: FluxE(:,:,:)
  end type Tfluxes

  type :: TTridiagonal
     real*8,allocatable :: a(:,:)
     real*8,allocatable :: b(:,:)
     real*8,allocatable :: c(:,:)
     real*8,allocatable :: S(:,:)
  end type TTridiagonal

  type :: TImplicit_coefficients
     real*8,allocatable :: west_density(:,:)
     real*8,allocatable :: east_density(:,:)
     real*8,allocatable :: west_velocity(:,:)
     real*8,allocatable :: east_velocity(:,:)
     real*8,allocatable :: west_temperature(:,:)
     real*8,allocatable :: east_temperature(:,:)
  end type TImplicit_coefficients

  type :: Tpenalisation_memories
     real*8,allocatable :: alpham(:,:)
     real*8,allocatable :: alphap(:,:)
  end type Tpenalisation_memories

  type :: Tspecies
     type(Telement) :: element
     integer*4 :: element_index
     integer*4 :: charge
     type(Ttransport_perp) :: transport_perp
     type(Ttransport_para) :: transport_para
     type(TVar_species) :: var(2) !1=OLD 2=NEW
     type(TCorners) :: corners    ! STEP_OLD
     type(TCorners) :: corners2   ! STEP_NEW
     type(TSources) :: sources
     type(TFluxes) :: fluxes
     type(TTridiagonal) :: tridiag
     type(TImplicit_coefficients) :: implicit_coefs
     type(Tpenalisation_memories) :: penalisation_memories
     logical :: compute_ionization
     logical :: compute_recombination
     type(TAm_vars) :: am_vars
     type(Tcoupling_terms) :: coupling_terms
     type(Tresiduals) :: residuals
     type(Tdrifts) :: drifts			!### Leybros Robin ###
  end type Tspecies

end module Mspecies
