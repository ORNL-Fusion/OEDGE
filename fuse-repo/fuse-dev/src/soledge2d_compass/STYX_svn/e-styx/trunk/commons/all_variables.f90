! Global variables inherited from SOLEDGE necessary for the coupling with EIRENE
module all_variables
  use Mreference_parameters  
  use Mphysics
  use MInterpolation_types
  use Melement_variables
  use Mglobal_variables

  implicit none

  Type :: Tglobals
     integer*4 :: N_Zones
     integer*4 :: N_Megazones
     integer*4 :: N_species
     integer*4 :: N_ions
     integer*4 :: N_iterations
     integer*4 :: N_solve_phi
     real*8 :: CFL
     Type(Telement),allocatable :: element_list(:)
     integer*4,allocatable :: ions_list(:,:) ! first column:  element number
                                                  ! second column: charge number
     integer*4,allocatable :: ind_ion_0(:) !ion index of single charge ions - size=n_species
     integer*4 :: N_psi_surfaces
     integer*4 :: N_phi_iterations
  end type Tglobals

  
  type(Treference_parameters) :: reference_parameters
  type(Tglobals) :: global_parameters
  type(Telement_variables),allocatable :: element_variables(:)
  type(TInterpolated_data2) :: Interp_Data2
  type(Tglobal_variables) :: global_variables
  type(TIntegrals) :: globals

end module all_variables
