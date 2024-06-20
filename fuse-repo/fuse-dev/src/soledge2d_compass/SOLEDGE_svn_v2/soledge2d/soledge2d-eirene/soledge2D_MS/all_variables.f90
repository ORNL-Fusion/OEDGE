module all_variables

#include "compile_opt.inc"

  use Mzone
  use Mreference_parameters
  use Minputs
  use Melement_variables
  use Mglobal_variables
  use Mflux_surface
  use MInterpolation_types
  use Mturbulence
  use Mplots
  
  implicit none
  type(Tzone),allocatable :: zones(:)
  type(Tmegazone),allocatable :: megazones(:) 
  type(Treference_parameters) :: reference_parameters
  type(Tflags) :: flags
  type(Tglobals) :: global_parameters
  type(Tdrift_flags) :: drift_flags
  type(TBoundary_Conditions) :: boundary_conditions
  type(TTransport_parameters) :: transport_parameters
  type(Tballooning) :: ballooning_parameters
  type(Telement_variables),allocatable :: element_variables(:)
  type(Tglobal_variables) :: global_variables
  type(Tpenalisation_parameters) :: penalisation_parameters
  type(Tflux_surface),allocatable :: flux_surfaces(:)
  type(TIntegrals) :: globals
  ! This interpolated data will be allocated 
  ! only if triangle mesh is provided
  type(TInterpolated_data2) :: Interp_Data2
  type(Tkepsilon_param) :: kepsilon_param
  integer*4 :: N_custom_plots
  Type(Tplot),allocatable :: custom_plots(:)

end module all_variables
