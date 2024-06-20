subroutine allocate_arrays_memory()
#include "compile_opt.inc"
  use all_variables
  implicit none
  integer*4 k
  do k=1,global_parameters%N_Zones
     zones(k)%number=k
     call allocate_metric_arrays(zones(k))
     allocate(zones(k)%species(0:global_parameters%N_ions))
     call broadcast_species_list(zones(k))
     call allocate_species_arrays(zones(k))
     call allocate_drift(zones(k)) 			!### Leybros modif ###
     call allocate_fluxes_and_sources(zones(k))
     call allocate_tridiagonal_systems(zones(k))
  end do
end subroutine allocate_arrays_memory



