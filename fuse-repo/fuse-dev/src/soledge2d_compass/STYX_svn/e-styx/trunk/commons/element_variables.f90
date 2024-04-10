! Module inherited from SOLEDGE necessary for the coupling with EIRENE
module Melement_variables

  Type :: Telement_variables
     real*8,allocatable :: core_outflux(:) ! in particle/s
     real*8 :: total_flux_core_nonionized  ! in particle/s
  end type Telement_variables

end module Melement_variables
