subroutine add_RHS_neutrals_sources(zone)
  use all_variables,only : global_variables
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  !what is gained for the plasma is lost for the neutrals
  zone%neutrals%RHS=zone%neutrals%RHS-zone%species(1)%sources%Sn_n&
       +zone%neutrals%Sn_nn
end subroutine add_RHS_neutrals_sources
