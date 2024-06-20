subroutine reset_source_terms(zone)
#include "compile_opt.inc"
  use all_variables, only : global_parameters, flags
  use Mzone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: n
  do n=0,global_parameters%N_ions
     zone%species(n)%fluxes%fluxn=0.D0
     zone%species(n)%fluxes%fluxG=0.D0
     zone%species(n)%fluxes%fluxE=0.D0
     zone%species(n)%sources%Sn=0.D0
     zone%species(n)%sources%SG=0.D0
     zone%species(n)%sources%SE=0.D0
     zone%species(n)%sources%Volumic_sources_n=0.D0
     zone%species(n)%sources%Volumic_sources_G=0.D0
     zone%species(n)%sources%Volumic_sources_E=0.D0
  end do
#if VORTICITY_PASTIX == 1
     zone%electric_fields(1)%j_perp=0.D0
     zone%electric_fields(1)%j_parallel=0.D0
     zone%electric_fields(1)%j_para_adv_W=0.D0
     zone%electric_fields(1)%j_diff_W=0.D0
#endif
  if(flags%turbulence_model.eq.1) then
     zone%kepsilon(1)%Sk=0.D0
     zone%kepsilon(1)%Sepsilon=0.D0
  end if
end subroutine reset_source_terms
