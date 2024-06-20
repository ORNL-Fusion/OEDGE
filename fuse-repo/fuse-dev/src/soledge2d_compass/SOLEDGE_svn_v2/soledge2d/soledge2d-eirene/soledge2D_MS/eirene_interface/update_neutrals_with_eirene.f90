subroutine update_neutrals_with_eirene(ntps)

#include "compile_opt.inc"

  use all_variables, only : global_parameters, Interp_data2, global_variables
  use Meirene_vars
  use Mfeedback_control
#if USE_EIRENE == 1
  use eirmod_cpes
  use styx2eirene
#endif
  implicit none
  integer*4,intent(in) :: ntps
  real*8 :: newPuff
  integer*4 :: n,nstrat0

  if (my_pe ==0) then
     do n=1,Npuffs
        Nstrat0=2*global_parameters%n_species
        Interp_data2%Puff(n)=Puffs(nstrat0+n)%rate ! all puffs
     end do
     ! update of puff number 1:
     if(eirene_vars%feedback.eq.1) then 
        call puff_feedback(Control_data,newPuff,error_data)
        Interp_data2%Puff(1)=newPuff
        call write_feedback(global_variables%tempus,newPuff,Control_data,ntps)
     end if
     if(eirene_vars%feedback.eq.2) then
        Interp_data2%Puff(1)=1.D15
     end if
  endif
  call calculate_sources(ntps)
  if (my_pe == 0) then
     call compute_total_ionization_source()
  endif
end subroutine update_neutrals_with_eirene
