function shall_we_run_eirene() result(run_eirene)
  use all_variables, only : flags
  implicit none
  logical :: run_eirene
  if(flags%neutral_model.eq.1) then
     run_eirene=.true.
  else
     run_eirene=.false.
  end if
end function shall_we_run_eirene
