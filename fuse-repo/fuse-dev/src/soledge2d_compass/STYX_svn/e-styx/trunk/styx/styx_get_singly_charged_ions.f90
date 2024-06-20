subroutine styx_get_singly_charged_ions
  use all_variables, only : global_parameters
  use eirmod_precision
  use styx2eirene
  implicit none
  integer*4 :: il,ipls

! get the list of singly charged ions

  allocate(singly_charged_ions(global_parameters%n_species))

  il=1
  do ipls=1,global_parameters%n_ions
    if (global_parameters%ions_list(ipls,2) == 1) then
      if (il>global_parameters%n_species) then
        write(*,*) 'Problem in infcop when looking for singly charged ions ...'
        stop
      endif
    singly_charged_ions(il)=ipls
    il=il+1
    endif
  enddo

end subroutine styx_get_singly_charged_ions
