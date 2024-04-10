subroutine styx_initialize_short_cycling_procedure
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use styx2eirene
  implicit none

  if (sc_level >= 2) then
    allocate(recflux_save(Nrecyc),puff_save(Npuffs))
    recflux_save=0._dp
    puff_save=0._dp
    ! nspami = natmi+nmoli+nioni
    allocate(scresc(nspami,5,NSTRATA))
    scresc=1._dp
  endif
  ! initialize short cycling procedure

  short_cycle=.false.
  enter_short_cycle=.false.  
  nscycles=0
  nrefresh=0

! initialization of global balances

  Nplasma_save=0._dp
  Eplasma_save=0._dp

end subroutine styx_initialize_short_cycling_procedure
