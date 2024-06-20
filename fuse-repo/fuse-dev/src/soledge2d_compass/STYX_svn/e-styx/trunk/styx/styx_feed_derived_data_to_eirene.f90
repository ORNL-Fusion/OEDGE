subroutine styx_feed_derived_data_to_eirene
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  implicit none

  call eirene_plasma_deriv(0)

  ! update rate coefficients
  ! hardwire AMJUEL1, ... and branch on setamd if non default is used

  !tst1 = omp_get_wtime()

  if (am_database == 1 .and. hardwired) then
    call styx_hardwired_amd
    !tend1 = omp_get_wtime()
    !write(*,*) ' time for hardwired am data (ms) = ', (tend1-tst1)*1e3_dp
  else
    call eirene_setamd(1)
  endif
   
end subroutine styx_feed_derived_data_to_eirene
