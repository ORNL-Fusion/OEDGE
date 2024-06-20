subroutine styx_reinit_eirene
 use eirmod_precision
 use eirmod_parmmod
 use eirmod_cestim
 use eirmod_cupd

 use styx2eirene
 implicit none
 integer :: j,istra
 
 ! reinitialize arrays used for estimation (all processess)

  call eirene_init_cestim(2)
  call eirene_init_cupd

  ! initialize incidence angles diagnostic arrays

  do j=1,Nsou
    sheath1D(j)%alpha_V=0._dp
    sheath1D(j)%beta_V=0._dp
    sheath1D(j)%hit_V=0
  enddo

   ! reinitialize all output structures (for ALL mpi process, otherwise bad things occur when reducing)

  do istra=0,nstrata
    eov(istra)%pdena = 0._dp
    eov(istra)%pdenm = 0._dp
    eov(istra)%pdeni = 0._dp
    eov(istra)%edena = 0._dp
    eov(istra)%edenm = 0._dp
    eov(istra)%edeni = 0._dp
    eov(istra)%vxdena = 0._dp
    eov(istra)%vxdenm = 0._dp
    eov(istra)%vxdeni = 0._dp
    eov(istra)%vydena = 0._dp
    eov(istra)%vydenm = 0._dp
    eov(istra)%vydeni = 0._dp
    eov(istra)%vzdena = 0._dp
    eov(istra)%vzdenm = 0._dp
    eov(istra)%vzdeni = 0._dp
  enddo

  if (direct_coupling) then
    do istra=0,nstrata
      eov(istra)%pael = 0._dp
      eov(istra)%pmel = 0._dp
      eov(istra)%piel = 0._dp
      eov(istra)%papl = 0._dp
      eov(istra)%pmpl = 0._dp
      eov(istra)%pipl = 0._dp
      eov(istra)%eael = 0._dp
      eov(istra)%emel = 0._dp
      eov(istra)%eiel = 0._dp
      eov(istra)%eapl = 0._dp
      eov(istra)%empl = 0._dp
      eov(istra)%eipl = 0._dp
      eov(istra)%eaplr = 0._dp
      eov(istra)%emplr = 0._dp
      eov(istra)%eiplr = 0._dp
      eov(istra)%mapl = 0._dp
      eov(istra)%mmpl = 0._dp
      eov(istra)%mipl = 0._dp
    enddo
  endif

end subroutine styx_reinit_eirene
