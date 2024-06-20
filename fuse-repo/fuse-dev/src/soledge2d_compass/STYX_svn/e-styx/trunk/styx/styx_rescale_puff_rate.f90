subroutine styx_rescale_puff_rate
  use all_variables, only : interp_data2, global_variables
  use Meirene_vars
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_ccona
  use eirmod_comsou
  use styx2eirene
  implicit none
  integer :: ipuff,istra,istra0
  real*8, save :: time
  real*8 :: tmax,puff_HWHM

  if(eirene_vars%feedback.ne.0) then
    total_puff=0._dp
    ! set puff rate (feedback scheme)
    istra0=Nrecyc+Nrecomb
    do ipuff=1,Npuffs
      istra=istra0+ipuff
      write(*,*) '############################################################'
      write(*,*) 'flux update', Interp_Data2%Puff(ipuff), Interp_Data2%Puff(ipuff)*elcha
      Write(*,*) '############################################################'
      FLUX(istra) = Interp_Data2%Puff(ipuff)*elcha ! flux in Amp, default= puff_rate(ISTRA)*elcha
      ! update total puff for particle balance diagnostics
      total_puff = total_puff + Interp_Data2%Puff(1) !main ion puff
      !    soreni(istra)   = ... ! if energy of puffed particles also modified, default= T0_puff(ISTRA) 
    enddo
  end if

  ! special treatment for SMBI injection
  ! t         : current time (adim)
  ! tmax      : time where injection is max (adim)
  ! puff_HWHM : HWHM of the time trace (adim)

  tmax=400
  puff_HWHM=100

  time=time+global_variables%dt*dfloat(n_short_cycles)

  ! force turning off recycling/recombination strata to improve load balancing

  !FLUX(1)=0._dp
  !FLUX(2)=0._dp

  !istra0=Nrecyc+Nrecomb
  !do ipuff=1,Npuffs
  !  istra=istra0+ipuff
  !  FLUX(istra)=Puffs(istra)%rate*elcha*&
  !    exp(-((time-tmax)/puff_HWHM)**2)
  !enddo

  write(*,*) '############################################################'
 ! write(*,*) 'updating SMBI puff rate ... REMINDER : recy/recom strata OFF'
 ! write(*,'(A7,f5.1,A6)') ' time = ',time,' tau_0'
 ! write(*,'(A7,es8.2,A7)') ' rate = ',Flux(istra)/elcha,' part/s'
  write(*,*) '############################################################'

 ! if (Npuffs>0) write(765,'(2e14.7)') time,flux(istra)/elcha

end subroutine styx_rescale_puff_rate
