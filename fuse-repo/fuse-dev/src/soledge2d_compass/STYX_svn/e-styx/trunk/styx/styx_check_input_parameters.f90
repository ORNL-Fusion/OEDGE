subroutine styx_check_input_parameters
  use all_variables, only : global_parameters
  use eirmod_precision
  use styx2eirene
  implicit none
  integer*4 :: last_active_strata,i,isp

  if (direct_coupling .and. sc_level == 3) then
    sc_level=2
    write(*,*) '************************************************************************'
    write(*,*) '*  WARNING: short cycling option 3 deactivated in direct coupling mode *'
    write(*,*) '*  sc_level = 2 enforced ...                                           *'
    write(*,*) '************************************************************************' 
  endif 


  ! check data for geometry plots (make sure last_active_strata does not appear anymore) 

  if (stratum_plot < 1) stratum_plot=1
  if (nhist_plot <1) nhist_plot=1

! calculate total puff rate (in terms of atoms = 2 x molecular flux)
  last_active_strata=Nrecyc+Nrecomb
  
  total_puff=0._dp
  do i=1,Npuffs
    total_puff=total_puff+puffs(Nrecyc+Nrecomb+i)%rate
    if (puffs(Nrecyc+Nrecomb+i)%rate > 0._dp) last_active_strata = Nrecyc+Nrecomb+i
  enddo

  if (stratum_plot > 1) then
    write(*,*) '************************************************************************'
    write(*,*) '*   Warning, geometry plot requested for stratum number > 1            *'  
    write(*,*) '*   this may deteriorate statistics for last iteration                 *'
    write(*,*) '*   because flux weighted allocation of CPU to stratum turned off      *'
    write(*,*) '************************************************************************'
    write(*,*)
  endif

  if ((stratum_plot > NSTRATA .and. (.not. timedep)) .or. & 
      (stratum_plot > NSTRATA+1 .and. timedep )) then
    write(*,*) '************************************************************************'
    write(*,*) '*   Warning, geometry plot requested for stratum number > NSTRATA      *'  
    write(*,*) '*   plot for last active gas puff stratum enforced ...                 *'
    write(*,*) '************************************************************************'
    write(*,*)  
    stratum_plot = last_active_strata
  endif

  if (stratum_plot > 2 .and. stratum_plot < NSTRATA+1) then
    if (puffs(stratum_plot)%rate == 0._dp) then
      write(*,*) '************************************************************************'
      write(*,*) '*   Warning, geometry plot requested for gas puff valve not in use     *'  
      write(*,*) '*   plot for last active gas puff stratum enforced ...                 *'
      write(*,*) '************************************************************************'
      write(*,*)
      stratum_plot = last_active_strata
    endif
  endif

  if (nhist_plot > min(Npart_Eirene(1),1000)) then
    write(*,*) '************************************************************************'
    write(*,*) '*   Warning, number of histories in geometry plot too large (> 1000)   *'     
    write(*,*) '*   Plotting 1000 histories ... 					   *'
    write(*,*) '************************************************************************'
  endif

  ! check that the puff locations are consistent with triangle properties

  do i=1,Npuffs
    if (IPROP(puffs(Nrecyc+Nrecomb+i)%iside,puffs(Nrecyc+Nrecomb+i)%itri) == 0 ) then
       write(*,*) ' Gas puff #',i,' on transparent surface, exit'
       write(*,*) ' triangle # = ',puffs(Nrecyc+Nrecomb+i)%itri,' side ',puffs(Nrecyc+Nrecomb+i)%iside
       call eirene_exit_own(1)
    elseif (IPROP(puffs(Nrecyc+Nrecomb+i)%iside,puffs(Nrecyc+Nrecomb+i)%itri) == 1 ) then
       write(*,*) ' Gas puff #',i,' on purely absorbing surface, exit'
       write(*,*) ' triangle # = ',puffs(Nrecyc+Nrecomb+i)%itri,' side ',puffs(Nrecyc+Nrecomb+i)%iside
       call eirene_exit_own(1)
    endif
    if (puffs(Nrecyc+Nrecomb+i)%itor <= 0 .or. puffs(Nrecyc+Nrecomb+i)%itor > Ntor_cells) then
       write(*,*) ' Toroidal location of Gas puff #',i,' out of range'
       write(*,*) ' Ntor_cells = ',Ntor_cells,' itor = ',puffs(Nrecyc+Nrecomb+i)%itor
    endif
  enddo

  ! check recycling coefficients

  do i=1,n_pumps
    if (any(pumps(i)%R > 1._dp)) then
  	write(*,*) ' Error in specification of recycling coefficient R (albedo) for pump #',i
        do isp=1,global_parameters%n_species
          write(*,*) 'R(',isp,')= ',pumps(i)%R(isp)
        enddo
  	call eirene_exit_own(1)
    endif
  enddo

  do i=1,n_pfc_types
    if (any(pfc_models(i)%R > 1._dp)) then
  	write(*,*) ' Error in specification of recycling coefficient R, for pfc type #',i
        do isp=1,global_parameters%n_species
          write(*,*) 'R(',isp,')= ',pfc_models(i)%R(isp)
        enddo
  	call eirene_exit_own(1)
    endif 
    if (pfc_models(i)%T > 3000._dp) then
      write(*,*) ' PFC are melting !!'
      write(*,*) ' Wall temperature set to T = ',pfc_models(i)%T,' K for pfc type #',i
      call eirene_exit_own(1)
    endif
  enddo

end subroutine styx_check_input_parameters
