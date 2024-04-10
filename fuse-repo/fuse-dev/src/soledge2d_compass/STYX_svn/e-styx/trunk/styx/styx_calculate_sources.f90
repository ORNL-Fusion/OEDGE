subroutine calculate_sources(ntps)
  use all_variables, only : global_parameters, interp_data2, global_variables, reference_parameters
  use eirmod_precision
  use styx2eirene
  use eirmod_cpes
  use eirmod_comprt, only : iunout
  use eirmod_coutau, only : FLUXT
#ifdef TK3X
  use tk3x_to_eirene
#endif

  implicit none

  integer, intent(in) :: ntps
  integer :: arg, ier, itor, itri, icell

  double precision :: t_bg_int,t_end_int
  double precision :: t_bg_int2,t_end_int2
  double precision :: omp_get_wtime

  real(dp) :: pi_
  real(dp) :: rhoL0

  logical :: refreshed

  real(dp), parameter :: elcha=1.6022e-19_dp
  real(dp), parameter :: amuakg=1.6606e-27_dp
  real(dp), parameter :: kB=1.3806e-23_dp     ! Boltzmann constant (J/K) 3807 or 3806 ? (check in EIRENE)

  include 'mpif.h'

  if (my_pe==0) then

     ! get iteration number for infcop
     curiter=ntps

     if (nscycles == n_short_cycles) short_cycle=.false.

     ! decide whether to call EIRENE or to short cycle
     ! all runs call EIRENE for the last iteration

     if (curiter == numiter) then

        arg=2

     elseif (n_short_cycles /= 0 .and. .not.short_cycle) then

        arg=2

     elseif (n_short_cycles /=0 .and. short_cycle) then

        arg=1

     elseif (n_short_cycles == 0) then

        arg =2

     endif

  endif !(my_pe==0)

  call mpi_barrier(mpi_comm_world,ier)


  call mpi_bcast(arg,1,mpi_integer,0,mpi_comm_world,ier)	

  select case(arg)

  case(1)

#ifdef S2D
     if (my_pe==0) then

        if (sc_level >= 2 .and. nrefresh == ns_refresh .and. nscycles< n_short_cycles) then

           nrefresh=0

           ! is the first call really needed in sc 2 ??
           call interpolate_plasma()
           call interpolate_wall_fluxes()

           call styx_short_cycle()

           refreshed=.true.

        else
           refreshed=.false.
        endif
        nrefresh=nrefresh+1
        nscycles=nscycles+1
     endif
     
#endif
#ifdef TK3X

        if (sc_level >= 2 .and. nrefresh == ns_refresh .and. nscycles< n_short_cycles) then

           nrefresh=0

           ! is the first call really needed in sc 2 ??
           call tk3xeirene_interpolate_plasma()
           call tk3xeirene_interpolate_fluxes()

           if (my_pe == 0) call styx_short_cycle()

           refreshed=.true.

        else
           refreshed=.false.
        endif
        nrefresh=nrefresh+1
        nscycles=nscycles+1
#endif

     t_bg_int=omp_get_wtime()
     t_end_int=omp_get_wtime()
     t_bg_int2=omp_get_wtime()
     t_end_int2=omp_get_wtime()

  case(2)

     ! full eirene run, first send styx2D data

     t_bg_int=omp_get_wtime()

     ! there we go ... all processes go there
     
#ifdef S2D
     if (my_pe == 0) then
        call interpolate_plasma()
        call interpolate_wall_fluxes()
     endif
#endif
#ifdef TK3X
        call tk3xeirene_interpolate_plasma()
        call tk3xeirene_interpolate_fluxes()
#endif

     t_bg_int=omp_get_wtime()

     call styx_run_eirene(ntps)

     if (my_pe == 0) then
        n_call_eir=n_call_eir+1
        short_cycle=.true.
     endif

     if (my_pe==0 .and. n_short_cycles > 0 .and. curiter < numiter) then
    	write(*,'(A80)') '--------------------------------------------------------------------------------'
     	write(*,'(A32,i2,A7,i5,A10)') '  Entering short cycling, level ',sc_level,' , for ',n_short_cycles,' steps ...'
        if (sc_level>1) write(*,'(A27,i4,A15)') '  Refreshing sources every ',ns_refresh,' time steps ...'                       
     	write(*,'(A80)') '--------------------------------------------------------------------------------'
     endif

     nscycles=1
     nrefresh=1
     refreshed=.true.

  end select
  
  ! recovers volumes of eirene triangles for TOKAM3X
#ifdef TK3X
  do itor=1,Ntor_cells
       do itri=1,Ntri_styx
         icell = itri + (itor-1)*(Ntri_styx+1)
         Interp_data2%tri_vol(itri,itor)      = vol_tri_eirene(icell)
       enddo
  enddo
#endif

  ! transform sources back to styx units


#ifdef S2D
  if (my_pe==0) then

     ! call this only if sources have been modified, i.e. after eirene call or after refreshing sources in sc_level>1

     if (refreshed) then
        pi_=4._dP*atan(1._dP)

   ! extraneous 1.5 removed from dat%T0*kb in energy terms 
        Interp_data2%tri_Sn=Interp_data2%tri_Sn/(reference_parameters%fields%n0/reference_parameters%fields%tau0)  
        Interp_data2%tri_SG=Interp_data2%tri_SG/(reference_parameters%fields%n0*reference_parameters%fields%c0/reference_parameters%fields%tau0)  
        Interp_data2%tri_SE=Interp_data2%tri_SE/(reference_parameters%fields%n0*reference_parameters%fields%T0*kb/reference_parameters%fields%tau0)

        Interp_data2%neutral_outflux(1)=FLUX_NEUTRAL_OUT/elcha

        t_bg_int2=omp_get_wtime()

        call interpolate_sources()
        call report_source_in_pen_cells_to_closest_plasma_cells()

        t_end_int2=omp_get_wtime()

        refreshed=.false. 

     endif
  endif
#endif
#ifdef TK3X
     ! call this only if sources have been modified, i.e. after eirene call or after refreshing sources in sc_level>1

     if (refreshed) then
        pi_=4._dP*atan(1._dP)

   ! extraneous 1.5 removed from dat%T0*kb in energy terms 
        Interp_data2%tri_Sn=Interp_data2%tri_Sn/(reference_parameters%fields%n0/reference_parameters%fields%tau0)  
        Interp_data2%tri_SG=Interp_data2%tri_SG/(reference_parameters%fields%n0*reference_parameters%fields%c0/reference_parameters%fields%tau0)  
        Interp_data2%tri_SE=Interp_data2%tri_SE/(reference_parameters%fields%n0*reference_parameters%fields%T0*kb/reference_parameters%fields%tau0)

!        Interp_data2%neutral_outflux(1)=FLUX_NEUTRAL_OUT/elcha
        
        rhoL0 = reference_parameters%fields%c0*reference_parameters%fields%tau0
            
        Interp_Data2%tri_vol = Interp_Data2%tri_vol / (rhoL0**3)

        Interp_data2%neutral_outflux = Interp_data2%neutral_outflux/(reference_parameters%fields%n0/reference_parameters%fields%tau0*rhoL0**3)

        t_bg_int2=omp_get_wtime()

        call tk3xeirene_interpolate_sources()

        t_end_int2=omp_get_wtime()

        refreshed=.false. 

     endif
#endif

  ! clean all the EIRENE stuff at last iteration
  if (curiter == numiter) call styx_deallocate_eirene_modules

end subroutine calculate_sources
