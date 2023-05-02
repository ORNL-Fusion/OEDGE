module mod_sol29

  ! This module is an interface for SOL 29, the "blobby" SOL. This SOL
  ! is fundamentally different from all the previous options in that
  ! it semi-ditches the fluid approximation. This SOL launches blobs
  ! from the separatrix and tracks them through the SOL, recording
  ! their locations and their Monte Carlo weights. Temperatures and
  ! densities are then the weighted averaged of the recorded blobs
  ! in each cell. A number of gross assumptions go into this model:
  !
  !   - The blobs are assumed to travel in the parallel direction at the 
  !       background plasma velocity, which we use the simple SOL 
  !       approximation for at first (zero halfway between the targets, 
  !       linearly increases to the sound speed at the targets).
  !   - The radial velocity of the blobs are chosen from a generalized 
  !       gamma distribution. The radial velocity does not change from
  !       that point onwards.
  !   - The parallel length of the blob extends from some starting 
  !       location along the separatrix to the nearest X-point.
  !   - The radial width of the blob is specified and the properties are 
  !       assumed to decay from the center according to a Gaussian 
  !       distribution.
  !   - The peak temperatures and densities of the blobs are assumed to 
  !       decay in time according to specified the characteristic (e.g. 
  !       1/e) times.
  
  ! Need these globally to 
  
  
contains

  subroutine sol29()
    use mod_sol29_input
    use mod_cgeom
    use mod_comtor
    use mod_pindata
    use nc_utils_generic
    use allocate_arrays
    use mod_cioniz
    !use omp_lib   ! Until I get smarter.
    implicit none
    
    integer :: ir, ik, attempts, i, j, ikmid, m, irold,  perc, flag
    integer :: it, retcode, iitersol, iparam, ierr, it_copy, nholes
    integer :: ikold, ifate, iz
    real :: vr, max_prob, seed, targ_cs, mi, targ_te, s1, s2, slope
    real :: min_dist, s, r, z, cross, theta, adjust, ran, blob_time
    real :: pintim, blob_volume, blob_tot_elec, ckkmin, k, smax
    real :: tmp_blob_ne, tmp_blob_te, tmp_blob_pe, imp_frac_charge
    real :: max_blob_counts, getranf, blob_absfac
    real, allocatable :: tizs(:,:,:), sol29_absfac(:)
    logical :: debug, warn_flag
    character(len=12) :: blob_char
    character(len=50) :: pin_command
    
    external getranf
    
    write(0,*) 'Begin SOL29'
    
    ! Mass of D atom in eV/c^2.
    mi = 931.49e6
    
    ! These variables are needed to pass into the crossfield step
    ! function but they aren't actually used for our purposes, so we
    ! just set them to some garbage value. I do this explicitly up front
    ! just so I know what I'm doing.
    adjust = 0.0
    flag = 0
    debug = .false.
    ckkmin = HI
    k = 0.0
    
    ! Open an output file to store results.
    !open(69, file='sol29.csv', status='new')
    !write(69,*) 'attempts,vr'
    
    ! Estimate the max value of the generalized gamma distribution.
    max_prob = max_gengamma(vr_gamma_loc, vr_gamma_scale, &
      vr_gamma_a, vr_gamma_c)
    
    ! Unfortunately a lot of subroutines would need to be modified to
    ! get seed from rundiv.f --> mod_sol29, so in order to be as least
    ! intrusive as possible we will just set seed to a random number.
    ! Later we can do the work to get seed from rundiv.f.
    !call surand2(4.0, 1, seed)
    
    ! If a previous DIVIMP run with impurities has been specified then
    ! load in the netCDF file. Open file, allocate array for TIZS, read.
    if (load_divimp.eq.1) then
      load_divimp_path = '/fusion/projects/codes/oedge/zamperinis/results/d3d-sic-wall-allgrap-009.nc'
      write(0,*) 'load_divimp_path = ',load_divimp_path
      ierr = open_nc_file(load_divimp_path, NC_READONLY)
      write(0,*) 'After loading DIVIMP run, ierr = ',ierr
      
      ! Get number of charge states CION to set dimension correctly.
      ! ...
      
      ! For some reason this gives a netcdf error (-38), so we have to 
      ! rely on the cion given in the input file.
      !ierr = read_nc('CION', cion)
      !write(0,*) 'After reading CION, ierr = ',ierr
      
      !write(0,*) 'cion = ', cion
      call allocate_array(tizs, cion+2, size(ktebs,1), size(ktebs, 2), 'tizs', ierr)
      !ierr = read_nc('TIZS', tizs)
      !write(0,*) 'After reading TIZS, ierr = ',ierr
      ierr = close_nc_file()
      
      do ir=1, nrs
        do ik=1, nks(ir)
          !write(0,*) tizs(4,ik,ir)
        enddo
      enddo
      
    endif
    
    ! Let's just assume a constant blob volume in the shape of a
    ! cylinder. Then calculate the total amount of electrons contained
    ! in each blob at launch.
    blob_volume = 3.1415 * blob_radius**2 * blob_length
    write(0,*) 'blob volume: ',blob_volume
    write(0,*) 'blob_ne: ', blob_ne
    blob_tot_elec = blob_volume * blob_ne
    write(0,*) 'blob_tot_elec: ',blob_tot_elec
    
    ! Calculate the absolute scaling factor. Value is assigned for each
    ! knot on the last core ring (where the blobs start). 
    allocate(sol29_absfac(nks(irsep-1)))
    call get_absfac(irsep-1, 0.001, blob_freq, sol29_absfac)
    
    ! I don't have any idea how the PFZ should be handled, so for now
    ! just fill them with a constant Te and ne. 
    do ir=irwall, nrs
      do ik=1, nks(ir)
        ktebs(ik,ir) = pfz_te
        ktibs(ik,ir) = ktebs(ik,ir)
        knbs(ik,ir) = pfz_ne
      enddo
    enddo
    
    ! Each iteration refines the SOL parameters to the blob-averaged 
    ! values.
    pinion = 0.0
    do i=1, niterations
    
      write(0,*) ' Iteration: ',i
      
      ! First iteration we set the targets to the seed values.
      if (i.eq.1) then
        !do ir=irsep, nrs
        !  ktebs(1,ir) = seed_targ_te
        !  ktebs(nks(ir),ir) = seed_targ_te
        !enddo
        
        do m=1, nds
          ktebs(ikds(m),irds(m)) = seed_targ_te
        enddo
        
      endif
        
      do ir=irsep, nrs
        ikmid = ikmids(ir) + 1
        
        do ik=1, nks(ir)
        
          ! We need to address if Te has ended up NaN (e.g., no blobs
          ! were found at this target element, and thus 0/0). If it is,
          ! then choose the nearest Te value as a fill-in for this
          ! iteration. This only happens after the first iteration.
          if (isnan(ktebs(1,ir)).or.ktebs(1,ir).eq.0) then
            write(0,*) 'NaN/0 at ik,ir = ',1,ir
            
            ! Using target index, map back to knot, ring and then
            ! assign the Te value. We continue to update until we
            ! try every option, making sure we got the closest 
            ! value.
            if (targ_te_fix.eq.0) then
              min_dist = 1e6
              do m=1,ndsin
                if (abs(sepdist2(m)-sepdist2(idds(ir,1))).lt.min_dist) then
                  if (.not.isnan(ktebs(ikds(m),irds(m))).and. &
                    ktebs(ikds(m),irds(m)).gt.0) then
                    !if (ikds(m).ne.1) then
                    !  write(0,*) ' Warning: ikds(m) should equal 1!'
                    !  write(0,*) ' m, ikds = ',m, ikds(m)
                    !endif
                    min_dist = abs(sepdist2(m)-sepdist2(idds(ir,1)))
                    ktebs(1,ir) = ktebs(ikds(m),irds(m))
                    ktibs(1,ir) = ktebs(1,ir)
                    !write(0,*) 'm, ir, ik = ',m,ir,ik,' -> ',ktebs(1,ir)
                  endif
                endif
              enddo
            
            ! Alternatively, set Te to the MAX value of this ring. 
            ! Practically, this is just to nudge the simulation towards
            ! higher speeds and thus encouraging deposition on this
            ! target element so we can get real values recorded.
            elseif (targ_te_fix.eq.1) then
              ktebs(1,ir) = maxval(ktebs(:,ir))
              ktibs(1,ir) = maxval(ktibs(:,ir))
            endif
            
          endif
          
          ! Repeat for other target.
          if (isnan(ktebs(nks(ir),ir)).or.ktebs(nks(ir),ir).eq.0) then
            write(0,*) 'NaN/0 at ik,ir = ',nks(ir),ir
            
            if (targ_te_fix.eq.0) then
              min_dist = 1e6
              do m=ndsin,nds
                if (abs(sepdist2(m)-sepdist2(idds(ir,2))).lt.min_dist) then
                  if (.not.isnan(ktebs(ikds(m),irds(m))).and. &
                    ktebs(ikds(m),irds(m)).gt.0) then
                    !if (ikds(m).ne.1) then
                    !  write(0,*) ' Warning: ikds(m) should equal nks(ir)!'
                    !  write(0,*) ' m, ikds, nks = ',m, ikds(m), nks(ir)
                    !endif
                    min_dist = abs(sepdist2(m)-sepdist2(idds(ir,1)))
                    ktebs(nks(ir),ir) = ktebs(ikds(m),irds(m))
                    ktibs(nks(ir),ir) = ktebs(nks(ir),ir)
                    !write(0,*) ' m, ir, ik, =',m,ir,ik,' ->',ktebs(nks(ir),ir)
                  endif
                endif
              enddo
            
            elseif (targ_te_fix.eq.1) then
              ktebs(nks(ir),ir) = maxval(ktebs(:,ir))
              ktibs(nks(ir),ir) = maxval(ktibs(:,ir))
            endif
          endif
        
          ! Setup the parallel flows to linearly increase to the sound
          ! speed at each target and zero at the midpoint. Point-slope.
          if (ik.lt.ikmid) then
            targ_te = ktebs(1,ir)
            s1 = 0.0
            s2 = kss(ikmid,ir)
            targ_cs = -sqrt(2.0 * targ_te / mi) * 3e8
          else
            targ_te = ktebs(nks(ir),ir)
            s2 = kss(ikmid,ir)
            s1 = kss(nks(ir),ir)
            targ_cs = sqrt(2.0 * targ_te / mi) * 3e8
          endif
          !write(0,*) 'ir,ik,ktebs,targ_cs,s1,s2 = ',ir,ik,ktebs(ik,ir),targ_cs,s1,s2
          slope = (0.0-targ_cs) / (s2-s1)
          kvhs(ik,ir) = slope * (kss(ik,ir)-s1) + targ_cs
          
        end do
      end do 
      
      ! Option to include a constant impurity fraction in the background
      ! calculation results. Requires IZTAU to have been ran for this 
      ! impurity to fill out the ionization rate (KFIZS) array. Needs to
      ! be run every iteration as we update KTEBS and KNBS (except the 
      ! first one because we haven't updated them yet).
      if ((include_imp.eq.1).and.(i.gt.1)) then
        write(0,*) 'Calling iztau... cion = ',cion
        
        ! Just set nizs = cion here.
        call iztau (crmi, crmb, cion, rizb, ciopta, cprint, cion)
        kfizs_sum = 0.0
        do iz=0, cion-1
          do ir=1, nrs
            do ik=1, nks(ir)
              kfizs_sum(ik,ir) = kfizs_sum(ik,ir) + kfizs(ik,ir,iz)
            enddo
          enddo
        enddo
      endif

      ! Initialize statistics arrays.
      blob_counts = 0
      blob_counts_time = 0
      ne_weights = 0.0
      te_weights = 0.0
      ne_neuts = 0.0
      pe_weights = 0.0
      ne_imps = 0.0
      

!!$OMP PARALLEL DO    

      ! Track one blob at a time.
      do j=1, nblobs
        
        ! Copy update line from div6.f.
        if ((nblobs/10).gt.0) then 
          if (mod(j,nblobs/10).eq.0) then 
            perc = int((j*10)/(nblobs/10))
            write (blob_char, '(i8)') j
            write(0,'(a,i3,a,a)') &
              'Following Blobs: ',perc,' % complete. Blob # = ', &
               adjustl(blob_char)
          endif
        endif
        
        ! For the starting location, we assume the radial starting
        ! location has a max probability at the separatrix and 
        ! exponentially decays out with R-Rsep@OMP. When tau < 0 we
        ! take that to mean the last core ring all the time.
        if (tau_rad_start.le.0) then
          ir = irsep - 1
        else
          write(0,*) 'tau_rad_start != -1 is obsolete!'
          write(0,*) 'Stopping'
          stop
          ir = find_start_ring(seed, tau_rad_start)
        endif
        
        ! Choose starting location along ring from between X-points 
        ! (e.g. the last core ring). For now, just uniform. Future 
        ! upgrade to specify a normal distribution centered at the 
        ! midplane.
        !call surand2(seed, 1, ran)
        ran = getranf()
        s = kss(nks(ir),ir) * ran
        do ik=1, nks(ir)
          if (s.ge.ksb(ik-1,ir).and.s.lt.ksb(ik,ir)) then
            exit
          endif
        enddo
        
        !ik = 0
        !write(0,*) 'before ir,ik,s:',ir,ik,s
        !call find_start_knot(ir, ik, s)
        !write(0,*) 'after ir,ik,s:',ir,ik,s
        
        !r = rs(ik,ir)
        !z = zs(ik,ir)
        cross = 0.0
        
        ! Pull random radial velocity from specified distribution.
        ! Gaussian distribution.
        if (blob_vr_type.eq.0) then
          vr = ran_gaussian(vr_gauss_loc, vr_gauss_scale)
          
        ! Generalized Gamma distribution.
        elseif (blob_vr_type.eq.1) then
          call ran_gengamma(seed, vr_gamma_loc, vr_gamma_scale, &
            vr_gamma_a, vr_gamma_c, vr, attempts, max_prob)
        endif
        
        ! Will probably remove this offset option...
        vr = vr + vr_offset
        
        ! Blob velocities should always be positive (experimentally 
        ! inward moving blobs have been detected though, so revisit 
        ! this as understanding improves).
        if (vr.le.0) then
          vr = 1.0
        endif
        !write(0,*) 'vr = ',vr
        
        ! TO-DO: Add option for ballooning transport approximation.
        
        ! Record starting point in statistics.
        blob_counts(ik,ir) = blob_counts(ik,ir) + 1
        te_weights(ik,ir) = te_weights(ik,ir) + blob_te
        ne_weights(ik,ir) = ne_weights(ik,ir) + blob_ne 
        pe_weights(ik,ir) = pe_weights(ik,ir) + (blob_te * blob_ne)
        
        ! Track particle until an exit condition.
        blob_absfac = sol29_absfac(ik)
        it = 1
        it_copy = 1
        cross = 0.0
        blob_time = 0.0
        warn_flag = .false.
        do
      
          ! Not needed anymore with update_parallel_sol29?
          ! Set theta value for non-orthogonal transport. Calculate 
          ! theta is in ion_parallel_transport.f. Theta only needed
          ! for these non-orthogonal options.
          !if (northopt.eq.1.or.northopt.eq.3) then
          !  call calculate_theta(ik, ir, s, theta)
          !endif
          
          ! Test if this blob is killed off (dispersed). Do this by
          ! rejection method from an exponential distribution.
          if (tau_life.gt.0) then
            ran = getranf()
            if (ran.gt.exp(-blob_time / tau_life)) then
              exit
            endif
          endif
          
          ! Perform crossfield step.
          cross = cross + (-vr * timestep)
          ! jk = ik
          irold = ir
          ikold = ik
          smax = ksmaxs(ir)
          !write(0,*) ' before: ir, ik, cross, s',ir,ik,cross,s         
          call update_crossfield_sol29(ik, ir, ikold, irold, s, &
            theta, cross, adjust, ckkmin, smax, k, debug, ifate)
          !write(0,*) ' after : ir, ik, cross, s',ir,ik,cross,s

          ! Perform parallel step. 
          s = s + kvhs(ik,ir) * timestep
          irold = ir
          ikold = ik
          call update_parallel_sol29(ikold, irold, ik, ir, s, theta, &
            adjust, cross, debug)
          !write(0,*) 'j, ir, ik, kvhs, s, ksmax= ', j, ir, ik, kvhs(ik,ir), s, ksmaxs(ir)
               
          ! Check if outside grid. I think this is how to check the 
          ! radial edges. The ends should just be if s > smax.
          if (ir.eq.irwall.or.ir.eq.irtrap.or.ir.eq.irwall2.or.&
            ir.eq.irtrap2.or.s.ge.ksmaxs(ir).or.s.le.0.0) then
            !write(0,*) '-gone-'
            if (warn_flag) then
              write(0,*) j,': final it = ',it_copy
            endif
            
            exit
          endif
          
          ! In time array, increment this cell by 1 if we haven't
          ! already done so for this spatial cell.
          !if ((it.ge.500).and.(.not.warn_flag)) then
          !  write(0,*) 'Error! Too many time steps. Increase timestep or number of time bins.'
            !write(0,*) 'blob_time = ', blob_time
          !  warn_flag = .true.
          !endif
          !blob_counts_time(it,ik,ir) = blob_counts_time(it,ik,ir) + 1
            
          ! Record counts and weights. See explanation after the loop
          ! for how ne and Te are calculated based off these weights.
          ! Do weights belong here? Not sure...
          tmp_blob_ne = exp(-blob_time / tau_ne) * blob_ne
          tmp_blob_te = exp(-blob_time / tau_te) * blob_te
          tmp_blob_pe = exp(-blob_time / tau_pe) * (blob_te * blob_ne)
          blob_counts(ik,ir) = blob_counts(ik,ir) + 1 
          ne_weights(ik,ir) = ne_weights(ik,ir) + tmp_blob_ne
          te_weights(ik,ir) = te_weights(ik,ir) + tmp_blob_te
          pe_weights(ik,ir) = pe_weights(ik,ir) + tmp_blob_pe
          
          ! The blobs have some parallel length, as specified above. For
          ! each location the length of blob add the weight to there too.
          ! To-do.
          
          ! The blobs have some radial width, here defined by the width
          ! variable which is used here as the standard deviation of
          ! a gaussian distribution in the radial direction.
          ! To-do.
          
          ! At each location, the blob will cause neutrals to ionize
          ! and thus increase the density at that location. PINION
          ! is returned from EIRENE and has units of ionizations/m3/s.
          ! So then during this time step this means an additional
          ! ne_neut = PINION(ik,ir) * timestep electrons/m3 have been
          ! freed. Keep track of that. First iteration everything is 0.
          !ne_neuts(ik,ir) = ne_neuts(ik,ir) + pinion(ik,ir) * timestep
          
          ! Trying something else... see after loop for how this gets
          ! converted to a neutral density (divide by blob_freq is all).
          ne_neuts(ik,ir) = ne_neuts(ik,ir) + pinion(ik,ir) * blob_absfac 
          
          ! Can perform a similar trick with a preloaded DIVIMP run
          ! of an impurity, e.g., carbon walls.
          if (load_divimp.eq.1) then
          
            ! Do once for each charge state.
            !ne_imps(ik,ir) = ne_imps(ik,ir) + tiz(ik,ir) * timestep
          endif
          
          ! Include additional contributions from a specified constant
          ! fraction of an impurity ion.
          if ((include_imp.eq.1).and.(i.gt.1)) then
            
            ! Go through for each charge state.
            do iz=0, cion-1
            
              ! KFIZS is in units of s-1 (it's the ionization rate 
              ! coefficient [m3 s-1] times KNBS [m-3]). The contribution
              ! to the electron density due to ionization is then
              ! KFIZS * (density of charge state). 
              ! A subtely is that imp_frac is nz/ne for all charge 
              ! states summed together, and thus each charge state has 
              ! its own fraction that should add up to imp_frac. As an
              ! approximation, we assume that imp_frac_charge is the
              ! imp_frac weighted by the relative ionization rate.
              imp_frac_charge = (kfizs(ik,ir,iz) / kfizs_sum(ik,ir)) &
                * imp_frac
                
              ! iz+1 because iz=0 is the 0 -> +1 ionization rate.
              ne_imps(ik,ir) = ne_imps(ik,ir) + knbs(ik,ir) &
                * imp_frac_charge * (iz + 1)
                
              if ((nblobs/10).gt.0) then 
                if (mod(j,nblobs/10).eq.0) then 
                !write(0,*) ir,ik,iz,imp_frac_charge,knbs(ik,ir), &
                !  knbs(ik,ir)*imp_frac_charge*(iz+1)
                endif
              endif
            enddo
            
          endif
          
          blob_time = blob_time + timestep
          if (blob_time.gt.3) then
            write(0,*) 'Warning! blob_time exceeded max. Skipping to next blob'
            write(0,*) ' vr = ',vr
            write(0,*) ' ik, ir = ',ik,ir
            write(0,*) ' s = ',s
            exit
          endif
          
          ! Debugging.
          if (.not.warn_flag) then
            it = it + 1
          endif
          it_copy = it_copy + 1
          !if (loop_count.gt.50) exit
        
        ! End of individual blob tracking loop.
        end do
      
      ! End of all blob tracking loop.
      end do
!!$OMP END PARALLEL DO

      ! Now repeat the above except for inward directed holes.
      nholes = int(nblobs * frac_holes)
      do j=1, nholes
        
        ! Not sure how or if I want to implement this for holes so
        ! don't allow the program to continue.
        if (pres_mode.ne.0) then
          write(0,*) 'Error! Pressure conservation mode not yet enabled for holes.'
          write(0,*) 'Stopping'
          stop
        endif
        
        write(0,*) 'Error! Blob code needs to be re-written!'
        write(0,*) 'Stopping'
        stop
        
        ! Same as above.
        if ((nholes/10).gt.0) then 
          if (mod(j,nholes/10).eq.0) then 
            perc = int((j*10)/(nholes/10))
            write (blob_char, '(i8)') j
            write(0,'(a,i3,a,a)') &
              'Following Holes: ',perc,' % complete. Hole # = ', &
               adjustl(blob_char)
          endif
        endif
        
        ! Same as above.
        ir = find_start_ring(seed, tau_rad_start)
        
        ! Same as above.
        !call surand2(seed, 1, ran)
        ran = getranf()
        s = kss(nks(ir),ir) * ran
        
        ! Same as above.
        ik = 1
        do while (ik.lt.nks(ir).and.s.gt.kss(ik,ir))
          ik = ik + 1
        end do
        
        r = rs(ik,ir)
        z = zs(ik,ir)
        cross = 0.0
        
        ! For holes, the velocity is going inwards. Right now it is not
        ! clear if they would have a different type of velocity
        ! distribution since that has never been measured, so we just
        ! assume it's the negative of the blobs. 
        if (blob_vr_type.eq.0) then
          vr = ran_gaussian(vr_gauss_loc, vr_gauss_scale)
        elseif (blob_vr_type.eq.1) then
          call ran_gengamma(seed, vr_gamma_loc, vr_gamma_scale, &
            vr_gamma_a, vr_gamma_c, vr, attempts, max_prob)
        endif
        vr = -(vr + vr_offset)
        if (vr.ge.0) then
          vr = -1.0
        endif
        !write(0,*) 'vr = ',vr
        
        ! Record starting point in statistics.
        blob_counts(ik,ir) = blob_counts(ik,ir) + 1
        te_weights(ik,ir) = te_weights(ik,ir) + hole_te
        ne_weights(ik,ir) = ne_weights(ik,ir) + hole_ne 
        
        ! Track particle until an exit condition.
        it = 1
        it_copy = 1
        cross = 0.0
        blob_time = 0.0
        warn_flag = .false.
        do
      
          ! Same as above.
          if (northopt.eq.1.or.northopt.eq.3) then
            call calculate_theta(ik, ir, s, theta)
          endif

          ! Same as above.
          cross = cross + (-vr * timestep)
          !jk = ik
          irold = ir
          !write(0,*) ' before: ir, ik, cross, s',ir,ik,cross,s
          ! call do_cfstep(jk,ik,ir,irold,cross,adjust,theta,flag,debug)
          !write(0,*) ' after : ir, ik, cross, s',ir,ik,cross,s
      
          ! Same as above.
          s = s + kvhs(ik,ir) * timestep
          !write(0,*) 'j, ir, ik, kvhs, s ksmax= ', j, ir, ik, kvhs(ik,ir), s, ksmaxs(ir)
                    
          ! Since holes are going inwards, we also check for crossing 
          ! into the core, at which point we stop following them.
          if (ir.eq.irwall.or.ir.eq.irtrap.or.ir.eq.irwall2.or. &
            ir.eq.irtrap2.or.s.ge.ksmaxs(ir).or.(s.le.0.0) &
            .or.ir.lt.irsep) then
            !write(0,*) '-gone-'
            if (warn_flag) then
              write(0,*) j,': final it = ',it_copy
            endif
            
            exit
          endif
          
          ! In time array, increment this cell by 1 if we haven't
          ! already done so for this spatial cell.
          !if ((it.ge.500).and.(.not.warn_flag)) then
          !  write(0,*) 'Error! Too many time steps. Increase timestep or number of time bins.'
            !write(0,*) 'blob_time = ', blob_time
          !  warn_flag = .true.
          !endif
          !blob_counts_time(it,ik,ir) = blob_counts_time(it,ik,ir) + 1
            
          ! Record counts and weights. See explanation after the loop
          ! for how ne and Te are calculated based off these weights.
          blob_counts(ik,ir) = blob_counts(ik,ir) + 1 
          ne_weights(ik,ir) = ne_weights(ik,ir) + &
            exp(-blob_time / hole_tau_ne) * hole_ne
          te_weights(ik,ir) = te_weights(ik,ir) + &
            exp(-blob_time / hole_tau_te) * hole_te
          
          ! The blobs have some parallel length, as specified above. For
          ! each location the length of blob add the weight to there too.
          ! To-do.
          
          ! The blobs have some radial width, here defined by the width
          ! variable which is used here as the standard deviation of
          ! a gaussian distribution in the radial direction.
          ! To-do.
          
          ! At each location, the blob will cause neutrals to ionize
          ! and thus increase the density at that location. PINION
          ! is returned from EIRENE and has units of ionizations/m3/s.
          ! So then during this time step this means an additional
          ! ne_neut = PINION(ik,ir) * timestep electrons/m3 have been
          ! freed. Keep track of that. First iteration everything is 0.
          ne_neuts(ik,ir) = ne_neuts(ik,ir) + pinion(ik,ir) * timestep
          
          
          ! Can perform a similar trick with a preloaded DIVIMP run
          ! of an impurity, e.g., carbon walls.
          if (load_divimp.eq.1) then
          
            ! Do once for each charge state.
            !ne_imps(ik,ir) = ne_imps(ik,ir) + tiz(ik,ir) * timestep
          endif
          
          blob_time = blob_time + timestep
          if (blob_time.gt.3) then
            write(0,*) 'Warning! blob_time exceeded max. Skipping to next blob'
            write(0,*) ' vr = ',vr
            write(0,*) ' ik, ir = ',ik,ir
            write(0,*) ' s = ',s
            exit
          endif
          
          ! Debugging.
          if (.not.warn_flag) then
            it = it + 1
          endif
          it_copy = it_copy + 1
          !if (loop_count.gt.50) exit
        
        ! End of individual blob tracking loop.
        end do
      
      ! End of all blob tracking loop.
      end do
      

      ! After we've followed all the blobs, convert the weights into
      ! average blob values per cell. For Te, we are just doing a 
      ! weighted average (weights summed up above, so multiply by
      ! by blob_te and divide by number of blobs to get average Te per 
      ! blob).
      ! Density is a little trickier, since it is not a property of 
      ! blobs in each cell like Te, but how many on average are in the 
      ! cell. 
      ! Note: Only doing the SOL here. PFZ was assigned constant values 
      ! above.
      max_blob_counts = maxval(blob_counts)
      do ir=irsep, irwall
        do ik=1, nks(ir)
          if (blob_counts(ik,ir).eq.0) then
            knbs(ik,ir) = 0.0
            ktebs(ik,ir) = 0.0
            ktibs(ik,ir) = 0.0
          else
            knbs(ik,ir) = ne_weights(ik,ir) / blob_counts(ik,ir)
            
            ! Since we are launching uniformly along irsep, the average number
            ! of blobs launched in a cell is nblobs / nks(irsep). Scale
            ! the densities in each cell so they equal blob_ne at irsep
            ! and scale according to count. Incorrect, need to scale to
            ! max counts at irsep.
            !knbs(ik,ir) = ne_weights(ik,ir) / (nblobs / nks(irsep)) * blob_ne
            
            !knbs(ik,ir) = blob_freq * timestep / (nblobs * &
            !  kvols(ik,ir)) * ne_weights(ik,ir) * blob_tot_elec
            !write(0,*) 'blob_freq, timestep, blob_tot_elec, kvols, ne_weights, blob_counts, knbs:',blob_freq, timestep, blob_tot_elec, kvols(ik,ir), ne_weights(ik,ir), blob_counts(ik,ir), knbs(ik,ir)

            ! We calculate ne the same way as DDLIMS. We calculate how
            ! many particles a blob represents (at weight=1), and then
            ! we multiply
            ! dabs = timestep / (nblobs * kareas)  s/m2
            ! absfac = nelec  electrons/m-tor
            ! knbs = ne_weights * dabs * absfac * fblob   electrons/m3
            !blob_area = 3.1415 * blob_radius * blob_radius 
            !nelec = blob_ne * blob_area  ! electrons / m-tor
            
            !knbs(ik,ir) = ne_weights(ik,ir)  &
            !  * timestep / (nblobs * kareas(ik,ir)) &
            !  * 3.1415 * blob_radius * blob_radius * blob_ne &
            !  * fblob
            
            ! Still need to divide by blob_counts first.
            ne_neuts(ik,ir) = ne_neuts(ik,ir) / blob_counts(ik,ir) &
              * timestep**2
            ne_imps(ik,ir) = ne_imps(ik,ir) / blob_counts(ik,ir)
            
            ! Not useful...
            ! The input blob frequency is assumed to be that at the
            ! separatrix. As we track the blobs, we know that will in
            ! theory need to decrease proportionally to the relative 
            ! amount of blobs in each cell. That is to say, we should 
            ! scale the fblob in each cell by the normalized blob_counts.
            !ne_neuts(ik,ir) = ne_neuts(ik,ir) / blob_counts(ik,ir) / &
            !  (blob_freq * blob_counts(ik,ir) / max_blob_counts)
            
            ! Density is the time averaged number of blobs times the
            ! density of the blobs. 
            ! To-do.

            ! Add on the neutral contribution to ne.
            knbs(ik,ir) = knbs(ik,ir) + ne_neuts(ik,ir) + ne_imps(ik,ir)
            
            if (pres_mode.eq.0) then
              ktebs(ik,ir) = te_weights(ik,ir) / blob_counts(ik,ir)
              
            ! If tracking the pressure of blobs, then Te = pe / ne.
            elseif (pres_mode.eq.1) then
              if (ne_weights(ik,ir).eq.0) then
                ktebs(ik,ir) = 0.0
              else
              
                ! Don't divide by blob_counts because it cancels out.
                ktebs(ik,ir) = pe_weights(ik,ir) / ne_weights(ik,ir)
              endif
            endif
            
            ktibs(ik,ir) = ktebs(ik,ir)
                      
          endif
          !write(0,*) 'ir,ik,counts,te_weights,ktebs:',ir,ik,blob_counts(ik,ir),te_weights(ik,ir),ktebs(ik,ir)
        end do
      end do
      
      ! After we have finished our first iteration, run EIRENE to get 
      ! neutral data. The second iteration onwards will have actual
      ! ionization data to add to the plasma density. Also don't run 
      ! after the last iteration since we're finished after that.
      ! actpin: Command to run EIRENE (Input file tag H04). Only testing
      !           with reire07
      ! pimtim: Time spent in EIRENE
      ! retcode: Return code.
      if (runeir29.eq.1.and.i.lt.niterations) then
        !if (pincode.ne.5) then
        !  write(0,*) 'Error! Stick to EIRENE07 for now. EIRENE not run.'
        !else
          iitersol = 1
          iparam = -1
          write(0,*) 'Writing Eirene files...'
          call writeeirenefiles_06(iitersol)
          write(0,*) 'Done.'
          call pr_trace('SOL29:PINEXE','AFTER WRTPIN 4-5')
          write(pin_command,'(A,I5.3)') trim(actpin), iparam
          write(0,*) 'PIN_COMMAND:', trim(pin_command)
          call invokepin(pin_command, pintim, retcode)
          call pr_trace('SOL29:PINEXE','AFTER PIN 4-5')
          write(0     ,'(a,i6,a)') ' Return from EIRENE after ', &
            nint(pintim),' s'
          !write(pinout,'(a,i6,a)') ' Return from EIRENE after ', &
          !  nint(pintim),' s'
          call readeireneresults_06(iitersol)
        
        ! After running, load the ionization data.
        !...
        
        !endif
      endif
      
      ! Too much? Run DIVIMP for C sputtering to get the ionization
      ! rate for carbon (TIZ), apply to the above similarly. Is it
      ! as simple as this? It is not.
      !call div (title, equil, nizs, nimps, nimps2, cpulim, iontim, &
      !  neutim, seed, nymfs, iter, nrand)   
      
      
    ! End of iteration loop.  
    end do
    
    ! Assign the target elements values.
    do m=1,nds
      if (isnan(ktebs(ikds(m),irds(m)))) then
        kteds(m) = 0.0
        knds(m) = 0.0
        blob_counts_targ(m) = 0
      endif
      kteds(m) = ktebs(ikds(m),irds(m))
      knds(m) = knbs(ikds(m),irds(m))
      blob_counts_targ(m) = blob_counts(ikds(m), irds(m))
    enddo
    
    ! Warning should be issued if the other "Run PIN" switch it on,
    ! because this will overwrite the last EIRENE values.
    if (cpinopt.gt.1) then
      write(0,*) 'Warning! PIN is set to run again after SOL29. This '
      write(0,*) 'may be okay, but SOL29 will not be run with it.'
      write(0,*) 'Turn off H03 to supress this warning.'
    endif
    
  if (allocated(tizs)) deallocate(tizs)
  deallocate(sol29_absfac)
  
  
  
  !close(69)
  write(0,*) 'End SOL29'
  
  end subroutine sol29
  
  
  subroutine ran_gengamma(seed, loc, scal, a, c, vr, attempts, &
    max_prob)
    
    ! Rejection method. Pretty inefficient, but it is too difficult to
    ! solve for the CDF for the generalized gamma distribution and then
    ! to solve it again for "x" (i.e. the direct method). It's clearly
    ! been done before, but I am not familair enough with it all to 
    ! implement it myself.
    implicit none
    real :: bound_vr, ran, vr, a, c, loc, scal, y, keep_prob, test_prob
    real :: fact, max_prob, seed, getranf
    integer :: attempts
    
    external :: getranf
    
     
    ! For a bounding value, let's just use something large but not too
    ! large or else it will severely hurt performance. 50,000 m/s seems
    ! reasonable, blobs haven't been seen going faster than that.
    bound_vr = 50000.0
    
    ! Group all the constants of the generalized gamma together to save 
    ! computational time below.
    fact = abs(c) / gamma(a) / scal
    
    y = 0.0
    attempts = 0
    keep_prob = 0.0
    test_prob = 1.0
    do while (test_prob.gt.keep_prob)
    
      ! Uniformly pick a random vr in our domain.
      !call surand2(seed, 1, ran)
      ran = getranf()
      vr = bound_vr * ran
      attempts = attempts + 1
    
      ! Pick random probability between 0-max_prob. If test_prob < f(vr), 
      ! where f(vr) = the generalized gamma distribution computed at vr, 
      ! keep vr. Else, try again.
      !call surand2(seed, 1, ran)
      ran = getranf()
      test_prob = ran * max_prob
      y = (vr - loc) / scal
      if (y.lt.0) cycle
      keep_prob = y**(c*a - 1.0) * exp(-y**c) * fact
      !write(0,*) 'attempt, test_prob, gamma, keep_prob, vr: ',attempts,test_prob,gamma(a),keep_prob,vr
      !write(0,*) '  ',loc,y,abs(c), y**(c*a - 1.0),exp(-y**c),scal
      !if (attempts.gt.10) then
      !  write(0,*) '--giving up--'
      !  exit
      !endif
    
    end do
    
  end subroutine ran_gengamma
  
  function max_gengamma(loc, scal, a, c)
    
    ! Estimate the max of the generalized gamma distribution. Doing this
    ! via brute force with 1 m/s resolution, but experimentally this is
    ! fine since 1 m/s is considered very fine resolution for blob 
    ! velocities.
    
    implicit none
    real :: loc, scal, a, c, max_gengamma, max_vr, vr
    integer :: i
    real, allocatable :: x(:), y(:)
    
    ! Generate a number of vr (x) bins with 1 m/s resolution and 
    ! calculate the generalized gamma distribution at each location (y).
    ! Hardcoding a maximum vr as 50,000 m/s since that probably captures
    ! the range of reasonable vr values. Bad coding practice though :(
    max_vr = 5e5  
    allocate(x(int(max_vr)))
    allocate(y(int(max_vr)))
    do i=1, size(x)
    
      ! i = vr since 1 m/s resolution
      x(i) = real(i)
      vr = (x(i) - loc) / scal 
      if (vr.lt.0) then
        y(i) = 0.0
      else
        y(i) = abs(c) * vr**(c*a - 1.0) * exp(-vr**c) / gamma(a) / scal
      endif

    end do

    ! Then simply just return the maximum.
    max_gengamma = maxval(y)
    
    deallocate(x)
    deallocate(y)
    
    return
  
  end function max_gengamma
  
  integer function find_start_ring(seed, tau)
  
    use mod_cgeom, only: middist, irsep, irwall
 
    real :: tau, ran, seed, start_r, tmp_dist, min_dist, getranf
    
    external :: getranf
    
    ! Since we are assuming the probability exponentially decays from
    ! the separatrix, we can just do direct inversion of an exponential
    ! (Monte Carlo thing) with 0 at the separatrix out to infinity. The
    ! result is:
    ! r = -tau * ln(1 - ran)
    !call surand2(seed, 1, ran)
    ran = getranf()
    start_r = -tau * log(1 - ran)
    
    ! start_r is the starting R-Rsep@OMP value. Find the nearest ring
    ! and return it.
    find_start_ring = irsep
    min_dist = 999.0
    do ir=irsep, irwall
      tmp_dist = abs(start_r - middist(ir, 2))
      !write(0,*) 'start_r, ir, middist(ir, 2) tmp_dist, min_dist = ',start_r, ir, middist(ir, 2), tmp_dist, min_dist
      if (tmp_dist.lt.min_dist) then
        find_start_ring = ir
        min_dist = tmp_dist
      endif
    enddo
    !write(0,*) 'start_r, find_start_ring = ', start_r, find_start_ring
    return
  
  end function find_start_ring
  
  function ran_gaussian(mean, stdev)
    
    ! Generate a random number from a specified Gaussian distribution.
    ! mean: Mean of the distribution.
    ! stdev: Standard deviation of the distribution. 
    
    real :: mean, stdev, u1, u2, pi, getranf
    external :: getranf
    
    pi = 3.14159265359
    
    ! Use Box-Muller algorithm. 
    u1 = getranf()
    u2 = getranf()
    ran_gaussian = stdev * sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2) &
      + mean
    return
  
  end function ran_gaussian
  
  
  ! I've copied this from ion_parallel_transport.f. Since it's a 
  ! subroutine and not a function, it uses the globally defines  ik, ir,
  ! etc. variables as used in DIVIMP. We have our own variables
  ! contained within SOL29 and thus need our own subroutine.
  !  - Removed core statistics keeping code
  !  - Replace old style if-loops with do-while loops
  subroutine update_parallel_sol29(ikold, irold, ik, ir, s, theta, &
    adjust, cross, debug)
   
    !use mod_params
    use mod_comtor
    use mod_cgeom
    !use mod_div1
    !use mod_div6
    !use mod_particle_specs
    implicit none
      
    integer, intent(inout) :: ikold, irold, ik, ir
    real, intent(inout) :: s, theta, adjust, cross
    logical, intent(in) :: debug
    
    
    ! Find nearest ik corresponding to distance s along contour ir
    ! Find theta value corresponding to distance s along contour ir
    ! Adjust cross field term for new distances between contours     
    ikold = ik
    irold = ir
    
    ! Find new cell    
    if (pdopt.eq.1) then
    
! 620     IF (IK.LT.NKS(IR).AND.S.GT.KSB(IK,IR)) THEN
!            IK = IK + 1
!            GOTO 620
!         ENDIF
         
      do while (ik.lt.nks(ir).and.s.gt.ksb(ik,ir))
          ik = ik + 1
      enddo
         
! 625     IF (IK.GT.1.AND.S.LT.KSB(IK-1,IR)) THEN
!            IK = IK - 1
!            GOTO 625
!         ENDIF
         
      do while (ik.gt.1.and.s.lt.ksb(ik-1,ir))
        ik = ik - 1
      enddo  
         
    elseif (pdopt.eq.0) then    
       
! 630     IF (IK.LT.NKS(IR).AND.S.GT.KSS(IK,IR)) THEN
!            IK = IK + 1
!            GOTO 630
!         ENDIF
         
      do while (ik.lt.nks(ir).and.s.gt.kss(ik,ir))
        ik = ik + 1
      enddo
            
! 635     IF (IK.GT.1.AND.S.LE.KSS(IK-1,IR)) THEN
!            IK = IK - 1
!            GOTO 635
!         ENDIF
         
      do while (ik.gt.1.and.s.le.kss(ik-1,ir))
        ik = ik - 1
      enddo
         
      if (ik.gt.1.and.s-kss(ik-1,ir).lt.kss(ik,ir)-s) ik = ik - 1   
      
    endif
     
    ! If using non-orthogonal transport - find theta value     
    call calculate_theta(ik, ir, s, theta) 

    ! Adjust cross-field term - if necessary     
    call adjust_cross(cross, adjust, ik, ir, ikold, irold, debug) 

  end subroutine update_parallel_sol29
  
  
  ! This subroutine has been copy/pasted from 
  ! ion_crossfield_transport.f to avoid potentially messing with any
  ! global DIVIMP variables, as well as some modifications removing
  ! code not relevant to SOL29. 
  !  - Removed core statistics keeping code
  !  - Removed mirror option
  !  - Removed target impact code, the subroutine does what it needs to
  !      by just setting s = 0 or s = smax. End of tracking conditions
  !      are handled int he main SOL29 loop.
  subroutine update_crossfield_sol29(ik, ir, ikold, irold, s, &
    theta, cross, adjust, ckkmin, smax, k, debug_all, ifate)
    
    !use mod_params
    use mod_comtor
    use mod_cgeom
    !use mod_crand
    implicit none

    integer, intent(inout) :: ik, ir, ikold, irold, ifate
    real, intent(inout) :: s, theta, adjust, cross, ckkmin, smax, k
    logical, intent(in) :: debug_all

    ! Local variables 
    real :: theta1
    integer :: flag, jk

    ! Move the particle across rings if the CROSS value
    ! has exceeded the distance to the next ring. Recalculate
    ! THETA if using NON-orthogonal transport
    call do_cfstep(jk, ik, ir, irold, cross, adjust, theta, flag, &
      debug_all)

    if (debug_all) write (6,1000) 'D3:', ik, ir, s, k, theta, smax, &
      cross, adjust, 0.0, 0.0, distin(ik,ir), -distout(ik,ir), &
      'UPDATED CROSS'


    ! Re-calculate S if the particle has changed rings.  
    
    ! There is a bug caused by adjusting S when entering 
    ! a virtual ring. The FP code does not properly deal with
    ! it - in order to fix this - try not adjusting S for 
    ! values in the virtual rings.
    
    if (ir.ne.irold.and.(ir.ne.irwall.and.ir.ne.irtrap)) then
      k = kks(ir)
      ckkmin = min (ckkmin, k)
      smax = ksmaxs(ir)
 
      ! ITER grid 
      if (cgridopt.eq.2) then
        if ((((ir.ge.irsep.and.ir.le.irwall2).or. &
          (ir.ge.irsep2.and.ir.le.irwall)).and. &
          (irold.ge.irtrap)).or. &
          (((irold.ge.irsep.and.irold.le.irwall2).or. &
          (irold.ge.irsep2.and.irold.le.irwall)).and. &
          (ir.ge.irtrap))) then
                  
          if (s.gt.kss(ikold,irold)) then
            s = kss(ik,ir) - (s-kss(ikold,irold)) * &
              (kbacds(ik,ir)/kfords(ikold,irold))
          else
            s = kss(ik,ir) + (kss(ikold,irold)-s) * &
              (kfords(ik,ir)/kbacds(ikold,irold))
          endif
          
        elseif (s.gt.kss(ikold,irold)) then
          s = kss(ik,ir) + (s-kss(ikold,irold)) * &
            (kfords(ik,ir)/kfords(ikold,irold))
        else
          s = kss(ik,ir) - (kss(ikold,irold)-s) * &
            (kbacds(ik,ir)/kbacds(ikold,irold))
        endif

      ! Orthogonal Transport 
      elseif (northopt.eq.0.or.northopt.eq.2) then
        if (s.gt.kss(ikold,irold)) then
          s = kss(ik,ir) + (s-kss(ikold,irold)) * &
            (kfords(ik,ir)/kfords(ikold,irold))
        else
          s = kss(ik,ir) - (kss(ikold,irold)-s) * &
            (kbacds(ik,ir)/kbacds(ikold,irold))
        endif


      ! Non-orthogonal Transport
      ! This block generates a new value for S after non-orthogonal 
      ! cross-field diffusion (that results in changing rings).
      elseif (northopt.eq.1.or.northopt.eq.3) then
        if (ir.lt.irsep.and.ik.eq.1.and.theta.lt.thetag(ik,ir)) then
          ik = nks(ir)

          ! Handle an error condition when IKOLD was also 1.
          ! The spacing from nks(ir) to nks(ir) -1 on core rings
          ! is the same as cell 1 to a mythical cell 0. since cell
          ! 1 and nks(ir) coincide. IKOLD should not be 1 for a 
          ! non-core 
          if (ikold.eq.1) then 
            theta = thetag(ik,ir)-(thetag(ik,ir)-thetag(ik-1,ir)) * &
              (thetag(ikold,irold) - theta) / &
              (thetag(nks(irold),irold) - thetag(nks(irold)-1,irold))

          ! Regular case
          else 
            theta = thetag(ik,ir) - &
              (thetag(ik,ir)       - thetag(ik-1,ir)      ) * &
              (thetag(ikold,irold) - theta                ) / &
              (thetag(ikold,irold) - thetag(ikold-1,irold)) 
          endif               


          if (debug_all) write (6,1000) 'D4:', ikold, irold, s, k, &
           theta, thetag(ik,ir), thetag(ik-1,ir), thetag(ikold,irold), &
           thetag(ikold-1,irold), 0.0, distin(ik,ir), -distout(ik,ir), &
           'UPDATED CROSS'

        endif


        if (debug_all) write (6,1000) 'D5:', ik, ir, s, k, theta, &
          smax, cross, adjust, 0.0, 0.0, distin(ik,ir), &
          -distout(ik,ir), 'UPDATED CROSS'

        ! Re-calculate THETA and S.

        ! First half of cell
        if (theta.lt.thetag(ik,ir)) then
          if (ik.eq.1) then
            if (ir.lt.irsep) then
              ik = nks(ir)             
              theta = theta + (thetag(ik,ir) - thetag(1,ir))

              theta1 = (thetag(ik,ir) - theta) / &
                (thetag(ik,ir) - thetag(ik-1,ir))
               
              s = kss(ik,ir) - kbacds(ik,ir) * theta1
            else
            
              ! Particle has struck target cross-field
              if (theta.le.thetat(idds(ir,2))) then 
                s = 0.0
                
              else
                s = kss(ik,ir) * &
                  (theta - thetat(idds(ir,2))) / &
                  (thetag(ik,ir) - thetat(idds(ir,2)))
              endif
            endif
          else
            theta1 = (thetag(ik,ir) - theta) / &
              (thetag(ik,ir) - thetag(ik-1,ir))
            
            s = kss(ik,ir) - kbacds(ik,ir) * theta1
          endif

        ! Particle in second half of cell. 
        else
          if( ik.eq.nks(ir) )then
            if (ir.lt.irsep) then
              ik = 1
              theta = theta - (thetag(nks(ir),ir) - thetag(1,ir))
  
               theta1 = (theta - thetag(ik,ir)) / &
                 (thetag(ik+1,ir) - thetag(ik,ir))
                  
                s = kss(ik,ir) + kfords(ik,ir) * theta1
            else
              if (theta.ge.thetat(idds(ir,1))) then
                s = smax
              else
                s = kss(ik,ir) + kfords(ik,ir) * &
                  (thetat(idds(ir,1)) - theta) / &
                  (thetat(idds(ir,1)) - thetag(ik,ir))
              endif
            endif
          else
            theta1 = (theta - thetag(ik,ir)) / &
              (thetag(ik+1,ir) - thetag(ik,ir))
               
            s = kss(ik,ir) + kfords(ik,ir) * theta1
          endif

        endif
    
      endif

      if (debug_all) write (6,1000) 'D6:', ik, ir, s, k, theta, &
        smax, cross, adjust, 0.0, 0.0, distin(ik,ir), &
        -distout(ik,ir), 'UPDATED CROSS'
      
    ! End of non-orthogonal transport.
    endif

      
    ! Adjust particle S value if required for the new ring 
    if (ir.lt.irsep) then 
      
      ! Leaving just in case SOL29 is modified to track blobs in core.
      ! Looping round main plasma contours
      if (s.lt.0.0) then      
        do while (s.lt.0.0)
          s = s + smax
        enddo      
      elseif (s.gt.smax) then      
        do while (s.gt.smax)
          s = s - smax
        enddo      
      endif

    ! Particle not in core 
    else

      ! Target impact code was here, but I removed it. I leave this 
      ! comment here just to make note of it.

    endif

    if (debug_all) write(6,'(a,i4,3(1x,g12.5))') 'UPD CROSS :', &
      ifate, s, theta, cross

    return

    ! Format statements
 1000 format(a,2i4,1p,10(g11.4),1x,a) 
  
      return
  end subroutine update_crossfield_sol29   
  
  subroutine find_start_knot(ir, ik, s)
    ! Find the starting knot and s value of a blob assuming a ballooning 
    ! nature.
  
    use mod_cgeom
    use mod_comtor
    implicit none
    
    integer, intent (in) :: ir
    integer, intent (out) :: ik
    real, intent (out) :: s
    
    real :: btotal, tmp_cdf_max, ran, getranf
    real, allocatable :: balloon_proxy(:), balloon_cdf(:)
    
    external :: getranf
    
    ! Allocate temporary storage and intialize.
    allocate(balloon_proxy(nks(ir)))
    allocate(balloon_cdf(nks(ir)))
    balloon_proxy = 0.0
    balloon_cdf = 0.0
    
    
    ! Fill in the ballooning proxy along the ring.
    do ik=1, nks(ir)
      btotal = sqrt(bts(ik,ir)**2 + (bts(ik,ir) * bratio(ik,ir))**2)
      balloon_proxy(ik) = (midplane_b(ir)**2) / (btotal**2)
    enddo
    !write(0,*) 'balloon_proxy: ',balloon_proxy
 
    ! Create CDF. For this ring, integrate the ballooning proxy 
    ! along s. 
    !write(0,*) 'ksb: ',ksb(:,ir)
    do ik=1, nks(ir)
      if (ik.eq.1) then
        balloon_cdf(ik) = 0.0
      else
      
        ! KSB starts indexing at 0. Ex: For ik=1 then the bin ranges are 
        ! ksb(ik-1, ir) to ksb(ik,ir)
        balloon_cdf(ik) = balloon_cdf(ik-1) + &
          (ksb(ik,ir) - ksb(ik-1,ir)) * &
          balloon_proxy(ik)
      endif
    enddo
    !write(0,*) 'cdf1: ',balloon_cdf
    
    ! Normalize to the total (the last CDF value) to create
    ! the normalized CDF.
    tmp_cdf_max = balloon_cdf(nks(ir))
    do ik=1, nks(ir)
      balloon_cdf(ik) = balloon_cdf(ik) / tmp_cdf_max
    enddo
    !write(0,*) 'cdf2: ',balloon_cdf
    
    ! Now sample a random number between 0-1, see which CDF bin it falls
    ! into to determine launch knot.
    ran = getranf()
    do ik=1, nks(ir)
    
      ! When condition is met we have our knot.
      if (ran.le.balloon_cdf(ik+1)) exit

    enddo
    
    ! Choose uniformally between this knot's s range.
    ran = getranf()
    s = ksb(ik,ir) + ran * (ksb(ik+1,ir) - ksb(ik,ir))
    
    ! Deallocate.
    deallocate(balloon_proxy)
    deallocate(balloon_cdf)
  
  end subroutine find_start_knot
  
  subroutine get_absfac(ir, probe_tip_radius, fblob_meas, &
    sol29_absfac)
    
    ! This function derives an "absolute factor" for the blobs similar 
    ! as to what is done for neutrals in DIVIMP. 
    
    use mod_cgeom
    use mod_comtor
    implicit none
  
    real, intent(in) :: probe_tip_radius, fblob_meas
    integer, intent(in) :: ir
    real, intent(out) :: sol29_absfac(:)
    
    integer :: ik
    real :: probe_area, pi, fblob_dens, max_proxy, btotal
    real, allocatable :: balloon_proxy(:), fblob_dens_prof(:)
    
    allocate(balloon_proxy(nks(ir)))
    !allocate(fblob_dens_prof(nks(ir)))
    balloon_proxy = 0.0
    !fblob_dens_prof = 0.0
    sol29_absfac = 0.0
    
    ! Calculate probe area.
    pi = 3.14159265359
    probe_area = pi * probe_tip_radius**2
    
    ! Calculate "blob frequency density" for lack of a better name.
    fblob_dens = fblob_meas / probe_area
    write(0,*) 'fblob_dens: ', fblob_dens
    
    ! Create a *poloidal* distribution of the frequency density
    ! according to a ballooning nature. First fill in the ballooning 
    ! proxy along the ring.
    do ik=1, nks(ir)
      btotal = sqrt(bts(ik,ir)**2 + (bts(ik,ir) * bratio(ik,ir))**2)
      balloon_proxy(ik) = (midplane_b(ir)**2) / (btotal**2)
    enddo
    
    ! Then we assume our calculate fblob_dens is the *maximum* value,
    ! (i.e., that it was measured at the OMP), and so the proxy above
    ! normalized and then multiplied by the frequency density.
    max_proxy = maxval(balloon_proxy)
    do ik=1, nks(ir)
      sol29_absfac(ik) = (balloon_proxy(ik) / max_proxy) * fblob_dens
    enddo
    write(0,*) 'sol29_absfac: ',sol29_absfac
    
    ! The absolute factor is then the total integral of the blob density
    ! profile. 
    !do ik=1, nks(ir)
    !  get_absfac = get_absfac + (kpb(ik,ir)-kpb(ik-1,ir)) &
    !    * fblob_dens_prof(ik)
    !enddo
    !write(0,*) 'get_absfac: ',get_absfac
    !return
  
  end subroutine get_absfac
  

end module mod_sol29
