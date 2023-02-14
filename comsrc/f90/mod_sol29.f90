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
  
contains

  subroutine sol29()
    use mod_solcommon
    use mod_cgeom
    use mod_comtor
    use mod_pindata
    use nc_utils_generic
    use allocate_arrays
    !use omp_lib   ! Until I get smarter.
    implicit none
    
    integer :: ir, ik, attempts, i, j, ikmid, m, irold,  perc, flag
    integer :: it, jk, retcode, iitersol, iparam, ierr, it_copy, nholes
    real :: vr, max_prob, seed, targ_cs, mi, targ_te, s1, s2, slope
    real :: min_dist, s, r, z, cross, theta, adjust, ran, blob_time
    real :: pintim, blob_volume, blob_tot_elec
    real, allocatable :: tiz(:,:)
    logical :: debug, warn_flag
    character(len=12) :: blob_char
    character(len=50) :: pin_command
    
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
    call surand2(4.0, 1, seed)
    
    ! If a previous DIVIMP run with impurities has been specified then
    ! load in the netCDF file. Open file, allocate array for TIZ, read.
    if (load_divimp.eq.1) then
      write(0,*) 'load_divimp_path = ',load_divimp_path
      ierr = open_nc_file(load_divimp_path)
      write(0,*) 'After loading DIVIMP run, ierr = ',ierr
      
      ! Get number of charge states CION to set dimension correctly.
      ! ...
      call allocate_array(tiz, size(ktebs,1), size(ktebs, 2), 'tiz', ierr)
      ierr = read_nc('TIZ', tiz)
      write(0,*) 'After reading TIZ, ierr = ',ierr
    endif
    
    ! Let's just assume a constant blob volume in the shape of a
    ! cylinder. Then calculate the total amount of electrons contained
    ! in each blob at launch.
    blob_volume = 3.1415 * blob_radius**2 * blob_length
    write(0,*) 'blob volume: ',blob_volume
    write(0,*) 'blob_ne: ', blob_ne
    blob_tot_elec = blob_volume * blob_ne
    write(0,*) 'blob_tot_elec: ',blob_tot_elec
    
    ! Each iteration refines the SOL parameters to the blob-averaged 
    ! values.
    pinion = 0.0
    do i=1, niterations
    
      write(0,*) ' Iteration: ',i
      
      ! First iteration we set the targets to the seed values.
      if (i.eq.1) then
        do ir=irsep, nrs
          ktebs(1,ir) = seed_targ_te
          ktebs(nks(ir),ir) = seed_targ_te
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
            min_dist = 1e6
            !do m=1,nds
            do m=1,ndsin
              if (abs(sepdist2(m)-sepdist2(idds(ir,1))).lt.min_dist) then
                if (.not.isnan(ktebs(ikds(m),irds(m))).and. &
                  ktebs(ikds(m),irds(m)).gt.0) then
                  if (ikds(m).ne.1) then
                    write(0,*) ' Warning: ikds(m) should equal 1!'
                    write(0,*) ' m, ikds = ',m, ikds(m)
                  endif
                  min_dist = abs(sepdist2(m)-sepdist2(idds(ir,1)))
                  ktebs(1,ir) = ktebs(ikds(m),irds(m))
                  ktibs(1,ir) = ktebs(1,ir)
                  !write(0,*) 'm, ir, ik = ',m,ir,ik,' -> ',ktebs(1,ir)
                endif
              endif
            enddo
          endif
          
          ! Repeat for other target.
          if (isnan(ktebs(nks(ir),ir)).or.ktebs(nks(ir),ir).eq.0) then
            write(0,*) 'NaN/0 at ik,ir = ',nks(ir),ir
            min_dist = 1e6
            !do m=1,nds
            do m=ndsin,nds
              if (abs(sepdist2(m)-sepdist2(idds(ir,2))).lt.min_dist) then
                if (.not.isnan(ktebs(ikds(m),irds(m))).and. &
                  ktebs(ikds(m),irds(m)).gt.0) then
                  if (ikds(m).ne.1) then
                    write(0,*) ' Warning: ikds(m) should equal nks(ir)!'
                    write(0,*) ' m, ikds, nks = ',m, ikds(m), nks(ir)
                  endif
                  min_dist = abs(sepdist2(m)-sepdist2(idds(ir,1)))
                  ktebs(nks(ir),ir) = ktebs(ikds(m),irds(m))
                  ktibs(nks(ir),ir) = ktebs(nks(ir),ir)
                  !write(0,*) ' m, ir, ik, =',m,ir,ik,' ->',ktebs(nks(ir),ir)
                endif
              endif
            enddo
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

      ! Initialize statistics arrays.
      blob_counts = 0
      blob_counts_time = 0
      ne_weights = 0.0
      te_weights = 0.0
      ne_neuts = 0.0

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
        ! exponentially decays out with R-Rsep@OMP.
        ir = find_start_ring(seed, tau_rad_start)
        
        ! Choose starting location along ring from between X-points 
        ! (e.g. the last core ring). For now, just uniform. Future 
        ! upgrade to specify a normal distribution centered at the 
        ! midplane.
        call surand2(seed, 1, ran)
        s = kss(nks(ir),ir) * ran
        
        ! Find nearest ik to starting s location.
        ik = 1
        do while (ik.lt.nks(ir).and.s.gt.kss(ik,ir))
          ik = ik + 1
        end do
        
        r = rs(ik,ir)
        z = zs(ik,ir)
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
        
        ! Record starting point in statistics.
        blob_counts(ik,ir) = blob_counts(ik,ir) + 1
        te_weights(ik,ir) = te_weights(ik,ir) + blob_te
        ne_weights(ik,ir) = ne_weights(ik,ir) + blob_ne 
        
        ! Track particle until an exit condition.
        it = 1
        it_copy = 1
        cross = 0.0
        blob_time = 0.0
        warn_flag = .false.
        do
      
          ! Set theta value for non-orthogonal transport. Calculate 
          ! theta is in ion_parallel_transport.f. Theta only needed
          ! for these non-orthogonal options.
          if (northopt.eq.1.or.northopt.eq.3) then
            call calculate_theta(ik, ir, s, theta)
          endif

          ! Current bug: This seems to fail to transport the blob into
          ! the extended grid region. Need this addressed first.
          ! Perform crossfield step. Luckily we can adapt the usual
          ! crossfield step modules used for ion transport.
          cross = cross + (-vr * timestep)
          jk = ik
          irold = ir
          !write(0,*) ' before: ir, ik, cross, s',ir,ik,cross,s
          call do_cfstep(jk,ik,ir,irold,cross,adjust,theta,flag,debug)
          !write(0,*) ' after : ir, ik, cross, s',ir,ik,cross,s
      
          ! Perform parallel step. I think the ik, ir update is handled
          ! when the next iteration comes around in do_cfstep.
          s = s + kvhs(ik,ir) * timestep
          !write(0,*) 'j, ir, ik, kvhs, s ksmax= ', j, ir, ik, kvhs(ik,ir), s, ksmaxs(ir)
                    
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
          blob_counts(ik,ir) = blob_counts(ik,ir) + 1 
          ne_weights(ik,ir) = ne_weights(ik,ir) + &
            exp(-blob_time / tau_ne) * blob_ne
          te_weights(ik,ir) = te_weights(ik,ir) + &
            exp(-blob_time / tau_te) * blob_te
          
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
!!$OMP END PARALLEL DO

      ! Now repeat the above except for inward directed holes.
      nholes = int(nblobs * frac_holes)
      do j=1, nholes
        
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
        call surand2(seed, 1, ran)
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
          jk = ik
          irold = ir
          !write(0,*) ' before: ir, ik, cross, s',ir,ik,cross,s
          call do_cfstep(jk,ik,ir,irold,cross,adjust,theta,flag,debug)
          !write(0,*) ' after : ir, ik, cross, s',ir,ik,cross,s
      
          ! Same as above.
          s = s + kvhs(ik,ir) * timestep
          !write(0,*) 'j, ir, ik, kvhs, s ksmax= ', j, ir, ik, kvhs(ik,ir), s, ksmaxs(ir)
                    
          ! Since holes are going inwards, we also check for crossing 
          ! into the core, at which point we stop following them.
          if (ir.eq.irwall.or.ir.eq.irtrap.or.ir.eq.irwall2.or.&
            ir.eq.irtrap2.or.s.ge.ksmaxs(ir).or.(s.le.0.0)&
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
      do ir=irsep, nrs
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

            ktebs(ik,ir) = te_weights(ik,ir) / blob_counts(ik,ir)
            ktibs(ik,ir) = ktebs(ik,ir)
            
            ! Still need to divide by blob_counts first.
            ne_neuts(ik,ir) = ne_neuts(ik,ir) / blob_counts(ik,ir)
            
            ! Density is the time averaged number of blobs times the
            ! density of the blobs. 
            ! To-do.

            ! Add on the neutral contribution to ne.
            knbs(ik,ir) = knbs(ik,ir) + ne_neuts(ik,ir)
            
            
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
    
  if (allocated(tiz)) deallocate(tiz)
  
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
    real :: fact, max_prob, seed
    integer :: attempts
    
     
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
      call surand2(seed, 1, ran)
      vr = bound_vr * ran
      attempts = attempts + 1
    
      ! Pick random probability between 0-max_prob. If test_prob < f(vr), 
      ! where f(vr) = the generalized gamma distribution computed at vr, 
      ! keep vr. Else, try again.
      call surand2(seed, 1, ran)
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
 
    real :: tau, ran, seed, start_r, tmp_dist, min_dist
    
    ! Since we are assuming the probability exponentially decays from
    ! the separatrix, we can just do direct inversion of an exponential
    ! (Monte Carlo thing) with 0 at the separatrix out to infinity. The
    ! result is:
    ! r = -tau * ln(1 - ran)
    call surand2(seed, 1, ran)
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
    
    real :: mean, stdev, u1, u2, pi
    
    pi = 3.14159265359
    
    ! Use Box-Muller algorithm. 
    u1 = getranf()
    u2 = getranf()
    ran_gaussian = stdev * sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2) &
      + mean
    return
  
  end function ran_gaussian

end module mod_sol29
