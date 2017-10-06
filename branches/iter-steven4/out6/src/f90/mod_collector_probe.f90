module mod_collector_probe
  use global_parameters

  private

  real :: impdens(maxizs+1,maxseg,3)
  real :: impflux(maxizs+1,maxseg,3)
  real :: slen(maxseg,3)
  real :: lcoll(maxseg)
  real :: tot_impdens(maxseg)

  real :: local_vals(maxseg,6)
  real :: local_outs(maxseg)
  real :: local_info(maxseg,8)

  public :: collector_probe,write_fp_main_density

  ! axis calculation
  real :: midplane_axis(maxnrs),rsep_out,rsep_in
  integer :: iouter,iinner

contains

  subroutine collector_probe(r1p,z1p,r2p,z2p,probe_diameter,dperp,axis_opt,iopt)
    use global_parameters
    use error_handling
    use debug_options
    use mod_cgeom
    use mod_comtor
    use mod_outcom
    use mod_fperiph
    use mod_fp_data
    implicit none
    real :: r1p,z1p,r2p,z2p,probe_diameter,dperp
    integer :: iopt,axis_opt

    !
    !Anyway, here is my design for a first order collector probe.
    !
    ! 1) Specify the probe location (R1,Z1) to (R2,Z2) and diameter d.
    ! 2) Calculate the intersections of this line with each relevant ring on the grid.
    ! 3) Use the ne and Te along the probe to estimate Cs, and Lcoll using the value of Dperp input for the impurity simulation ... along the probe line for each ring.
    ! 4) Using the Lcoll for each ring ... calculate the average upstream and downstream impurity densities within Lcoll of the probe by charge state.
    ! 5) Calculate nu and alpha for each ring and charge state.
    ! 6) Calculate Gamma_imp onto the upstream and downstream probe sides using the average nimp(Z) (summed over charge states) calculated in (4). This will give the impurity flux/sec for each ring of the grid along the probe length.
    ! 7) Leaving out all the sputtering effects the deposition would then be Cm= Gamma_imp * t where t is the exposure time. 
    !
    !
    ! Lcoll= d**2 * Cs / (8 * Dperp)
    !
    ! nu = 1.4e-12 * (ne * Z**2 / (T**3/2) * (mD/mIMP)
    !
    ! alpha = Z * kTe / (mIMP * Lcoll**2) + ne * Cs * Lcoll
    !
    ! Gamma_IMP = 1/4 * nIMP * nu * Lcoll * [ -1 + (1 + 8 *alpha / nu**2)]
    !
    !
    ! Axis options:
    ! 1 - PSIN
    ! 2 - R-Rsep_OMP (if available ... for rings that do not reach the OMP this is not available)
    ! 3 - R-Rsep_probe (R-Rsep at the probe location - measured along the probe assuming perpendicular insertion)
    !

    !
    !     Local Variables
    !
    integer :: in,ip,ik,ir,icnt,irref
    real :: sint,psin,pint
    real :: tmp,rsect,zsect
    real :: cs

    real :: ne, te,ti,vb,ef,tgrade,tgradi

    integer :: outunit
    logical :: exist
    character*1024 :: fname
    
    integer :: fp_reg,fp_in,fp_ik

    !
    call pr_trace('COLLECTOR_PROBE','START')

    call print_debug_data
    !
    call pr_trace('COLLECTOR_PROBE','AFTER PRINT_DEBUG_DATA')

    ! Initialization
    impdens = 0.0
    impflux = 0.0
    slen = 0.0

    !
    ! Calculate midplane axis
    !
    call calc_axis(midplane_axis,rsep_out,rsep_in,r1p,z1p,r2p,z2p,axis_opt)

    call pr_trace('COLLECTOR_PROBE','AFTER CALC_MIDPLANE_AXIS')

    !
    !     Loop through the main SOL rings starting at the separatrix and the
    !     outer target - find the first intersection with the specified LOS
    !     and extract/interpolate the values of Ne, Te, Ti and pressure.
    !
    !     Counter for number of intersections found
    !
    icnt = 0

    !
    ! Turn off fp calculations ... turn on later if fp was active
    !
    fp_reg = 0
    fp_in  = 0
    fp_ik  = 0
    !
    !     jdemod - allow for collection of data from core as well - might help with plots
    !
    do ir = 1,irwall-1
       !
       !      do ir = irsep,irwall-1
       !
       do ik = 1,nks(ir)

          call find_intsect(ik,ir,r1p,z1p,r2p,z2p,rsect,zsect,sint,pint,psin)

          if (sint.gt.0.0.and.midplane_axis(ir).ne.0.0) then
             !
             !              Increment counter
             !
             icnt = icnt + 1

             !write(0,'(a,3i4,1p,12(1x,g12.5))') 'OSM1:',ik,ir,icnt,rsect,zsect,sint,kss(ik,ir),kss2(ik,ir),ksb(ik-1,ir),ksb(ik,ir),midplane_axis(ir)

             !
             !              Assign axis value depending on option
             !
             local_info(icnt,1) = ik
             local_info(icnt,2) = ir
             local_info(icnt,3) = sint
             local_info(icnt,4) = rsect
             local_info(icnt,5) = zsect

             ! just use outer target as reference
             local_info(icnt,6) = psitarg(ir,outer_targid)
             !local_info(icnt,7) = middist(ir,outer_targid)
             local_info(icnt,7) = midplane_axis(ir)
             local_info(icnt,8) = pint

             ! distance from end of probe
             ! Default axis but could be changed using an option
             local_outs(icnt) = sqrt((rsect-r1p)**2+(zsect-z1p)**2)

             !
             !              Intersection greater than cell center
             !
             if (sint.ge.kss(ik,ir)) then
                !
                !                 Grid point at target
                !
                if (ik.eq.nks(ir)) then
                   !
                   !                    Ne
                   !
                   local_vals(icnt,1) = knbs(ik,ir)     +&
                        (knds(idds(ir,1))-knbs(ik,ir) ) *&
                        (sint-kss(ik,ir))/(ksmaxs(ir)-kss(ik,ir))
                   !
                   !                    Te
                   !

                   local_vals(icnt,2) = ktebs(ik,ir)     +&
                        (kteds(idds(ir,1))-ktebs(ik,ir)) *&
                        (sint-kss(ik,ir))/(ksmaxs(ir)-kss(ik,ir))
                   !
                   !                    Ti
                   !

                   local_vals(icnt,3) = ktibs(ik,ir)     +&
                        (ktids(idds(ir,1))-ktibs(ik,ir)) *&
                        (sint-kss(ik,ir))/(ksmaxs(ir)-kss(ik,ir))
                   !
                   !                    Vb
                   !
                   local_vals(icnt,5) = (kvhs(ik,ir)     +&
                        (kvds(idds(ir,1))-kvhs(ik,ir)) *&
                        (sint-kss(ik,ir))/(ksmaxs(ir)-kss(ik,ir)))&
                        /qtim
                   !
                   !                    Pressure
                   !
                   local_vals(icnt,4) = local_vals(icnt,1) * ( &
                        (local_vals(icnt,2)+local_vals(icnt,3)) *ech +&
                        crmb * amu * local_vals(icnt,5)**2)
                   !
                   !                 Grid points away from target
                   !
                else
                   !
                   !                    Ne
                   !
                   local_vals(icnt,1) = knbs(ik,ir)     +&
                        (knbs(ik+1,ir)-knbs(ik,ir)) *&
                        (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
                   !
                   !                    Te
                   !

                   local_vals(icnt,2) = ktebs(ik,ir)     +&
                        (ktebs(ik+1,ir)-ktebs(ik,ir)) *&
                        (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
                   !
                   !                    Ti
                   !

                   local_vals(icnt,3) = ktibs(ik,ir)     +&
                        (ktibs(ik+1,ir)-ktibs(ik,ir)) *&
                        (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
                   !
                   !                    Vb
                   !
                   local_vals(icnt,5) = (kvhs(ik,ir)     +&
                        (kvhs(ik+1,ir)-kvhs(ik,ir)) *&
                        (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir)))&
                        /qtim
                   !
                   !                    Pressure
                   !
                   local_vals(icnt,4) = local_vals(icnt,1) * (&
                        (local_vals(icnt,2)+local_vals(icnt,3)) *ech +&
                        crmb * amu * local_vals(icnt,5)**2)
                   !
                endif
                !
                !              Intersection less than cell center
                !

             else
                !
                !                 Grid points at the target
                !
                if (ik.eq.1) then

                   !
                   !                    Ne
                   !
                   local_vals(icnt,1) = knbs(ik,ir)     +&
                        (knds(idds(ir,2))-knbs(ik,ir)) *&
                        (kss(ik,ir)-sint)/(kss(ik,ir))
                   !
                   !                    Te
                   !

                   local_vals(icnt,2) = ktebs(ik,ir)     +&
                        (kteds(idds(ir,2))-ktebs(ik,ir)) *&
                        (kss(ik,ir)-sint)/(kss(ik,ir))
                   !
                   !                    Ti
                   !

                   local_vals(icnt,3) = ktibs(ik,ir)     +&
                        (ktids(idds(ir,2))-ktibs(ik,ir)) *&
                        (kss(ik,ir)-sint)/(kss(ik,ir))
                   !
                   !                    Vb
                   !
                   local_vals(icnt,5) = (kvhs(ik,ir)     +&
                        (kvds(idds(ir,2))-kvhs(ik,ir)) *&
                        (kss(ik,ir)-sint)/(kss(ik,ir)))&
                        /qtim
                   !
                   !                    Pressure
                   !
                   local_vals(icnt,4) = local_vals(icnt,1) * (&
                        (local_vals(icnt,2)+local_vals(icnt,3)) *ech +&
                        crmb * amu * local_vals(icnt,5)**2)
                   !
                   !                 Grid points away from the target
                   !
                else

                   !
                   !                    Ne
                   !
                   local_vals(icnt,1) = knbs(ik,ir)     +&
                        (knbs(ik-1,ir)-knbs(ik,ir)) *&
                        (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
                   !
                   !                    Te
                   !

                   local_vals(icnt,2) = ktebs(ik,ir)     +&
                        (ktebs(ik-1,ir)-ktebs(ik,ir)) *&
                        (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
                   !
                   !                    Ti
                   !

                   local_vals(icnt,3) = ktibs(ik,ir)     +&
                        (ktibs(ik-1,ir)-ktibs(ik,ir)) *&
                        (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
                   !
                   !                    Vb
                   !
                   local_vals(icnt,5) = (kvhs(ik,ir)     +&
                        (kvhs(ik-1,ir)-kvhs(ik,ir)) *&
                        (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir)))&
                        /qtim
                   !
                   !                    Pressure
                   !
                   local_vals(icnt,4) = local_vals(icnt,1) * (&
                        (local_vals(icnt,2)+local_vals(icnt,3)) *ech +&
                        crmb * amu * local_vals(icnt,5)**2)
                   !
                endif

             endif
             !
             !              With all the plasma data and the s_intersection ... calculate the 
             !              averaged local, upstream and downstream impurity densities
             !              upper is S+Lcoll, lower is S-Lcoll and local is +/- Lcoll/4
             !
             !
             ! should we use OEDGE calculated Ti when finding local Cs? Assume so for now
             !
             !              Save sound speed 
             !
             local_vals(icnt,6) = SQRT((local_vals(icnt,2)+local_vals(icnt,3)) *ech/(CRMB*amu))
             cs = local_vals(icnt,6)

             lcoll(icnt) = probe_diameter**2 * cs / (8.0 * dperp)

             !write(0,'(a,i8,1p,10(1x,g12.5))') 'DATA:',icnt,local_vals(icnt,1),local_vals(icnt,2),local_vals(icnt,3),local_vals(icnt,6),lcoll(icnt)

             ! 
             ! Loop through charge states and calculate the impurity density for the collector probe analysis
             ! Leave out absfac scaling for now
             !
             call calc_impdens(ir,sint,lcoll(icnt),icnt,fp_reg,fp_in,fp_ik)

             ! Calculate the collected impurity flux onto the upper and lower sides by charge state and total

             call calc_collector_probe_flux(icnt)

             !
             !              Exit inner loop - proceed with next ring
             !
             exit
             !
          endif

       end do

    end do


    call pr_trace('COLLECTOR_PROBE','AFTER REGULAR SOL')

    !
    ! Code to add the points from the FP 
    !

    !write(0,*) 'fpopt:',fpopt,fp_n_bins

    if (fpopt.eq.5.or.fpopt.eq.6) then 
       ! far periphery option 5 was in use so impurity data beyond the grid should be availble
       ! Add the points from the far periphery mesh to the profile. 
       ! This code will not work for double null grids

       fp_reg = fp_main

       !ir = irwall-1

       !
       ! Find the intersection with the irwall-1 ring at the grid edge. This is the reference
       ! ring for the initial fp implementation.
       !
       irref = irwall
       
       fp_ik = 0
       do ik = 1,nks(irref)
          ir = irins(ik,irref)

          call find_intsect(ik,r1p,z1p,r2p,z2p,rsect,zsect,sint,pint,psin)
          if (sint.gt.0) then 
             fp_ik = ik
             exit
          endif

       end do

       write(6,*) 'fp_ik:', fp_ik,fp_reg,fp_n_bins

       if (fp_ik.eq.0) then 
          call errmsg('Collector probe:','FP DATA: No intersection of probe with outermost grid ring: FP data not included')
       else

          do in = 1,fp_n_bins

          ! Check to see if the bin with the probe is flagged as being beyond the wall ... if it is then we are done here
          
             if (fp_grid_flag(fp_ik,in,fp_reg).ne.0) then
                write(6,*) 'FP_GRID_FLAG:',fp_ik,in,fp_reg,fp_grid_flag(fp_ik,in,fp_reg),fp_grid_dist(fp_ik,fp_reg)
                exit
             endif

             icnt = icnt + 1

             ! list local - use in instead of ir ... though all fp use the reference ring at the moment
             local_info(icnt,1) = fp_ik
             local_info(icnt,2) = in
             local_info(icnt,3) = sint

             ! offset both the rsect value and the Romp values by the distance into the fp from fp_grid_dist
             local_info(icnt,4) = rsect + fp_grid_dist(in,fp_reg)
             local_info(icnt,5) = zsect

             ! just use outer target as reference
             local_info(icnt,6) = 0.0
             !local_info(icnt,7) = middist(ir,outer_targid) + fp_grid_dist(in,fp_reg)
             local_info(icnt,7) = midplane_axis(ir) + fp_grid_dist(in,fp_reg)
             local_info(icnt,8) = pint

             call fp_get_plasma(fp_ik,in,fp_reg,ne,te,ti,vb,ef,tgrade,tgradi)

             ! plasma is not interpolated in fp ... just take cell values
             !
             !                    Ne
             !
             local_vals(icnt,1) = ne

             !
             !                    Te
             !

             local_vals(icnt,2) = te

             !
             !                    Ti
             !

             local_vals(icnt,3) = ti

             !
             !                    Vb
             !
             local_vals(icnt,5) = vb

             !
             !                    Pressure
             !
             local_vals(icnt,4) = local_vals(icnt,1) * ( &
                  (local_vals(icnt,2)+local_vals(icnt,3)) *ech +&
                  crmb * amu * local_vals(icnt,5)**2)


             ! use the modified rsect value
             local_outs(icnt) = sqrt((local_info(icnt,4)-r1p)**2+(zsect-z1p)**2)


             !
             !              With all the plasma data and the s_intersection ... calculate the 
             !              averaged local, upstream and downstream impurity densities
             !              upper is S+Lcoll, lower is S-Lcoll and local is +/- Lcoll/4
             !
             !
             ! should we use OEDGE calculated Ti when finding local Cs? Assume so for now
             !
             !              Save sound speed 
             !
             local_vals(icnt,6) = SQRT((local_vals(icnt,2)+local_vals(icnt,3)) *ech/(CRMB*amu))
             cs = local_vals(icnt,6)

             lcoll(icnt) = probe_diameter**2 * cs / (8.0 * dperp)

             write(6,'(a,i8,1p,10(1x,g12.5))') 'FP DATA:',icnt,local_vals(icnt,1),local_vals(icnt,2),local_vals(icnt,3),local_vals(icnt,6),lcoll(icnt),fp_grid_dist(in,fp_reg)

             ! 
             ! Loop through charge states and calculate the impurity density for the collector probe analysis
             ! Leave out absfac scaling for now
             !
             call calc_impdens(ir,sint,lcoll(icnt),icnt,fp_reg,in,fp_ik)

             ! Calculate the collected impurity flux onto the upper and lower sides by charge state and total

             call calc_collector_probe_flux(icnt)


          end do


       endif

    endif

    call pr_trace('COLLECTOR_PROBE','AFTER PERIPHERY')


    !
    !  now print out the data and plot it 
    !
    numthe = icnt

    if (outer_targid .eq.1) then 
       ELABS(1)='IDF IDF  '
       ELABS(2)='ODF ODF '
       ELABS(3)='CENTCENT'
    else
       ELABS(1)='ODF ODF '
       ELABS(2)='IDF  IDF  '
       ELABS(3)='CENTCENT'
    endif
    !
    !       Select normalized profiles no matter what value of IOPT was
    !       input. 
    !
    ngs=3

    if (axis_opt.eq.1) then 
       XLAB = '   PSIN'
    elseif (axis_opt.eq.2) then
       XLAB = '   R-Rsep (OMP-GRID [M])'
    elseif (axis_opt.eq.3) then 
       XLAB = '   R-Rsep (Along probe [M])'
    endif

    YLAB = '   IMPURITY FLUX (PART/M2/S)'
    REF  = 'COLLECTOR PROBE IMPURITY FLUX'
    write(anly,'(a,f9.5)')  'DIAM(M) =',probe_diameter
    write(plane,'(a,f9.5)') 'DP(M2/S)=',dperp
    nview = ' '

    do in = 1,numthe
       !touts(in) = local_outs(in)
       ! use the Romp estimate for the plot axis
       touts(in) = local_info(in,7)
       tvals(in,1) = impflux(maxizs+1,in,1)*absfac
       tvals(in,2) = impflux(maxizs+1,in,2)*absfac
       tvals(in,3) = impflux(maxizs+1,in,3)*absfac
       ! set twids to 1.0 for now ... better would be estimate of radial width of ring but no easy way to get that at the moment
       twids(in) = 1.0

       !write(0,'(a,i5,6(1x,g12.3))') 'CP:',in,touts(in),tvals(in,1),tvals(in,2),tvals(in,3)

    end do

    themin = touts(1)
    themax = touts(numthe)
    !
    !       Generate plot
    !
    CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs, &
         ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS, &
         JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE, &
         TABLE,IOPT,2,1.0,0)

    call pr_trace('COLLECTOR_PROBE','AFTER PLOT')

    !
    ! Write out all of the collector probe results to a separate file
    ! 


    call find_free_unit_number(outunit)

    fname = 'collector_probe.out'

    inquire(file=trim(fname),exist=exist)

    if (exist) then 
       open(unit=outunit,file=trim(fname),status='old',form='formatted',position='append')
    else
       open(unit=outunit,file=trim(fname),status='new',form='formatted')
    endif

    write(outunit,'(a)') 
    write(outunit,'(a)') ' COLLECTOR PROBE SUMMARY INFORMATION:'
    write(outunit,'(a,4(1x,f12.5))') ' LOCATION:', r1p,z1p,r2p,z2p
    write(outunit,'(a,4(1x,f12.5))') ' DIAMETER,DPERP:', probe_diameter,dperp
    write(outunit,'(a)') 

    write(outunit,'(a)') ' TOTAL FLUX AND DENSITY SUMMARY:'
    write(outunit,'(a,10(1x,g18.8))') ' ABSFAC:',absfac_neut,absfac
    write(outunit,'(31(1x,a))') 'INDEX','IK','IR','SINT','RSECT','ZSECT','PSI','ROMP','NE','TE','TI','CS','LCOLL','LEN1','LEN2','LEN3','DIST',&
               'IMPDENS_'//trim(elabs(1)(1:4)),'IMPDENS_'//trim(elabs(2)(1:4)),'IMPDENS_CENT','IMPFLUX_'//trim(elabs(1)(1:4)),'IMPFLUX_'//trim(elabs(2)(1:4)),'IMPFLUX_CENT','PINT','PMAX'
    do in = 1,numthe
       write(outunit,'(i8,30(1x,g12.5))') in,local_info(in,1),local_info(in,2),local_info(in,3),local_info(in,4),local_info(in,5),local_info(in,6),local_info(in,7),&
            local_vals(in,1),local_vals(in,2),local_vals(in,3),local_vals(in,6),&
            lcoll(in),slen(in,1),slen(in,2),slen(in,3),&
            local_outs(in),impdens(maxizs+1,in,1),impdens(maxizs+1,in,2),impdens(maxizs+1,in,3),&
            impflux(maxizs+1,in,1),impflux(maxizs+1,in,2),impflux(maxizs+1,in,3),local_info(in,8),kpmaxs(int(local_info(in,2)))
    end do
    write(outunit,'(a)') 

    close(outunit)

    call pr_trace('COLLECTOR_PROBE','END')



  end subroutine collector_probe




  subroutine calc_impdens(ir,sint,slcoll,icnt,fp_reg,fp_in,fp_ik)
    use global_parameters
    use mod_outcom
    use mod_dynam2
    use mod_cgeom
    use mod_comtor
    use mod_fp_data
    implicit none
    integer :: icnt,ir,fp_reg,fp_in,fp_ik
    real :: slcoll,sint

    ! Locals

    real :: smin,smax
    real :: smin_cent,smax_cent
    real :: slen_tmp
    real :: swall_start,swall_end
    integer :: ik

    !
    ! Need 3 average densities
    ! 
    ! 1)  Smin to sint ... lower density
    ! 2)  sint to smax ... upper density
    ! 3)  smin_cent to smax_cent ... central density
    !

    !
    ! Modify this so that the smin and smax extend at most 1/2 way
    ! to the nearest surface
    ! 


    !smin = max(sint-slcoll,0.0)
    !smax = min(sint+slcoll,ksmaxs(ir))

    if (fp_reg.eq.0) then 

       smin = max(sint-slcoll,sint/2.0)
       smax = min(sint+slcoll,(ksmaxs(ir)+sint)/2.0)
       ! arbitrarily specify the local region as +/- LCOLL/4
       smin_cent = max(sint-slcoll/4.0,sint/2.0)
       smax_cent = min(sint+slcoll/4.0,(ksmaxs(ir)+sint)/2.0)

    else
       ! for fp need to check fp_grid_flag for cells that may impinge on wall and then 
       ! set the min and max values to be either sint+/-slcoll or 1/2 the distance to the wall intersection
       ! Need the ik index of the intersecting cell for this process


       ! scan down
       swall_start = 0.0
       do ik = fp_ik,1,-1 
          ! looking for wall incursion onto grid going down
          if (fp_grid_flag(ik,fp_in,fp_reg).ne.0) then 
             swall_start = kss(ik,ir)
             exit
          endif
       end do

       ! scan up
       
       swall_end = ksmaxs(ir)
       do ik = fp_ik,nks(ir)
          ! looking for wall incursion onto grid going down
          if (fp_grid_flag(ik,fp_in,fp_reg).ne.0) then 
             swall_end = kss(ik,ir)
             exit
          endif
       end do

       smin = max(sint-slcoll,(sint+swall_start)/2.0)
       smax = min(sint+slcoll,(swall_end+sint)/2.0)
       ! arbitrarily specify the local region as +/- LCOLL/4
       smin_cent = max(sint-slcoll/4.0,(sint+swall_start)/2.0)
       smax_cent = min(sint+slcoll/4.0,(swall_end+sint)/2.0)

       write(6,'(a,10(1x,g12.5))') 'FP SMIN SMAX:',smin,smax,smin_cent,smax_cent,swall_start,swall_end

    endif


    call update_impdens(ir,icnt,smin,sint,1,fp_reg,fp_in)
    call update_impdens(ir,icnt,sint,smax,2,fp_reg,fp_in)
    call update_impdens(ir,icnt,smin_cent,smax_cent,3,fp_reg,fp_in)

  end subroutine calc_impdens


  subroutine update_impdens(ir,icnt,smin,smax,ireg,fp_reg,fp_in)
    use global_parameters
    use mod_cgeom
    use mod_outcom
    implicit none

    integer :: ireg,ir,icnt,ik,iz,fp_reg,fp_in
    real :: smin,smax
    real :: slen_tmp


    do ik = 1,nks(ir)
       ! lower average region
       if (.not.(smin.ge.ksb(ik,ir).or.smax.le.ksb(ik-1,ir))) then 

          ! both ends of test range in the same cell ... would require very small value of Lcoll or very big grid
          if ((smin.ge.ksb(ik-1,ir).and.smin.le.ksb(ik,ir)).and.(smax.ge.ksb(ik-1,ir).and.smax.le.ksb(ik,ir))) then 

             slen_tmp = smax - smin
             call update_impdens_data(icnt,ik,ir,ireg,slen_tmp,fp_reg,fp_in)

             ! Deal with end contributions - first half cell
          elseif (smin.ge.ksb(ik-1,ir).and.smin.le.ksb(ik,ir)) then 
             slen_tmp = ksb(ik,ir) - smin
             call update_impdens_data(icnt,ik,ir,ireg,slen_tmp,fp_reg,fp_in)

             ! last half cell
          elseif (smax.ge.ksb(ik-1,ir).and.smax.le.ksb(ik,ir)) then
             slen_tmp = smax - ksb(ik-1,ir)
             call update_impdens_data(icnt,ik,ir,ireg,slen_tmp,fp_reg,fp_in)

             ! full cell
          else
             slen_tmp = ksb(ik,ir) - ksb(ik-1,ir)
             call update_impdens_data(icnt,ik,ir,ireg,slen_tmp,fp_reg,fp_in)
          endif

       endif
    end do

    ! re-normalize to get average impurity density (w/o absfac)

    if (slen(icnt,ireg).gt.0.0) then 
       do iz = 1,nizs
          impdens(iz,icnt,ireg) = impdens(iz,icnt,ireg) / slen(icnt,ireg)
       enddo
       impdens(maxizs+1,icnt,ireg) = impdens(maxizs+1,icnt,ireg) / slen(icnt,ireg)
    endif


  end subroutine update_impdens



  subroutine update_impdens_data(icnt,ik,ir,ireg,slen_tmp,fp_reg,fp_in)
    use global_parameters
    use mod_outcom
    use mod_dynam2
    use mod_fp_data
    implicit none

    integer :: icnt,ik,ir,iz,ireg,fp_reg,fp_in
    real :: slen_tmp

    slen(icnt,ireg) = slen(icnt,ireg) + slen_tmp

    if (fp_reg.eq.0) then 
       ! use regular grid

       do iz = 1,nizs
 
          impdens(iz,icnt,ireg) = impdens(iz,icnt,ireg) + sdlims(ik,ir,iz) * slen_tmp
          impdens(maxizs+1,icnt,ireg) = impdens(maxizs+1,icnt,ireg) + sdlims(ik,ir,iz) * slen_tmp

       end do

    else
       ! add content from far periphery grid

       do iz = 1,nizs
 
          impdens(iz,icnt,ireg) = impdens(iz,icnt,ireg) + fp_density(ik,fp_in,iz,fp_reg) * slen_tmp
          impdens(maxizs+1,icnt,ireg) = impdens(maxizs+1,icnt,ireg) + fp_density(ik,fp_in,iz,fp_reg) * slen_tmp

       end do

   endif


  end subroutine update_impdens_data



  subroutine calc_collector_probe_flux(icnt)
    use global_parameters
    use mod_outcom
    use mod_dynam2
    use mod_cgeom
    use mod_comtor
    implicit none
    integer :: icnt
    real :: slcoll

    ! Locals
    integer :: ik, ir, iz
    real :: ne,te,ti,cs,alpha,nu,riz

    ne = local_vals(icnt,1)
    te = local_vals(icnt,2)
    ti = local_vals(icnt,3)
    cs = local_vals(icnt,6)
    slcoll = lcoll(icnt)

    do iz = 1,nizs

       riz = real(iz)
       nu = 1.4e-12 * (ne * riz**2/ te**(3/2))*(crmb/crmi)
       alpha = riz * ech * te /(crmi*amu * slcoll**2) + nu * cs /slcoll
       impflux(iz,icnt,1) = 0.25 * impdens(iz,icnt,1) * nu * slcoll * ( -1.0 + sqrt(1.0 + 8.0 * alpha/nu**2))
       impflux(iz,icnt,2) = 0.25 * impdens(iz,icnt,2) * nu * slcoll * ( -1.0 + sqrt(1.0 + 8.0 * alpha/nu**2))
       impflux(iz,icnt,3) = 0.25 * impdens(iz,icnt,3) * nu * slcoll * ( -1.0 + sqrt(1.0 + 8.0 * alpha/nu**2))

       impflux(maxizs+1,icnt,1) = impflux(maxizs+1,icnt,1) + impflux(iz,icnt,1)
       impflux(maxizs+1,icnt,2) = impflux(maxizs+1,icnt,2) + impflux(iz,icnt,2)
       impflux(maxizs+1,icnt,3) = impflux(maxizs+1,icnt,3) + impflux(iz,icnt,3)

    end do


  end subroutine calc_collector_probe_flux


  subroutine print_debug_data
    use global_parameters
    use mod_cgeom
    use mod_dynam2
    use mod_fp_data
    use mod_fperiph
    implicit none
    integer :: ik,ir,in,iz
    integer :: fp_reg
    integer :: fp_in

    fp_reg = fp_main
    fp_in = 1
    ir = irwall-1


    if (fpopt.eq.5.or.fpopt.eq.6.and.allocated(fp_density).and.allocated(fp_grid_area)) then 

       do ik = 1,nks(ir)
          !write(0,'(a,4i6,200(1x,i5,1x,2(1x,g12.5)))') 'Density:',ik,ir,fp_in,fp_reg,(iz,sdlims(ik,ir,iz),fp_density(ik,fp_in,iz,fp_reg),iz=1,10)
          !write(0,'(a,4i6,200(1x,i5,1x,2(1x,g12.5)))') 'Density:',ik,ir,fp_in,fp_reg,(iz,sdlims(ik,ir,iz),fp_density(ik,fp_in,iz,fp_reg),iz=11,20)

          write(6,'(a,4i6,3(1x,g12.5),200(1x,i5,1x,2(1x,g12.5)))') 'Density:',ik,ir,fp_in,fp_reg,kareas(ik,ir),fp_grid_area(ik,fp_reg),kareas(ik,ir)/fp_grid_area(ik,fp_reg),&
                   &(iz,sdlims(ik,ir,iz),fp_density(ik,fp_in,iz,fp_reg),iz=1,10)
          write(6,'(a,4i6,3(1x,g12.5),200(1x,i5,1x,2(1x,g12.5)))') 'Density:',ik,ir,fp_in,fp_reg,kareas(ik,ir),fp_grid_area(ik,fp_reg),kareas(ik,ir)/fp_grid_area(ik,fp_reg),&
                   &(iz,kareas(ik,ir)*sdlims(ik,ir,iz),fp_density(ik,fp_in,iz,fp_reg)*fp_grid_area(ik,fp_reg),iz=1,10)

          !write(6,'(a,4i6,200(1x,i5,1x,2(1x,g12.5)))') 'Density:',ik,ir,fp_in,fp_reg,(iz,sdlims(ik,ir,iz),fp_density(ik,fp_in,iz,fp_reg),iz=11,20)

       end do
    endif


  end subroutine print_debug_data


  subroutine write_fp_main_density(ounit,maxnk,nizs,sol_axis,sol_impdens,ring_axis,ring_data)
    use global_parameters
    use error_handling
    use mod_fperiph
    use mod_fp_data
    use mod_cgeom
    use mod_comtor
    implicit none
    integer :: ounit,maxnk,nizs
    integer :: fp_reg
    real,allocatable :: sol_impdens(:),sol_axis(:),ring_axis(:),ring_data(:)
    integer :: ir,fp_in,ik,iz,in
    integer,external :: ipos
    real :: fact

    if ((.not.allocated(sol_impdens)).or.(.not.allocated(sol_axis))) then
       call errmsg('MOD_COLLECTOR_PROBE:WRITE_FP_MAIN_DENSITY','ROUTINE CALLED BUT sol_impdens or sol_axis not allocated')
       return
    endif

    if ((.not.allocated(ring_data)).or.(.not.allocated(ring_axis))) then
       call errmsg('MOD_COLLECTOR_PROBE:WRITE_FP_MAIN_DENSITY','ROUTINE CALLED BUT ring_data or ring_axis not allocated')
       return
    endif


    fp_reg = fp_main
    ir = irwall-1

    !write(0,*) 'fps:',fp_reg,ir,fp_n_bins,nks(ir),nizs

    do fp_in = 1,fp_n_bins

       sol_impdens = 0.0

       ring_axis = 0.0
       ring_data = 0.0         
       ring_axis(1) = 0.0
       ring_axis(nks(ir)+2) = 1.0

       do ik = 1,nks(ir)
          ring_axis(ik+1) = kps(ik,ir)/kpmaxs(ir)
          do iz = 0,nizs
             ring_data(ik+1) = ring_data(ik+1) +  fp_density(ik,fp_in,iz,fp_reg) * absfac
          enddo
       end do
       ring_data(1) = ring_data(2)
       ring_data(nks(ir)+2) = ring_data(nks(ir)+1)


       ! interpolate

       do ik = 1, maxnk

          in = ipos(sol_axis(ik),ring_axis,nks(ir)+2)

          if (in.le.1) then 
             sol_impdens(ik) = ring_data(1)
          elseif (in.ge.nks(ir)+2) then 
             sol_impdens(ik) = ring_data(nks(ir)+2)
          else
             fact = (sol_axis(ik)-ring_axis(in-1)) /(ring_axis(in)-ring_axis(in-1))
             sol_impdens(ik) = fact * (ring_data(in)-ring_data(in-1))   + ring_data(in-1)
          endif

       end do

       !write(ounit,'(201(1x,g18.8))') 0.0, middist(ir,outer_targid) + fp_grid_dist(fp_in,fp_reg), (sol_impdens(in),in=1,maxnk)
       write(ounit,'(201(1x,g18.8))') 0.0, midplane_axis(ir) + fp_grid_dist(fp_in,fp_reg), (sol_impdens(in),in=1,maxnk)

    end do



  end subroutine write_fp_main_density


  subroutine calc_axis(midplane_axis,rsep_out,rsep_in,r1p,z1p,r2p,z2p,axis_opt)
    use global_parameters
    use mod_cgeom
    implicit none

    real :: r1p,z1p,r2p,z2p
    integer:: axis_opt
    real :: midplane_axis(maxnrs)
    real :: rsep_out,rsep_in
    !
    !     Note: on extended grids the rings may not be ordered consecutively
    !           so the usual routine to find midplane coordinates may not be 
    !           adequate.
    !
    !   Revise: -Intersection of the line Z=Z0 with the grid to obtain the
    !             value of Rmid for any rings intersecting this line.
    !             -Calculate the value of R at the separatrix
    !             -Estimate the R-Rsep distance based on these intersection 
    !              values. 
    !             -Linearly interpolate in PSI across the separatrix to find
    !              the value for Rsep_in and Rsep_out
    !
    !
    !     First find Rsep_in and Rsep_out
    !
    !     Line segments Z=Z0, R=[RMIN,R0] and R=[R0,RMAX]       
    !     
    !
    !     Check every cell on every ring for an intersection with either 
    !     the inner or outer midplane line. There should be only one 
    !     intersection of this type on 
    !
    !     Irwall is a boundary ring so do not include it
    !     if irwall2 > irwall then use it instead 
    !     Not checking for intersections in the PFZ
    !

    !
    !     locals
    !
    real*8 :: rcell1,zcell1,rcell2,zcell2
    real*8 :: rin1,zin1,rin2,zin2
    real*8 :: rout1,zout1,rout2,zout2
    real*8 :: rint,zint
    real :: midplane_psi(maxnrs)
    integer :: flag,sect
    integer :: ik,ir,in,ikmid


    midplane_axis = 0.0
    rsep_out = 0.0

    do ir = 1,nrs

       ikmid = nks(ir)/2

       !
       !        Use grid psi values if available
       !        Note: the two should be about the same
       !
       if (psifl(ikmid,ir).ne.0.0) then 
          midplane_psi(ir) = psifl(ikmid,ir)
       else
          midplane_psi(ir) = psitarg(ir,1)
       endif

    end do




    ! PSI - does not work for peripheral mesh
    if (axis_opt.eq.1) then 

       midplane_axis = midplane_psi


       ! 2: R-Rsep OMP - rings that do not intersect the outer midplane don't have a value - this can be a problem
       ! when there are shorter rings below the OMP where the collector probe samples
       ! 3: R-Rsep along probe location radial to plasma
    elseif (axis_opt.eq.2.or.axis_opt.eq.3) then 
       !

       if (axis_opt.eq.2) then 

          rout1 = r0
          zout1 = z0
          rout2 = rmax
          zout2 = z0

          
       elseif (axis_opt.eq.3) then 

          if (z1p.eq.z2p) then 
             rout1 = r0
             zout1 = z1p
             rout2 = rmax
             zout2 = z1p
          elseif (r1p.eq.r2p) then 
             rout1 = r1p
             zout1 = z0
             rout2 = r2p
             zout2 = zmax
          else
             rout1 = r1p
             zout1 = z1p
             rout2 = r2p
             zout2 = z2p
          endif

       endif

       !
       !     Loop over entire grid
       !     midplane coordinates will only exist for rings with intersections
       !
       !     Initialize to zero
       !

       do ir = 1,nrs

          do ik = 1,nks(ir)
             !
             !           Define center line of cell
             !         
             rcell1= krb(ik-1,ir)
             zcell1= kzb(ik-1,ir)
             rcell2= krb(ik,ir)
             zcell2= kzb(ik,ir)

             !
             !           Check for intersection with outer midplane
             !
             call intsect2dp(rcell1,zcell1,rcell2,zcell2,&
                  rout1,zout1,rout2,zout2,rint,zint,sect,flag)
             !
             !           Intersection between both sets of end points AND non-colinear
             !           Give outer Rsep value for ring
             !
             if (sect.eq.1.and.flag.eq.0) then 
                midplane_axis(ir) = rint
             endif

          end do

       end do

       rsep_out = (1.0 - midplane_psi(irsep-1))&
            /(midplane_psi(irsep)-midplane_psi(irsep-1))&
            *(midplane_axis(irsep)-midplane_axis(irsep-1))&
            + midplane_axis(irsep-1)

       do ir = 1,nrs
          if (midplane_axis(ir).ne.0.0) then 
             midplane_axis(ir) = midplane_axis(ir) - rsep_out
          endif
       end do

    endif



    do ir = 1,nrs
       write(6,'(a,3i8,10(1x,g18.8))') 'MIDPLANE_AXIS1:',ir,irsep,axis_opt,midplane_axis(ir),rsep_out
    end do






    return
  end subroutine calc_axis





end module mod_collector_probe
