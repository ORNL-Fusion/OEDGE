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
  real :: local_info(maxseg,6)

  public :: collector_probe


contains

  subroutine collector_probe(r1p,z1p,r2p,z2p,probe_diameter,dperp,iopt)
    use global_parameters
    use mod_cgeom
    use mod_comtor
    use mod_outcom
    real :: r1p,z1p,r2p,z2p,probe_diameter
    integer :: iopt

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
    !     Local Variables
    !
    integer :: in,ip,ik,ir,icnt
    real :: sint,psin
    real :: tmp,rsect,zsect
    real :: cs

    integer :: outunit
    logical :: exist
    character*1024 :: fname

    !

    ! Initialization
    impdens = 0.0
    impflux = 0.0
    slen = 0.0


    !
    !     Loop through the main SOL rings starting at the separatrix and the
    !     outer target - find the first intersection with the specified LOS
    !     and extract/interpolate the values of Ne, Te, Ti and pressure.
    !
    !     Counter for number of intersections found
    !
    icnt = 0
    !
    !     jdemod - allow for collection of data from core as well - might help with plots
    !
    do ir = 1,irwall-1
       !
       !      do ir = irsep,irwall-1
       !
       do ik = 1,nks(ir)

          call find_intsect(ik,ir,r1p,z1p,r2p,z2p,rsect,zsect,sint,psin)

          if (sint.gt.0.0) then
             !
             !              Increment counter
             !
             icnt = icnt + 1

             !write(0,'(a,3i4,1p,12(1x,g12.5))') 'OSM1:',ik,ir,icnt,rsect,zsect,sint,kss(ik,ir),kss2(ik,ir),ksb(ik-1,ir),ksb(ik,ir)

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
             local_info(icnt,7) = middist(ir,outer_targid)

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
             call calc_impdens(ir,sint,lcoll(icnt),icnt)

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

    XLAB = '   PROBE DISTANCE (M)'
    YLAB = '   IMPURITY FLUX (PART/M2/S)'
    REF  = 'COLLECTOR PROBE IMPURITY FLUX'
    write(anly,'(a,f9.5)')  'DIAM(M) =',probe_diameter
    write(plane,'(a,f9.5)') 'DP(M2/S)=',dperp
    nview = ' '

    do in = 1,numthe
       !touts(in) = local_outs(in)
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
    write(outunit,'(30(1x,a))') 'INDEX','IK','IR','SINT','RSECT','ZSECT','PSI','ROMP','NE','TE','TI','CS','LCOLL','LEN1','LEN2','LEN3','DIST',&
               'IMPDENS_'//trim(elabs(1)(1:4)),'IMPDENS_'//trim(elabs(2)(1:4)),'IMPDENS_CENT','IMPFLUX_'//trim(elabs(1)(1:4)),'IMPFLUX_'//trim(elabs(2)(1:4)),'IMPFLUX_CENT'
    do in = 1,numthe
       write(outunit,'(i8,30(1x,g12.5))') in,local_info(in,1),local_info(in,2),local_info(in,3),local_info(in,4),local_info(in,5),local_info(in,6),local_info(in,7),&
            local_vals(in,1),local_vals(in,2),local_vals(in,3),local_vals(in,6),&
            lcoll(in),slen(in,1),slen(in,2),slen(in,3),&
            local_outs(in),impdens(maxizs+1,in,1),impdens(maxizs+1,in,2),impdens(maxizs+1,in,3),&
            impflux(maxizs+1,in,1),impflux(maxizs+1,in,2),impflux(maxizs+1,in,3)
    end do
    write(outunit,'(a)') 

    close(outunit)

  end subroutine collector_probe




  subroutine calc_impdens(ir,sint,slcoll,icnt)
    use global_parameters
    use mod_outcom
    use mod_dynam2
    use mod_cgeom
    use mod_comtor
    implicit none
    integer :: icnt,ir
    real :: slcoll,sint

    ! Locals

    real :: smin,smax
    real :: smin_cent,smax_cent
    real :: slen_tmp

    !
    ! Need 3 average densities
    ! 
    ! 1)  Smin to sint ... lower density
    ! 2)  sint to smax ... upper density
    ! 3)  smin_cent to smax_cent ... central density
    !

    smin = max(sint-slcoll,0.0)
    smax = min(sint+slcoll,ksmaxs(ir))

    ! arbitrarily specify the local region as +/- LCOLL/4
    smin_cent = max(sint-slcoll/4.0,0.0)
    smax_cent = min(sint+slcoll/4.0,ksmaxs(ir))

    call update_impdens(ir,icnt,smin,sint,1)
    call update_impdens(ir,icnt,sint,smax,2)
    call update_impdens(ir,icnt,smin_cent,smax_cent,3)

  end subroutine calc_impdens


  subroutine update_impdens(ir,icnt,smin,smax,ireg)
    use global_parameters
    use mod_cgeom
    use mod_outcom
    implicit none

    integer :: ireg,ir,icnt,ik,iz
    real :: smin,smax
    real :: slen_tmp

    do ik = 1,nks(ir)
       ! lower average region
       if (.not.(smin.ge.ksb(ik,ir).or.smax.le.ksb(ik-1,ir))) then 

          ! both ends of test range in the same cell ... would require very small value of Lcoll or very big grid
          if ((smin.ge.ksb(ik-1,ir).and.smin.le.ksb(ik,ir)).and.(smax.ge.ksb(ik-1,ir).and.smax.le.ksb(ik,ir))) then 

             slen_tmp = smax - smin
             call update_impdens_data(icnt,ik,ir,ireg,slen_tmp)

             ! Deal with end contributions - first half cell
          elseif (smin.ge.ksb(ik-1,ir).and.smin.le.ksb(ik,ir)) then 
             slen_tmp = ksb(ik,ir) - smin
             call update_impdens_data(icnt,ik,ir,ireg,slen_tmp)

             ! last half cell
          elseif (smax.ge.ksb(ik-1,ir).and.smax.le.ksb(ik,ir)) then
             slen_tmp = smax - ksb(ik-1,ir)
             call update_impdens_data(icnt,ik,ir,ireg,slen_tmp)

             ! full cell
          else
             slen_tmp = ksb(ik,ir) - ksb(ik-1,ir)
             call update_impdens_data(icnt,ik,ir,ireg,slen_tmp)
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



  subroutine update_impdens_data(icnt,ik,ir,ireg,slen_tmp)
    use global_parameters
    use mod_outcom
    use mod_dynam2
    implicit none

    integer :: icnt,ik,ir,iz,ireg
    real :: slen_tmp

    slen(icnt,ireg) = slen(icnt,ireg) + slen_tmp

    do iz = 1,nizs

       impdens(iz,icnt,ireg) = impdens(iz,icnt,ireg) + sdlims(ik,ir,iz) * slen_tmp
       impdens(maxizs+1,icnt,ireg) = impdens(maxizs+1,icnt,ireg) + sdlims(ik,ir,iz) * slen_tmp

    end do

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



end module mod_collector_probe
