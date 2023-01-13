module mod_sol22_interface


  
  use mod_assignpp
  use mod_sol22
  use mod_sol22_support
  use mod_sol22_divimp
  implicit none




contains



  subroutine calcsol_interface (irlim1,irlim2,ikopt)
    use error_handling
    use debug_options
    use mod_sol22_input
    use sol22_debug
    use mod_params
    use mod_comtor
    use mod_cgeom
    use mod_pindata
    use mod_cedge2d
    use mod_slcom
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    use allocate_arrays
    implicit none

    !     The subroutine calcsol uses numerical methods (Runge-Kutta) to
    !     solve the fluid equations along the field lines. This routine
    !     acts as an interface between DIVIMP and the calcsol subroutine.
    !     This approach was chosen so that the variables and components
    !     of the solascv module could remain as independent of DIVIMP as
    !     possible - thus facilitating changes and corrections to the code
    !     based on developments in the stand-alone program module solascv.f
    !     ... in particular ... at this time the calcsol series of
    !     subroutines calculate all values in extended precision for
    !     accuracy - which is not generally wanted or needed inside DIVIMP.
    !     In addition, the calcsol routines have their own variable names
    !     for many quantities already included in DIVIMP - in order to avoid
    !     re-writing the code at this time - the interface routine assures
    !     that all the quantities required in the calcsol subroutine series
    !     are properly loaded.

    !     David Elder,  Jan 24, 1995



    !     include 'params'
    !     include 'comtor'
    !     include 'cgeom'
    !     include 'pindata'
    !     include 'cedge2d'
    ! slmod begin - new
    !     INCLUDE 'slcom'

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'



    integer irlim1, irlim2,ikopt
    !COMMON /POWERFLOW/ cve        ,cvi        ,cde        ,cdi
    REAL               cve(MXSPTS),cvi(MXSPTS),cde(MXSPTS),cdi(MXSPTS)

    REAL Clock2
    ! slmod end

    REAL deltat

    integer ir,ircor,midnks,ik,ikstart,ikend
    real gperpa(maxnks,maxnrs)
    real oldknbs(maxnks,maxnrs)
    real grad_oldknbs(maxnks,maxnrs)
    real oldktibs(maxnks,maxnrs)
    real oldktebs(maxnks,maxnrs)

    !     Radiation loss term from previous DIVIMP run
    !     Impurity Ionization/recombination potential energy loss to e-

    real oldkvds(maxnds)
    real div_tpowls(maxnks,maxnrs),div_tcooliz(maxnks,maxnrs)

    real div_cool(maxnks,maxnrs) 

    real,allocatable :: ext_epowsrc(:,:)  ! external electron and ion power terms
    real,allocatable :: ext_ipowsrc(:,:)

    real gradi,grado,flux,len1,len2
    real quant2grad,polysidelen
    external quant2grad,polysidelen

    integer rc
    integer errcode,seterror,new_errlevel

    !     For calling calcsol

    real*8 serr
    integer npts,id
    real*8 spts(mxspts),fact,slope,gam
    real*8 te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts),exp_press(mxspts),act_press(mxspts),prad(mxspts)
    real*8 ttarg

    !     Assigning ACTFFRIC

    real*8 int_powrat(3)
    !real*8 find_ffric 
    !external find_ffric

    integer :: ierr

    !     Local variables

    real ds1,dp1,dt1,nb1,ds2,dp2,dt2,nb2
    integer applydef(maxnrs,2),in,pplasma,ikfirst,iklast
    integer sol22_iter

    !     Temporary local variables until these are loaded in input file

    data sol22_iter /0/
    integer:: pfz_dist_opt

    !     TEMPORARY MULTIPLIER

    !      real qesrc_mult
    !      parameter (qesrc_mult=1.0)

    !     Locals for smoothing

    real*8 :: pfz_dist_param(2)
    real temid,timid,nmid,smax,asmexp,tmp

    ! Initialize some output options in SOL22 using values from slcom
    call init_solcommon(osm_mode,outmode)


    !     Make a call to initialize debugging if it is on

    call pr_trace('MOD_SOL22_INTERFACE','START CALCSOL_INTERFACE')
    if (debug_sol22.ne.0) then 
       ! debug_sol22 is specified in the input file as unstructured input option 284
       ! parameters to the call are the ring number and ikopt for the half ring for which high res
       ! debugging data is required. 
       call init_sol22_debug(debug_sol22_ir,debug_sol22_ikopt) 

       !     Set a flag internal to the SOLASCV code that indicates that
       !     PIN data is available if the ionization option for it has been
       !     specified.

    endif


    !     If Switch(swe2d).lt.0 - use this to indicate using the
    !     Edge2D background plasma as the seed to PIN for the
    !     SOL option 22 iteration - flip the switch "+" for
    !     subsequent iterations.

    pinavail = piniter

    if (.not.pinavail) then

       e2dstart = 0

       if (switch(swe2d).lt.0.0.or.switch(swioni).eq.12.or.switch(swioni).eq.13.or.switch(swioni).eq.14) then

          sol22_iter = sol22_iter + 1
          switch(swe2d) = abs(switch(swe2d))

          e2dstart = 1
          if (cre2d.eq.0) then
             !               cre2d = 1
             call redge2d(0)

          endif
          if (switch(swioni).eq.14) then
             do ir = 1,nrs
                do ik = 1,nks(ir)
                   oldknbs(ik,ir) = e2dnbs(ik,ir)
                   oldktibs(ik,ir) = e2dtibs(ik,ir)
                   oldktebs(ik,ir) = e2dtebs(ik,ir)
                   knbs(ik,ir) = e2dnbs(ik,ir)
                   ktebs(ik,ir)= e2dtebs(ik,ir)
                   ktibs(ik,ir)= e2dtibs(ik,ir)
                   kvhs(ik,ir) = e2dvhs(ik,ir)
                   kes(ik,ir)  = e2des(ik,ir)

                end do
                if (ir.lt.irsep) then
                   oldknbs(nks(ir),ir) = e2dnbs(1,ir)
                   oldktibs(nks(ir),ir) = e2dtibs(1,ir)
                   oldktebs(nks(ir),ir) = e2dtebs(1,ir)
                   knbs(nks(ir),ir) = e2dnbs(1,ir)
                   ktebs(nks(ir),ir)= e2dtebs(1,ir)
                   ktibs(nks(ir),ir)= e2dtibs(1,ir)
                   kvhs(nks(ir),ir) = e2dvhs(1,ir)
                   kes(nks(ir),ir)  = e2des(1,ir)

                endif

             end do

          else

             do ir = irlim1,irlim2

                call set_ikvals(ir,ikstart,ikend,ikopt)
                do ik = ikstart,ikend
                   oldknbs(ik,ir) = e2dnbs(ik,ir)
                   oldktibs(ik,ir) = e2dtibs(ik,ir)
                   oldktebs(ik,ir) = e2dtebs(ik,ir)
                   knbs(ik,ir) = e2dnbs(ik,ir)
                   ktebs(ik,ir)= e2dtebs(ik,ir)
                   ktibs(ik,ir)= e2dtibs(ik,ir)
                   kvhs(ik,ir) = e2dvhs(ik,ir)
                   kes(ik,ir)  = e2des(ik,ir)
                end do

             end do

          endif

          return

          !           Set-up for directly read Edge2D ionization

       elseif (switch(swioni).eq.11) then
          if (cre2d.eq.0) then
             !               cre2d = 1
             call redge2d(0)

          endif
          do ir = irlim1,irlim2
             do ik = 1,nks(ir)
                oldknbs(ik,ir) =  e2dnbs(ik,ir)
                oldktibs(ik,ir) = e2dtibs(ik,ir)
                oldktebs(ik,ir) = e2dtebs(ik,ir)
             end do

          end do

          !           Set-up for directly read Edge2D ionization

       elseif (switch(swioni).eq.15) then
          if (cre2d.eq.0) then
             !               cre2d = 1
             call redge2d(0)

          endif
          do ir = 1,nrs
             do ik = 1,nks(ir)
                oldknbs(ik,ir) =  e2dnbs(ik,ir)
                oldktibs(ik,ir) = e2dtibs(ik,ir)
                oldktebs(ik,ir) = e2dtebs(ik,ir)
                knbs(ik,ir) = e2dnbs(ik,ir)
                ktebs(ik,ir)= e2dtebs(ik,ir)
                ktibs(ik,ir)= e2dtibs(ik,ir)
                kvhs(ik,ir) = e2dvhs(ik,ir)
                kes(ik,ir)  = e2des(ik,ir)

             end do

             if (ir.lt.irsep) then
                oldknbs(nks(ir),ir) = e2dnbs(1,ir)
                oldktibs(nks(ir),ir) = e2dtibs(1,ir)
                oldktebs(nks(ir),ir) = e2dtebs(1,ir)
                knbs(nks(ir),ir) = e2dnbs(1,ir)
                ktebs(nks(ir),ir)= e2dtebs(1,ir)
                ktibs(nks(ir),ir)= e2dtibs(1,ir)
                kvhs(nks(ir),ir) = e2dvhs(1,ir)

                kes(nks(ir),ir)  = e2des(1,ir)

             endif

             !        Load previously calculated DIVIMP background for the first
             !        iteration of SOL 22.

          end do

          !           Set-up for directly read DIV background plasma

       elseif (switch(swioni).eq.16.or.switch(swioni).eq.17.or.switch(swioni).eq.18) then

          !           Increment SOL 22 iteration count

          call readdivbg

          sol22_iter = sol22_iter + 1
          do ir = 1,nrs
             do ik = 1,nks(ir)
                oldknbs(ik,ir) =  knbs(ik,ir)
                oldktibs(ik,ir) = ktibs(ik,ir)
                oldktebs(ik,ir) = ktebs(ik,ir)
             end do

          end do

          do id = 1,nds

             oldkvds(id) = kvds(id)

             !           Exit

          end do

          return



       endif

    elseif ((switch(swioni).eq.13.or.switch(swioni).eq.14).and.sol22_iter.eq.1) then

       sol22_iter = sol22_iter + 1
       if (switch(swioni).eq.14) then
          do ir = 1,nrs
             do ik = 1,nks(ir)
                oldknbs(ik,ir) = e2dnbs(ik,ir)
                oldktibs(ik,ir) = e2dtibs(ik,ir)
                oldktebs(ik,ir) = e2dtebs(ik,ir)
                knbs(ik,ir) = e2dnbs(ik,ir)
                ktebs(ik,ir)= e2dtebs(ik,ir)
                ktibs(ik,ir)= e2dtibs(ik,ir)
                kvhs(ik,ir) = e2dvhs(ik,ir)
                kes(ik,ir)  = e2des(ik,ir)

             end do
             if (ir.lt.irsep) then
                knbs(nks(ir),ir) = e2dnbs(1,ir)
                ktebs(nks(ir),ir)= e2dtebs(1,ir)
                ktibs(nks(ir),ir)= e2dtibs(1,ir)
                kvhs(nks(ir),ir) = e2dvhs(1,ir)
                kes(nks(ir),ir)  = e2des(1,ir)

             endif

          end do

       else
          do ir = irlim1,irlim2
             do ik = 1,nks(ir)
                oldknbs(ik,ir) = e2dnbs(ik,ir)
                oldktibs(ik,ir) = e2dtibs(ik,ir)
                oldktebs(ik,ir) = e2dtebs(ik,ir)
                knbs(ik,ir) = e2dnbs(ik,ir)
                ktebs(ik,ir)= e2dtebs(ik,ir)
                ktibs(ik,ir)= e2dtibs(ik,ir)
                kvhs(ik,ir) = e2dvhs(ik,ir)
                kes(ik,ir)  = e2des(ik,ir)
             end do

          end do

       endif

       !     ENDIF for .not.pinavail

       return

       !     For initial ionization option 17 - re-load the core from the
       !     original plasma solution on every iteration. If the Tgrad option
       !     is set to 0 then also load the previous private plasma
       !     solution.

    endif

    !        Copy over target velocities so that mach 1 is not forced.

    if (switch(swioni).eq.16.or.switch(swioni).eq.17) then
       do ir = irsep,nrs
          kvds(idds(ir,1)) = oldkvds(idds(ir,1))
          kvds(idds(ir,2)) = oldkvds(idds(ir,2))

       end do

    endif

    !        Reload core from original plasma solution - uses e2d variables
    !        but will be loading a previous DIVIMP solution.

    if (switch(swioni).eq.17) then
       do ir = 1,irsep-1
          do ik = 1,nks(ir)
             knbs(ik,ir) = e2dnbs(ik,ir)
             ktebs(ik,ir)= e2dtebs(ik,ir)
             ktibs(ik,ir)= e2dtibs(ik,ir)
             kvhs(ik,ir) = e2dvhs(ik,ir)
             kes(ik,ir)  = e2des(ik,ir)
          end do

          !        If private plasma solution is also off - then reload original
          !        for private plasma as well.

       end do

       if (ciopto.eq.0) then
          do ir = irtrap,nrs
             do ik = 1,nks(ir)
                knbs(ik,ir) = e2dnbs(ik,ir)
                ktebs(ik,ir)= e2dtebs(ik,ir)
                ktibs(ik,ir)= e2dtibs(ik,ir)
                kvhs(ik,ir) = e2dvhs(ik,ir)
                kes(ik,ir)  = e2des(ik,ir)
             end do

          end do

       endif



       !     FOR Radiation option 4 - read in the radiation data from
       !     a previous DIVIMP run.

    endif
    call rzero(div_tpowls,maxnks*maxnrs) 
    call rzero(div_tcooliz,maxnks*maxnrs) 

    call rzero(div_cool,maxnks*maxnrs) 
    !
    if (switch(swprad).eq.4) then 
       !        Sum components of impurity electron cooling terms together

       call readdivaux('TPOWLS:',div_tpowls,maxnks,maxnrs,1,1,ierr)
       call readdivaux('TCOOLIZ:',div_tcooliz,maxnks,maxnrs,1,1,ierr)
       do ir = 1,nrs
          do ik = 1,nks(ir)
             div_cool(ik,ir) = div_tpowls(ik,ir)+div_tcooliz(ik,ir)
          end do
       end do
    endif

    !
    ! Load external electron power terms if the option is turned on
    !
    if (switch(swepow).ne.0.0) then 
       ! Load Epower from divimp auxiliary input file
       ! allocate storage for data 
       call allocate_array(ext_epowsrc,maxnks,maxnrs,'ext_epowsrc',ierr)
       if (ierr.ne.0) then
          call errmsg('MOD_SOL22_INTERFACE: ERROR ALLOCATING STORAGE FOR ext_epowsrc : setting switch(swepow)=0.0 : IERR=', ierr)
          switch(swepow) = 0.0
       else
          if (switch(swepow).eq.1.0) then 
             call readdivaux('EXTEPOW:',ext_epowsrc,maxnks,maxnrs,1,1,ierr)
             if (ierr.ne.0) then
                call errmsg('MOD_SOL22_INTERFACE: ERROR loading ext_epowsrc from readdivaux: setting switch(swepow)=0.0 : IERR=', ierr)
                switch(swepow) = 0.0
             endif
          elseif (switch(swepow).eq.2.0) then
             call load_extpowsrc(ext_epowsrc,ext_epow_fn,maxnks,maxnrs,nrs,nks,rs,zs,ierr)
             if (ierr.ne.0) then
                call errmsg('MOD_SOL22_INTERFACE: ERROR loading ext_epowsrc from load_extpowsrc: setting switch(swepow)=0.0 : IERR=', ierr)
                switch(swepow) = 0.0
             endif
          endif
       endif
    endif

    !
    ! Load external ion power terms if the option is turned on
    !

    if (switch(swipow).ne.0.0) then 
       ! Load Epower from divimp auxiliary input file
       ! allocate storage for data 
       call allocate_array(ext_ipowsrc,maxnks,maxnrs,'ext_ipowsrc',ierr)
       if (ierr.ne.0) then
          call errmsg('MOD_SOL22_INTERFACE: ERROR ALLOCATING STORAGE FOR ext_ipowsrc : setting switch(swipow)=0.0 : IERR=', ierr)
          switch(swipow) = 0.0
       else
          if (switch(swipow).eq.1.0) then 
             call readdivaux('EXTIPOW:',ext_ipowsrc,maxnks,maxnrs,1,1,ierr)
             if (ierr.ne.0) then
                call errmsg('MOD_SOL22_INTERFACE: ERROR loading ext_ipowsrc from readdivaux: setting switch(swipow)=0.0 : IERR=', ierr)
                switch(swipow) = 0.0
             endif
          elseif (switch(swipow).eq.2.0) then
             call load_extpowsrc(ext_ipowsrc,ext_ipow_fn,maxnks,maxnrs,nrs,nks,rs,zs,ierr)
             if (ierr.ne.0) then
                call errmsg('MOD_SOL22_INTERFACE: ERROR loading ext_ipowsrc from load_extpowsrc: setting switch(swipow)=0.0 : IERR=', ierr)
                switch(swipow) = 0.0
             endif
          endif
       endif
    endif


    !     For radiation option 5 calculate the distribution of radiation over the grid
    !     given total radiation from each region.
    !     1) Inner divertor
    !     2) Outer divertor
    !     3) Inner and outer PFZ
    !     4) Inner SOL to top
    !     5) Outer SOL to top

    !     Distribute Pin_region proportional to targ_flux(ir) * kvols(ik,ir) and
    !     integrate over region to get the appropriate total prad
    !     Note: Pinqe and Pinqi options should be off when this option is used.



    !     Increase iteration count

    call pr_trace('MOD_CALCSOL_INTERFACE','BEFORE INIT')

    !     Initialization

    sol22_iter = sol22_iter + 1
    call izero(cerr,maxnrs*2)
    call izero(cdeferr,maxnrs*2)
    call rzero(cserr,maxnrs*2)
    call rzero(cdefserr,maxnrs*2)
    call rzero(cdeferropt,maxnrs*2)

    !     Zero the correction factor array.

    call rzero(gperpa,maxnks*maxnrs)

    !     Zero the pressure array

    !      call rzero(kpress,maxnks*maxnrs*2)

    call qzero(rconst,mxspts)

    call iinit(applydef,maxnrs*2,-1)
    if (ndef.gt.0) then
       do in = 1,ndef
          applydef(INT(deflist(in,1)),INT(deflist(in,2))) =deflist(in,3)
       end do

       !     Read in EDGE 2D target data if Edge 2D compatibility is ON

    end if
    if (switch(swe2d).ne.0.0)then
       call redge2d(1)


       !     Set up fixed input quantities

    endif
    title = 'DIVIMP'
    mb = crmb
    zb = rizb
    k0e = ck0

    !     The following graph option will generate a series of plots based
    !     on the SOL Option 22 data.They are currently turned off and in the
    !     original code were read in as options.These are left for debugging
    !     purposes and can be turned on by setting the options equal to 1.
    !     This does mean however that the GHOST libraries would need to be
    !     bound to DIVIMP for these plots. The grahing code can be easily
    !     exised by commenting out the graph specific calls in each of the
    !     subroutines. This will preserve the ability to print out the
    !     relevant tables of information without needing to bind the GHOST
    !     libraries.

    !     David Elder                1995, May 16

    k0i = ck0i
    graph    = 1
    graphaux = 1
    graphvel = 1

    !     Adjust Ionization switch settings if appropriate

    !     Check to see if option set to use PIN data to
    !     normalize the analytic sources.

    graphran = 0.0

    pinnorm = 0
    if (switch(swion).eq.8.0.and.pinavail) then
       switch(swion) = switch(swioni)
       switch(swioni) = 8.0
       pinnorm = 1

       !     Calculate the second gradient taken from the OLD or
       !     previous solution - uses values stored in KNBS -

    endif

    call rzero(grad_oldknbs,maxnks*maxnrs)

    if (switch(swgperp).eq.7.or.switch(swgperpp).eq.7.or.switch(swgperp).eq.8.or.switch(swgperpp).eq.8) then
       do ir = irlim1,irlim2
          do ik = 1,nks(ir)

             call quantgrad(ir,kss(ik,ir),oldknbs,knds,gradi,grado,2)
             len1 = polysidelen(ik,ir,INWARD41,rc)
             len2 = polysidelen(ik,ir,OUTWARD23,rc)

             flux = gradi*len1 -  grado*len2
             if (kareas(ik,ir).ne.0.0) then
                grad_oldknbs(ik,ir)= flux/kareas(ik,ir)
             else
                grad_oldknbs(ik,ir)= 1.0


                !               write (6,'(a,2i4,15(1x,g12.5))') 'GP:',ik,ir,
                !     >                 cdperp*grad_oldknbs(ik,ir),
                !     >                 gradi,grado,flux,karea2(ik,ir),
                !     >                 len1,len2,len1/len2

             endif
          end do

       end  do


       !      write (6,*) 'Swion:',switch(swion),switch(swioni)


       !     Loop through rings

       !     All quantities that are used by the calcsol routine must be
       !     copied into appropriate variables. There are two reasons for
       !     this - first to allow the calcsol routines to use a different
       !     level of precision than the rest of the code. Second to make the
       !     contents of calcsol and its subroutines as independent of DIVIMP
       !     as possible - so that changes and additions may be easily made.

    endif
    if (ctestsol.gt.0.0) then
       irlim1 = int(ctestsol)
       irlim2 = int(ctestsol)

       !     Set debugging ring

    endif

    !      write (6,*) 'IR debug:',irdebug

    !     Set the starting cell for all rings - if the E2D compatibility option
    !     is set to 9 - set the starting cell not equal to 1 - may be an
    !     input parameter in later revisions.

    irdebug = irlim1

    if (switch(swe2d).eq.9) then

       ike2d_start = ike2d

    else

       ike2d_start = 1

       !     Calculate the target particle and power outfluxes.

       !     Temporarily set pfz_dist_opt and pfz_dist_param until these
       !     are added to the input file

    endif
    pfz_dist_opt = 1
    pfz_dist_param(1) = sepdist2(idds(irsep+4,1))

    pfz_dist_param(2) = sepdist2(idds(irsep+4,2))

    call pr_trace('MOD_CALCSOL_INTERFACE','BEFORE CALCFLUXES')

    !     If pin is available
    !     Call routine to calculate GPERP CORection factors

    call calcfluxes(gtarg,ionptarg,elecptarg,e2dgtarg,presstarg,gamcor,gamecor,ike2d_start,g_pfzsol,pe_pfzsol,&
         pi_pfzsol,pr_pfzsol,pfz_dist_opt,pfz_dist_param)
    call pr_trace('MOD_CALCSOL_INTERFACE','BEFORE CALCSOLIZ')

    !        Call routine to print comparison of sources.

    if (pinavail) then

       if (cprint.eq.9) call ioniz_comp

       call calcsoliz(rconst,recfrac,gtarg,areasum,gperpa,oldknbs,grad_oldknbs,pinion,pinrec,gperprat,ike2d_start)
    elseif ((switch(swioni).eq.11.or.switch(swioni).eq.15).and.(.not.pinavail)) then
       if (switch(swrecom).eq.2.0) then
          call calcsoliz(rconst,recfrac,gtarg,areasum,gperpa,oldknbs,grad_oldknbs,e2dion,e2dhrec,gperprat,ike2d_start)
       else
          call calcsoliz(rconst,recfrac,gtarg,areasum,gperpa,oldknbs,grad_oldknbs,e2dion,pinrec,gperprat,ike2d_start)
       endif

       !     Use the EDGE2D background through ALL iterations
       !     - essentially bypasses SOL 22 COMPLETELY -
       !     Used as a method of iterative PIN testing while still
       !     obtaining the SOL 22 flux debugging information from calcsoliz.

    endif
    call pr_trace('MOD_CALCSOL_INTERFACE','BEFORE EDGE2D SAVE')
    if (switch(swe2d).eq.4.or.switch(swe2d).eq.5.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7) then
       do ir = irlim1,irlim2
          do ik = 1,nks(ir)
             oldknbs(ik,ir) = e2dnbs(ik,ir)
             oldktibs(ik,ir) = e2dtibs(ik,ir)
             oldktebs(ik,ir) = e2dtebs(ik,ir)
             knbs(ik,ir) = e2dnbs(ik,ir)
             ktebs(ik,ir)= e2dtebs(ik,ir)
             ktibs(ik,ir)= e2dtibs(ik,ir)
             kvhs(ik,ir) = e2dvhs(ik,ir)
             kes(ik,ir)  = e2des(ik,ir)
          end do
       end do
       return

    endif
    call pr_trace('MOD_CALCSOL_INTERFACE','BEFORE RING LOOP')

    do ir = irlim1, irlim2


       !        Bypass ALL of the code for ir = irwall and ir = irtrap
       !        since these are rings with only virtual boundary cells.
       !        They should not be used in SOL option 22 - the solutions
       !        for rings irwall-1 and irtrap + 1 will eb copied in to
       !        ensure valid values for ne,te,ti just in case they are
       !        required.

       call set_ikvals(ir,ikstart,ikend,ikopt)

       !        Set the values of lengths to a maximum of 1/2 of the field
       !        line length. (Same for lenr)

       if (ir.eq.irwall.or.ir.eq.irtrap) goto 1000
       ringlen = ksmaxs2(ir)

       !        Set private plasma ring indicator

       ringnum = ir
       if (ringnum.gt.irwall) then
          pplasma = 1
       else
          pplasma = 0

          !        If the private plasma ionization option specifies the
          !        ANALYTIC option - then call this for Private Plasma
          !        rings.

       endif

       if (switch(swionp).eq.-2.0.and.pplasma.eq.1) then

          call specplas(ir,ir,ikopt)
          do ik = 1,nks(ir)
             oldknbs(ik,ir) = knbs(ik,ir)
             oldktibs(ik,ir) = ktibs(ik,ir)
             oldktebs(ik,ir) = ktebs(ik,ir)

          end do

          !            cycle

          !...  Assign PP plasma using Thomson data:

          goto 1000

          ! slmod begin - new
          !...LITTLE FIX:
       elseif ((switch(swionp).eq.-3.0.or.switch(swionp).eq.-4.0).and.pplasma.eq.1) then

          !            call thompp(ir,ir,ikopt,int(abs(swionp)))
          ! slmod end

          call thompp(ir,ir,ikopt,int(abs(switch(swionp))))
          do ik = 1, nks(ir)
             oldknbs (ik,ir) = knbs (ik,ir)
             oldktibs(ik,ir) = ktibs(ik,ir)
             oldktebs(ik,ir) = ktebs(ik,ir)

          end do
          ! slmod begin - new
          goto 1000
          !...        Assign private plasma from a listing of T and n in the input
          !           file:
       elseif ((switch(swionp).eq.-5.0.or.switch(swionp).eq.-6.0).and.pplasma.eq.1) then
          call assignpp(ir,ir,ikopt,int(abs(switch(swionp))))
          do ik = 1, nks(ir)
             oldknbs (ik,ir) = knbs (ik,ir)
             oldktibs(ik,ir) = ktibs(ik,ir)
             oldktebs(ik,ir) = ktebs(ik,ir)
          end do
          ! slmod end

          !            cycle



          goto 1000

          !        Set up the effective GPERP/core correction term.

       end if


       !         if ((pinavail
       !     >       .and.((pinnorm.eq.1.or.actswion.eq.2.0))
       !     >       .and.(actswgperp.eq.2.0.or.actswgperp.eq.6.0
       !     >             .or.actswgperp.eq.7.or.actswgperp.eq.8.0))
       !     >       .or.
       !     >       ((actswioni.eq.11.or.actswioni.eq.15).and.
       !     >        (actswgperp.eq.2.0.or.actswgperp.eq.6.0
       !     >         .or.actswgperp.eq.7.or.actswgperp.eq.8.0)
       !     >       .and.(.not.pinavail))
       !     >       ) then
       !            gnet = rconst(ir)


       !         else
       !            gnet = 0.0
       !         endif



       !        Set up the power distribution length for the
       !        power input options 7 and 8

       gnet = rconst(ir)
       spow = spowbeg * ringlen

       !         write (6,*) 'SPOW:',spow,spow2

       !        Set up gperp lengths for particle source compensation

       spow2= spowlen * ringlen
       sgperpbeg = gperpbegf * ringlen

       sgperpend = gperpendf * ringlen
       !         write(0,*) 'solascv:midnks:',midnks,ir,sol22_halfringlen_opt
       midnks = ikmids(ir)

       !-----------------------------------------------------------------
       !        OUTER TARGET
       !-----------------------------------------------------------------


       !        Only execute for the selected half-rings - ikopt


       if (ikopt.eq.1.or.ikopt.eq.3) then

          !        Check for sol22 debug start

          !           The parameters are current ring snd local IKOPT
          call pr_trace('MOD_CALCSOL_INTERFACE','START OUTER TARGET')

          !        Set the values for a call - do 1/2 of the ring at a time.


          !        Outer target first - i.e. elements 1 to 1/2 nks(ir)
          !                                  target designated as 2


          !        Outer target


          !        Reset switches if necessary

          if (debug_sol22.ne.0) call check_init_record_data(ir,1)

          !        Set active switches

          seterror = 0
          if (applydef(ir,2).eq.-1) then

             !           Reset Gperp switch for option 3 in case of under ionization

             call setsw(-1,pplasma,new_errlevel)
             if (( (pinavail.or.((.not.pinavail.and.(actswioni.eq.11.or.actswioni.eq.15)))).and.actswgperp.eq.3.0.and.&
                  rconst(ir).gt.0.0)) then
                actswgperp = 2.0
                write(6,*) 'Swgperp:',rconst(ir),actswgperp
             elseif ((.not.pinavail).and.(actswioni.ne.11.and.actswioni.ne.15).and.&
                  (actswgperp.eq.2.0.or.actswgperp.eq.3.0.or.actswgperp.eq.7.0.or.actswgperp.eq.8.0)) then
                actswgperp = 0.0
             elseif (actswgperp.eq.7.0.and.gperprat(ir).gt.5.0) then
                actswgperp = 2.0

             endif
          elseif (applydef(ir,2).ge.0) then
             call setsw(applydef(ir,2),pplasma,new_errlevel)

             actswerror = applydef(ir,2)

             if (actswerror.eq.0.0) seterror = 2
             cdeferr(ir,2) =  8
             cdefserr(ir,2) = 0.0
             cdeferropt(ir,2) = applydef(ir,2)

             !        Initialize HALF ring length to be used

             !        halfringlen = ksb(midnks,ir)

          endif
          if (sol22_halfringlen_opt.eq.0) then 
             halfringlen  = 0.5d0*ringlen
          elseif (sol22_halfringlen_opt.eq.1) then 
             ! ring midpoint is upper boundary of middle cell just below midpoint
             !            write(6,'(a,2i8,10(1x,g12.5))') 'Halfringlen1:',
             !     >           midnks,ir,halfringlen,ringlen/2.0
             !            write(0,'(a,2i8,10(1x,g12.5))') 'Halfringlen1:',
             !     >           midnks,ir,halfringlen,ringlen/2.0
             halfringlen = ksb(midnks+1,ir)

             !        Initialize the ionization and radiation source lengths

          endif

          !        Initialize the friction parameter for momentum loss options

          call initlen
          actffric = find_ffric(ir,2,actlenmom,actlammom)
          !
          !        Initialize major radius and target condition data

          call assign_radiation_parameters(ir,2)

          targfact = 1.0
          if (actswmajr.eq.1.0) then
             if (actswe2d.eq.0.0) then
                id = idds(ir,2)
                targfact = rp(id) / r0
             elseif (actswe2d.ne.0.0) then
                targfact = rs(1,ir) / r0
             endif

          endif
          te0 = kteds(idds(ir,2))
          ti0 = ktids(idds(ir,2))
          n0 = knds(idds(ir,2))

          v0 = kvds(idds(ir,2))
          if (forcet.eq.1) then
             ttarg = (((5.0* te0 + 3.5 * ti0)*sqrt(te0+ti0))/ (8.5*sqrt(2.0)))**(2.0/3.0)
             te0 = ttarg
             ti0 = ttarg
             kteds(idds(ir,2)) = te0
             ktids(idds(ir,2)) = ti0

          endif
          if (actswe2d.eq.9.0) then
             n1 = e2dnbs(ike2d_start,ir)
             te1= e2dtebs(ike2d_start,ir)
             ti1= e2dtibs(ike2d_start,ir)
             v1e2d = e2dvhs(ike2d_start,ir)

             vpe2d = v1e2d

             vpg = (e2dgpara(ike2d_start,ir)+e2dgpara(ike2d_start+1,ir))/ (2.0 * n1)
             if (forcet.eq.1) then
                ttarg = (((5.0* te1 + 3.5 * ti1) * sqrt(te1+ti1))/ (8.5 * sqrt(2.0)))**(2.0/3.0)
                te1 = ttarg
                ti1 = ttarg


             endif

             write(6,'(a,i4,g16.8,2f9.4,3g14.6)')'E2d-9:SOL 22:First:',ir,n1,te1,ti1,v1e2d,vpe2d,vpg
          elseif (actswe2d.ne.0.0) then
             n1 = cellvals(ir,1,2)
             te1= cellvals(ir,2,2)
             ti1= cellvals(ir,3,2)
             v1e2d = cellvals(ir,4,2)
             vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*mb))

             vpg = gtarg(ir,2) / n1
             if (forcet.eq.1) then
                ttarg = (((5.0* te1 + 3.5 * ti1) * sqrt(te1+ti1))/ (8.5 * sqrt(2.0)))**(2.0/3.0)
                te1 = ttarg
                ti1 = ttarg

             endif

             write(6,'(a,i4,g16.8,2f9.4,3g14.6)')'E2d:SOL 22:First:',ir,n1,te1,ti1,v1e2d,vpe2d,vpg

          endif
          errcode = 0

          serr = 0.0

          !        This loop would also be used to load other
          !        quantities - like ionization source or background
          !        neutral densities if these were to be used in the
          !        calcsol subroutine.

          npts = midnks

          if (npts.le.0) then
             call errmsg('MOD_SOL22_INTERFACE:CALCSOL_INTERFACE: NO POINTS ON HALF RING (1ST): NPTS = ',npts)
             return
          endif


          rbnd(0) = krb(0,ir)

          sbnd(0) = ksb(0,ir)

          do ik = 1,midnks

             !            write(6,*) 'in1:',ir,ik,spts(ik)

             spts(ik) = kss2(ik,ir)

             fact = recfrac
             if (actswmajr.eq.2.0) then
                fact = rs(ik,ir)/r0 * fact
             elseif (actswmajr.eq.3.0) then
                fact = r0/rs(ik,ir) * fact

             endif
             if ((actswioni.eq.11.or.actswioni.eq.15).and.(.not.pinavail)) then
                ionsrc(ik) = e2dion(ik,ir) * fact
             else
                ionsrc(ik) = pinion(ik,ir) * fact

             endif
             if (actswrecom.eq.2.and.(.not.pinavail)) then
                recsrc(ik) = e2dhrec(ik,ir)
             else
                recsrc(ik) = pinrec(ik,ir)

             endif

             gperp(ik)  = gperpa(ik,ir)
             if (oldktebs(ik,ir).lt.tcutqe) then
                qesrc(ik)  = 0.0
                write(6,'(a5,i4,2(1x,g13.5))') 'QE:',ik,oldktebs(ik,ir),tcutqe
             else
                qesrc(ik)  = -pinqe(ik,ir) * fact* qesrc_mult

             endif
             if (oldktibs(ik,ir).lt.tcutcx) then
                qisrc(ik)  = 0.0
                write(6,'(a5,i4,2(1x,g13.5))') 'QI:',ik,oldktibs(ik,ir),tcutcx
             else
                qisrc(ik)  = -pinqi(ik,ir) * fact

             endif

             radsrc(ik) = div_cool(ik,ir)
             nhs(ik) = PINATOM(IK,IR)
             nh2s(ik) = PINMOL(ik,ir)
             ths(ik) = pinena(ik,ir) * 0.6666

             ! Assign external electron power term
             if (switch(swepow).ne.0.0) then
                epowsrc(ik) = ext_epowsrc(ik,ir)
             else
                epowsrc = 0.0
             endif

             ! Assign external ion power term
             if (switch(swipow).ne.0.0) then
                ipowsrc(ik) = ext_ipowsrc(ik,ir)
             else
                ipowsrc = 0.0
             endif
             
             !           Load values from last iteration for PINQID calculations

             nhs0(ik) = pinvdist(1,1,ik,ir) + pinvdist(2,1,ik,ir)+ pinvdist(3,1,ik,ir)
             oldne(ik) = oldknbs(ik,ir)
             oldte(ik) = oldktebs(ik,ir)

             oldti(ik) = oldktibs(ik,ir)
             rbnd(ik) = krb(ik,ir)

             sbnd(ik) = ksb(ik,ir)

             momsrc(ik) = pinmp(ik,ir) * fact

             !        The end of the Spts array can either be the cell boundary
             !        or the mid-point of the ring. The cell boundary will result
             !        in more correct integrations since the contribution of the
             !        cells near the mid-point will be calculated more accurately.
             !        (Though the difference would be expected to be small.)
             !        However, the ring midpoint (0.5 * ringlen) makes it easier to
             !        scale or extend sources over the same spatial region by allowing
             !        input to be specified proportional to the ring length.

             !         spts(npts+1) = sbnd(npts)

          end do

          !        Find an estimated S-value of the X-point

          spts(npts+1) = halfringlen

          sxp = 0.0
          do ik = 1,midnks
             if (zs(1,ir).lt.zxp) then
                if (zs(ik,ir).ge.zxp) then
                   sxp = (kss2(ik,ir) - kss2(ik-1,ir)) *(zxp-zs(ik-1,ir))/(zs(ik,ir)-zs(ik-1,ir))+ kss2(ik-1,ir)
                   goto 100
                endif
             elseif (zs(1,ir).gt.zxp) then
                if (zs(ik,ir).lt.zxp) then
                   sxp = (kss2(ik,ir) - kss2(ik-1,ir)) *(zxp-zs(ik-1,ir))/(zs(ik,ir)-zs(ik-1,ir))+ kss2(ik-1,ir)
                   goto 100
                endif
             endif
          end do

100       continue

          !        Set both target power values -

          !        Note:  pae_end and pai_end have
          !        a "-" ve because a "-" sign is used in the embedded power
          !        equations.

          write (6,*) 'Sxp:',ir,nks(ir),kss2(midnks,ir),sxp,ringlen
          pae_start = -elecptarg(ir,2)

          pai_start = -ionptarg(ir,2)
          pae_end = -elecptarg(ir,1)

          !        Private Plasma target power loss re-distribution to main SOL

          !        Turned off in private plasma. Electron then ion.

          pai_end = -ionptarg(ir,1)
          if (actswppelec.eq.0.0.or.pplasma.eq.1) then
             ppelecpow = 0.0

             !           Determine corresponding PP ring to current ring.

          elseif (actswppelec.eq.1.0.or.actswppelec.eq.2.0) then

             ircor = nrs - (ir-irsep)

             !              Turn option off and set power to zero

             if (ircor.le.irtrap) then
                actswppelec = 0.0

                ppelecpow   = 0.0

             else

                ppelecpow = -elecptarg(ircor,2)
             endif

             !           see comment above regarding "-" sign

          elseif (actswppelec.eq.3.0) then

             ppelecpow = -pe_pfzsol(ir,2)

             !        PP Ion power loss term

          endif
          if (actswppion.eq.0.0.or.pplasma.eq.1) then
             ppionpow = 0.0

             !           Determine corresponding PP ring to current ring.

          elseif (actswppion.eq.1.0.or.actswppion.eq.2.0) then

             ircor = nrs - (ir-irsep)

             !              Turn option off and set power to zero

             if (ircor.le.irtrap) then
                actswppion = 0.0

                ppionpow   = 0.0

             else

                ppionpow = -ionptarg(ircor,2)

             endif

             !           see comment above regarding "-" sign

          elseif (actswppelec.eq.3.0) then

             ppionpow = -pi_pfzsol(ir,2)

             !        PP pressure loss term

          endif
          if (actswppress.eq.0.0.or.pplasma.eq.1) then
             pp_press = 0.0

             !           Determine corresponding PP ring to current ring.

          elseif (actswppress.eq.1.0.or.actswppress.eq.2.0) then


             ircor = nrs - (ir-irsep)

             !              Turn option off and set power to zero

             if (ircor.le.irtrap) then
                actswppress = 0.0

                pp_press   = 0.0

             else

                pp_press = presstarg(ircor,2)

             endif

             write(0,*) 'Press2:',actswppress,ir,ircor,pp_press

             !           see comment above regarding "-" sign

          elseif (actswppress.eq.3.0) then

             pp_press = pr_pfzsol(ir,2)

             !        Set up extra perpenicular source parameters

          endif

          !           Calculate source start and stop positions - for
          !           first 1/2 ring use the parameters as specified.

          !           Note: These stop and start positions may be
          !           on either half of the flux tube or even overlap
          !           the top. The integrations procedd outward from
          !           each target so the values of the start and
          !           stop positions will have to be adjusted but
          !           the rest will remain the same.

          if (switch(swextra).gt.0.0) then
             start_gextra_src = gextra_src_start * ringlen
             stop_gextra_src  = gextra_src_stop  * ringlen
             start_gextra_sink= gextra_sink_start * ringlen

             !           Calculate extra source strength


             stop_gextra_sink = gextra_sink_stop * ringlen
             if (start_gextra_src.ne.stop_gextra_src) then
                gextra_src = abs(gtarg(ir,3)) * gextra_mult /abs(start_gextra_src-stop_gextra_src)
             else
                gextra_src = abs(gtarg(ir,3)) * gextra_mult

                !           Calculate extra sink strength

             endif
             if (start_gextra_sink.ne.stop_gextra_sink) then
                gextra_sink = -abs(gtarg(ir,3)) * gextra_mult /abs(start_gextra_sink-stop_gextra_sink)
             else
                gextra_sink = -abs(gtarg(ir,3)) * gextra_mult

             endif

             write (6,*) 'Gextra 1:',ir,gextra_src,gextra_sink,gtarg(ir,3),ringlen,start_gextra_src,stop_gextra_src,&
                  start_gextra_sink,stop_gextra_sink

             !        Error correcting branch point

          endif

          !        Algorithmic ionization options need to come after the error check point so that they are set correctly


          !        Set up the ionization source - if choice option is selected.

150       continue

          !            Analyse target conditions and select appropriate analytic
          !            ionization option

          !            From TN1369

          if (actswion.eq.5.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.5.0.and.&
               (.not.pinavail))) then

             if (n0.gt.1.0e19) then
                if (actswion.eq.5) then
                   actswion = 3.0
                else
                   actswioni = 3.0

                endif
                if (te0.le.1.3) then
                   ssrcst = 13.0 - 10.0 * te0
                   ssrcfi = ssrcst + min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = ssrcst + min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0

                   !                write(0,'(a,i8,20(1x,g12.5))') 'Ion 5a:',ir,n0,te0,
                   !     >                ssrcst, ssrcfi,ssrclen,alg_ion_src_len,halfringlen
                   !
                endif

             elseif (n0.le.1.0e19) then
                if (actswion.eq.5) then
                   actswion = 4.0
                else
                   actswioni= 4.0

                endif
                if (te0.le.10.0) then
                   !                   ssrcfi = min(13.0 - te0,halfringlen)
                   ssrcst = 0.0
                   ssrcfi = min(13.0 - te0,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif

             endif

             write(6,'(a,6g13.5)') 'Ion5O:',actswion,actswioni,ssrcst,ssrcfi,ssrcmid,ssrclen

             !           Adjust parameters for S**5 Gaussian ionization option.

             !           This ionization source always starts at zero.
             !           Thus the quantity in the input for the start position
             !           is used instead as the width factor for the distribution.
             !           If the Start point is listed as 0.0 then this is adjusted
             !           to 1.0 for the width factor for the gaussian.

             !           This code only needs to be executed once on each ring.

          elseif ((actswion.eq.6.0.or.actswion.eq.9.0).or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.&
               (actswioni.eq.6.0.or.actswioni.eq.9).and.(.not.pinavail))) then
             if (ssrcdecay.eq.0.0) then
                s5gausslen = 1.0
             else
                s5gausslen = ssrcdecay

             endif

             ssrcst = 0.0

             ssrclen = ssrcfi - ssrcst

             !            Analyse target conditions and select appropriate analytic
             !            ionization option

             !            From TN1372

          elseif (actswion.eq.7.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.7.0.and.&
               (.not.pinavail))) then

             if (n0.gt.1.0e19) then
                if (actswion.eq.7) then
                   actswion = 6.0
                else
                   actswioni= 6.0

                endif
                if (te0.le.1.3) then
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 14.0 - 10.0 * te0
                else
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.0

                endif

             elseif (n0.le.1.0e19) then
                if (actswion.eq.7) then
                   actswion = 4.0
                else
                   actswioni= 4.0

                endif
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif

             endif

             write(6,'(a,6g13.5)') 'Ion7O:',actswion,actswioni,ssrcst,ssrcfi,ssrcmid,ssrclen

             !            Analyse target conditions and select appropriate analytic
             !            ionization option

             !            From TN1372

          elseif (actswion.eq.10.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.10.0.and.&
               (.not.pinavail))) then

             if (n0.gt.1.0e19) then
                if (actswion.eq.10.0) then
                   actswion = 9.0
                else
                   actswioni= 9.0

                endif
                if (te0.le.1.3) then
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.5*(14.0 - 10.0 * te0)
                else
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.5

                endif

             elseif (n0.le.1.0e19) then
                if (actswion.eq.10.0) then
                   actswion = 4.0
                else
                   actswioni= 4.0

                endif
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif

             endif

             write(6,'(a,6g13.5)') 'Ion10O:',actswion,actswioni,ssrcst,ssrcfi,ssrcmid,ssrclen
          endif
          !

          !        Initialize output arrays to default error values


          call dinit(te,mxspts,10.0d0)
          call dinit(ti,mxspts,10.0d0)
          call dinit(ne,mxspts,1.0d19)

          call dinit(vb,mxspts,0.0d0)
          write (6,'(a,i4,a)') '************   START FIRST HALF ',ir,'  ***********'

          !        Zero power ratio information in case calcsol is not called.

          write (6,*) 'First:',ir,te0,ti0,n0,midnks

          !        Override with detached plasma solution for the outer

          call qzero (int_powrat,3)

          if (switch(swdetach).eq.1.0) then

             !        Solve using SOL22

             call detached_plasma(spts,npts,errcode,serr,te,ti,ne,vb,ir,switch(swdetach),te0,ti0,n0,v0,act_press,mb)

             ! slmod begin - new
             !...Need a system routine for the other system files.  Rename
             !   CLOCK2:
          else
             IF (osm_mode.GE.2) THEN
                deltat = Clock2()
                CALL OpenStorageFiles(ir,IKLO,'tmp.dat')
                ! slmod end
             ENDIF
             ! slmod begin - new

             ! jdemod - move initialization of ierror out of mod_sol22 (calcsol) to the calling routine
             ierror = MAXNKS

             call calcsol (spts,npts,errcode,serr,te,ti,ne,vb,exp_press,act_press,prad,ir,irlim1,int_powrat,cprint,&
                  cve,cvi,cde,cdi)

             IF (osm_mode.GE.2) THEN
                CALL SOL22Status(IKLO,ir,deltat,serr,spts,npts,errcode)
                ! slmod end

             ENDIF

             !        Handle Error Condition if Error Switch is set

          endif

          !           If a negative N error or NaNQ occurs without error correction -
          !           Turn error correction ON.

          if ((   actswerror.ge.1.0.and.(errcode.eq.3.or.errcode.eq.4.or.errcode.eq.5.or.errcode.eq.6.or.errcode.eq.7)).or.&
               (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7))) then
             if (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7).and.seterror.eq.0) then
                actswerror = 10
                seterror = 1
             elseif (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7).and.seterror.eq.1) then
                seterror = 2

                !               ERROR - negative N has been found EVEN with highest level
                !                       of error correction - issue error messages and stop.

             elseif (seterror.eq.2) then

                !                write (7,*) 'SOLASCV: SOL22:'//
                !     >               ' Unsolvable Negative N error encountered'
                !                write (7,*) 'Program will STOP'

                call errmsg('SOLASCV:SOL22',' Unsolvable Negative N error encountered.'//' Program Stopping')

                stop 'SOL22:NEG N'

                !           Record error

             endif
             cdeferr(ir,2) = errcode

             !           Set switches

             cdefserr(ir,2) = serr

             !           Record error level of next attempt

             !            cdeferropt(ir,2) = actswerror

             call setsw(-2,pplasma,new_errlevel)

             cdeferropt(ir,2) = new_errlevel

             !            actffric = find_ffric(ir,2,actlenmom)
             !            call assign_radiation_parameters(ir,2)

             call initlen

             !           Re-do Outer 1/2 ring

             write (6,*) 'Error Handler: SET OUTER:',cdeferr(ir,2),cdefserr(ir,2),cdeferropt(ir,2)

             goto 150

             !        Save background density and temperature
             !        Target first in case it has changed

          endif
          knds(idds(ir,2)) = n0

          !        If E2D option 9 was in use - fill in the missing
          !        plasma with a linear fit.

          kvds(idds(ir,2)) = v0

          !           Fill in missing values - depending on option specified.

          if (actswe2d.eq.9.0.and.switch(swdetach).ne.1.0) then

             if (fillopt.eq.0) then
                do ik = 1,ike2d_start-1
                   te(ik) = te0 + (te(ike2d_start)-te0) *spts(ik)/spts(ike2d_start)
                   ti(ik) = ti0 + (ti(ike2d_start)-ti0) *spts(ik)/spts(ike2d_start)
                   ne(ik) = n0 + (ne(ike2d_start)-n0) *spts(ik)/spts(ike2d_start)
                   vb(ik) = v0 + (vb(ike2d_start)-v0) *spts(ik)/spts(ike2d_start)

                end do

             elseif (fillopt.eq.1) then

                !              Extrapolate target electron temperature

                gam = n0 * v0

                slope = (te(ike2d_start+1) - te(ike2d_start)) /(spts(ike2d_start+1) - spts(ike2d_start))

                te0 = - spts(ike2d_start) * slope + te(ike2d_start)

                !              Extrapolate target ion temperature

                if (te0.le.0.0) te0 = te(ike2d_start)

                slope = (ti(ike2d_start+1) - ti(ike2d_start)) /(spts(ike2d_start+1) - spts(ike2d_start))

                ti0 = - spts(ike2d_start) * slope + ti(ike2d_start)

                if (ti0.le.0.0) ti0 = ti(ike2d_start)

                v0 = -sqrt((te0+ti0)/crmb * econv/mconv)

                n0 = gam/v0
                do ik = 1,ike2d_start-1
                   te(ik) = te0 + (te(ike2d_start)-te0) *spts(ik)/spts(ike2d_start)
                   ti(ik) = ti0 + (ti(ike2d_start)-ti0) *spts(ik)/spts(ike2d_start)
                   ne(ik) = n0 + (ne(ike2d_start)-n0) *spts(ik)/spts(ike2d_start)
                   vb(ik) = v0 + (vb(ike2d_start)-v0) *spts(ik)/spts(ike2d_start)

                   !              Reset target values

                end do
                knds(idds(ir,2)) = n0
                kvds(idds(ir,2)) = v0
                kteds(idds(ir,2)) = te0

                ktids(idds(ir,2)) = ti0

                !              Fill with constant values from last point

                !              Target values are NOT affected

             elseif (fillopt.eq.2) then
                do ik = 1,ike2d_start-1
                   te(ik) = te(ike2d_start)
                   ti(ik) = ti(ike2d_start)
                   ne(ik) = ne(ike2d_start)
                   vb(ik) = vb(ike2d_start)

                end do

                !              Fill with values held constant at target conditions

             elseif (fillopt.eq.3) then
                do ik = 1,ike2d_start-1
                   te(ik) = te0
                   ti(ik) = ti0
                   ne(ik) = n0
                   vb(ik) = v0

                end do

             endif

             !        Assign background ...

          endif

          do ik = 1,midnks
             ktebs(ik,ir) = te(ik)
             ktibs(ik,ir) = ti(ik)
             knbs(ik,ir) = ne(ik)
             oldknbs(ik,ir) = ne(ik)
             oldktibs(ik,ir) = ti(ik)
             oldktebs(ik,ir) = te(ik)
             ! slmod begin - new
             !            oldkvhs(ik,ir) = vb(ik)
             kvhs(ik,ir) = vb(ik)
             osmcde(ik,ir) = cde(ik)
             osmcdi(ik,ir) = cdi(ik)
             osmcve(ik,ir) = cve(ik)
             ! slmod end

             osmcvi(ik,ir) = cvi(ik)
             kpress(ik,ir,1) = exp_press(ik)
             kpress(ik,ir,2) = act_press(ik)

             kprad(ik,ir)    = prad(ik)
             ! slmod begin - new
             !... .TRUE. temporary:
          end do
          ! slmod end

          !        Save mach numbers and error codes

          IF (.NOT.pinavail.AND.(rel_opt.EQ.1.OR.rel_opt.EQ.3))CALL CalcInitSrc(IKLO,ir)
          cmachno(ir,2) = m0
          cerr(ir,2)  = errcode

          !        Save Power Ratios

          cserr(ir,2) = serr
          sol22_power_ratio(ir,2,1) = int_powrat(1)
          sol22_power_ratio(ir,2,2) = int_powrat(2)

          sol22_power_ratio(ir,2,3) = int_powrat(3)

          !        End of 1/2 ring selection

          !           Check to print debugging data
          !           The parameters are current ring snd local IKOPT
          !           data will be printed if it has been collected
          write(6,'(a,i4,3(1x,g12.5))') 'PR2:',ir,int_powrat(1),int_powrat(2),int_powrat(3)
          if (debug_sol22.ne.0) call check_print_data(ir,1)


          !-----------------------------------------------------------------
          !        INNER TARGET
          !-----------------------------------------------------------------


          !        Only execute for the selected half-rings - ikopt

       endif
       if (ikopt.eq.2.or.ikopt.eq.3) then


          !        Check for sol22 debug start

          !           The parameters are current ring snd local IKOPT
          call pr_trace('MOD_CALCSOL_INTERFACE','START INNER TARGET')

          !        Set the values for a call - do 1/2 of the ring at a time.

          !        Do second half of the ring.

          !        Inner target second - i.e. elements 1/2 nks(ir) to nks(ir)
          !                                   target designated as 1

          !        Inner target

          !        Reset switches


          if (debug_sol22.ne.0) call check_init_record_data(ir,2)

          seterror = 0
          if (applydef(ir,1).eq.-1) then

             !           Reset Gperp switch for option 3 in case of under ionization

             call setsw(-1,pplasma,new_errlevel)
             if (( (pinavail.or.((.not.pinavail.and.(actswioni.eq.11.or.actswioni.eq.15)))).and.actswgperp.eq.3.0.and.&
                  rconst(ir).gt.0.0)) then
                actswgperp = 2.0
                write(6,*) 'Swgperp:',rconst(ir),actswgperp
             elseif ((.not.pinavail).and.(actswioni.ne.11.and.actswioni.ne.15).and. (actswgperp.eq.2.0.or.actswgperp.eq.3.0.or.&
                  actswgperp.eq.7.or.actswgperp.eq.8.0)) then
                actswgperp = 0.0
             elseif (actswgperp.eq.7.0.and.gperprat(ir).gt.5.0) then
                actswgperp = 2.0

             endif
          elseif (applydef(ir,1).ge.0) then
             call setsw(applydef(ir,1),pplasma,new_errlevel)

             actswerror = applydef(ir,1)

             if (actswerror.eq.0.0) seterror = 2
             cdeferr(ir,1) =  8
             cdefserr(ir,1) = 0.0
             cdeferropt(ir,1) = applydef(ir,1)

             !        Initialize HALF ring length to be used
             !        - keep in mind that the second half of the ring works down
             !          from KSMAXS2(IR).
             !
             !        halfringlen = ksmaxs2(ir) - ksb(midnks,ir)

          endif
          if (sol22_halfringlen_opt.eq.0) then 
             halfringlen  = 0.5d0*ringlen
          elseif (sol22_halfringlen_opt.eq.1) then 
             ! ring midpoint is upper boundary of middle cell just below midpoint
             !            write(6,'(a,2i8,10(1x,g12.5))') 'Halfringlen2:',
             !     >           midnks,ir,halfringlen,ringlen/2.0
             !            write(0,'(a,2i8,10(1x,g12.5))') 'Halfringlen2:',
             !     >           midnks,ir,halfringlen,ringlen/2.0
             halfringlen = ksmaxs2(ir) - ksb(midnks+1,ir)

             !        Initialize the ionization and radiation source lengths

          endif

          !        Initialize the friction parameter for momentum loss options

          call initlen
          actffric = find_ffric(ir,1,actlenmom,actlammom)

          !        Initialize major radius and target condition data

          call assign_radiation_parameters(ir,1)

          targfact = 1.0
          if (actswmajr.eq.1.0) then
             if (actswe2d.eq.0.0) then
                id = idds(ir,1)
                targfact = rp(id) / r0
             elseif (actswe2d.ne.0.0) then
                targfact = rs(nks(ir),ir)/r0
             endif


             !        The solvers always deal with 1/2 a ring at a time and need a
             !        target velocity < 0 for all targets.


          endif
          te0 = kteds(idds(ir,1))
          ti0 = ktids(idds(ir,1))
          n0 = knds(idds(ir,1))

          v0 = -kvds(idds(ir,1))
          if (forcet.eq.1) then
             ttarg = (((5.0* te0 + 3.5 * ti0) * sqrt(te0+ti0))/ (8.5 * sqrt(2.0)))**(2.0/3.0)
             te0 = ttarg
             ti0 = ttarg
             kteds(idds(ir,1)) = te0
             ktids(idds(ir,1)) = ti0

          endif
          if (actswe2d.eq.9.0) then
             n1 = e2dnbs(nks(ir)-ike2d_start+1,ir)
             te1= e2dtebs(nks(ir)-ike2d_start+1,ir)
             ti1= e2dtibs(nks(ir)-ike2d_start+1,ir)
             v1e2d = -e2dvhs(nks(ir)-ike2d_start+1,ir)

             vpe2d = v1e2d

             vpg = -(e2dgpara(nks(ir)-ike2d_start+1,ir)+e2dgpara(nks(ir)-ike2d_start+1+1,ir))/ (2.0 * n1)
             if (forcet.eq.1) then
                ttarg = (((5.0* te1 + 3.5 * ti1) * sqrt(te1+ti1))/ (8.5 * sqrt(2.0)))**(2.0/3.0)
                te1 = ttarg
                ti1 = ttarg

             endif

             write(6,'(a,i4,g16.8,2f9.4,3g14.6)')'E2d-9:SOL 22:First:',ir,n1,te1,ti1,v1e2d,vpe2d,vpg
          elseif (actswe2d.ne.0.0) then
             n1 = cellvals(ir,1,1)
             te1= cellvals(ir,2,1)
             ti1= cellvals(ir,3,1)
             v1e2d = -cellvals(ir,4,1)
             vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*mb))

             vpg = gtarg(ir,1) / n1
             if (forcet.eq.1) then
                ttarg = (((5.0* te1 + 3.5 * ti1) * sqrt(te1+ti1))/ (8.5 * sqrt(2.0)))**(2.0/3.0)
                te1 = ttarg
                ti1 = ttarg

             endif

             write(6,'(a,i4,g16.8,2f9.4,3g14.6)')'E2d:SOL 22:First:',ir,n1,te1,ti1,v1e2d,vpe2d,vpg

          endif
          errcode = 0

          serr    = 0.0

          npts = nks(ir) - midnks

          if (npts.le.0) then
             call errmsg('MOD_SOL22_INTERFACE:CALCSOL_INTERFACE: NO POINTS ON HALF RING (2ND): NPTS = ',npts)
             return
          endif



          do ik = nks(ir), midnks + 1 , -1

             !           write(6,*) 'in2:',ir,ik,nks(ir)-ik+1,
             !     >                   spts(nks(ir)-ik+1)

             spts(nks(ir) - ik + 1) =  ksmaxs2(ir) - kss2(ik,ir)

             fact = recfrac
             if (actswmajr.eq.2.0) then
                fact = rs(ik,ir)/r0 * fact
             elseif (actswmajr.eq.3.0) then
                fact = r0/rs(ik,ir) * fact

             endif
             if ((actswioni.eq.11.or.actswioni.eq.15).and.(.not.pinavail)) then
                ionsrc(nks(ir)-ik+1) = e2dion(ik,ir) * fact
             else
                ionsrc(nks(ir)-ik+1) = pinion(ik,ir) * fact

             endif
             if (actswrecom.eq.2.and.(.not.pinavail)) then
                recsrc(nks(ir)-ik+1) = e2dhrec(ik,ir)
             else
                recsrc(nks(ir)-ik+1) = pinrec(ik,ir)

             endif

             gperp(nks(ir)-ik+1)  = gperpa(ik,ir)
             if (oldktebs(ik,ir).lt.tcutqe) then
                qesrc(nks(ir)-ik+1)  = 0.0
                write(6,'(a5,i4,2(1x,g13.5))') 'QE:',ik,oldktebs(ik,ir),tcutqe
             else
                qesrc(nks(ir)-ik+1)  = -pinqe(ik,ir) * fact* qesrc_mult

             endif
             if (oldktibs(ik,ir).lt.tcutcx) then
                qisrc(nks(ir)-ik+1)  = 0.0
                write(6,'(a5,i4,2(1x,g13.5))') 'QI:',ik,oldktibs(ik,ir),tcutcx
             else
                qisrc(nks(ir)-ik+1)  = -pinqi(ik,ir) * fact

             endif

             radsrc(ik) = div_cool(ik,ir)
             nhs(nks(ir)-ik+1) = PINATOM(IK,IR)
             nh2s(nks(ir)-ik+1) = PINMOL(ik,ir)
             ths(nks(ir)-ik+1) = pinena(ik,ir) * 0.6666

             ! Assign external electron power term
             if (switch(swepow).ne.0.0) then
                epowsrc(ik) = ext_epowsrc(ik,ir)
             else
                epowsrc = 0.0
             endif

             ! Assign external ion power term
             if (switch(swipow).ne.0.0) then
                ipowsrc(ik) = ext_ipowsrc(ik,ir)
             else
                ipowsrc = 0.0
             endif
                          
             !           Load values from last iteration for PINQID calculations

             nhs0(nks(ir)-ik+1) = pinvdist(1,1,ik,ir)+ pinvdist(2,1,ik,ir) + pinvdist(3,1,ik,ir)
             oldne(nks(ir)-ik+1) = oldknbs(ik,ir)
             oldte(nks(ir)-ik+1) = oldktebs(ik,ir)


             !           Note that the rbnd array is shifted by one - it starts
             !           with a 0 index.

             oldti(nks(ir)-ik+1) = oldktibs(ik,ir)
             rbnd(nks(ir)-ik) = krb(ik,ir)

             sbnd(nks(ir)-ik) = ksmaxs2(ir) - ksb(ik,ir)

             momsrc(nks(ir)-ik+1) = - pinmp(ik,ir) * fact

             !        Add last rbnd and sbnd elements.


          end do
          rbnd(nks(ir)-midnks) = krb(midnks,ir)

          !        The end of the Spts array can either be the cell boundary
          !        or the mid-point of the ring. The cell boundary will result
          !        in more correct integrations since the contribution of the
          !        cells near the mid-point will be calculated more accurately.
          !        (Though the difference would be expected to be small.)
          !        However, the ring midpoint (0.5 * ringlen) makes it easier to
          !        scale or extend sources over the same spatial region by allowing
          !        input to be specified proportional to the ring length.

          sbnd(nks(ir)-midnks) = ksmaxs2(ir) - ksb(midnks,ir)

          !        Find an estimated S-value of the X-point

          spts(npts+1) = halfringlen

          sxp = 0.0
          do ik = nks(ir), midnks + 1 , -1
             if (zs(nks(ir),ir).lt.zxp) then
                if (zs(ik,ir).ge.zxp) then
                   sxp = (kss2(ik+1,ir) - kss2(ik,ir)) *(zxp-zs(ik+1,ir))/(zs(ik,ir)-zs(ik+1,ir))+ (ksmaxs2(ir) - kss2(ik+1,ir))
                   goto 200
                endif
             elseif (zs(nks(ir),ir).gt.zxp) then
                if (zs(ik,ir).lt.zxp) then
                   sxp = (kss2(ik+1,ir) - kss2(ik,ir)) *(zxp-zs(ik+1,ir))/(zs(ik,ir)-zs(ik+1,ir))+ (ksmaxs2(ir)-kss2(ik+1,ir))
                   goto 200
                endif
             endif
          end do

200       continue


          !        Set both target power values

          !        Note:  pae_end and pai_end have
          !        a "-" ve because a "-" sign is used in the embedded power
          !        equations.

          write (6,*) 'Sxp:',ir,nks(ir),kss2(midnks,ir),sxp,ringlen
          pae_start = -elecptarg(ir,1)

          pai_start = -ionptarg(ir,1)
          pae_end = -elecptarg(ir,2)


          !        Private Plasma target power loss re-distribution to main SOL

          !        Turned off in private plasma. Electron then ion.

          pai_end = -ionptarg(ir,2)
          if (actswppelec.eq.0.0.or.pplasma.eq.1) then
             ppelecpow = 0.0

             !           Determine corresponding PP ring to current ring.

          elseif (actswppelec.eq.1.0.or.actswppelec.eq.2.0) then

             ircor = nrs - (ir-irsep)

             !              Turn option off and set power to zero

             if (ircor.le.irtrap) then
                actswppelec = 0.0

                ppelecpow   = 0.0

             else

                ppelecpow = -elecptarg(ircor,1)

             endif

             !        PP Ion power loss term

          endif
          if (actswppion.eq.0.0.or.pplasma.eq.1) then
             ppionpow = 0.0

             !           Determine corresponding PP ring to current ring.

          elseif (actswppion.eq.1.0.or.actswppion.eq.2.0) then

             ircor = nrs - (ir-irsep)

             !              Turn option off and set power to zero

             if (ircor.le.irtrap) then
                actswppion = 0.0

                ppionpow   = 0.0

             else

                ppionpow = -ionptarg(ircor,1)

             endif

             !        PP pressure loss term

          endif
          if (actswppress.eq.0.0.or.pplasma.eq.1) then
             pp_press = 0.0

             !           Determine corresponding PP ring to current ring.

          elseif (actswppress.eq.1.0.or.actswppress.eq.2.0) then

             ircor = nrs - (ir-irsep)

             !              Turn option off and set power to zero

             if (ircor.le.irtrap) then
                actswppress = 0.0

                pp_press   = 0.0

             else

                pp_press = presstarg(ircor,1)
             endif

             write(0,*) 'Press3:',actswppress,ir,ircor,pp_press


          endif

          !        Set up extra perpenicular source parameters


          !           Calculate source start and stop positions - for
          !           first 1/2 ring use the parameters as specified.

          !           Note: These stop and start positions may be
          !           on either half of the flux tube or even overlap
          !           the top. The integrations procedd outward from
          !           each target so the values of the start and
          !           stop positions will have to be adjusted but
          !           the rest will remain the same.

          !           For the second half ring - looking from the other
          !           target - the positions of "start" and "stop" are
          !           reversed.

          !           e.g.               SOURCE
          !                  start  stop
          !           0        |      |                          SMAX
          !           |-------------------------------------------|
          !                                   |         |
          !                                 start     stop
          !                               SINK

          !           But for S measured from the SMAX end of the
          !           ring  SINK_stop becomes the start of the SINK
          !           region and so on for the rest of the positions.


          if (switch(swextra).gt.0.0) then
             start_gextra_src = (1.0-gextra_src_stop) * ringlen
             stop_gextra_src  = (1.0-gextra_src_start)  * ringlen
             start_gextra_sink= (1.0-gextra_sink_stop)* ringlen

             !           Calculate extra source strength

             stop_gextra_sink = (1.0-gextra_sink_start) * ringlen
             if (start_gextra_src.ne.stop_gextra_src) then
                gextra_src = abs(gtarg(ir,3)) * gextra_mult /abs(start_gextra_src-stop_gextra_src)
             else
                gextra_src = abs(gtarg(ir,3)) * gextra_mult

                !           Calculate extra sink strength

             endif
             if (start_gextra_sink.ne.stop_gextra_sink) then
                gextra_sink = -abs(gtarg(ir,3)) * gextra_mult /abs(start_gextra_sink-stop_gextra_sink)
             else
                gextra_sink = -abs(gtarg(ir,3)) * gextra_mult

             endif

             write (6,*) 'Gextra 1:',ir,gextra_src,gextra_sink,gtarg(ir,3),ringlen,start_gextra_src,stop_gextra_src,&
                  start_gextra_sink,stop_gextra_sink

             !        Error correcting branch point

          endif

          !        Algorithmic ionization options need to come after the error check point so that they are set correctly


          !        Set up the ionization source - if choice option is selected.

250       continue

          !            Analyse target conditions and select appropriate analytic
          !            ionization option

          !            From TN1369

          if (actswion.eq.5.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.5.0.and.&
               (.not.pinavail))) then

             if (n0.gt.1.0e19) then
                if (actswion.eq.5) then
                   actswion = 3.0
                else
                   actswioni= 3.0

                endif
                if (te0.le.1.3) then
                   ssrcst = 13.0 - 10.0 * te0
                   ssrcfi = ssrcst + min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = ssrcst + min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0

                   !                write(0,'(a,i8,20(1x,g12.5))') 'Ion 5b:',ir,n0,te0,
                   !     >                ssrcst, ssrcfi,ssrclen,alg_ion_src_len,halfringlen
                   !

                endif

             elseif (n0.le.1.0e19) then
                if (actswion.eq.5) then
                   actswion = 4.0
                else
                   actswioni= 4.0

                endif
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = min(13.0 - te0,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif

             endif

             write(6,'(a,6g13.5)') 'Ion5I:',actswion,actswioni,ssrcst,ssrcfi,ssrcmid,ssrclen

             !           Adjust parameters for S**5 Gaussian ionization option.

             !           This ionization source always starts at zero.
             !           Thus the quantity in the input for the start position
             !           is used instead as the width factor for the distribution.
             !           If the Start point is listed as 0.0 then this is adjusted
             !           to 1.0 for the width factor for the gaussian.

             !           This code only needs to be executed once on each ring.

          elseif ((actswion.eq.6.0.or.actswion.eq.9.0).or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.&
               (actswioni.eq.6.0.or.actswioni.eq.9).and.(.not.pinavail))) then
             if (ssrcdecay.eq.0.0) then
                s5gausslen = 1.0
             else
                s5gausslen = ssrcdecay

             endif

             ssrcst = 0.0

             ssrclen = ssrcfi - ssrcst

             !            Analyse target conditions and select appropriate analytic
             !            ionization option

             !            From TN1372

          elseif (actswion.eq.7.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.7.0.and.&
               (.not.pinavail))) then

             if (n0.gt.1.0e19) then
                if (actswion.eq.7) then
                   actswion = 6.0
                else
                   actswioni= 6.0

                endif
                if (te0.le.1.3) then
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 14.0 - 10.0 * te0
                else
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.0

                endif

             elseif (n0.le.1.0e19) then
                if (actswion.eq.7) then
                   actswion = 4.0
                else
                   actswioni= 4.0

                endif
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif

             endif
             write(6,'(a,6g13.5)') 'Ion7I:',actswion,actswioni,ssrcst,ssrcfi,ssrcmid,ssrclen

             !            Analyse target conditions and select appropriate analytic
             !            ionization option

             !            From TN1372

          elseif (actswion.eq.10.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.10.0.and.&
               (.not.pinavail))) then

             if (n0.gt.1.0e19) then
                if (actswion.eq.10.0) then
                   actswion = 9.0
                else
                   actswioni= 9.0

                endif
                if (te0.le.1.3) then
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.5* (14.0 - 10.0 * te0)
                else
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.5

                endif

             elseif (n0.le.1.0e19) then
                if (actswion.eq.10.0) then
                   actswion = 4.0
                else
                   actswioni= 4.0

                endif
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = min(alg_ion_src_len,halfringlen)
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif

             endif

             write(6,'(a,6g13.5)') 'Ion10I:',actswion,actswioni,ssrcst,ssrcfi,ssrcmid,ssrclen
             !
             !        Initialize output arrays to default error values

          endif
          call dinit(te,mxspts,10.0d0)
          call dinit(ti,mxspts,10.0d0)
          call dinit(ne,mxspts,1.0d19)

          call dinit(vb,mxspts,0.0d0)
          write (6,'(a,i4,a)') '************   START SECOND HALF ',ir,'  ***********'

          !        Zero power ratio information in case calcsol is not called.

          write (6,*) 'Next :',ir,te0,ti0,n0,npts,nks(ir)

          !        Override with detached plasma solution for the inner

          call qzero (int_powrat,3)

          if (switch(swdetach).eq.2.0) then

             !        Solve using SOL22

             call detached_plasma(spts,npts,errcode,serr,te,ti,ne,vb,ir,switch(swdetach),te0,ti0,n0,v0,act_press,mb)

             ! slmod begin - new
             !...Need a system routine for the other system files.  Rename
             !   CLOCK2:
          else
             IF (osm_mode.EQ.2) THEN
                CALL OpenStorageFiles(ir,IKLO,'tmp.dat')
                deltat = Clock2()
                ! slmod end
             ENDIF
             ! slmod begin - new

             ! jdemod - move initialization of ierror out of mod_sol22 (calcsol) to the calling routine
             ierror = MAXNKS

             call calcsol (spts,npts,errcode,serr,te,ti,ne,vb,exp_press,act_press,prad,ir,irlim1,int_powrat,cprint,&
                  cve,cvi,cde,cdi)

             IF (osm_mode.EQ.2) THEN
                CALL SOL22Status(IKHI,ir,deltat,serr,spts,npts,errcode)
                ! slmod end

             ENDIF

             !        Handle Error Condition if Error Switch is set

          endif

          !           If a negative N error occurs without error correction -
          !           Turn error correction ON.

          if ((   actswerror.ge.1.0.and.(errcode.eq.3.or.errcode.eq.4.or.errcode.eq.5.or.errcode.eq.6.or.errcode.eq.7)).or.&
               (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7))) then
             if (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7).and.seterror.eq.0) then
                actswerror = 5.0
                seterror = 1
             elseif (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7).and.seterror.eq.1) then
                seterror = 2

                !               ERROR - negative N has been found EVEN with highest level
                !                       of error correction - issue error messages and stop.

             elseif (seterror.eq.2) then
                write (6,*) 'SOLASCV: SOL22:'//' Unsolvable Negative N error encountered'

                write (6,*) 'Program will STOP'
                write (7,*) 'SOLASCV: SOL22:'//' Unsolvable Negative N error encountered'

                write (7,*) 'Program will STOP'

                stop

                !           Record error

             endif
             cdeferr(ir,1) = errcode

             !           Set switches

             cdefserr(ir,1) = serr

             !           Record error level of next attempt

             !            cdeferropt(ir,1) = actswerror

             call setsw(-2,pplasma,new_errlevel)


             cdeferropt(ir,1) = new_errlevel

             !            actffric = find_ffric(ir,1,actlenmom)
             !            call assign_radiation_parameters(ir,1)

             call initlen

             !           Re-do Inner 1/2 ring

             write (6,*) 'Error Handler: SET INNER:',cdeferr(ir,1),cdefserr(ir,1),cdeferropt(ir,1)

             goto 250

             !        Save background density and temperature

          endif
          knds(idds(ir,1)) =  n0

          !        If E2D option 9 was in use - fill in the missing
          !        plasma with a linear fit.

          kvds(idds(ir,1)) = -v0

          !           Fill in missing values - depending on option specified.

          if (actswe2d.eq.9.0.and.switch(swdetach).ne.2.0) then

             if (fillopt.eq.0) then
                do ik = 1,ike2d_start-1
                   te(ik) = te0 + (te(ike2d_start)-te0) *spts(ik)/spts(ike2d_start)
                   ti(ik) = ti0 + (ti(ike2d_start)-ti0) *spts(ik)/spts(ike2d_start)
                   ne(ik) = n0 + (ne(ike2d_start)-n0) *spts(ik)/spts(ike2d_start)
                   vb(ik) = v0 + (vb(ike2d_start)-v0) *spts(ik)/spts(ike2d_start)

                end do

             elseif (fillopt.eq.1) then

                !              Extrapolate target electron temperature

                gam = n0 * v0

                slope = (te(ike2d_start+1) - te(ike2d_start)) /(spts(ike2d_start+1) - spts(ike2d_start))

                te0 = - spts(ike2d_start) * slope + te(ike2d_start)

                !              Extrapolate target ion temperature

                if (te0.le.0.0) te0 = te(ike2d_start)

                slope = (ti(ike2d_start+1) - ti(ike2d_start)) /(spts(ike2d_start+1) - spts(ike2d_start))

                ti0 = - spts(ike2d_start) * slope + ti(ike2d_start)

                if (ti0.le.0.0) ti0 = ti(ike2d_start)

                v0 = -sqrt((te0+ti0)/crmb * econv/mconv)

                n0 = gam/v0
                do ik = 1,ike2d_start-1
                   te(ik) = te0 + (te(ike2d_start)-te0) *spts(ik)/spts(ike2d_start)
                   ti(ik) = ti0 + (ti(ike2d_start)-ti0) *spts(ik)/spts(ike2d_start)
                   ne(ik) = n0 + (ne(ike2d_start)-n0) *spts(ik)/spts(ike2d_start)
                   vb(ik) = v0 + (vb(ike2d_start)-v0) *spts(ik)/spts(ike2d_start)

                   !              Reset target values

                end do
                knds(idds(ir,1)) = n0
                kvds(idds(ir,1)) = -v0
                kteds(idds(ir,1)) = te0

                ktids(idds(ir,1)) = ti0

                !              Fill with constant values from last point

                !              Target values are NOT affected

             elseif (fillopt.eq.2) then
                do ik = 1,ike2d_start-1
                   te(ik) = te(ike2d_start)
                   ti(ik) = ti(ike2d_start)
                   ne(ik) = ne(ike2d_start)
                   vb(ik) = vb(ike2d_start)

                end do

                !              Fill with values held constant at target conditions

             elseif (fillopt.eq.3) then
                do ik = 1,ike2d_start-1
                   te(ik) = te0
                   ti(ik) = ti0
                   ne(ik) = n0
                   vb(ik) = v0

                end do

             endif

             !        Assign Background

          endif

          do ik = nks(ir), midnks + 1 , -1
             ktebs(ik,ir) = te(nks(ir) - ik + 1)
             ktibs(ik,ir) = ti(nks(ir) - ik + 1)
             knbs(ik,ir) = ne(nks(ir) - ik + 1)
             oldknbs(ik,ir) = ne(nks(ir) - ik + 1)
             oldktibs(ik,ir) = ti(nks(ir) - ik + 1)
             oldktebs(ik,ir) = te(nks(ir) - ik + 1)
             ! slmod begin - new
             !            oldkvhs(ik,ir) = -vb(nks(ir) - ik + 1)
             kvhs(ik,ir) = -vb(nks(ir) - ik + 1)
             osmcde(ik,ir) = -cde(nks(ir)-ik+1)
             osmcdi(ik,ir) = -cdi(nks(ir)-ik+1)
             osmcve(ik,ir) = -cve(nks(ir)-ik+1)
             ! slmod end

             osmcvi(ik,ir) = -cvi(nks(ir)-ik+1)
             kpress(ik,ir,1) = exp_press(nks(ir) - ik + 1)
             kpress(ik,ir,2) = act_press(nks(ir) - ik + 1)

             kprad(ik,ir)    = prad(nks(ir)-ik+1)
             ! slmod begin - new
             !... .TRUE. temporary:
          end do
          ! slmod end

          !        Save mach numbers and error codes

          IF (.NOT.pinavail.AND.(rel_opt.EQ.1.OR.rel_opt.EQ.3))CALL CalcInitSrc(IKHI,ir)
          cmachno(ir,1) = m0
          cerr(ir,1) = errcode

          !        Save Power Ratios

          cserr(ir,1) = serr
          sol22_power_ratio(ir,1,1) = int_powrat(1)
          sol22_power_ratio(ir,1,2) = int_powrat(2)

          sol22_power_ratio(ir,1,3) = int_powrat(3)

          !        End of 1/2 ring selection

          !           The parameters are current ring snd local IKOPT
          !           data will be printed if it has been collected
          write(6,'(a,i4,3(1x,g12.5))') 'PR1:',ir,int_powrat(1),int_powrat(2),int_powrat(3)
          if (debug_sol22.ne.0) call check_print_data(ir,2)
       endif

       !       CALCULATE ELECTRIC FIELD


       !       IN THE FOLLOWING EQUATIONS THE FACTOR E CANCELS WITH THE
       !       SAME FACTOR USED IN CONVERTING T IN EV TO KT.

       !       If OFIELD is turned ON then set the electric field to
       !       zero for this case. Note: the electric field is
       !       initialized to zero in the plasma.d3a module.

       !       For OFIELD ... 0=off   1=on

       call pr_trace('MOD_CALCSOL_INTERFACE','CALCULATE EFIELD')


       if (ofield.eq.0) then
          if (ikopt.eq.1.or.ikopt.eq.3) then
             if (kss2(1,ir).eq.0.0) then
                DS1 = KSS2(2,IR) - KSS2(1,IR)
                DP1 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
                DT1 = (KTEBS(2,IR)-KTEBS(1,IR))

                NB1 = 0.5*(KNBS(2,IR)+KNBS(1,IR))

                KES(1,IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
                KEDS(idds(ir,2))= kes(1,ir)
             else
                DS2 = KSS2(2,IR) - KSS2(1,IR)
                DP2 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
                DT2 = (KTEBS(2,IR)-KTEBS(1,IR))

                NB2 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
                DS1 = KSS2(1,IR)
                DP1 = KNBS(1,IR)*KTEBS(1,IR)-KNDS(idds(ir,2))*KTEDS(idds(ir,2))
                DT1 = KTEBS(1,IR)-KTEDS(idds(ir,2))

                NB1 = 0.5*(KNBS(1,IR)+KNDS(idds(ir,2)))

                KES(1,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)+ (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
                KEDS(idds(ir,2)) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1

             endif


          endif

          if (ikopt.eq.2.or.ikopt.eq.3) then
             if (kss2(nks(ir),ir).eq.ksmaxs2(ir)) then
                DS1 = KSS2(NKS(IR),IR) - KSS2(NKS(IR)-1,IR)
                DP1 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)-KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
                DT1 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))

                NB1 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))

                KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
                KEDS(idds(ir,1))= kes(nks(ir),ir)

             else
                DS2 = KSS2(NKS(IR),IR) - KSS2(NKS(IR)-1,IR)
                DP2 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)-KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
                DT2 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))

                NB2 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))
                DS1 = ksmaxs2(ir) - KSS2(nks(ir),IR)
                DP1 = KNDS(idds(ir,1))*KTEDS(idds(ir,1))-KNBS(nks(ir),IR)*KTEBS(nks(ir),IR)
                DT1 = KTEDS(idds(ir,1))-KTEBS(nks(ir),IR)

                NB1 = 0.5*(KNBS(nks(ir),IR)+KNDS(idds(ir,1)))

                KES(nks(ir),IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)+ (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
                KEDS(idds(ir,1)) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1

             endif

          endif
          if (ikopt.eq.1) then
             ikfirst = ikstart +1
             iklast = ikmids(ir)
          elseif (ikopt.eq.2) then
             ikfirst = ikmids(ir) + 1
             iklast = ikend - 1
          elseif (ikopt.eq.3) then
             ikfirst = ikstart +1
             iklast = ikend -1

          endif
          DO 500 IK = ikfirst,iklast
             DS1 = KSS2(IK,IR) - KSS2(IK-1,IR)
             DP1 = KNBS(IK,IR)*KTEBS(IK,IR)-KNBS(IK-1,IR)*KTEBS(IK-1,IR)
             DT1 = (KTEBS(IK,IR)-KTEBS(IK-1,IR))
             NB1 = 0.5*(KNBS(IK,IR)+KNBS(IK-1,IR))
             DS2 = KSS2(IK+1,IR) - KSS2(IK,IR)
             DP2 = KNBS(IK+1,IR)*KTEBS(IK+1,IR)-KNBS(IK,IR)*KTEBS(IK,IR)
             DT2 = (KTEBS(IK+1,IR)-KTEBS(IK,IR))
             NB2 = 0.5*(KNBS(IK+1,IR)+KNBS(IK,IR))

             !            WRITE(6,*) 'KES:',IK,IR,KES(IK,IR)

             KES(IK,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)+ (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))

             !         End of test on OFIELD

500          CONTINUE
             write(6,*) 'KES:',keds(idds(ir,2)),keds(idds(ir,1)),idds(ir,2),idds(ir,1),ir
             do ik =1,nks(ir)
                write (6,*) ik,kes(ik,ir)

             end do

             !       If smoothing is turned ON

          endif

          !        WF'95:    CORRECTING(?) FOR PLATE ASSYMETRY

          !        profile scaling of ne,Te,Ti to match midpoint values,
          !        exponent asmexp determines smoothing range

          if (switch(swsmooth).eq.1.0) then
             TEMID = 0.5*(KTEBS(MIDNKS+1,IR) + KTEBS(MIDNKS,IR))
             TIMID = 0.5*(KTIBS(MIDNKS+1,IR) + KTIBS(MIDNKS,IR))
             !         VMID = 0.5*(KVHS(MIDNKS+1,IR) + KVHS(MIDNKS,IR))
             NMID = 0.5*(KNBS(MIDNKS+1,IR) + KNBS(MIDNKS,IR))
             SMAX = KSMAXS2(IR)
             !         write(6,*) 'assymetry', ir,ktibs(midnks,ir),
             !     >   ktibs(midnks+1,ir), timid

             ASMEXP = 1.0
             DO IK = ikstart , ikend
                IF (IK.LE.MIDNKS) THEN
                   TMP = (KSS2(IK,IR) / (SMAX / 2.0))**ASMEXP
                   KTEBS(IK,IR) = KTEBS(IK,IR)*(1.0 +TMP*(TEMID/KTEBS(MIDNKS,IR) - 1.0))
                   KTIBS(IK,IR) = KTIBS(IK,IR)*(1.0 +TMP*(TIMID/KTIBS(MIDNKS,IR) - 1.0))
                   !           KVHS(IK,IR) = KVHS(IK,IR)*(1.0 +
                   !     >        TMP*(VMID/KVHS(MIDNKS,IR) - 1.0))
                   KNBS(IK,IR) = KNBS(IK,IR)*(1.0 +TMP*(NMID/KNBS(MIDNKS,IR) - 1.0))
                ELSE
                   TMP = ((SMAX-KSS2(IK,IR)) / (SMAX / 2.0))**ASMEXP
                   KTEBS(IK,IR) = KTEBS(IK,IR)*(1.0 +TMP*(TEMID/KTEBS(MIDNKS+1,IR) - 1.0))
                   KTIBS(IK,IR) = KTIBS(IK,IR)*(1.0 +TMP*(TIMID/KTIBS(MIDNKS+1,IR) - 1.0))
                   !           KVHS(IK,IR) = KVHS(IK,IR)*(1.0 +
                   !     >        TMP*(VMID/KVHS(MIDNKS+1,IR) - 1.0))
                   KNBS(IK,IR) = KNBS(IK,IR)*(1.0 +TMP*(NMID/KNBS(MIDNKS+1,IR) - 1.0))
                ENDIF

                !      End of Smoothing IF

             ENDDO

             !      Branch location for virtual rings

          endif

          !     End of DO loop for IR

1000      continue

          !     Reset the ionization options for swion=8.0

       end do
       if (switch(swioni).eq.8.0.and.pinavail) then
          switch(swioni) = switch(swion)
          switch(swion) = 8.0

          !     End interface routine

       endif

       !     Call routine to print out table of values to file for Excel plotting
       !     (req. by Wojciech Fundamenski)

       call pr_trace('MOD_CALCSOL_INTERFACE','BEFORE EXCEL PRINT')

       call print_sol_excel
       return



     end subroutine calcsol_interface






     subroutine calcfluxes(gtarg,ionptarg,elecptarg,e2dgtarg,presstarg,gamcor,gamecor,ike2d_start,g_pfzsol,pe_pfzsol,&
          pi_pfzsol,pr_pfzsol,pfz_dist_opt,pfz_dist_param)
       use debug_options
       use mod_params
       use mod_solparams
       use mod_solswitch
       use mod_comtor
       use mod_cgeom
       use mod_pindata
       use mod_cedge2d
       !     include 'params'
       !     include 'solparams'
       !     include 'solswitch'

       !     include 'comtor'
       !     include 'cgeom'
       !     include 'pindata'
       !     include 'cedge2d'

       implicit none
       real*8 gtarg(mxspts,3)
       real*8 e2dgtarg(mxspts,3)
       real*8 ionptarg(mxspts,3)
       real*8 elecptarg(mxspts,3)
       real*8 presstarg(mxspts,3)
       real*8 gamcor,gamecor
       integer ike2d_start
       ! record total actual flux to sol and pfz target regions - inner and outer
       real*8 :: solpfz_fluxes(2,2,4) ! note: this code only works properly for single null/non-extended geometries
       ! distribute extra sources to sol from pfz
       real*8 :: g_pfzsol(maxpts,3)
       real*8 :: pe_pfzsol(maxpts,3)
       real*8 :: pi_pfzsol(maxpts,3)
       real*8 :: pr_pfzsol(maxpts,3)
       integer :: pfz_dist_opt
       real*8 :: pfz_dist_param(2)

       !     This subroutine calculates the target partcle and
       !     power fluxes based on the target data and the
       !     SOL 22 options selected.


       !     Local Variables


       real*8 :: dist_fact(maxpts,2)
       real*8 :: maxpress(2)
       integer :: maxpress_ir(2)
       integer ik,ir,in,id,ierr,ringno
       external ringno
       real rfact,v0
       real n1,te1,ti1,vpe2d,v1e2d,te0
       real mach0o,mach0i,rmeano,rmeani

       real gae,gaio,gaii

       !     Initialization

       call pr_trace('MOD_SOL22_INTERFACE','START CALCFLUXES')

       solpfz_fluxes = 0.0

       !        Calculate mach numbers -

       do ir = irsep,nrs

          !           Outer

          if (switch(swmach).eq.3.0) then

             v0 = -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb* econv/mconv)
             if (v0.ne.0.0) then 
                mach0o = kvds(idds(ir,2)) / v0
             else
                mach0o = 1.0

                !           Inner

             endif

             v0 = -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb* econv/mconv)
             if (v0.ne.0.0) then 
                mach0i = kvds(idds(ir,1)) / v0
             else
                mach0i = 1.0

             endif
          else
             mach0o = 1.0
             mach0i = 1.0

             !        Calculate heat transmission

          endif
          gae = 5.0 + gamecor
          gaio = 2.5 +  0.5* mach0o**2* (1.0 + kteds(idds(ir,2))/ktids(idds(ir,2)))+ gamcor

          !        Calculate target fluxes (both particles and heat)
          !        Assume target MACH number = 1.0

          !        Outer

          gaii = 2.5 +  0.5* mach0i**2* (1.0 + kteds(idds(ir,1))/ktids(idds(ir,1)))+ gamcor

          if (switch(swe2d).eq.0.0.or.switch(swe2d).eq.5) then
             if (switch(swmajr).eq.4.0) then
                rmeano = rp(idds(ir,2))
                rmeani = rp(idds(ir,1))
             else
                rmeano = 1.0
                rmeani = 1.0

                !           v0 = -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb
                !     >             * econv/mconv)

             endif

             v0 = -abs(kvds(idds(ir,2)))

             gtarg(ir,2)= knds(idds(ir,2)) * v0 * rmeano
             if (e2dtargopt.eq.5) then
                elecptarg(ir,2)= e2dpepara(1,ir)
                ionptarg(ir,2) = e2dpipara(1,ir)
             else
                elecptarg(ir,2)=gae*kteds(idds(ir,2))*econv*gtarg(ir,2)
                ionptarg(ir,2)=gaio*ktids(idds(ir,2))*econv*gtarg(ir,2)
             endif
             ! calculate target pressrure

             !          Inner

             !           v0 = -sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))/crmb
             !     >             * econv/mconv)

             presstarg(ir,2) = knds(idds(ir,2)) *(kteds(idds(ir,2))+ktids(idds(ir,2)))* (1.0 + mach0o**2) * econv
             v0 = -abs(kvds(idds(ir,1)))

             gtarg(ir,1)= knds(idds(ir,1)) * v0 * rmeani
             if (e2dtargopt.eq.5) then
                elecptarg(ir,1)= -e2dpepara(nks(ir)+1,ir)

                ionptarg(ir,1) = -e2dpipara(nks(ir)+1,ir)

             else
                elecptarg(ir,1)=gae*kteds(idds(ir,1))*econv*gtarg(ir,1)

                ionptarg(ir,1)=gaii*ktids(idds(ir,1))*econv*gtarg(ir,1)
             endif
             ! calculate target pressrure

             !          Totals

             presstarg(ir,1) = knds(idds(ir,1)) *(kteds(idds(ir,1))+ktids(idds(ir,1)))* (1.0 + mach0i**2) * econv
             gtarg(ir,3) = gtarg(ir,1) + gtarg(ir,2)
             elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
             ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)

             presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 

             write (6,'(a,i4,3g18.7)') 'DIVIMP  FLUXES:',ir,gtarg(ir,1),gtarg(ir,2),gtarg(ir,3)

             write (6,'(a,i4,3g18.7)') 'DIVIMP  ELEC P:',ir,elecptarg(ir,1),elecptarg(ir,2),elecptarg(ir,3)
             write (6,'(a,i4,3g18.7)') 'DIVIMP  ION  P:',ir,ionptarg(ir,1),ionptarg(ir,2),ionptarg(ir,3)

             write (6,'(a,i4,3g18.7)') 'DIVIMP  PRESS :',ir,presstarg(ir,1),presstarg(ir,2),presstarg(ir,3)

             !            Find ring for EDGE2D data

             if (fluxpts.gt.0.0) then

                !            Inner target

                in = ringno(ir,fluxinfo,fluxpts,maxins,4,ierr)

                !            Outer target

                e2dgtarg(ir,1) = -abs(fluxinfo(in,2))  /(dds(idds(ir,1))* 2.0 * PI * rp(idds(ir,1)))

                e2dgtarg(ir,2) = -abs(fluxinfo(in,3))  /(dds(idds(ir,2))* 2.0 * PI * rp(idds(ir,2)))

                e2dgtarg(ir,3) = e2dgtarg(ir,1)+e2dgtarg(ir,2)

                !            Modify for non-orthogonality and magnetic field

                if (cprint.eq.9)write (6,'(a,i4,3g18.7)') 'E2D RAW FLUXES:',ir,e2dgtarg(ir,1),e2dgtarg(ir,2),e2dgtarg(ir,3)

                e2dgtarg(ir,1)=e2dgtarg(ir,1)*kbfst(ir,1)/costet(idds(ir,1))

                e2dgtarg(ir,2)=e2dgtarg(ir,2)*kbfst(ir,2)/costet(idds(ir,2))

                !            Print outs

                e2dgtarg(ir,3) = e2dgtarg(ir,1)+e2dgtarg(ir,2)
                if (cprint.eq.9) then
                   write (6,'(a,i4,3g18.7)') 'E2D COR FLUXES:',ir,e2dgtarg(ir,1),e2dgtarg(ir,2),e2dgtarg(ir,3)
                   write (6,'(a,i4,3g18.7)') 'RATIOS        :',ir,e2dgtarg(ir,1)/gtarg(ir,1),e2dgtarg(ir,2)/gtarg(ir,2),&
                        e2dgtarg(ir,3)/gtarg(ir,3)

                   write (6,*)
                   write (6,'(a,i4,3g18.7)') 'E2D TARG FLUXES-E2D   :',ir,e2dtarg(ir,6,1),e2dtarg(ir,6,2)
                   write (6,'(a,i4,3g18.7)') 'E2D TARG FLUXES-OLDDIV:',ir,e2dtarg(ir,5,1),e2dtarg(ir,5,2)
                endif


             endif

          elseif (switch(swe2d).eq.1.0.or.switch(swe2d).eq.6.0) then
             if (switch(swmajr).eq.4.0) then
                rmeano = rs(1,ir)
                rmeani = rs(nks(ir),ir)
             else
                rmeano = 1.0
                rmeani = 1.0

             endif
             n1 = cellvals(ir,1,2)
             te0= kteds(idds(ir,2))
             te1= cellvals(ir,2,2)
             ti1= cellvals(ir,3,2)
             v1e2d = cellvals(ir,4,2)

             vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
             gtarg(ir,2)= n1 * vpe2d * rmeano
             elecptarg(ir,2)=gae*te1*econv*gtarg(ir,2)
             ionptarg(ir,2)=gaio*ti1*econv*gtarg(ir,2)

             presstarg(ir,2) = n1 * (te1+ti1) * econv+ n1 * vpe2d**2 *crmb * mconv
             n1 = cellvals(ir,1,1)
             te0= kteds(idds(ir,1))
             te1= cellvals(ir,2,1)
             ti1= cellvals(ir,3,1)
             v1e2d = -cellvals(ir,4,1)

             vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
             gtarg(ir,1)= n1 * vpe2d * rmeani
             elecptarg(ir,1)=gae*te1*econv*gtarg(ir,1)
             ionptarg(ir,1) =gaii*ti1*econv*gtarg(ir,1)

             presstarg(ir,1) = n1 * (te1+ti1) * econv+ n1 * vpe2d**2 *crmb * mconv
             gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
             elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
             ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)

             presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 

          elseif (switch(swe2d).eq.2.0.or.switch(swe2d).eq.4.0) then
             if (switch(swmajr).eq.4.0) then
                rmeano = rs(1,ir)
                rmeani = rs(nks(ir),ir)
             else
                rmeano = 1.0
                rmeani = 1.0

             endif
             n1 = cellvals(ir,1,2)
             te0= kteds(idds(ir,2))
             te1= cellvals(ir,2,2)
             ti1= cellvals(ir,3,2)
             v1e2d = cellvals(ir,4,2)

             vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
             gtarg(ir,2)= n1 * v1e2d * rmeano
             elecptarg(ir,2)=gae*te1*econv*gtarg(ir,2)
             ionptarg(ir,2)=gaio*ti1*econv*gtarg(ir,2)

             presstarg(ir,2) = n1 * (te1+ti1) * econv+ n1 * v1e2d**2 *crmb * mconv
             n1 = cellvals(ir,1,1)
             te0= kteds(idds(ir,1))
             te1= cellvals(ir,2,1)
             ti1= cellvals(ir,3,1)
             v1e2d = -cellvals(ir,4,1)

             vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
             gtarg(ir,1)= n1 * v1e2d * rmeani
             elecptarg(ir,1)=gae*te1*econv*gtarg(ir,1)
             ionptarg(ir,1) =gaii*ti1*econv*gtarg(ir,1)

             presstarg(ir,1) = n1 * (te1+ti1) * econv+ n1 * v1e2d**2 *crmb * mconv
             gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
             elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
             ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)

             presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 

          elseif (switch(swe2d).eq.3.0.or.switch(swe2d).eq.7.0) then
             if (switch(swmajr).eq.4.0) then
                rmeano = rs(1,ir)
                rmeani = rs(nks(ir),ir)
             else
                rmeano = 1.0
                rmeani = 1.0

             endif
             n1 = cellvals(ir,1,2)
             te0= kteds(idds(ir,2))
             te1= cellvals(ir,2,2)
             ti1= cellvals(ir,3,2)
             v1e2d = cellvals(ir,4,2)

             vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
             gtarg(ir,2)= -0.5 * (abs(e2dflux(1,ir))+ abs(e2dflux(2,ir))) * rmeano
             elecptarg(ir,2)=gae*te1*econv*gtarg(ir,2)
             ionptarg(ir,2)=gaio*ti1*econv*gtarg(ir,2)

             !            write (6,*) 'E2d1:',ir,nks(ir),e2dflux(1,ir),
             !     >                            e2dflux(2,ir),
             !     >                            gtarg(ir,2)


             presstarg(ir,2) = n1 * (te1+ti1) * econv+ n1 * v1e2d**2 *crmb * mconv
             n1 = cellvals(ir,1,1)
             te0= kteds(idds(ir,1))
             te1= cellvals(ir,2,1)
             ti1= cellvals(ir,3,1)
             v1e2d = -cellvals(ir,4,1)

             vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
             gtarg(ir,1)= -0.5 * (abs(e2dflux(nks(ir),ir))+ abs(e2dflux(nks(ir)+1,ir)))* rmeani
             elecptarg(ir,1)=gae*te1*econv*gtarg(ir,1)
             ionptarg(ir,1) =gaii*ti1*econv*gtarg(ir,1)

             !            write (6,*) 'E2d1:',ir,nks(ir),e2dflux(nks(ir),ir),
             !     >                            e2dflux(nks(ir)+1,ir),
             !     >                            gtarg(ir,1)



             presstarg(ir,1) = n1 * (te1+ti1) * econv+ n1 * v1e2d**2 *crmb * mconv
             gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
             elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
             ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)

             presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 

          elseif (switch(swe2d).eq.8.0) then
             if (switch(swmajr).eq.4.0) then
                rmeano = rs(1,ir)
                rmeani = rs(nks(ir),ir)
             else
                rmeano = 1.0
                rmeani = 1.0

             endif
             te1= cellvals(ir,2,2)

             ti1= cellvals(ir,3,2)

             gtarg(ir,2)= -abs((e2dgpara(1,ir)+e2dgpara(2,ir))/2.0)
             if (e2dtargopt.eq.5) then
                elecptarg(ir,2)= (e2dpepara(1,ir)+e2dpepara(2,ir))/2.0
                ionptarg(ir,2) = (e2dpipara(1,ir)+e2dpipara(2,ir))/2.0
             else
                elecptarg(ir,2)=gae*te1*econv*gtarg(ir,2)
                ionptarg(ir,2) =gaio*ti1*econv*gtarg(ir,2)

             endif
             te1= cellvals(ir,2,1)

             ti1= cellvals(ir,3,1)

             gtarg(ir,1)= -abs((e2dgpara(nks(ir)+1,ir)+e2dgpara(nks(ir),ir))/2.0)

             if (e2dtargopt.eq.5) then
                elecptarg(ir,2)= -(e2dpepara(nks(ir),ir)+e2dpepara(nks(ir)+1,ir))/2.0

                ionptarg(ir,2) = -(e2dpipara(nks(ir),ir)+e2dpipara(nks(ir)+1,ir))/2.0

             else
                elecptarg(ir,1)=gae*te1*econv*gtarg(ir,1)

                ionptarg(ir,1) =gaii*ti1*econv*gtarg(ir,1)

             endif
             gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
             elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
             ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)
             ! jdemod - turn off target pressure compensation for these options
             presstarg(ir,1) = 0.0
             presstarg(ir,2) = 0.0

             presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 

          elseif (switch(swe2d).eq.9.0) then
             if (switch(swmajr).eq.4.0) then
                rmeano = rs(1,ir)
                rmeani = rs(nks(ir),ir)
             else
                rmeano = 1.0
                rmeani = 1.0

             endif

             gtarg(ir,2)=      (e2dgpara(ike2d_start,ir)+e2dgpara(ike2d_start+1,ir))/2.0

             elecptarg(ir,2)=  (e2dpepara(ike2d_start,ir)+e2dpepara(ike2d_start+1,ir))/2.0

             ionptarg(ir,2) =  (e2dpipara(ike2d_start,ir)+e2dpipara(ike2d_start+1,ir))/2.0

             gtarg(ir,1)= -(e2dgpara(nks(ir)-ike2d_start+1,ir)+e2dgpara(nks(ir)-ike2d_start+1+1,ir))/2.0

             elecptarg(ir,1)=-(e2dpepara(nks(ir)-ike2d_start+1,ir)+e2dpepara(nks(ir)-ike2d_start+1+1,ir))/2.0

             ionptarg(ir,1) =-(e2dpipara(nks(ir)-ike2d_start+1,ir)+e2dpipara(nks(ir)-ike2d_start+1+1,ir))/2.0
             gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
             elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
             ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)
             ! jdemod - turn off target pressure compensation for these options
             presstarg(ir,1) = 0.0
             presstarg(ir,2) = 0.0

             presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 

          endif
          write(6,'(a,i8,10(1x,g12.5))') 'Calcfluxes:',ir,presstarg(ir,1),presstarg(ir,2),presstarg(ir,3),elecptarg(ir,1),&
               elecptarg(ir,2),elecptarg(ir,3),ionptarg(ir,1),ionptarg(ir,2),ionptarg(ir,3)
       end do
       !
       ! Accumulate data for main SOL and PFZ
       !
       call pr_trace('MOD_SOL22_INTERFACE:CALCFLUXES','BEFORE SOL/PFZ')
       do id = 1,2
          do ir = irsep,nrs
             if (ir.eq.irwall.or.ir.eq.irtrap) cycle
             if (ir.le.irwall) then 
                in=1
             else
                in=2
             endif
             solpfz_fluxes(in,id,1)=solpfz_fluxes(in,id,1)+ gtarg(ir,id) * dds(idds(ir,id))
             solpfz_fluxes(in,id,2) = solpfz_fluxes(in,id,2)+ elecptarg(ir,id)  * dds(idds(ir,id))
             solpfz_fluxes(in,id,3) = solpfz_fluxes(in,id,3)+ionptarg(ir,id)   * dds(idds(ir,id))
             solpfz_fluxes(in,id,4) = solpfz_fluxes(in,id,4)+presstarg(ir,id)  * dds(idds(ir,id))
          end do
          in = 1
          write(6,'(a,2i8,6(1x,g12.5))') 'SOL FLUX:',in,id,solpfz_fluxes(in,id,1),solpfz_fluxes(in,id,2),solpfz_fluxes(in,id,3),&
               solpfz_fluxes(in,id,4)
          in = 2
          write(6,'(a,2i8,6(1x,g12.5))') 'PFZ FLUX:',in,id,solpfz_fluxes(in,id,1),solpfz_fluxes(in,id,2),solpfz_fluxes(in,id,3),&
               solpfz_fluxes(in,id,4)
       end do
       !
       ! Calculate assignment of pfz fluxes to rings
       ! - use some decay profile ... linear or exponential? or something else
       !
       ! May not use them all but need to be calculated where the target
       ! geometry information is still available
       ! 
       ! sepdist2 contains the distance along the target to the target element 
       !

       ! find maximum pressure inner and outer
       maxpress=0.0
       maxpress_ir = 0
       call pr_trace('MOD_SOL22_INTERFACE:CALCFLUXES','BEFORE PRESSURE')
       do ir = irsep,nrs
          do id = 1,2
             if (presstarg(ir,id).gt.maxpress(id)) then 
                maxpress(id) = presstarg(ir,id)
                maxpress_ir(id) = ir
             endif
          end do
       end do
       dist_fact = 0
       ! flat distribution over distance pfz_dist_param
       if (pfz_dist_opt.eq.1) then 
          do ir = irsep,irwall
             do id = 1,2
                if (sepdist2(idds(ir,id)).le.pfz_dist_param(id)) then 
                   dist_fact (ir,id) = dds(idds(ir,id))
                   dist_fact(irwall+1,id) = dist_fact(irwall+1,id)+ dist_fact(ir,id)
                endif
             end do
          end do
          ! linear decay distribution over distance pfz_dist_param (largest at separatrix)
       elseif (pfz_dist_opt.eq.2) then

          do ir = irsep,irwall
             do id = 1,2
                if (sepdist2(idds(ir,id)).le.pfz_dist_param(id)) then 
                   dist_fact (ir,id) = dds(idds(ir,id)) *(pfz_dist_param(id)-sepdist2(idds(ir,id)))/pfz_dist_param(id)
                   dist_fact(irwall+1,id) = dist_fact(irwall+1,id)+ dist_fact(ir,id)
                endif
             end do
          end do
          ! exponential decay distribution with lambda = pfz_dist_param (largest at separatrix)
       elseif (pfz_dist_opt.eq.3) then
          do ir = irsep,irwall
             do id = 1,2
                dist_fact (ir,id) = dds(idds(ir,id))* exp(-sepdist2(idds(ir,id))/pfz_dist_param(id))
                dist_fact(irwall+1,id) = dist_fact(irwall+1,id)+ dist_fact(ir,id)
             end do
          end do

          ! linear decay to peak of pressure profile on target
       elseif (pfz_dist_opt.eq.4) then
          do ir = irsep,irwall
             do id = 1,2
                if (ir.le.maxpress_ir(id)) then 
                   dist_fact (ir,id) = dds(idds(ir,id)) *(sepdist2(idds(maxpress_ir(id),id))-&
                        sepdist2(idds(ir,id)))/sepdist2(idds(maxpress_ir(id),id))
                   dist_fact(irwall+1,id) = dist_fact(irwall+1,id)+ dist_fact(ir,id)
                endif
             end do
          end do
          ! exponential decay to peak of pressure profile on target
       elseif (pfz_dist_opt.eq.5) then
          do ir = irsep,irwall
             do id = 1,2
                if (ir.le.maxpress_ir(id)) then 
                   dist_fact (ir,id) = dds(idds(ir,id))* exp(-sepdist2(idds(ir,id))/sepdist2(idds(maxpress_ir(id),id)))
                   dist_fact(irwall+1,id) = dist_fact(irwall+1,id)+ dist_fact(ir,id)
                endif
             end do
          end do
       endif
       call pr_trace('MOD_SOL22_INTERFACE:CALCFLUXES','BEFORE NORMALIZE')
       ! normalize distribution 
       do ir = irsep,irwall
          do id = 1,2
             if (dist_fact(irwall+1,id).gt.0.0) then   
                dist_fact(ir,id) = dist_fact(ir,id)/dist_fact(irwall+1,id)
             else
                dist_fact(ir,id) = 0.0
             endif
          end do
       end do
       ! distribute pfz fluxes
       do ir = irsep,irwall
          do id = 1,2
             if (dds(idds(ir,id)).gt.0.0) then
                g_pfzsol(ir,id) = solpfz_fluxes(2,id,1)*dist_fact(ir,id)/dds(idds(ir,id))
                pe_pfzsol(ir,id) = solpfz_fluxes(2,id,2)*dist_fact(ir,id)/dds(idds(ir,id))
                pi_pfzsol(ir,id) = solpfz_fluxes(2,id,3)*dist_fact(ir,id)/dds(idds(ir,id))
                pr_pfzsol(ir,id) = solpfz_fluxes(2,id,3)*dist_fact(ir,id)/dds(idds(ir,id))
             else
                g_pfzsol(ir,id) = 0.0
                pe_pfzsol(ir,id) = 0.0
                pi_pfzsol(ir,id) = 0.0
                pr_pfzsol(ir,id) = 0.0
             endif
             write(6,'(a,2i8,6(1x,g12.5))') 'ADDITIONAL FLUX:',ir,id,g_pfzsol(ir,id),pe_pfzsol(ir,id),pi_pfzsol(ir,id),&
                  pr_pfzsol(ir,id)
          end do


       end do
       !      lambda = 0.01
       !      dist_opt = 1

       ! Calculate distribution factors
       !do ir = irsep,irwall
       !         elseif (dist_opt.eq.2) then   ! linear decay
       !  if (dist_opt.eq.1) then     ! flat
       !    elseif (dist_opt.eq.3) then ! exponential decay
       !       endif
       ! normalize the distribution
       !      dist_fact(ir,1) 

       call pr_trace('MOD_SOL22_INTERFACE:CALCFLUXES','END CALCFLUXES')
       return
     end subroutine calcfluxes

     
     subroutine calcsoliz(rconst,recfrac,gtarg,areasum,gperpa,oldknbs,grad_oldknbs,ioniz_src,rec_src,gperprat,ike2d_start)
       use debug_options
       use mod_params
       use mod_solparams
       use mod_solswitch
       use mod_printopt
       use mod_comtor
       use mod_cgeom
       use mod_pindata
       use mod_cedge2d
       !     include 'params'
       !     include 'solparams'
       !     include 'solswitch'
       !     include 'printopt'
       !      include 'solcommon'
       implicit none
       real    gperpa(maxnks,maxnrs),oldknbs(maxnks,maxnrs)
       real    grad_oldknbs(maxnks,maxnrs)
       real    ioniz_src(maxnks,maxnrs),rec_src(maxnks,maxnrs)
       real*8 rconst(mxspts),gperprat(mxspts)
       real*8 gtarg(mxspts,3),areasum(mxspts)
       real*8 gamcor,gamecor,recfrac


       !     include 'comtor'
       !     include 'cgeom'
       !     include 'pindata'
       !     include 'cedge2d'

       !     CALCSOLIZ:

       !     This routine calculates the total
       !     ionization on each ring in order
       !     to evaluate the net total excess or deficit of flux
       !     onto the ring. These values for each ring are then
       !     returned and used in the solver if gperp option 2 has
       !     been specified.



       !     Local variables

       integer ike2d_start
       integer ik,ir,in,id,startik,endik
       real rfact,v0,rmean1,rmean2
       real srcinteg (maxnks+3,maxnrs)
       real recinteg (maxnks+3,maxnrs)
       real srcsum,recsum,nesum,gradsum
       real gradsumneg,gradsumpos

       !     New variables for summaries

       real n1,te1,ti1,vpe2d,v1e2d,te0
       real intgrad,intgradpos,intgradneg,inte2dsrc,intdist
       real inte2dhrec,intgrad2,ds,dcell,mulfact
       integer ikin,irin,ikout,irout
       real    influx,outflux,netflux,intnflux,flux1,flux2,flux3
       real    flux4
       real fluxst,fluxend,ioniz,rec,gperpd,gperpn,gdiv,sider
       real intgperpd,intgnet,intgdiv
       real initflux,dp,dp1,dp2,brat1,brat2,startflux

       real flux_const,endflux

       logical debug

       !     Recylcling flux factor

       call pr_trace('MOD_SOL22_INTERFACE','START CALCSOLIZ')
       rfact = recfrac

       !     Calculate the net fluxes on each ring after summing over
       !     the ionization source.


       debug = .false.


       !        Now to sum up ionization and recombination (if on).

       !        For this specific case ONLY the integral will
       !        be done simply on a cell by cell basis - and
       !        so the R-value used can be the mean between the
       !        two end-points.


       do ir = irsep,nrs
          srcsum = 0.0
          recsum = 0.0
          areasum(ir) = 0.0
          nesum = 0.0
          gradsum = 0.0
          gradsumpos = 0.0

          !        Set IK loop limits based on E2D option

          gradsumneg = 0.0
          if (switch(swe2d).eq.0.0) then
             startik = 1
             endik   = nks(ir)
          elseif (switch(swe2d).gt.0.0) then
             startik = ike2d_start
             endik   = nks(ir) - ike2d_start + 1

             !         write (6,*) 'E2d switch:',switch(swe2d),ike2d_start,
             !     >                startik,endik


          endif

          do ik = startik,endik
             if (switch(swmajr).eq.4.0) then
                rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
             else
                rmean1 = 1.0
                rmean2 = 1.0

             endif

             if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
                  switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
                  switch(swe2d).eq.9))) then
                srcsum = srcsum+ ioniz_src(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))*rmean1
                recsum = recsum+ rec_src(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))*rmean1

                areasum(ir) = areasum(ir)+ (kss2(ik,ir)-ksb(ik-1,ir)) * rmean1
                nesum = nesum+ oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))*rmean1

                gradsum = gradsum+ grad_oldknbs(ik,ir)* (kss2(ik,ir)-ksb(ik-1,ir))*rmean1
                if (grad_oldknbs(ik,ir).lt.0.0) then
                   gradsumneg = gradsumneg+ grad_oldknbs(ik,ir)* (kss2(ik,ir)-ksb(ik-1,ir))*rmean1
                else
                   gradsumpos = gradsumpos+ grad_oldknbs(ik,ir)* (kss2(ik,ir)-ksb(ik-1,ir))*rmean1

                endif


                !           Do in two parts to get integral at cell centre.

             endif
             srcinteg(ik,ir) = srcsum

             recinteg(ik,ir) = recsum

             if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
                  switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
                  switch(swe2d).eq.9))) then
                srcsum = srcsum+ ioniz_src(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))*rmean2
                recsum = recsum+ rec_src(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))*rmean2

                areasum(ir) = areasum(ir)+ (ksb(ik,ir)-kss2(ik,ir)) * rmean2

                nesum = nesum+ oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))*rmean2

                gradsum = gradsum+ grad_oldknbs(ik,ir)* (ksb(ik,ir)-kss2(ik,ir))*rmean2
                if (grad_oldknbs(ik,ir).lt.0.0) then
                   gradsumneg = gradsumneg+ grad_oldknbs(ik,ir)* (ksb(ik,ir)-kss2(ik,ir))*rmean2
                else
                   gradsumpos = gradsumpos+ grad_oldknbs(ik,ir)* (ksb(ik,ir)-kss2(ik,ir))*rmean2

                endif

             endif

             !              Assign amount of outer ionization source

             if ( ( (ksmaxs2(ir)/2.0).ge.ksb(ik-1,ir) ).and.((ksmaxs2(ir)/2.0).lt.ksb(ik,ir))) then
                srcinteg(nks(ir)+2,ir) = srcsum
                recinteg(nks(ir)+2,ir) = recsum

             endif

          end do
          srcinteg(nks(ir)+1,ir) = srcsum

          recinteg(nks(ir)+1,ir) = recsum
          srcinteg(nks(ir)+3,ir) = srcsum - srcinteg(nks(ir)+2,ir)

          !        Calculate Rconst

          recinteg(nks(ir)+3,ir) = recsum - recinteg(nks(ir)+2,ir)
          if (switch(swrecom).eq.0.0) then
             if (areasum(ir).ne.0.0) then 
                rconst(ir)= - (gtarg(ir,3)+rfact*srcinteg(nks(ir)+1,ir))/areasum(ir)
             else
                rconst(ir)= 0.0
             endif
          elseif (switch(swrecom).eq.1.0.or.switch(swrecom).eq.2) then
             if (areasum(ir).ne.0.0) then 
                rconst(ir)= - (gtarg(ir,3)+rfact*srcinteg(nks(ir)+1,ir)- recinteg(nks(ir)+1,ir))/areasum(ir)
             else
                rconst(ir)= 0.0
             endif

             !        Calculate Gperp fraction ratios - for use in Gperp option 7.

          endif

          if (switch(swgperp).eq.7.0.or.switch(swgperpp).eq.7.0) then
             if (gradsum.eq.0.0) then
                gperprat(ir)= 100.0
             else
                if (gradsum.ne.0.0) then 
                   gperprat(ir)= max(abs(gradsumpos/gradsum),abs(gradsumneg/gradsum))
                else
                   gperprat(ir)= 0.0
                endif

             endif

             !         write (6,'(a,i4,1p,10(1x,g12.5))') 'Rconst:',ir,rconst(ir),
             !     >                        areasum(ir),
             !     >         gradsum,gradsumpos,gradsumneg,gradsumpos/gradsum,
             !     >                    gradsumneg/gradsum

             !        Calculate a distributed gperp source if that is
             !        required.


             !        Distributed with d2ne/dr2 component and constant.

          endif

          if ((switch(swgperp).eq.8.and.ir.le.irwall).or.(switch(swgperpp).eq.8.and.ir.ge.irtrap)) then

             !           Modify RCONST for GPERP term at given DPERP

             if (cprint.eq.9)write (6,'(a,i4,5g13.5)') 'GPERPA:',ir,nesum,gradsum,rconst(ir),areasum(ir),rconst(ir)*areasum(ir)
             if (areasum(ir).ne.0.0) then 
                rconst(ir) = - (-rconst(ir)*areasum(ir)+cdperp*gradsum)/areasum(ir)
             else
                rconst(ir) = 0.0

             endif
             intgrad = 0.0
             intgradpos = 0.0
             intgradneg = 0.0

             inte2dsrc = 0.0

             do ik = startik,endik
                if (switch(swmajr).eq.4.0) then
                   rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                   rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
                else
                   rmean1 = 1.0
                   rmean2 = 1.0

                endif

                gperpa(ik,ir) = cdperp * grad_oldknbs(ik,ir)+ rconst(ir)


                if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
                     switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
                     switch(swe2d).eq.9))) then

                   intgrad = intgrad+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
                   if (gperpa(ik,ir).le.0.0) then
                      intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
                   else
                      intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))

                   endif

                endif


                if (cprint.eq.9)write(6,'(2i4,7(1x,g13.5))') ir,ik,gperpa(ik,ir),oldknbs(ik,ir),cdperp*grad_oldknbs(ik,ir),&
                     intgrad,intgradpos,intgradneg


                if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
                     switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
                     switch(swe2d).eq.9))) then

                   intgrad = intgrad+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
                   if (gperpa(ik,ir).le.0.0) then
                      intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
                   else
                      intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))

                   endif
                endif
             end do

             if (cprint.eq.9)write (6,'(a,i4,7g13.5)') 'GPERPA-END:',ir,gradsum,rconst(ir),areasum(ir),&
                  rconst(ir)*areasum(ir),intgrad,intgradpos,intgradneg



             !        Distributed proportional to ne or d2ne/dr2

          endif

          if ( ((switch(swgperp).eq.3.0.or.switch(swgperp).eq.7.0).and.ir.le.irwall).or.&
               ((switch(swgperpp).eq.3.or.switch(swgperpp).eq.7.0).and.ir.gt.irwall)) then

             if (cprint.eq.9)write (6,'(a,i4,5g13.5)') 'GPERPA:',ir,nesum,gradsum,rconst(ir),areasum(ir),rconst(ir)*areasum(ir)
             intgrad = 0.0
             intgradpos = 0.0
             intgradneg = 0.0

             inte2dsrc = 0.0

             do ik = startik,endik
                if (switch(swmajr).eq.4.0) then
                   rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                   rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
                else
                   rmean1 = 1.0
                   rmean2 = 1.0

                endif

                if (switch(swgperp).eq.3.0.or.switch(swgperpp).eq.3.0)    then
                   if (nesum.ne.0.0) then 
                      gperpa(ik,ir)=oldknbs(ik,ir)* rconst(ir)*areasum(ir)/ nesum
                   else
                      gperpa(ik,ir)=0.0

                   endif

                elseif (switch(swgperp).eq.7.0.or.switch(swgperpp).eq.7.0)   then
                   if (gradsum.ne.0.0) then 
                      gperpa(ik,ir)=grad_oldknbs(ik,ir)* rconst(ir)*areasum(ir)/ gradsum
                   else
                      gperpa(ik,ir)=0.0

                   endif

                else

                   gperpa(ik,ir)=rconst(ir)

                endif


                if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.&
                     switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.&
                     switch(swe2d).eq.8))) then

                   intgrad = intgrad+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
                   if (gperpa(ik,ir).le.0.0) then
                      intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
                   else
                      intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))

                   endif

                endif


                if (cprint.eq.9)write(6,'(2i4,7(1x,g13.5))') ir,ik,gperpa(ik,ir),oldknbs(ik,ir),grad_oldknbs(ik,ir),&
                     intgrad,intgradpos,intgradneg


                if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
                     switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8))) then

                   intgrad = intgrad+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
                   if (gperpa(ik,ir).le.0.0) then
                      intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
                   else
                      intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))

                   endif
                endif
             end do
             if (cprint.eq.9)write (6,'(a,i4,7g13.5)') 'GPERPA-END:',ir,gradsum,rconst(ir),areasum(ir),&
                  rconst(ir)*areasum(ir),intgrad,intgradpos,intgradneg

             !     ------------------------------------------------------------------


             !        The following are a series of DIAGNOSTIC print outs related
             !        mostly to the EDGE2D sources that are read in from the GHOST file.

             !        They are only printed for a complete print out situation.


          endif

          if (cprint.eq.10) then
             call print_edge2d_flux_analysis(rconst,gtarg,areasum,ir,gperpa,oldknbs,grad_oldknbs,srcinteg,recinteg,ike2d_start)
             !     ------------------------------------------------------------------

             !        Close IR loop

          endif




          !     Print out the results


       end do
       write (6,*) 'Ring Ionization and Flux Integration:'
       write (6,*) 'Recycling factor (rfact) = ',rfact
       write (6,*)
       if (switch(swrecom).eq.0.0) then
          write (6,*) 'Recombination is OFF'
       elseif (switch(swrecom).eq.1.0.or.switch(swrecom).eq.2) then
          write (6,*) 'Recombination is ON'
       endif
       write (6,*)
       write (6,100) outer,inner,outer,inner,outer,inner
       write (6,*)
       do ir = irsep,nrs
          write (6,200) ir,gtarg(ir,3),srcinteg(nks(ir)+1,ir),rfact*srcinteg(nks(ir)+1,ir),rconst(ir),areasum(ir),&
               srcinteg(nks(ir)+2,ir),srcinteg(nks(ir)+3,ir),gtarg(ir,2),gtarg(ir,1),recinteg(nks(ir)+1,ir),&
               recinteg(nks(ir)+2,ir),recinteg(nks(ir)+3,ir)

          !     If the debug option is compiled ON then calculate
          !     all the potential fluxes - all of the integrals and
          !     their respective balances and corresponding Rconst
          !     values for BOTH the Edge2d ionization AND the
          !     DIVIMP calculated Pinion.

       end do
       if (debug) then
          call debugsoliz(rfact)
       endif

       !     Format statements

       call pr_trace('MOD_SOL22_INTERFACE:CALCSOLIZ','END CALCSOLIZ')
100    format(10x,'Total Flux',4x,'Total Ioniz',4x,'Rfact*Ioniz',5x,'Net Rconst',4x,'R(s) Integral',3x,a5,' Ioniz',3x,a5,&
            ' Ioniz',4x,a5,' Flux',4x,a5,' Flux',4x,'Total Recom',4x, a5,' Recom',4x, a5,' Recom' )



200    format(i4,2x,12(e14.6,1x))
       return

     end subroutine calcsoliz


     subroutine debugsoliz(rfact)
       use mod_params
       use mod_solparams
       use mod_solswitch
       use mod_comtor
       use mod_cgeom
       use mod_pindata
       use mod_cedge2d
       !     include 'params'
       !     include 'solparams'
       !     include 'solswitch'
       implicit none

       !     include 'comtor'
       !     include 'cgeom'
       !     include 'pindata'
       !     include 'cedge2d'

       !     This routine analyses the Flux/ionization balance
       !     for the BG plasma using various methods of calculating
       !     the target fluxes and ionization integrals including
       !     ds integration and major radius corrected integration.




       !     Local variables

       real rfact
       integer ik,ir,in,id
       real v0,rmean
       real srcinteg (maxnks+3,maxnrs)
       real srcsum

       real n1,te1,ti1,vpe2d,v1e2d,te0
       real*8 rconst(mxspts)
       real*8 gtarg(mxspts,3),areasum(mxspts)
       real gammas(maxnrs,8,2,3)
       real intsrc(maxnks+1,maxnrs,16)
       real srccomp(maxnks,maxnrs,24)
       real intarea(maxnks+1,maxnrs,16)

       !      real fluxes(maxnks+1,maxnrs,16)

       real newrconst(maxnrs,16)


       real srcsume,srcsumre,srcsumr,srchalfe,srchalfre,srchalf,srchalfr,rmean1,rmean2,ds1,ds2,ds,ds1h,ds2h,&
            dsh,areasuma,areasumr,areahalf,areahalfr,srcsuma,srchalfa,da,da1,da2,dah1,dah2,srcsumra,srchalfra,&
            asuma,asumra,ahalfa,ahalfra,dah,e2dsrc,fluxtot,fluxout,fluxin,fluxtotr,fluxinr,fluxoutr,e2dsrcr
       e2dsrc = 0.0

       fluxtot = 0.0
       do ir = 1,nrs
          if (ir.ge.irsep) then
             fluxout =-e2dflux(1,ir) * dds2(idds(ir,2))* costet(idds(ir,2)) / kbfst(ir,2)
             fluxin =e2dflux(nks(ir)+1,ir) * dds2(idds(ir,1))* costet(idds(ir,1)) / kbfst(ir,1)

             fluxtot = fluxtot + fluxout + fluxin
             fluxoutr = -e2dflux(1,ir)* rp(idds(ir,2)) * dds2(idds(ir,2))* costet(idds(ir,2)) / kbfst(ir,2)

             fluxinr = e2dflux(nks(ir)+1,ir)* rp(idds(ir,1)) * dds2(idds(ir,1))* costet(idds(ir,1)) / kbfst(ir,1)
             fluxtotr = fluxtotr + fluxoutr + fluxinr
             write(6,*) 'E2d Flux: ',ir,e2dflux(1,ir),e2dbvel(1,ir),e2dnbs(1,ir),e2dflux(nks(ir)+1,ir),&
                  e2dbvel(nks(ir)+1,ir),e2dnbs(nks(ir),ir),dds2(idds(ir,2)),dds2(idds(ir,1))
          endif
          do ik = 1,nks(ir)
             e2dsrc = e2dsrc + e2dion(ik,ir) * karea2(ik,ir)
             e2dsrcr = e2dsrcr + e2dion(ik,ir) * karea2(ik,ir)* rs(ik,ir)
          end do



       end do
       write (6,*) 'EDGE2D: Total Target Flux    = ',fluxtot
       write (6,*) 'EDGE2D: Total Ionization     = ',e2dsrc

       write (6,*) 'EDGE2D: Ratio Ioniz/Flux     = ',e2dsrc/fluxtot
       write (6,*) 'EDGE2D: Total Target Flux (R)= ',fluxtotr
       write (6,*) 'EDGE2D: Total Ionization  (R)= ',e2dsrcr

       write (6,*) 'EDGE2D: Ratio Ioniz/Flux  (R)= ',e2dsrcr/fluxtotr

       !        Calculate target fluxes

       !        Outer


       !        DIVIMP  G0

       do ir = irsep,nrs
          v0 = -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb* econv/mconv)
          gammas(ir,1,1,2)= knds(idds(ir,2)) * v0

          !          Inner

          gammas(ir,1,2,2)= knds(idds(ir,2)) * v0 *rp(idds(ir,2))
          v0 = -sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))/crmb* econv/mconv)
          gammas(ir,1,1,1)= knds(idds(ir,1)) * v0

          gammas(ir,1,2,1)= knds(idds(ir,1)) * v0 *rp(idds(ir,1))
          gammas(ir,1,1,3)= gammas(ir,1,1,1)+gammas(ir,1,1,2)

          !        EDGE2D G0

          !          Outer

          gammas(ir,1,2,3)= gammas(ir,1,2,1)+gammas(ir,1,2,2)
          gammas(ir,2,1,2)= e2dflux(1,ir)

          !          Inner

          gammas(ir,2,2,2)= e2dflux(1,ir) * rp(idds(ir,2))
          gammas(ir,2,1,1)= -e2dflux(nks(ir)+1,ir)

          gammas(ir,2,2,1)= -e2dflux(nks(ir)+1,ir) * rp(idds(ir,1))
          gammas(ir,2,1,3)= gammas(ir,2,1,1)+gammas(ir,2,1,2)

          !        EDGE2D G1


          !          Outer

          gammas(ir,2,2,3)= gammas(ir,2,2,1)+gammas(ir,2,2,2)
          gammas(ir,3,1,2)= e2dflux(2,ir)

          !          Inner

          gammas(ir,3,2,2)= e2dflux(2,ir) * krb(1,ir)
          gammas(ir,3,1,1)= -e2dflux(nks(ir),ir)

          gammas(ir,3,2,1)= -e2dflux(nks(ir),ir) * krb(nks(ir)-1,ir)
          gammas(ir,3,1,3)= gammas(ir,3,1,1)+gammas(ir,3,1,2)

          !        EGDE2D G1/2 - A

          !          Outer

          gammas(ir,3,2,3)= gammas(ir,3,2,1)+gammas(ir,3,2,2)
          gammas(ir,4,1,2)= 0.5*(gammas(ir,2,1,2)+gammas(ir,3,1,2))

          !          Inner

          gammas(ir,4,2,2)= 0.5*(gammas(ir,2,2,2)+gammas(ir,3,2,2))
          gammas(ir,4,1,1)= 0.5*(gammas(ir,2,1,1)+gammas(ir,3,1,1))

          gammas(ir,4,2,1)= 0.5*(gammas(ir,2,2,1)+gammas(ir,3,2,1))
          gammas(ir,4,1,3)= gammas(ir,4,1,1)+gammas(ir,4,1,2)

          !        EDGE2D G1/2 - B

          gammas(ir,4,2,3)= gammas(ir,4,2,1)+gammas(ir,4,2,2)
          n1 = cellvals(ir,1,2)
          te0= kteds(idds(ir,2))
          te1= cellvals(ir,2,2)
          ti1= cellvals(ir,3,2)

          vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
          gammas(ir,5,1,2)= n1 * vpe2d

          gammas(ir,5,2,2)= n1 * vpe2d * rs(1,ir)
          n1 = cellvals(ir,1,1)
          te0= kteds(idds(ir,1))
          te1= cellvals(ir,2,1)
          ti1= cellvals(ir,3,1)

          vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
          gammas(ir,5,1,1)= n1 * vpe2d

          gammas(ir,5,2,1)= n1 * vpe2d * rs(nks(ir),ir)
          gammas(ir,5,1,3)= gammas(ir,5,1,1)+gammas(ir,5,1,2)


          !        EDGE2D G1/2 - C

          gammas(ir,5,2,3)= gammas(ir,5,2,1)+gammas(ir,5,2,2)
          n1 = cellvals(ir,1,2)

          v1e2d = cellvals(ir,4,2)
          gammas(ir,6,1,2)= n1 * v1e2d

          gammas(ir,6,2,2)= n1 * v1e2d * rs(1,ir)
          n1 = cellvals(ir,1,1)

          v1e2d = -cellvals(ir,4,1)
          gammas(ir,6,1,1)= n1 * v1e2d

          gammas(ir,6,2,1)= n1 * v1e2d * rs(nks(ir),ir)
          gammas(ir,6,1,3)= gammas(ir,6,1,1)+gammas(ir,6,1,2)



          !        EDGE2D G0 - KAREA

          !          Outer

          gammas(ir,6,2,3)= gammas(ir,6,2,1)+gammas(ir,6,2,2)
          id = idds(ir,2)
          gammas(ir,7,1,2)= e2dflux(1,ir) * dds2(id) * costet(id)/ kbfst(ir,2)

          !          Inner

          gammas(ir,7,2,2)= e2dflux(1,ir) * rp(idds(ir,2))* dds2(id) * costet(id) / kbfst(ir,2)
          id = idds(ir,1)
          gammas(ir,7,1,1)= -e2dflux(nks(ir)+1,ir)* dds2(id) * costet(id)/ kbfst(ir,1)

          gammas(ir,7,2,1)= -e2dflux(nks(ir)+1,ir) * rp(idds(ir,1))* dds2(id) * costet(id) / kbfst(ir,1)
          gammas(ir,7,1,3)= gammas(ir,7,1,1)+gammas(ir,7,1,2)

          !        EGDE2D G1/2 - A - KAREA

          !          Outer

          !          Outer

          gammas(ir,7,2,3)= gammas(ir,7,2,1)+gammas(ir,7,2,2)
          id = idds(ir,2)
          gammas(ir,8,1,2)= 0.5*(e2dflux(1,ir)+e2dflux(2,ir))* dds2(id) * costet(id)/ kbfs(1,ir)

          !          Inner

          gammas(ir,8,2,2)= 0.5*(e2dflux(1,ir)+e2dflux(2,ir))* rs(1,ir)* dds2(id) * costet(id) / kbfs(1,ir)
          id = idds(ir,1)
          gammas(ir,8,1,1)= -0.5*(e2dflux(nks(ir)+1,ir)+e2dflux(nks(ir),ir))* dds2(id) * costet(id)/ kbfs(nks(ir),ir)

          gammas(ir,8,2,1)= -0.5*(e2dflux(nks(ir)+1,ir)+e2dflux(nks(ir),ir))* rs(nks(ir),ir)* dds2(id) * &
               costet(id) / kbfs(nks(ir),ir)
          gammas(ir,8,1,3)= gammas(ir,8,1,1)+gammas(ir,8,1,2)

          !          That should be all of the gammas


          !          Calculate the ionization source terms and all of
          !          its components.


          !        Now to sum up ionization.

          !        For this specific case ONLY the integral will
          !        be done simply on a cell by cell basis - and
          !        so the R-value used can be the mean between the
          !        two end-points.


          gammas(ir,8,2,3)= gammas(ir,8,2,1)+gammas(ir,8,2,2)
          srcsume  = 0.0
          srcsumre = 0.0
          srcsum  = 0.0
          srcsumr = 0.0
          srchalfe = 0.0
          srchalfre = 0.0
          srchalf = 0.0

          srchalfr = 0.0
          srcsuma = 0.0
          srchalfa= 0.0
          srcsumra = 0.0

          srchalfra= 0.0
          areasuma = 0.0
          areahalf = 0.0
          areasumr = 0.0

          areahalfr= 0.0
          asuma = 0.0
          asumra = 0.0
          ahalfa = 0.0


          ahalfra = 0.0

          do ik = 1,nks(ir)

             rmean = (krb(ik-1,ir) + krb(ik,ir)) /2.0

             rmean1 = (krb(ik-1,ir) + rs(ik,ir)) /2.0

             rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
             ds1 = kss2(ik,ir) - ksb(ik-1,ir)
             ds2 = ksb(ik,ir)  - kss2(ik,ir)
             ds = ds1 + ds2
             ds1h = ds1
             ds2h = ds2
             da = karea2(ik,ir)
             da1 = da /2.0
             da2 = da /2.0
             dah1 = da / 2.0
             dah2 = da / 2.0
             if (ik.eq.1) then
                ds1h = 0.0
                dah1 = 0.0
             endif
             if (ik.eq.nks(ir)) then
                ds2h = 0.0
                dah2 = 0.0
             endif
             dsh = ds1h + ds2h
             dah = dah1 + dah2
             srcsumr = srcsumr + pinion(ik,ir)*ds1*rmean1
             srcsumre = srcsumre + e2dion(ik,ir)*ds1*rmean1
             srcsum = srcsum  + pinion(ik,ir)*ds1

             srcsume = srcsume + e2dion(ik,ir)*ds1
             srcsuma = srcsuma + e2dion(ik,ir) * da1
             srcsumra = srcsumra + e2dion(ik,ir) * da1 * rmean1
             srchalfa = srchalfa + e2dion(ik,ir) * dah1

             srchalfra = srchalfra + e2dion(ik,ir) * dah1 * rmean1
             srchalfr = srchalfr + pinion(ik,ir)*ds1h*rmean1
             srchalfre = srchalfre + e2dion(ik,ir)*ds1h*rmean1
             srchalf = srchalf + pinion(ik,ir)*ds1h

             srchalfe = srchalfe + e2dion(ik,ir)*ds1h
             areasuma = areasuma + ds1
             areahalf = areahalf + ds1h
             areasumr = areasumr + rmean1 * ds1

             areahalfr = areahalfr + rmean1 * ds1h
             asuma = asuma + da1
             asumra = asumra + da1 * rmean1
             ahalfa = ahalfa + dah1


             !           Do in two parts to get integral at cell centre.

             ahalfra = ahalfra + dah1 * rmean1
             intsrc(ik,ir,1) = srcsum
             intsrc(ik,ir,2) = srcsumr
             intsrc(ik,ir,3) = srcsume
             intsrc(ik,ir,4) = srcsumre
             intsrc(ik,ir,5) = srchalf
             intsrc(ik,ir,6) = srchalfr
             intsrc(ik,ir,7) = srchalfe

             intsrc(ik,ir,8) = srchalfre
             intsrc(ik,ir,9) =  srcsuma
             intsrc(ik,ir,10) = srchalfa
             intsrc(ik,ir,11) = srcsumra

             intsrc(ik,ir,12) = srchalfra
             intarea(ik,ir,1) = areasuma
             intarea(ik,ir,2) = areahalf
             intarea(ik,ir,3) = areasumr

             intarea(ik,ir,4) = areahalfr
             intarea(ik,ir,5) = asuma
             intarea(ik,ir,6) = ahalfa
             intarea(ik,ir,7) = asumra


             intarea(ik,ir,8) = ahalfra
             srcsumr = srcsumr + pinion(ik,ir)*ds2*rmean2
             srcsumre = srcsumre + e2dion(ik,ir)*ds2*rmean2
             srcsum = srcsum  + pinion(ik,ir)*ds2

             srcsume = srcsume + e2dion(ik,ir)*ds2
             srcsuma = srcsuma + e2dion(ik,ir) * da2
             srcsumra = srcsumra + e2dion(ik,ir) * da2 * rmean2
             srchalfa = srchalfa + e2dion(ik,ir) * dah2

             srchalfra = srchalfra + e2dion(ik,ir) * dah2 * rmean2
             srchalfr = srchalfr + pinion(ik,ir)*ds2h*rmean2
             srchalfre = srchalfre + e2dion(ik,ir)*ds2h*rmean2
             srchalf = srchalf + pinion(ik,ir)*ds2h

             srchalfe = srchalfe + e2dion(ik,ir)*ds2h
             areasuma = areasuma + ds2
             areahalf = areahalf + ds2h
             areasumr = areasumr + rmean2 * ds2

             areahalfr = areahalfr + rmean2 * ds2h
             asuma = asuma + da2
             asumra = asumra + da2 * rmean2
             ahalfa = ahalfa + dah2

             !           Calculate Areas and other components.

             ahalfra = ahalfra + dah2 * rmean2
             srccomp(ik,ir,1) = pinion(ik,ir)
             srccomp(ik,ir,2) = e2dion(ik,ir)
             srccomp(ik,ir,3) = pinion(ik,ir) * ds
             srccomp(ik,ir,4) = e2dion(ik,ir) * ds
             srccomp(ik,ir,5) = pinion(ik,ir) * dsh
             srccomp(ik,ir,6) = e2dion(ik,ir) * dsh
             srccomp(ik,ir,7) = pinion(ik,ir) * rmean
             srccomp(ik,ir,8) = e2dion(ik,ir) * rmean
             srccomp(ik,ir,9) = pinion(ik,ir) * rmean * ds
             srccomp(ik,ir,10) = e2dion(ik,ir) * rmean * ds
             srccomp(ik,ir,11) = pinion(ik,ir) * rmean * dsh

             srccomp(ik,ir,12) = e2dion(ik,ir) * rmean * dsh
             srccomp(ik,ir,17) = e2dion(ik,ir) * da
             srccomp(ik,ir,18) = e2dion(ik,ir) * dah
             srccomp(ik,ir,19) = e2dion(ik,ir) * rmean * da

             srccomp(ik,ir,20) = e2dion(ik,ir) * rmean * dah
             srccomp(ik,ir,13) = rmean * ds
             srccomp(ik,ir,14) = rmean * dsh
             srccomp(ik,ir,15) = rmean * da

             srccomp(ik,ir,16) = rmean * dah
             srccomp(ik,ir,21) = ds



             srccomp(ik,ir,22) = dsh
             intsrc(ik,ir,13) = srcsuma
             intsrc(ik,ir,14) = srchalfa
             intsrc(ik,ir,15) = srcsumra

             intsrc(ik,ir,16) = srchalfra
             intarea(ik,ir,9) = asuma
             intarea(ik,ir,10) = ahalfa
             intarea(ik,ir,11) = asumra

             intarea(ik,ir,12) = ahalfra

          end do
          intsrc(nks(ir)+1,ir,1) = srcsum
          intsrc(nks(ir)+1,ir,2) = srcsumr
          intsrc(nks(ir)+1,ir,3) = srcsume
          intsrc(nks(ir)+1,ir,4) = srcsumre
          intsrc(nks(ir)+1,ir,5) = srchalf
          intsrc(nks(ir)+1,ir,6) = srchalfr
          intsrc(nks(ir)+1,ir,7) = srchalfe

          intsrc(nks(ir)+1,ir,8) = srchalfre
          intsrc(nks(ir)+1,ir,9) =  srcsuma
          intsrc(nks(ir)+1,ir,10) = srchalfa
          intsrc(nks(ir)+1,ir,11) = srcsumra

          intsrc(nks(ir)+1,ir,12) = srchalfra
          intsrc(nks(ir)+1,ir,13) = srcsuma
          intsrc(nks(ir)+1,ir,14) = srchalfa
          intsrc(nks(ir)+1,ir,15) = srcsumra

          intsrc(nks(ir)+1,ir,16) = srchalfra
          intarea(nks(ir)+1,ir,1) = areasuma
          intarea(nks(ir)+1,ir,2) = areahalf
          intarea(nks(ir)+1,ir,3) = areasumr

          intarea(nks(ir)+1,ir,4) = areahalfr
          intarea(nks(ir)+1,ir,5) = asuma
          intarea(nks(ir)+1,ir,6) = ahalfa
          intarea(nks(ir)+1,ir,7) = asumra

          intarea(nks(ir)+1,ir,8) = ahalfra
          intarea(nks(ir)+1,ir,9) = asuma
          intarea(nks(ir)+1,ir,10) = ahalfa
          intarea(nks(ir)+1,ir,11) = asumra

          intarea(nks(ir)+1,ir,12) = ahalfra
          newrconst(ir,1) = - (gammas(ir,1,1,3)+ rfact*intsrc(nks(ir)+1,ir,1))/ intarea(nks(ir)+1,ir,1)

          newrconst(ir,2) = - (gammas(ir,1,2,3)+ rfact*intsrc(nks(ir)+1,ir,2))/ intarea(nks(ir)+1,ir,3)
          newrconst(ir,3) = - (gammas(ir,2,1,3)+ intsrc(nks(ir)+1,ir,3))/ intarea(nks(ir)+1,ir,1)

          newrconst(ir,4) = - (gammas(ir,2,2,3)+ intsrc(nks(ir)+1,ir,4))/ intarea(nks(ir)+1,ir,3)
          newrconst(ir,5) = - (gammas(ir,4,1,3)+ intsrc(nks(ir)+1,ir,7))/ intarea(nks(ir)+1,ir,2)

          newrconst(ir,6) = - (gammas(ir,4,2,3)+ intsrc(nks(ir)+1,ir,8))/ intarea(nks(ir)+1,ir,4)
          newrconst(ir,7) = - (gammas(ir,4,1,3)+ rfact*intsrc(nks(ir)+1,ir,5))/ intarea(nks(ir)+1,ir,2)

          newrconst(ir,8) = - (gammas(ir,4,2,3)+ rfact*intsrc(nks(ir)+1,ir,6))/ intarea(nks(ir)+1,ir,4)
          newrconst(ir,9) = - (gammas(ir,5,1,3)+ rfact*intsrc(nks(ir)+1,ir,7))/ intarea(nks(ir)+1,ir,2)

          newrconst(ir,10) = - (gammas(ir,5,2,3)+ rfact*intsrc(nks(ir)+1,ir,8))/ intarea(nks(ir)+1,ir,4)
          newrconst(ir,11) = - (gammas(ir,6,1,3)+ rfact*intsrc(nks(ir)+1,ir,7))/ intarea(nks(ir)+1,ir,2)

          newrconst(ir,12) = - (gammas(ir,6,2,3)+ rfact*intsrc(nks(ir)+1,ir,8))/ intarea(nks(ir)+1,ir,4)
          newrconst(ir,13) = - (gammas(ir,7,1,3)+ intsrc(nks(ir)+1,ir,9))/ intarea(nks(ir)+1,ir,5)

          newrconst(ir,14) = - (gammas(ir,7,2,3)+ intsrc(nks(ir)+1,ir,11))/ intarea(nks(ir)+1,ir,7)
          newrconst(ir,15) = - (gammas(ir,8,1,3)+ intsrc(nks(ir)+1,ir,10))/ intarea(nks(ir)+1,ir,6)

          !        Calculate the Fluxes at each grid point

          newrconst(ir,16) = - (gammas(ir,8,2,3)+ intsrc(nks(ir)+1,ir,12))/ intarea(nks(ir)+1,ir,8)
          do ik = 1,nks(ir)
             fluxes(ik,ir,1) = (gammas(ir,7,1,2) + intsrc(ik,ir,9)+  newrconst(ir,13) * intarea(ik,ir,5))* &
                  srccomp(ik,ir,21)/karea2(ik,ir)
             fluxes(ik,ir,2) = ((gammas(ir,7,2,2) + intsrc(ik,ir,11)+  newrconst(ir,14) * intarea(ik,ir,7))*&
                  srccomp(ik,ir,21)/karea2(ik,ir))/  rs(ik,ir)
             fluxes(ik,ir,3) = (gammas(ir,8,1,2) + intsrc(ik,ir,10)+  newrconst(ir,15) * intarea(ik,ir,6))*  &
                  srccomp(ik,ir,22)/karea2(ik,ir)/rs(ik,ir)

             !           Fluxes for regular Edge2d full and 1/2 cell cases
             !           not using kareas - normal and R-corrected. Cell centred.

             fluxes(ik,ir,4) = ((gammas(ir,8,2,2) + intsrc(ik,ir,12)+  newrconst(ir,16) * intarea(ik,ir,8))* &
                  srccomp(ik,ir,22)/karea2(ik,ir)/rs(ik,ir))/  rs(ik,ir)
             fluxes(ik,ir,9) = (gammas(ir,2,1,2) + intsrc(ik,ir,3)+  newrconst(ir,3) * intarea(ik,ir,1))
             fluxes(ik,ir,10) = (gammas(ir,2,2,2) + intsrc(ik,ir,4)+  newrconst(ir,4) * intarea(ik,ir,3))/  rs(ik,ir)
             fluxes(ik,ir,11) = (gammas(ir,4,1,2) + intsrc(ik,ir,7)+  newrconst(ir,5) * intarea(ik,ir,2))

             !           Fluxes for regular DIVIMP/PIN full and 1/2 cell cases
             !           not using kareas - normal and R-corrected. Cell centred.

             fluxes(ik,ir,12) = (gammas(ir,4,2,2) + intsrc(ik,ir,8)+  newrconst(ir,6) * intarea(ik,ir,4))/  rs(ik,ir)
             fluxes(ik,ir,13) = (gammas(ir,1,1,2) + intsrc(ik,ir,1)+  newrconst(ir,1) * intarea(ik,ir,1))
             fluxes(ik,ir,14) = (gammas(ir,1,2,2) + intsrc(ik,ir,2)+  newrconst(ir,2) * intarea(ik,ir,3))/  rs(ik,ir)
             fluxes(ik,ir,15) = (gammas(ir,4,1,2) + intsrc(ik,ir,5)+  newrconst(ir,7) * intarea(ik,ir,2))
             fluxes(ik,ir,16) = (gammas(ir,4,2,2) + intsrc(ik,ir,6)+  newrconst(ir,8) * intarea(ik,ir,4))/  rs(ik,ir)

          end do
          do ik = 1,nks(ir)+1
             if (ik.eq.1) then
                fluxes(ik,ir,5) = gammas(ir,7,1,2)*srccomp(1,ir,21)/karea2(1,ir)
                fluxes(ik,ir,6) = gammas(ir,7,2,2)*srccomp(1,ir,21)/karea2(1,ir)/ rp(idds(ir,2))
                fluxes(ik,ir,7) = gammas(ir,8,1,2)
                fluxes(ik,ir,8) = gammas(ir,8,2,2) / rs(ik,ir)
             else
                fluxes(ik,ir,5) = (gammas(ir,7,1,2) + intsrc(ik-1,ir,13)+  newrconst(ir,13) * intarea(ik-1,ir,9))* &
                     srccomp(ik-1,ir,21)/karea2(ik-1,ir)
                fluxes(ik,ir,6) = (gammas(ir,7,2,2) + intsrc(ik-1,ir,15)+  newrconst(ir,14) * intarea(ik-1,ir,11))* &
                     srccomp(ik-1,ir,21)/karea2(ik-1,ir)/rs(ik-1,ir)
                fluxes(ik,ir,7) = (gammas(ir,8,1,2) + intsrc(ik-1,ir,14)+  newrconst(ir,15) * intarea(ik-1,ir,10))* &
                     srccomp(ik-1,ir,22)/karea2(ik-1,ir)/rs(ik-1,ir)
                fluxes(ik,ir,8) = (gammas(ir,8,2,2) + intsrc(ik-1,ir,16)+  newrconst(ir,16) * intarea(ik-1,ir,12))* &
                     srccomp(ik-1,ir,22)/karea2(ik-1,ir)/rs(ik-1,ir)

             endif

             !     And NOW to somehow print ALL of this in a usable fashion!!!


             !     Print only for the separatrix ring for now.

          end do

          if (ir.eq.irsep.or.ir.eq.irsep+1.or.ir.eq.irsep+2.or.ir.eq.irsep+5.or.ir.eq.irsep+7.or.ir.eq.irsep+11) then
             write (6,*) 'ANALYSIS for ring:', ir
             write (6,*) ' The Recycling Fraction is applied to'//' the ionization source calculated below - THEN '
             write (6,*) ' Used in the formulae for RCONST'

             !        Target Fluxes

             write (6,*) ' RECYCLING FRACTION = ',rfact
             write (6,*) ' ALL Fluxes and Integrals are per meter'//' toroidally and not for the whole torus.'


             write (6,*) ' The units of the target flux is: ions/m2/s'
             write (6,*) 'Target Fluxes - calculated various ways:'

             write (6,*) 'Flux name:           Inner           Outer'//'         Total         R-Inner         R-outer'//&
                  '         R-Total'

             write (6,300) 'DIVIMP 0',gammas(ir,1,1,1),gammas(ir,1,1,2),gammas(ir,1,1,3),gammas(ir,1,2,1),gammas(ir,1,2,2),&
                  gammas(ir,1,2,3)

             write (6,300) 'EDGE2D 0',gammas(ir,2,1,1),gammas(ir,2,1,2),gammas(ir,2,1,3),gammas(ir,2,2,1),gammas(ir,2,2,2),&
                  gammas(ir,2,2,3)

             write (6,300) 'EDGE2D G1/2',gammas(ir,4,1,1),gammas(ir,4,1,2),gammas(ir,4,1,3),gammas(ir,4,2,1),gammas(ir,4,2,2),&
                  gammas(ir,4,2,3)

             write (6,300) 'EDGE2D VP',gammas(ir,5,1,1),gammas(ir,5,1,2),gammas(ir,5,1,3),gammas(ir,5,2,1),gammas(ir,5,2,2),&
                  gammas(ir,5,2,3)

             write (6,300) 'EDGE2D VCAV',gammas(ir,6,1,1),gammas(ir,6,1,2),gammas(ir,6,1,3),gammas(ir,6,2,1),gammas(ir,6,2,2),&
                  gammas(ir,6,2,3)

             write (6,300) 'EDGE2D POL 0',gammas(ir,7,1,1),gammas(ir,7,1,2),gammas(ir,7,1,3),gammas(ir,7,2,1),gammas(ir,7,2,2),&
                  gammas(ir,7,2,3)

             !        Total Ionization values and integrated areas

             write (6,300) 'EDGE2D POL.5',gammas(ir,8,1,1),gammas(ir,8,1,2),gammas(ir,8,1,3),gammas(ir,8,2,1),gammas(ir,8,2,2),&
                  gammas(ir,8,2,3)
             write (6,*) 'Ionization and AREA Totals for each type:'

             write (6,*) 'Name:            Total Src       Total R-Src'//'      Total Area      Total R-Area'

             write (6,300) 'DIVIMP FULL',intsrc(nks(ir)+1,ir,1),intsrc(nks(ir)+1,ir,2),intarea(nks(ir)+1,ir,1),&
                  intarea(nks(ir)+1,ir,3)

             write (6,300) 'EDGE2D FULL',intsrc(nks(ir)+1,ir,3),intsrc(nks(ir)+1,ir,4),intarea(nks(ir)+1,ir,1),&
                  intarea(nks(ir)+1,ir,3)

             write (6,300) 'DIVIMP 1/2',intsrc(nks(ir)+1,ir,5),intsrc(nks(ir)+1,ir,6),intarea(nks(ir)+1,ir,2),&
                  intarea(nks(ir)+1,ir,4)

             write (6,300) 'EDGE2D 1/2',intsrc(nks(ir)+1,ir,7),intsrc(nks(ir)+1,ir,8),intarea(nks(ir)+1,ir,2),&
                  intarea(nks(ir)+1,ir,4)

             write (6,300) 'EDGE2D POL F',intsrc(nks(ir)+1,ir,9),intsrc(nks(ir)+1,ir,11),intarea(nks(ir)+1,ir,5),&
                  intarea(nks(ir)+1,ir,7)


             !        NET RCONST values

             write (6,300) 'EDGE2D POL.5',intsrc(nks(ir)+1,ir,10),intsrc(nks(ir)+1,ir,12),intarea(nks(ir)+1,ir,6),&
                  intarea(nks(ir)+1,ir,8)
             write (6,*)  ' RCONST values:'
             write (6,*)  ' DIVIMP ionization sources * recycling factor'

             write (6,*)  'Name:           Rconst          R-Rconst'
             write (6,300) 'DIVIMP FULL',newrconst(ir,1),newrconst(ir,2)
             write (6,300) 'EDGE2D FULL',newrconst(ir,3),newrconst(ir,4)
             write (6,300) 'EDGE2D G1/2',newrconst(ir,5),newrconst(ir,6)
             write (6,300) 'DIVIMP G1/2',newrconst(ir,7),newrconst(ir,8)
             write (6,300) 'EDGE2D VP',newrconst(ir,9),newrconst(ir,10)
             write (6,300) 'EDGE2D VCAV',newrconst(ir,11),newrconst(ir,12)
             write (6,300) 'EDGE2D POL F',newrconst(ir,13),newrconst(ir,14)


             !        Now for the Ionization Source details!

             write (6,300) 'EDGE2D POL.5',newrconst(ir,15),newrconst(ir,16)
             write (6,*)

             write (6,*) ' And NOW for the Ionization Source Details:'
             write (6,*) ' EDGE2D DATA - FULL RING -'//' REGULAR and R-corrected - POL - KAREA'
             write (6,*) ' Integrals are to cell centres:'
             write (6,700)
             do ik = 1,nks(ir)
                write(6,600) ik,karea2(ik,ir),srccomp(ik,ir,21),intsrc(ik,ir,9),intsrc(ik,ir,11),e2dion(ik,ir),&
                     srccomp(ik,ir,8),srccomp(ik,ir,15),srccomp(ik,ir,19),srccomp(ik,ir,17),intarea(ik,ir,5),&
                     intarea(ik,ir,7),fluxes(ik,ir,1),fluxes(ik,ir,2)

             end do
             write (6,*) ' EDGE2D DATA - RING - 1/2 CELLS'//' -  REGULAR and R-corrected - POL - KAREA'
             write (6,*) ' Integrals are to cell centres:'
             write (6,700)
             do ik = 1,nks(ir)
                write(6,600) ik,karea2(ik,ir),srccomp(ik,ir,22),intsrc(ik,ir,10),intsrc(ik,ir,12),e2dion(ik,ir),&
                     srccomp(ik,ir,8),srccomp(ik,ir,16),srccomp(ik,ir,20),srccomp(ik,ir,18),intarea(ik,ir,6),&
                     intarea(ik,ir,8),fluxes(ik,ir,3),fluxes(ik,ir,4)

             end do
             write (6,*) ' EDGE2D DATA - FULL RING -'//' REGULAR and R-corrected - POL - KAREA -'//' AT CELL BOUNDARIES'
             write (6,*) ' Integrals are to cell boundaries:'
             write (6,750)
             do ik = 1,nks(ir)+1
                write(6,650) ik,karea2(ik,ir),srccomp(ik,ir,21),intsrc(ik,ir,13),intsrc(ik,ir,15),e2dion(ik,ir),&
                     srccomp(ik,ir,8),srccomp(ik,ir,15),srccomp(ik,ir,19),srccomp(ik,ir,17),intarea(ik,ir,9),&
                     intarea(ik,ir,11),fluxes(ik,ir,5),fluxes(ik,ir,6),e2dflux(ik,ir),ksb(ik,ir)
             end do
             write (6,*) ' EDGE2D DATA - RING - 1/2 CELLS'//' -  REGULAR and R-corrected - POL - KAREA'//' AT CELL BOUNDARIES'
             write (6,*) ' Integrals are to cell boundaries:'
             write (6,750)
             do ik = 1,nks(ir)+1
                write(6,650) ik,karea2(ik,ir),srccomp(ik,ir,22),intsrc(ik,ir,14),intsrc(ik,ir,16),e2dion(ik,ir),&
                     srccomp(ik,ir,8),srccomp(ik,ir,16),srccomp(ik,ir,20),srccomp(ik,ir,18),intarea(ik,ir,10),&
                     intarea(ik,ir,12),fluxes(ik,ir,7),fluxes(ik,ir,8),e2dflux(ik,ir),ksb(ik,ir)



             end do
             write (6,*) ' EDGE2D DATA - FULL RING -'//' REGULAR and R-corrected'
             write (6,*) ' Integrals are to cell centres:'
             write (6,450)
             do ik = 1,nks(ir)
                write(6,550) ik,kss2(ik,ir),intsrc(ik,ir,3),intsrc(ik,ir,4),e2dion(ik,ir),srccomp(ik,ir,8),&
                     srccomp(ik,ir,13),srccomp(ik,ir,10),srccomp(ik,ir,4),fluxes(ik,ir,9),fluxes(ik,ir,10)

             end do
             write (6,*) ' EDGE2D DATA - RING - 1/2 CELLS'//' -  REGULAR and R-corrected'
             write (6,*) ' Integrals are to cell centres:'
             write (6,450)
             do ik = 1,nks(ir)
                write(6,550) ik,kss2(ik,ir),intsrc(ik,ir,7),intsrc(ik,ir,8),e2dion(ik,ir),srccomp(ik,ir,8),&
                     srccomp(ik,ir,14),srccomp(ik,ir,12),srccomp(ik,ir,6),fluxes(ik,ir,11),fluxes(ik,ir,12)



             end do
             write (6,*) ' DIVIMP/PINION DATA - FULL RING'//' - REGULAR and R-corrected'
             write (6,*) ' Integrals are to cell centres:'
             write (6,450)
             do ik = 1,nks(ir)
                write(6,550) ik,kss2(ik,ir),intsrc(ik,ir,1),intsrc(ik,ir,2),srccomp(ik,ir,1),srccomp(ik,ir,7),&
                     srccomp(ik,ir,13),srccomp(ik,ir,9),srccomp(ik,ir,3),fluxes(ik,ir,15),fluxes(ik,ir,16)

             end do
             write (6,*) ' DIVIMP/PINION DATA - RING - 1/2 CELLS'//' -  REGULAR and R-corrected'
             write (6,*) ' Integrals are to cell centres:'
             write (6,450)
             do ik = 1,nks(ir)
                write(6,550) ik,kss2(ik,ir),intsrc(ik,ir,5),intsrc(ik,ir,6),srccomp(ik,ir,1),srccomp(ik,ir,7),&
                     srccomp(ik,ir,14),srccomp(ik,ir,11),srccomp(ik,ir,5),fluxes(ik,ir,15),fluxes(ik,ir,16)

             end do

             write (6,*) ' Cell AREA calculations for ring:',ir
             write(6,800)
             do ik = 1,nks(ir)
                in = korpg(ik,ir)
                write(6,900) ik,karea2(ik,ir),areap(in),hro(ik,ir)*drho(ik,ir)*hteta(ik,ir)*dthetag(ik,ir),&
                     hro(ik,ir),drho(ik,ir),hro(ik,ir)*drho(ik,ir),hteta(ik,ir),dthetag(ik,ir),hteta(ik,ir)*dthetag(ik,ir),&
                     karea2(ik,ir)/srccomp(ik,ir,21),hro(ik,ir)*drho(ik,ir)/kbfs(ik,ir)

             end do
800          format(2x,'IK',7x,'KAREA2',8x,'AREAP', 5x,'HrdrHtdt',9x,'Hrho',9x,'dRHO',4x,'Hrho*dRHO',8x,'Hteta',6x,'dTHETAG',3x,&
                  'Hteta*dTHE',4x,'KAREA2/dS',2x,'HroDro*Bt/B')

900          format(i4,11(1x,e12.5))

300          format(a12,6(2x,e14.6))
400          format(2x,'IK',3x,'S (m)',5x,'Intsrc',4x,'Intsrc*R',9x,'SRC',5x,'SRC * R',6x,'R * ds',2x,'SRC * R * ds',4x,'SRC * ds')

500          format(i4,f8.4,7e12.5)
450          format(2x,'IK',3x,'S (m)',5x,'Intsrc',4x,'Intsrc*R',9x,'SRC',5x,'SRC * R',6x,'R * ds',1x,'SRC * R * ds',4x,&
                  'SRC * ds',8x,'G(s)',4x,'G(s,R)/R')

550          format(i4,f8.4,9e12.5)
600          format(i4,e12.5,f8.4,11e12.5)

700          format(2x,'IK',6x,'Area ',4x,' S (m) ',3x,'Intsrc',4x,'Intsrc*R',9x,'SRC',5x,'SRC * R',6x,'R * ds',2x,&
                  'SRC * R * ds',4x,'SRC * ds',4x,'Int Area',4x,'Int A (R)',6x,'G(s)',4x,'G(s,R)/R')
650          format(i4,e12.5,f8.4,13e12.5)

750          format(2x,'IK',6x,'Area ',4x,' S (m) ',3x,'Intsrc',4x,'Intsrc*R',9x,'SRC',5x,'SRC * R',6x,'R * ds',2x,&
                  'SRC * R * ds',4x,'SRC * ds',4x,'Int Area',4x,'Int A (R)',6x,'G(s) ',3x,'G(s,R)/R',5x,'E2DFlux',5x,'S bound')

             !        Close IR loop

          endif
       end do
       return



     end subroutine debugsoliz


     subroutine load_extpowsrc(ext_powsrc,filename,maxnks,maxnrs,nrs,nks,rs,zs,ierr)
       use error_handling
       use common_utilities
       use plasma_overlay  ! shared code functionality from this module
       implicit none
       integer :: maxnrs,maxnks
       integer :: nrs
       integer :: nks(maxnrs)
       real :: ext_powsrc(maxnks,maxnrs),rs(maxnks,maxnrs),zs(maxnks,maxnrs)
       character*(*) :: filename
       integer :: ierr
       integer :: infile
       integer :: ik,ir
       integer :: nr,nz
       character*256 :: line
       logical :: data_read
       !
       ! This reuses the code for loading external plasma background overlays from DTS data
       ! Overlay routines were modified to support loading ti
       !


       ierr = 0
       call find_free_unit_number(infile)

       open(infile,file=trim(filename), status='old',form='formatted',iostat=ierr)

       if (ierr.ne.0) then
          call errmsg('ERROR: MOD_SOL22_INTERFACE : LOAD_POWSRC : ERROR OPENING FILE='//trim(filename)//': IERR=',ierr)
          return
       endif
          
       !
       ! Read in file data - everything is ignored until the DATA: line
       !
       ! DATA:  nrows  ncols 
       !   R(iy) Z(ix) DATA(ix,iy)     - one data point/line 
       !       ...
       !
       ! Multiple data sets can be in the file. Interpolation happens for
       ! each data set in order so if R,Z space for a data set overlaps with
       ! a previous data set only the most recent in the file is used. 
       !
       ! Data is expected on an evenly spaced regular mesh
       ! Both Raxis and zaxis should be in ascending order - these
       ! will be used to identify the cell for interpolation
       ! Loop iterates IR on outside, IZ on inside - data listed by columns
       !
       ! Data is read for each region
       ! minr, maxr, minz, maxz are determined
       !
       ! Read size parameters.
       ! Load data
       ! 
       ! Loop through rs,zs for each cell 
       !  - interpolate R,Z into tabulated data to obtain power(r,z)
       !
       data_read=.false.
       
       do while (ierr.eq.0)
          read(infile,'(a256)',iostat=ierr) line
          ! Note ierr < 0 from this read is treated as an EOF and will exit the loop without issuing an error message
          if (ierr.gt.0) then
             call errmsg('ERROR: MOD_SOL22_INTERFACE : LOAD_EXTPOWSRC : ERROR READING FROM FILE:',filename)
          elseif (ierr.eq.0) then 
             if (ucase(  trim(   line(1:len('DATA:'))     )  ).eq.'DATA:') then                
                read(line(len('DATA:')+1:),*,iostat=ierr) nr,nz
                if (ierr.ne.0) then
                   call errmsg('ERROR: MOD_SOL22_INTERFACE : LOAD_POWSRC : ERROR READING ROW/COLUMNS FROM DATA LINE: IERR=',ierr)
                else

                   ! The load data routine handles storage allocation within the overlay/interpolation module. 
                   ! Read in the data in format R(ir)  Z(iz)   DATA(ir,iz) - R,Z must be in ascending order, listed for each data point and form a regular mesh. 
                   call load_extdata(infile,nr,nz,1,1,ierr)

                   if (ierr.eq.0) then 
                      ! Loop through and interpolate the data onto the SOL22 power array
                      ! Unfortunately, have to loop the entire mesh for each region but the interpolation code first checks [rmin,rmax], [zmin,zmax] to make the process more efficient
                      do ir = 1,nrs
                         do ik = 1,nks(ir)
                            call interpolate_extdata(rs(ik,ir),zs(ik,ir),ext_powsrc(ik,ir),1)
                         enddo
                      enddo
                      data_read = .true.
                   endif
                   ! free up any storage used in the data load/interpolate routine
                   call po_deallocate_storage

                endif
             endif
          endif
       end do

       if (data_read) then ! normal end of execution at end of file - data was read - no error to flag
          ierr = 0
       endif
       
       close(infile)

     end subroutine load_extpowsrc




     
   end module mod_sol22_interface
