module mod_calcsol_interface


  implicit none


contains



  subroutine calcsol_interface (n0_in,te0_in,ti0_in,ringlen_in,npts_in,spts,ne,te,ti,vb,crmb,rizb)
    use error_handling
    use debug_options
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    use mod_sol22
    use mod_sol22_solver
    use mod_sol22_sources
    use mod_sol22_input
    use sol22_debug
    use mod_io_units
    implicit none

    real*8 :: n0_in, te0_in, ti0_in,ringlen_in
    real*8 :: crmb,rizb
    integer :: npts_in
    !real*8 :: spts(*),te(*),ti(*),ne(*),vb(*)
    real*8 :: spts(mxspts),te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts)
    integer :: cprint
    
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



    !COMMON /POWERFLOW/ cve        ,cvi        ,cde        ,cdi
    REAL               cve(MXSPTS),cvi(MXSPTS),cde(MXSPTS),cdi(MXSPTS)

    !REAL Clock2
    ! slmod end

    !REAL deltat

    !integer ir,ircor,midnks,ik,ikstart,ikend
    !real gperpa(maxnks,maxnrs)
    !real oldknbs(maxnks,maxnrs)
    !real grad_oldknbs(maxnks,maxnrs)
    !real oldktibs(maxnks,maxnrs)
    !real oldktebs(maxnks,maxnrs)

    !     Radiation loss term from previous DIVIMP run
    !     Impurity Ionization/recombination potential energy loss to e-

    !real oldkvds(maxnds)
    !real div_tpowls(maxnks,maxnrs),div_tcooliz(maxnks,maxnrs)

    !real div_cool(maxnks,maxnrs) 
    !real gradi,grado,flux,len1,len2
    !real quant2grad,polysidelen
    !external quant2grad,polysidelen

    !integer rc
    !integer errcode,seterror,new_errlevel

    !     For calling calcsol

    real*8 serr
    !integer id
    ! these all have to be passed in to calcsol_inteface since these are the outputs

    real*8 exp_press(mxspts),act_press(mxspts),prad(mxspts)
    real*8 ttarg

    !     Assigning ACTFFRIC

    real*8 int_powrat(3)
    !real*8 find_ffric 
    !external find_ffric

    !      integer ierr

    !     Local variables

    !real ds1,dp1,dt1,nb1,ds2,dp2,dt2,nb2
    !integer applydef(maxnrs,2),
    integer in,pplasma
    !integer in,pplasma,ikfirst,iklast
    integer sol22_iter

    !     Temporary local variables until these are loaded in input file

    data sol22_iter /0/
    !integer:: pfz_dist_opt

    !     TEMPORARY MULTIPLIER

    !      real qesrc_mult
    !      parameter (qesrc_mult=1.0)

    !     Locals for smoothing

    !real*8 :: pfz_dist_param(2)
    !real temid,timid,nmid,smax,asmexp,tmp
    real*8 :: gtarget

    integer :: npts
    integer :: ierr, errcode, seterror, new_errlevel
    !real,external :: getcs
    real*8 :: ck0,ck0i

    integer :: new_unit
    

    
    call pr_trace('MOD_CALCSOL_INTERFACE','START CALCSOL_INTERFACE')

    ! set debugging flag - leave this off - it generates an immense amount of output
    debug_s22 = .false.
    !debug_s22 = .true.


    cprint = 0
    !cprint = 9   ! debugging - produces a lot of output - leave off
    ! this code could be updated by passing through the LIM input value for cprint - but not sure that it is worthwhile

    ck0  = 2000.0
    ck0i =   58.9
    ! turn off pplasma switches - this might be another mechanism to use in LIM for changing options on collector probe flux tubes for example 
    pplasma = 0
    
    call pr_trace('MOD_CALCSOL_INTERFACE','BEFORE READSOL')

    !     jdemod - supplemented by sol22_debug module functionality
    if (debug_sol22.ne.0) call init_sol22_debug(0,0,12)

    ! since only 1/2 ring - turn on debugging immediately
    if (debug_sol22.ne.0) call check_init_record_data(0,0)   
    

    !
    ! Assign some default input values that could change
    !
    
    !     Make a call to initialize debugging if it is on

    call pr_trace('MOD_CALCSOL_INTERFACE','START CALCSOL_INTERFACE')

    pinavail = .false. ! pin data is only avaialble in divimp

    !     Initialization

    sol22_iter = sol22_iter + 1

    rconst = 0.0d0  ! rconst is an array
    !call qzero(rconst,mxspts)


    title = 'DIVIMP'
    mb = crmb
    zb = rizb
    k0e = ck0
    k0i = ck0i

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

    graph    = 1
    graphaux = 1
    graphvel = 1

    !     Adjust Ionization switch settings if appropriate

    !     Check to see if option set to use PIN data to
    !     normalize the analytic sources.

    graphran = 0.0

    pinnorm = 0


    ! Initialization for solver



    te0 = te0_in
    ti0 = ti0_in
    n0  = n0_in
    v0  = -getcs_sol22(sngl(te0),sngl(ti0))
    npts = npts_in
    ringlen = ringlen_in
    halfringlen = ringlen/2.0

    
    actffric = find_ffric(0,1,actlenmom,actlammom)

    !write(0,*) 'FFRIC:',actffric,actlenmom,actlammom
    
    call assign_radiation_parameters(0,0)


    if (forcet.eq.1) then
       ttarg = (((5.0* te0 + 3.5 * ti0)*sqrt(te0+ti0))/ (8.5*sqrt(2.0)))**(2.0/3.0)
       te0 = ttarg
       ti0 = ttarg
    endif


    !
    ! Set cell boundaries
    !
       
    sbnd(1) = 0.0
    do in = 1,npts-1
       !spts(in) = (in-1) * (halfringlen/npts)
       sbnd(in+1) = (spts(in+1)-spts(in))/2.0
    end do
    sbnd(npts+1) = halfringlen
    !sbnd(npts+1) = halfringlen

    sxp = 0.0


    ! Initialize sources
    
    !gperp(ik) = 
    gperp = 0.0

    !qesrc(ik) = 
    qesrc = 0.0

    !qisrc(ik) = 
    qisrc = 0.0

    !radsrc(ik) = div_cool(ik,ir)
    radsrc = 0.0

    !nhs(ik) = PINATOM(IK,IR)
    nhs = 0.0

    !nh2s(ik) = PINMOL(ik,ir)
    nh2s = 0.0

    !ths(ik) = pinena(ik,ir) * 0.6666
    ths = 1.0   ! just set non-zero in case it causes an issue having a zero temperature somewhere - shouldn't be used anyway

    !nhs0(ik) = pinvdist(1,1,ik,ir) + pinvdist(2,1,ik,ir)+ pinvdist(3,1,ik,ir)
    nhs0 = 0.0

    !           Load values from last iteration for PINQID calculations
    !oldne(ik) = oldknbs(ik,ir)
    oldne = 0.0

    !oldte(ik) = oldktebs(ik,ir)
    oldte = 1.0

    !oldti(ik) = oldktibs(ik,ir)
    oldti = 1.0

    !rbnd(ik) = krb(ik,ir)
    rbnd     = 1.0
    
    !momsrc(ik) = pinmp(ik,ir) * fact
    momsrc = 0.0




    ! initialization
    prad = 0.0
    exp_press = 0.0
    act_press = 0.0
    te = 0.0
    ti = 0.0
    ne = 0.0
    vb = 0.0

    ! Initialize power terms to zero to force calculation in mod_sol22_sources:initval
    pae_start = 0.0
    pae_end   = 0.0
    pai_start = 0.0
    pai_end   = 0.0
    
    ! need these for calculating power
    !gammae = 5.0 + gamecor
    !gammai = 2.5 +  0.5* m0**2 * (1.0 + te0/ti0) + gamcor

    ! target power flux
    !pae_start = gammae * te0 * econv * n0 * abs(v0)
    !pai_start = gammai * ti0 * econv * n0 * abs(v0)
    !pae_end = pae_start
    !pai_end = pai_start

    !write(6,*) 'PAE:', pae_start,te0,ti0,n0,abs(v0),econv,gammae
    !write(6,*) 'PAI:', pai_start,te0,ti0,n0,abs(v0),econv,gammai
    
    gtarget = n0 * v0


    spow = spowbeg * ringlen

    !        Set up gperp lengths for particle source compensation

    spow2= spowlen * ringlen
    sgperpbeg = gperpbegf * ringlen
    sgperpend = gperpendf * ringlen


    seterror = 0

    !        Set active switches
    call setsw(-1,pplasma,new_errlevel)

    !        Initialize the ionization and radiation source lengths

    call initlen

    targfact = 1.0

    errcode = 0

    serr = 0.0



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

       stop_gextra_sink = gextra_sink_stop * ringlen

       !           Calculate extra source strength

       if (start_gextra_src.ne.stop_gextra_src) then
          !gextra_src = abs(gtarg(ir,3)) * gextra_mult /abs(start_gextra_src-stop_gextra_src)
          gextra_src = abs(gtarget) * gextra_mult /abs(start_gextra_src-stop_gextra_src)
       else
          !gextra_src = abs(gtarg(ir,3)) * gextra_mult
          gextra_src = abs(gtarget) * gextra_mult
       endif

       !           Calculate extra sink strength
       if (start_gextra_sink.ne.stop_gextra_sink) then
          !gextra_sink = -abs(gtarg(ir,3)) * gextra_mult /abs(start_gextra_sink-stop_gextra_sink)
          gextra_sink = -abs(gtarget) * gextra_mult /abs(start_gextra_sink-stop_gextra_sink)
       else
          !gextra_sink = -abs(gtarg(ir,3)) * gextra_mult
          gextra_sink = -abs(gtarget) * gextra_mult
       endif

    endif


    !        Error correcting branch point
    !        Algorithmic ionization options need to come after the error check point so that they are set correctly
    !        Set up the ionization source - if choice option is selected.

150 continue

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
             ssrcfi = ssrcst + min(dble(alg_ion_src_len),halfringlen)
             ssrclen = ssrcfi - ssrcst
             ssrcmid = (ssrcfi + ssrcst) / 2.0
          else
             ssrcst = 0.0
             ssrcfi = ssrcst + min(dble(alg_ion_src_len),halfringlen)
             ssrclen = ssrcfi - ssrcst
             ssrcmid = (ssrcfi + ssrcst) / 2.0
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
             ssrcfi = min(dble(13.0 - te0),halfringlen)
             ssrclen = ssrcfi - ssrcst
             ssrcmid = (ssrcfi + ssrcst) / 2.0
          else
             ssrcst = 0.0
             ssrcfi = min(dble(alg_ion_src_len),halfringlen)
             ssrclen = ssrcfi - ssrcst
             ssrcmid = (ssrcfi + ssrcst) / 2.0
          endif

       endif

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
             ssrcfi = min(dble(alg_ion_src_len),halfringlen)
             ssrclen = ssrcfi - ssrcst
             ssrcmid = (ssrcfi + ssrcst) / 2.0
          endif

       endif

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
             ssrcfi = min(dble(alg_ion_src_len),halfringlen)
             ssrclen = ssrcfi - ssrcst
             ssrcmid = (ssrcfi + ssrcst) / 2.0
          endif

       endif

    endif
    !

    !        Initialize output arrays to default error values

    te = 10.0
    ti = 10.0
    ne = 1.0d19
    vb = 0.0
    !call dinit(te,mxspts,10.0d0)
    !call dinit(ti,mxspts,10.0d0)
    !call dinit(ne,mxspts,1.0d19)
    !call dinit(vb,mxspts,0.0d0)


    !write (0,'(a,i4,a)') '************   STARTING  ***********'


    !write (0,*) 'First:',te0,ti0,n0,npts


    !        Zero power ratio information in case calcsol is not called.
    int_powrat = 0.0

    ! jdemod - move initialization of ierror out of mod_sol22 (calcsol) to the calling routine
    ierror = 0


    !call calcsol (spts,npts,errcode,serr,te,ti,ne,vb,exp_press,act_press,prad,ir,irlim1,int_powrat,cprint,&
    !                  cve,cvi,cde,cdi) 
    call calcsol (spts,npts,errcode,serr,te,ti,ne,vb,exp_press,act_press,prad,0,0,int_powrat,cprint,&
         cve,cvi,cde,cdi) 


    !           If a negative N error or NaNQ occurs without error correction -
    !           Turn error correction ON.

    if ((   actswerror.ge.1.0.and.(errcode.eq.3.or.errcode.eq.4.or.errcode.eq.5.or.errcode.eq.6.or.errcode.eq.7)).or.&
         (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7))) then

       write(0,*) 'ERROR:',actswerror,errcode

       if (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7).and.seterror.eq.0) then
          actswerror = 10
          seterror = 1
          !               ERROR - negative N has been found EVEN with highest level
          !                       of error correction - issue error messages and stop.

       elseif (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7).and.seterror.eq.1) then
          seterror = 2

       elseif (seterror.eq.2) then


          call errmsg('SOLASCV:SOL22',' Unsolvable Negative N error encountered. Program Stopping')

          stop 'SOL22:NEG N'

          !           Record error

       endif

       call setsw(-2,pplasma,new_errlevel)


       call initlen

       !write(0,*) 'ERROR loop back'
       goto 150

       !        Save background density and temperature
       !        Target first in case it has changed

    endif


    if (debug_sol22.ne.0) call check_print_data(0,0)


    !
    ! Add code to calculate electric field 
    !


    !DO 500 IK = ikfirst,iklast
    !   DS1 = KSS2(IK,IR) - KSS2(IK-1,IR)
    !   DP1 = KNBS(IK,IR)*KTEBS(IK,IR)-KNBS(IK-1,IR)*KTEBS(IK-1,IR)
    !   DT1 = (KTEBS(IK,IR)-KTEBS(IK-1,IR))
    !   NB1 = 0.5*(KNBS(IK,IR)+KNBS(IK-1,IR))
    !   DS2 = KSS2(IK+1,IR) - KSS2(IK,IR)
    !   DP2 = KNBS(IK+1,IR)*KTEBS(IK+1,IR)-KNBS(IK,IR)*KTEBS(IK,IR)
    !   DT2 = (KTEBS(IK+1,IR)-KTEBS(IK,IR))
    !   NB2 = 0.5*(KNBS(IK+1,IR)+KNBS(IK,IR))

    !            WRITE(6,*) 'KES:',IK,IR,KES(IK,IR)

    !   KES(IK,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)+ (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))

    !         End of test on OFIELD

    !500          CONTINUE

    !       If smoothing is turned ON

    !        WF'95:    CORRECTING(?) FOR PLATE ASSYMETRY

    !        profile scaling of ne,Te,Ti to match midpoint values,
    !        exponent asmexp determines smoothing range

    !if (switch(swsmooth).eq.1.0) then
    !   TEMID = 0.5*(KTEBS(MIDNKS+1,IR) + KTEBS(MIDNKS,IR))
    !   TIMID = 0.5*(KTIBS(MIDNKS+1,IR) + KTIBS(MIDNKS,IR))
    !         VMID = 0.5*(KVHS(MIDNKS+1,IR) + KVHS(MIDNKS,IR))
    !   NMID = 0.5*(KNBS(MIDNKS+1,IR) + KNBS(MIDNKS,IR))
    !   SMAX = KSMAXS2(IR)
    !         write(6,*) 'assymetry', ir,ktibs(midnks,ir),
    !     >   ktibs(midnks+1,ir), timid

    !   ASMEXP = 1.0
    !   DO IK = ikstart , ikend
    !      IF (IK.LE.MIDNKS) THEN
    !         TMP = (KSS2(IK,IR) / (SMAX / 2.0))**ASMEXP
    !         KTEBS(IK,IR) = KTEBS(IK,IR)*(1.0 +TMP*(TEMID/KTEBS(MIDNKS,IR) - 1.0))
    !         KTIBS(IK,IR) = KTIBS(IK,IR)*(1.0 +TMP*(TIMID/KTIBS(MIDNKS,IR) - 1.0))
    !           KVHS(IK,IR) = KVHS(IK,IR)*(1.0 +
    !     >        TMP*(VMID/KVHS(MIDNKS,IR) - 1.0))
    !         KNBS(IK,IR) = KNBS(IK,IR)*(1.0 +TMP*(NMID/KNBS(MIDNKS,IR) - 1.0))
    !      ELSE
    !         TMP = ((SMAX-KSS2(IK,IR)) / (SMAX / 2.0))**ASMEXP
    !         KTEBS(IK,IR) = KTEBS(IK,IR)*(1.0 +TMP*(TEMID/KTEBS(MIDNKS+1,IR) - 1.0))
    !         KTIBS(IK,IR) = KTIBS(IK,IR)*(1.0 +TMP*(TIMID/KTIBS(MIDNKS+1,IR) - 1.0))
    !           KVHS(IK,IR) = KVHS(IK,IR)*(1.0 +
    !     >        TMP*(VMID/KVHS(MIDNKS+1,IR) - 1.0))
    !         KNBS(IK,IR) = KNBS(IK,IR)*(1.0 +TMP*(NMID/KNBS(MIDNKS+1,IR) - 1.0))
    !      ENDIF

    !      End of Smoothing IF

    !   ENDDO

    !      Branch location for virtual rings

    !endif

    return

  end subroutine calcsol_interface







  
end module mod_calcsol_interface
