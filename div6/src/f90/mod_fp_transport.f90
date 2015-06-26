module mod_fp_transport

  !
  ! This module implements the rules for ion transport in the peripheral plasma region. The options 
  ! implemented are a subset of the transport options supported for the main grid. This module is designed
  ! with a certain amount of flexibility - the long term objective is to replace the DIVIMP ion transport code
  ! with a self-contained module - this will be done when time and needs permit.This will make the code more maintainable
  ! in the long run and should allow for code reuse in other contexts - like Yarongs code or sharing an ion transport 
  ! module between the HC code and the core DIVIMP routines. 
  !
  !
  !
  ! transport coefficients
  !
  ! Update_s
  !
  ! Update_cross
  !
  ! Forces
  !
  ! Exit conditions - 1) parallel loss - struck wall at end of FP - determine approximate location and wall/target index struck
  !                   2) crossfield loss - struck wall awsy from target - at S != 0
  !                   3) return to grid - possibly in a new cell - new ik and new S
  !
  ! Use specifications for IRWALL and IRTRAP for special conditions within FP? NO - use FP specific arrays 
  !
  ! To do and possible enhancements:
  !
  ! Add change of state code - use a callback? or check for HC vs. normal?
  !                   - handle transitions to neutral state within FP
  !                   - handle initial ionization of particle within the FP (?)
  !                   - introduce neutral FP transport (?)
  !
  !                   - change of state requires calls to set the particle characteristics including mass if it has changed
  !                     followed by a call to init_taus and eval_taus to set the characteristic times appropriately - it also requires
  !                     a call to the fp_set_tgrad_coeffs routine if the mass or charge has changed. 

  type particle_characteristics

     real :: s,cross,vel
     real :: temi
     real :: slast
     !
     ! No use in this code at the present time for r,z data
     ! real :: r,z
     !
     integer :: ik,ir,iz,istate
     real :: mass
     integer :: fp_reg

  end type particle_characteristics

  type bg_plasma

     real :: ne,te,ti,ef,vb,tge,tgi

  end type bg_plasma

  !
  ! The geometry checking code is handled outside of the module for the time being
  ! because none of the geometry data is presently available in a module format 
  !
  !type cell_characteristics
  !
  !   real :: inner_cross_bnd,outer_cross_bnd
  !   real :: lower_s_bnd,upper_s_bnd
  !
  !end type cell_characteristics

  !type transport_coefs
  !   real :: lfps,lllfps,lfss,lfts
  !   real :: alphai, betai
  !end type transport_coefs


  ! make the local module variables static

  save

  type (particle_characteristics),private :: part
  type (bg_plasma),private :: plasma 

  !type (cell_characteristics),private :: cell


  !
  ! Local options and variables
  !


  !
  ! DIVIMP options - hc, collision, heating, stopping, temperature gradient, drifts
  !
  integer,private :: fp_hc_follow_opt
  integer,private :: fp_optb,fp_optc,fp_optd,fp_optm,fp_optn 
  integer,private :: fp_cpdrft 
  integer,private :: fp_pinchopt

  real,private :: mb,fp_rizb,mi,mfp ! masses

  real,private :: fp_timestep ! timestep


  real,private :: fp_dperp,fp_kperps ! Cross field diffusion and step size

  !
  ! legacy options and inputs 
  !
  real,private ::   fp_czo,fp_chzo,fp_ctemav,fp_czenh
  real,private :: fp_cstgrad
  integer,private :: fp_cizeff
  integer,private :: fp_irspec 

  !
  ! Temperature gradient range
  !
  real,private :: zero_tgrad_smin,zero_tgrad_smax

  !
  ! Drift velocity range
  ! 
  !real,private :: fp_cdrftv_start,fp_cdrftv_end  ! factors for start and end of drift region
  real,private :: fp_drfts_start, fp_drfts_end ! actual S values for start and end of drift region
  real,private :: fp_drftv                     ! value of drift velocity in peripheral region
  !
  ! Count of random numbers used in fp for current particle
  !
  integer,private :: ran_used 

  !
  ! Transport coefficients
  !

  real,private :: lfps,lllfps,lfss,lfts
  real,private :: alphai, betai


  logical,private :: mod_debug_fp

contains


  subroutine fp_set_options(cioptb,cioptc,cioptd,cioptn,cioptm,cpdrft,pinchopt,crmb,crmi,rizb,&
                           &hc_follow_option,cdperpfp,&
                           &cstgrad,czo,chzo,czenh,cizeff,ctemav,irspec,qtim,debug_fp)

    implicit none
    integer :: cioptb,cioptc,cioptd,cioptn,cioptm,hc_follow_option,cpdrft,pinchopt
    real :: crmb,crmi,rizb,cdperpfp
    real :: cstgrad

    integer :: cizeff,irspec,czo
    real :: chzo,czenh,ctemav,qtim

    logical :: debug_fp

    fp_optb = cioptb   ! friction option
    fp_optc = cioptc   ! stopping option
    fp_optd = cioptd   ! heating  option
    fp_hc_follow_opt = hc_follow_option ! hc flag indicating whether hc electric field mod is active

    fp_optm = cioptm   ! electron temperature gradient force option
    fp_optn = cioptn   ! ion temperature gradient force option

    fp_cpdrft = cpdrft ! parallel drift option 
    fp_pinchopt = pinchopt

    mb = crmb ! background plasma mass
    fp_rizb= rizb ! background plasma charge
    mi = crmi ! base impurity ion mass

    fp_irspec = irspec ! special ring designation not really needed for fp - but some legacy code requires it

    fp_timestep = qtim  ! fp time step
    fp_dperp = cdperpfp ! fp diffusion rate
    fp_kperps = sqrt(2.0*fp_dperp*fp_timestep) ! fp cross field step size 

    ! The following are all inputs required to calculate the parallel collision times used for the velocity diffusion options
    ! These are loaded once since they should not change during the run. 
    !
    ! CTEMAV is supposed to be the average neutral temperature - however the application of the code that uses this is so outdated 
    ! that it is useless. It only uses one value for CTEMAV applied over the entire grid which is likely not a very useful 
    ! approximation. In addition, the options which utilize it are not used at the present time.However. since existing legacy 
    ! code requires a value be assigned - a default value of 2.0eV is assigned if CTEMAV does not contain reasonable data (e.g. 0.0) 
    !

    fp_czo    = czo
    fp_chzo   = chzo
    fp_czenh  = czenh
    fp_cizeff = cizeff
    fp_cstgrad= cstgrad
    !
    ! Assign default value of CTEMAV - CTEMAV is the average temperature of the neutrals at ionization - calculated by neut
    !                                  but it is not available in the HC module or when the fp module is first initialized
    !                                  so a default value is assigned for now. 
    !

    if (ctemav.le.0.0) then
       fp_ctemav = 2.0
    else
       fp_ctemav = ctemav
    endif

    mod_debug_fp = debug_fp
    !mod_debug_fp = .true.
    !write(0,*) 'SETTING MOD_DEBUG_FP:', mod_debug_fp

  end subroutine fp_set_options

  subroutine fp_set_ctemav(ctemav)
    implicit none
    real ctemav
    ! 
    ! called to set CTEMAV to the value returned from neut
    !
    fp_ctemav = ctemav
  end subroutine fp_set_ctemav




  subroutine fp_init_particle(s,cross,vel,temi,ik,ir,iz,istate,fp_reg,crmfp,smax,drftv,sdrftv_start,sdrftv_end)

    use taus

    implicit none
    real,intent(in) :: s,cross,vel,temi,crmfp,smax
    integer,intent(in) :: ik,ir,iz,istate,fp_reg
    real,intent(in) :: drftv,sdrftv_start,sdrftv_end ! drift velocity for specific far periphery region

    !
    ! Initialize transport coefficient routine
    !
    mfp = crmfp
    !
    ! Initialize the particle characteristics
    !
    call fp_set_particle_data(s,cross,vel,temi,ik,ir,iz,istate,mfp,fp_reg)

    !
    ! Initialize transport coefficient routine
    ! 
    if (mod_debug_fp) then 
       !write(6,'(a,5i6,10g18.10)') 'FP INIT TAUS:',fp_optb,fp_optc,fp_optd,fp_irspec,fp_cizeff,fp_rizb,mb,mfp,&
       !                                           &fp_czenh,fp_ctemav,fp_timestep
       !write(0,'(a,5i6,10g18.10)') 'FP INIT TAUS:',fp_optb,fp_optc,fp_optd,fp_irspec,fp_cizeff,fp_rizb,mb,mfp,&
       !                                           &fp_czenh,fp_ctemav,fp_timestep
    endif

    call init_taus(mb,mfp,fp_rizb,fp_optb,fp_optc,fp_optd,fp_czenh,fp_cizeff,fp_ctemav,fp_irspec,fp_timestep)

    !
    ! Initialize the plasma characteristics - this includes setting the transport coefficients with a call to eval_taus
    !
    ! All calls to set_fp_plasma as the particle changes cells will update the transport coefficients
    !
    call fp_set_plasma(ik,fp_reg)

    !
    ! Initialize the local conditions
    !

    !call set_cell_characteristics(ik,fp_reg)

    !
    ! Initialize temperature gradient coefficients - needs to be called again if charge or mass of fp particle changes
    ! The particle initialization routine must be called before this so that IZ (the particle charge state) is properly set
    !
    
    call fp_set_tgrad_coefs

    !
    ! Set S limits for temperature gradient forces
    !
    zero_tgrad_smin = fp_cstgrad * smax
    zero_tgrad_smax = (1.0-fp_cstgrad)* smax

    fp_drfts_start  = sdrftv_start
    fp_drfts_end    = sdrftv_end
    fp_drftv        = drftv


  end subroutine fp_init_particle


  subroutine fp_follow_particle(rc,imp,cist,cistfp,cstmax,rsect,zsect,nrands)
    implicit none
    real*8 :: cistfp,cist
    real :: cstmax
    integer rc
    integer nrands,imp
    real :: rsect,zsect

    integer is_followed,target_impact,wall_impact,return_to_grid,time_limit_exceeded

    parameter( is_followed         = 0,&
               return_to_grid      = 1,&
               time_limit_exceeded = 2,&
               wall_impact         = 3,&
               target_impact       = 4)



    !write(6,'(a,4i6,4(1x,g12.5))') 'FP:DEBUG1:',imp,rc,part%ik,part%ir,part%s,part%cross,cistfp

    !
    !     Initialize
    !
    !     time in far periphery is initialized to 0.0
    !
    cistfp = 0.0

    !
    ! Initialize random number count
    !
    ran_used = 0
    !
    !     Initialize the coordinates of wall impact to 0
    !
    rsect = 0.0
    zsect = 0.0
    !
    !     Initialize particle status
    !
    rc=Is_followed
    !
    !     Loop while particle remains in the periphery
    !
    do 
       !
       !        Execute parallel step - recalculate S 
       !
       call fp_parallel_step
       !
       !        Check for target impact and reassign cell association
       !
       call fp_check_s(rc)
       !
       !        Exit if particle has reached the target
       !
       if (rc.ne.0) exit
       !
       !        Execute cross-field step
       !
       call fp_cross_step(imp,cist+cistfp)
       !
       !        Check cross-field location for grid re-entry or wall loss
       !
       call fp_check_cross(rsect,zsect,rc)
       !
       !        Exit if particle has reached wall or returned to grid
       !
       if (rc.ne.0) exit
       !
       !        Update particle temperature
       !
       call fp_update_ion_temperature

       !
       !        Check for particle change of state - could be either atomic or molecular - need flag
       !
       !         call fp_check_change_state
       !

       !
       !        Update time step count
       !
       cistfp = cistfp + 1.0
       !
       !        check total time 
       !
       if ((cist+cistfp).gt.cstmax) then
          rc = time_limit_exceeded
          !
          !           Add loop exit statement in case other exit conditions are later added after this one
          !     
          exit
       endif

       !write(6,'(a,4i6,4(1x,g12.5))') 'FP:DEBUG2:',imp,rc,part%ik,part%ir,part%s,part%cross,cistfp
 
    end do

    nrands = ran_used


  end subroutine fp_follow_particle



  subroutine fp_set_plasma(ik,fp_reg)
    use taus
    implicit none
    integer :: ik,fp_reg
    real :: ne,te,ti,ef,vb,tgrade,tgradi

    call fp_get_plasma(ik,fp_reg,ne,te,ti,vb,ef,tgrade,tgradi)

    if (ti.eq.0.0) then 
       write(0,'(a,2i6,10(1x,g12.5))') 'FP_SET_PLASMA Ti=0.0',ik,fp_reg,ne,te,ti,vb,ef,tgrade,tgradi
    endif

    plasma%ne = ne
    plasma%te = te
    plasma%ti = ti
    plasma%vb = vb
    plasma%ef = ef
    plasma%tge = tgrade
    plasma%tgi = tgradi


    !
    ! Update transport coefficients
    !

    call eval_taus(part%ik,part%ir,part%iz,plasma%ne,plasma%ti,lfps,lllfps,lfss,lfts)

    if (mod_debug_fp) then 
       !write(6,'(a,2i6,12g18.10)') 'FP TAUS:',part%ik,part%ir,part%iz,plasma%ne,plasma%ti,lfps,lllfps,lfss,lfts
       !write(0,'(a,2i6,12g18.10)') 'FP TAUS:',part%ik,part%ir,part%iz,plasma%ne,plasma%ti,lfps,lllfps,lfss,lfts
    endif

  end subroutine fp_set_plasma



  subroutine fp_get_plasma(ik,fp_reg,ne,te,ti,vb,ef,tgrade,tgradi)
    use mod_fperiph
    implicit none
    integer ik,fp_reg,id
    real,intent(out) :: ne,te,ti,ef,vb,tgrade,tgradi
    !
    !      include 'params'
    !      include 'fperiph_com'
    !
    !     Return the fp plasma for the specified cell
    !
    ne = fp_plasma(ik,fp_reg,1)
    te = fp_plasma(ik,fp_reg,2)
    ti = fp_plasma(ik,fp_reg,3)
    vb = fp_plasma(ik,fp_reg,4)
    ef = fp_plasma(ik,fp_reg,5)
    tgrade = fp_plasma(ik,fp_reg,6)
    tgradi = fp_plasma(ik,fp_reg,7)

    if (mod_debug_fp) then 
       !write(6,'(a,2i6,7g12.5)') 'FP PLASMA:',fp_reg,ik,(fp_plasma(ik,fp_reg,id),id=1,7)
       !write(0,'(a,2i6,7g12.5)') 'FP PLASMA:',fp_reg,ik,(fp_plasma(ik,fp_reg,id),id=1,7)
    endif 

  end subroutine fp_get_plasma




  subroutine fp_set_tgrad_coefs
    use global_parameters
    implicit none
    real :: mu,fact

    !
    ! Only recalculate Alpha and beta if the mass changes
    !
    ! Note: the temperature gradient data itself is scaled by qtim * qtim * emi /crmi
    !       In order to correct for changing masses properly - the coefficients of the forces
    !       must be set to compensate for any difference between the current particle mass 
    !       and the base impurity mass - in addition to the explicit dependence in the coefficients
    !
    ! This routine only needs to be called again if the mass or charge of the particle in the fp changes
    !
    !

    IF     (FP_OPTM.EQ.0) THEN
       ALPHAI = 0.0
    ELSEIF (FP_OPTM.EQ.1.or.fp_optm.eq.3) THEN
       ALPHAI = 0.71 * REAL(PART%IZ*PART%IZ)
    ELSEIF (FP_OPTM.EQ.2) THEN
       ALPHAI = 1.5*(1.0-0.6934*(1.3167**(-PART%IZ)))*REAL(PART%IZ*PART%IZ)
    ENDIF

    IF     (FP_OPTN.EQ.0) THEN
       BETAI = 0.0
    ELSEIF (FP_OPTN.EQ.1.or.fp_optn.eq.3) THEN
       MU = MFP / (MFP+MB)
       BETAI=-3.0*(1.0-MU-5.0*ROOT2*(1.1*MU**2.5-0.35*MU**1.5)*  REAL(PART%IZ*PART%IZ)) / (2.6 - 2.0*MU + 5.4*MU*MU)
    ELSEIF (FP_OPTN.EQ.2) THEN
       BETAI= fp_CHZO * REAL(PART%IZ*PART%IZ) / (fp_CZO + SQRT(0.5*(1.0+MB/MFP)))
    ENDIF


    ! At the moment this correction is now included in the force calculation routines - if this code
    ! ends up being too slow it can be moved here to save a few CPU cycles/particle step
    !
    !fact = mi/mfp
    !
    !alphai = alphai * fact
    !betai  = betai * fact

  end subroutine fp_set_tgrad_coefs



  subroutine fp_set_particle_data(s,cross,vel,temi,ik,ir,iz,istate,mass,fp_reg)
    implicit none
    real,intent(in) :: s,cross,vel,temi,mass
    integer,intent(in) :: ik,ir,iz,istate,fp_reg

    part%s = s
    part%cross = cross
    part%vel = vel
    part%temi = temi
    part%slast = s
    part%ik = ik
    part%ir = ir
    part%iz = iz
    part%istate=istate
    part%mass = mass
    part%fp_reg = fp_reg

    
    if (mod_debug_fp) then 
       write(6,'(a,5i6,10g12.5)') 'FP SET PART:',ik,ir,iz,istate,fp_reg,s,cross,vel,temi,mass
    endif

  end subroutine fp_set_particle_data


  subroutine fp_get_particle_data(s,cross,vel,temi,ik,ir,iz,istate)
    implicit none
    real,intent(out) :: s,cross,vel,temi
    integer,intent(out) :: ik,ir,iz,istate

    s = part%s
    cross = part%cross
    vel = part%vel
    temi=part%temi
    ik = part%ik
    ir = part%ir
    iz = part%iz
    istate = part%istate

    if (mod_debug_fp) then 
       write(6,'(a,5i6,10g12.5)') 'FP GET PART:',ik,ir,iz,istate,part%fp_reg,s,cross,vel,temi,part%mass
    endif

  end subroutine fp_get_particle_data


  real function fp_force_fe()

    implicit none

    real local_efield

    if (fp_hc_follow_opt.ne.0.and.part%fp_reg.eq.0) then 

       call hc_electric_field_mod(part%ik,part%ir,part%iz,part%s,local_efield)

    else

       local_efield =  plasma%ef

    endif

    fp_force_FE    = real(part%iz) * local_efield * mi/mfp

    if (mod_debug_fp) then 
       !write(6,'(a,i6,12g18.10)') 'FP FE :',part%iz,local_efield,mi,mfp,fp_force_fe
    endif

  end function fp_force_fe





  real function fp_force_ff(fvh)
    use hc_get
    implicit none
    integer ik,ir
    real fvh

    real kfssmod

    if (part%fp_reg.gt.0) then 
       kfssmod = 1.0
    else
       kfssmod = gkfssmod(part%ik,part%ir)
    endif

    fp_force_FF    = KFSSMOD  * lfss * (FVH-part%vel)

    if (mod_debug_fp) then 
       !write(6,'(a,12g18.10)') 'FP FF :',part%fp_reg,kfssmod,lfss,fvh,part%vel,fp_force_ff
    endif

  end function fp_force_ff


  real function fp_force_fig()
    implicit none

    if (fp_optn.eq.3.and.part%s.gt.zero_tgrad_smin.and.part%s.lt.zero_tgrad_smax) then
       fp_force_fig = 0.0
    else
       fp_force_FIG   = BETAI * plasma%tgi * mi/mfp
    endif

    if (mod_debug_fp) then 
       !write(6,'(a,12g18.10)') 'FP FIG:',betai,plasma%tgi,mi,mfp,fp_force_fig
    endif

  end function fp_force_fig



  real function fp_force_feg()
    implicit none


    !     Calculate modifications to forces if any

    if (fp_optm.eq.3.and.part%s.gt.zero_tgrad_smin.and.part%s.lt.zero_tgrad_smax) then

       fp_force_feg = 0.0

    else 

       fp_force_FEG   = ALPHAI * plasma%tge * mi/mfp

    endif

    if (mod_debug_fp) then 
       !write(6,'(a,12g18.10)') 'FP FEG:',alphai,plasma%tge,mi,mfp,fp_force_feg
    endif

  end function fp_force_feg



  real function fp_force_fvg()
    implicit none
    integer ik,ir

    fp_force_FVG   = 0.0

  end function fp_force_fvg



  subroutine fp_force_col(dvpara,dspara,vpara,spara)
    implicit none

    real dvpara,dspara,vpara,spara


    ! 
    real ran


    ! 
    !     Only need one random number since VPARA and SPARA based methods of 
    !     diffusion based collisonal transport are mutually exclusive.
    !

    ran_used = ran_used+1 

    call getran(ran)

    dspara= SIGN(SPARA,ran-0.5)
    dvpara= sign(vpara,ran-0.5)

    if (mod_debug_fp) then 
       !write(6,'(a,12g18.10)') 'FP COL:',spara,dspara,vpara,dvpara
    endif

  end subroutine fp_force_col






  subroutine fp_parallel_step
    ! Calculate the new S value based on forces in FP calculated from the
    ! specified fp plasma. 
    ! - reiser force formulation is not supported at the present time
    ! - changing mass dependence on forces and coefficients is supported

    implicit none

    !
    ! local
    ! 
    real bg_drftvel,imp_drftvel

    real spara,dspara,vpara,dvpara

    real fvh
    real ff,feg,fig,fvg,fe
    real quant

    real ds_dperpz
    real fp_delta_s_dperpz
    external fp_delta_s_dperpz

    !
    ! jdemod - put in support for parallel pinch effect but it is not activated - need to consider
    !          whether this type of drift would exist in fp.
    !
    ! real ds_pinch,ds_kpinchs
    ! external ds_kpinchs
    ! 
    ! ds_pinch = ds_kpinchs(part%ik,part%ir)
    !

    !if (mod_debug_fp) then 
    !   write(6,'(a,12g18.10)') 'FP UPDATE S1:',part%s,part%cross,part%vel,quant,dvpara,dspara,ff,fe,feg,fig,fvg
    !   write(0,'(a,12g18.10)') 'FP UPDATE S1:',part%s,part%cross,part%vel,quant,dvpara,dspara,ff,fe,feg,fig,fvg
    !endif

    call fp_set_drift_velocity(bg_drftvel,imp_drftvel)


    call fp_set_collisional_step(spara,vpara,lfps,lllfps)


    ! plasma flow velocity

    fvh = plasma%vb + bg_drftvel

    !
    ! Reiser not implemented in the periphery
    !

    ! forces

    FF    = fp_force_ff(fvh)
! slmod begin
!...ifort compiler needs brackets:
    FIG   = fp_force_fig()
    FEG   = fp_force_feg()
    FVG   = fp_force_fvg()
    FE    = fp_force_fe()
!
!    FIG   = fp_force_fig
!    FEG   = fp_force_feg
!    FVG   = fp_force_fvg
!    FE    = fp_force_fe
! slmod end
    !
    ! diffusion in 3D
    !     
    DS_DPERPZ = fp_delta_s_dperpz(part%ik,part%ir,ran_used)

    !
    ! Collisional term
    !
    call fp_force_col(dvpara,dspara,vpara,spara)       

    !     
    !     Calculate change in velocity due to forces
    !     
    QUANT = FF + FE + FEG + FIG + FVG
    !     
    !     
    !     changed order of s and vel calculation... Krieger IPP 12/94
    !     
    !     
    !     Record last S value in order to detect parallel reflection situations  
    !     
    part%slast = part%s
    !     
    part%S = part%S + part%VEL + 0.5 * (QUANT+dvpara) + dspara + IMP_DRFTVEL + ds_dperpz
    !
    !     Update particle parallel velocity
    !     

    part%vel   = part%vel + QUANT + dvpara

    if (mod_debug_fp) then 
       write(6,'(a,12g18.10)') 'FP UPDATE S:',part%slast,part%s,part%cross,part%vel,quant,dvpara,dspara,&
                                             &ff,fe,feg,fig,fvg
    endif

  end subroutine fp_parallel_step




  subroutine fp_cross_step(imp,cist)
    use mod_fperiph
    implicit none
    integer ik,ir,imp
    real*8,intent(in) :: cist

    ! NOTE: cist is total time and is not a separate value - it should not be changed in this routine

    !
    ! In/OUT diffusion probability of 0.5 is used. The algorithms in DIVIMP 
    ! are based on cell geometry and the changing ratio of the inner cell length
    ! to outer cell length over a distance of one cross field step. This concept does
    ! not really apply to the fp as constituted here. 
    !

    logical vr_assigned
    integer ierr
    real kprob,pinchvel

    real getranf
    external getranf

    !
    ! Set IN/OUT probability to 0.5
    !
    kprob = 0.5

    ! 
    ! FP needs to consider the sign conventions for pinch velocity. 
    ! In DIVIMP +ve Pinchvel is inward while -ve is outward. 
    ! However, in the FP - CROSS starts at 0.0 and INCREASES to the value
    ! at the wall. Thus an outward pinch would need to be +ve in the fp. 
    ! This means that the sign of the pinch velocity has to be changed in 
    ! the FP so that the interpretation is the same. 
    ! This is dependent on the FP region - for the main FP a sign change is 
    ! required - this is not the case for the PFZ since the PFZ wall is 
    ! the same topologically as the inner edge of the confined plasma. 
    !
    call set_pinch_velocity(part%ik,part%ir,ran_used,imp,cist,pinchvel,vr_assigned,ierr)

    !
    ! Reverse sign of pinchvel for main FP to remain consistent with DIVIMP sign convention
    !
    if (part%fp_reg.eq.fp_main) pinchvel = -pinchvel

    if (fp_pinchopt.eq.4.and.vr_assigned) then 
       part%CROSS = part%CROSS + PINCHVEL
    else
       ran_used = ran_used + 1
       part%CROSS = part%CROSS + SIGN (fp_kperps, kprob-getranf()) + PINCHVEL
    endif

    if (mod_debug_fp) then
       write(6,'(a,12g18.10)') 'FP UPDATE C:',part%cross,fp_kperps,pinchvel,part%temi
    endif

  end subroutine fp_cross_step

 
  subroutine fp_set_drift_velocity(bg_drftvel,imp_drftvel)
    implicit none
    real bg_drftvel,imp_drftvel

    real s
    integer ir
    !     
    !-----------------------------------------------------------------------
    !     
    !     Set the value of the drift velocity displacement 
    !     depending on the particles position on the ring.
    !     
    !     Drifts only apply in the main SOL or PFZ at the moment
    !     
    !     Assign zero values as default
    !     
    imp_drftvel = 0.0
    bg_drftvel = 0.0
    !     
    !     Replace with assigned values (which may be zero) depending
    !     on options - option 1,3 are direct to the impurity particle -
    !     option 2 is assigned as a background flow component and coupled
    !     to the impurity through friction
    !     
    !     Region limitations are enforced in the setup_drftv routine
    !     
    if ((fp_cpdrft.eq.1..or.fp_cpdrft.eq.3).and.(part%s.ge.fp_drfts_start.and.part%s.le.fp_drfts_end)) then 
       imp_drftvel = fp_drftv
       bg_drftvel = 0.0
    elseif (fp_cpdrft.eq.2.and.(s.ge.fp_drfts_start.and.s.le.fp_drfts_end)) then 
       imp_drftvel= 0.0
       bg_drftvel = fp_drftv
    endif

  end subroutine fp_set_drift_velocity




  subroutine fp_update_ion_temperature
    use global_parameters
    implicit none

    part%TEMI = MAX (LO, part%TEMI + (plasma%ti-part%TEMI) * LFTS)


  end subroutine fp_update_ion_temperature



  subroutine fp_set_collisional_step(spara,vpara,lfps,lllfps)
    use global_parameters
    implicit none
    !
    ! NOTE: this only supports a subset of the main code options - in particular - it is setup to 
    !       use the velocity diffusion options and not the spatial diffusion ones. 
    !


    real spara,vpara
    real lfps,lllfps

    logical :: warning_issued = .false.
    real :: ran1, ran2, rgauss
    real :: getranf
    external getranf



    !
    ! Collision option 1 turns off parallel diffusive transport
    !
    if (fp_optb.eq.1) then
       spara = 0.0
       vpara = 0.0
    elseif (fp_optb.eq.6.or.fp_optb.eq.11) then
       !     
       !     Collision option 6 ... changes on each time-step
       !     
       spara  = 0.0
       vpara = lllfps
       !     
       !     Collision option 14 - set to zero for regions of low 
       !     collisionality - specified in input as a fraction of SMAX - 
       !     ideally less than 0.5.
       !     
    elseif (fp_optb.eq.14.and.(part%s.gt.zero_tgrad_smin.and.part%s.lt.zero_tgrad_smax)) then
       spara  = 0.0
       vpara  = 0.0
       !     
       !     Collision options 12, 13 and 14 (when 14 is not set to zero)
       !     
    elseif (fp_optb.eq.12.or.fp_optb.eq.13.or.fp_optb.eq.14) then
       !     
       !     Collision option 12 ... changes on each time-step
       !     Collision option 13 ... tau para / (1+Mb/Mi)
       !     
       spara  = 0.0
       vpara = lllfps
7702   ran_used=ran_used + 1
       call getran(ran1)
       if (ran1.eq.0.0) goto 7702
       ran_used=ran_used+1
       call getran(ran2)
       rgauss = sqrt(-2.0* log(ran1))*cos(2.0*PI*ran2)
       vpara = vpara * rgauss

       !
       ! IF another option has been specified - the code issues a warning and turns off parallel diffusion
       !

    else
       spara = 0.0
       vpara = 0.0
       if (.not.warning_issued) then 
          write (0,*) 'FP_SET_COLLISIONAL_STEP:',' WARNING:',' INVALID PARALLEL DIFFUSION OPTION&
               & SPECIFIED FOR PERIPHERY TRANSPORT - PARALLEL DIFFUSION TURNED OFF IN PERIPHERY' 
          write (6,*) 'FP_SET_COLLISIONAL_STEP:',' WARNING:',' INVALID PARALLEL DIFFUSION OPTION&
               & SPECIFIED FOR PERIPHERY TRANSPORT - PARALLEL DIFFUSION TURNED OFF IN PERIPHERY' 
          warning_issued = .true.
       endif
    endif


  end subroutine fp_set_collisional_step




  subroutine fp_check_s(rc)
    use mod_fperiph
    implicit none
    integer rc
    !
    !      include 'param'
    !      include 'fperiph_com'
    !
    !     RC = 0 - cell index assigned
    !     RC = 4 - struck target
    !
    !        Status codes set to match those in FPERIPH
    !
    !        Is_Followed   : particle_status=0   (Particle is still being followed) 
    !        Return to grid: particle_status=1   (Cross < 0 - assuming Cross is set to 0 when entering FP) 
    !        Maximum time  : particle_status=2   (cistfp+cist > cstmax)
    !        Wall Impact   : particle_status=3   (Cross > Wall_dist(S))
    !        Target impact : particle_status=4   (S<0, S>Smax)

    !
    !     Scan through the fp_s array to find which cell the particle occupies
    !      
    !     Keep in mind that fp_s runs from 0 to 2*nks(ir) with the center s for 
    !     each cell at the coordinate 2*ik-1 and the bounds at +/- 1 from this 
    !
    !     If flag is set to 0 then this code will search the entire fp_s array to find the 
    !     proper index - rc is set to 0 if the cell is found - 1 if S < 0 and 2 if S > Smax
    !
    !     Smax is stored in fp_s(2*fp_cells,fp_reg)
    !
    integer direction,test_ik
    integer ik,fp_reg
    real    s
    !
    !     Initialize return code
    !
    rc = 0 

    !
    ! Set local values from module variable to make the code more readable
    !
    s = part%s
    ik = part%ik
    fp_reg = part%fp_reg
    !
    !     Particle still in current cell - no changes required
    !
    if (s.ge.fp_s(2*ik-2,fp_reg).and.s.le.fp_s(2*ik,fp_reg)) then 

       rc = 0

    else
       !
       ! Particle not in current cell - look for new one
       !
       if (s.lt.fp_s(2*ik-2,fp_reg)) then 
          !
          !        Scan for new cell - scan down
          !
          direction = -1 
       else 
          !
          !        scan up - since it is not in the current cell and is not lower - it must be higher
          !
          direction =  1
       endif
       !
       !        Scan to find cell 
       !
       test_ik = ik + direction
       !
       !        As long as test_ik is in range and a new value has not been assigned
       !
       do while (test_ik.ge.1.and.test_ik.le.fp_cells(fp_reg).and.ik.ne.test_ik) 

          if (s.ge.fp_s(2*test_ik-2,fp_reg).and.s.le.fp_s(2*test_ik,fp_reg)) then 
             ik = test_ik
             exit
          endif

          test_ik = test_ik + direction

       end do
       !
       !        Either < 0 or > SMAX 
       !
       if (ik.ne.test_ik) then

          if (s.le.fp_s(0,fp_reg)) then 
             rc = 4
             ik = 1
             s = 0.0
          elseif (s.ge.fp_s(2*fp_cells(fp_reg),fp_reg)) then 
             rc = 4
             ik = fp_cells(fp_reg)
             s = fp_s(2*fp_cells(fp_reg),fp_reg)
          endif
       else
          !        
          !  Update the periphery plasma and transport conditions in the transport module
          !

          call fp_set_plasma(ik,fp_reg)

       endif


       !
       ! Assign particle values from local variables since particle has changed cells 
       ! S is set because it may have been set to 0.0 or SMAX
       !
       part%s = s
       part%ik = ik

    endif

    if (mod_debug_fp) then 
       write(6,'(a,3i6,12g18.10)') 'FP CHECK  S:',rc,part%ik,test_ik,part%s,part%cross
       if (rc.eq.4) then 
          write(6,'(a,3i6,12g18.10)') 'FP TARGET  :'
       endif
    endif

  end subroutine fp_check_s


  subroutine fp_check_cross(rsect,zsect,rc)
    use mod_fperiph
    implicit none
    real rsect,zsect
    integer rc
    !
    !      include 'params'
    !      include 'fperiph_com'
    !
    !     RC = 0 - no change - particle still in cell
    !     RC = 1 - return to grid
    !     RC = 3 - struck wall 
    !
    !        Status codes set to match those in FPERIPH
    !
    !        Is_Followed   : particle_status=0   (Particle is still being followed) 
    !        Return to grid: particle_status=1   (Cross < 0 - assuming Cross is set to 0 when entering FP) 
    !        Maximum time  : particle_status=2   (cistfp+cist > cstmax)
    !        Wall Impact   : particle_status=3   (Cross > Wall_dist(S))
    !        Target impact : particle_status=4   (S<0, S>Smax)
    !
    !     Determine whether the particle has returned to the grid or struck the wall
    !
    real walldist,frac
    integer direction

    real s,cross
    integer ik,fp_reg

    !
    !     Set default return code
    !
    rc = 0

    !
    ! Assign local variables from module variables to make code more legible
    ! NOTE: this routine does not change any of these values - only calculates whether
    ! the particle is inside the walls or back on the grid
    !
    s      = part%s
    cross  = part%cross
    ik     = part%ik
    fp_reg = part%fp_reg

    if (cross.lt.0.0) then
    !
    !     Cross less than 0.0 - particle returns to grid - later code handles assignment of on grid cross value
    !
       rc = 1
    elseif (cross.gt.min_fp_walldist(ik,fp_reg)) then 
       !
       !     Particle is close to wall - need to check for wall impact
       !     Calculate walldist at specific S value
       !     
       !
       !        Calculate walldist at current S
       !
       !        Check S against cell center value
       !
       if (s.le.fp_s(2*ik-1,fp_reg)) then 
          direction = -1
       else
          direction = 1
       endif
       !
       frac =   (s-fp_s(2*ik-1,fp_reg))  / (fp_s(2*ik-1+direction,fp_reg) - fp_s(2*ik-1,fp_reg))

       walldist = fp_walldist(2*ik-1,fp_reg) + frac * ( fp_walldist(2*ik-1+direction,fp_reg) -fp_walldist(2*ik-1,fp_reg)) 
       !
       !        Check for wall impact
       !
       if (cross.gt.walldist) then

          rc = 3
          !
          !           Calculate approximate R,Z of wall intersection
          !
          rsect = fp_wallcoords(2*ik-1,fp_reg,1) +  &
               & frac * ( fp_wallcoords(2*ik-1+direction,fp_reg,1)-fp_wallcoords(2*ik-1,fp_reg,1)) 

          zsect = fp_wallcoords(2*ik-1,fp_reg,2) +  &
               & frac * ( fp_wallcoords(2*ik-1+direction,fp_reg,2)-fp_wallcoords(2*ik-1,fp_reg,2)) 


       endif

    endif
    !
    !     If cross is within the appropriate region then exit
    !

    if (mod_debug_fp) then 
       write(6,'(a,3i6,12g18.10)') 'FP CHECK  C:',rc,ik,fp_reg,cross,walldist,min_fp_walldist(ik,fp_reg),rsect,zsect
       if (rc.eq.3) then 
          write(6,'(a,2i6,12g18.10)') 'FP WALL    :'
       endif   
    endif

  end subroutine fp_check_cross



end module mod_fp_transport
