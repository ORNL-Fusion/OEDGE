module velocity_dist


  implicit none

  private
  save

  !
  ! This module is used to debug the neutral velocity distributions in 
  ! both DIVIMP and the HC module. It can record both the initial 
  ! distributions as well as the integrated distributions over time - these
  ! could be different if different velocity results give different ionization
  ! or loss times. They could also be different based on lying within the 
  ! viewing cone of a specific line of sight. Particles with high lateral 
  ! velocities will leave the cone faster. 
  !


  !
  !     Add velocity and angular distribution diagnostics
  !     - vin
  !     - anglan
  !     - tanlan
  !     - anglan+tanlan
  !     - psi, beta - where appropriate
  !
  integer,parameter :: maxangs = 360

  real    :: vin_scale,eiv
  real,allocatable :: vin_dist(:),psi_dist(:),beta_dist(:),&
       &  anglan_dist(:),tanlan_dist(:),angtot_dist(:)  
  real,allocatable :: vlp_dist(:)
  real,allocatable :: vtor_dist(:),vpol_dist(:),vr_dist(:),vz_dist(:)

  integer :: debug_neutv_mod
  integer :: debug_neutv_nbins
  real    :: debug_neutv_einmax


  public :: init_velocity_dist,record_vdist,print_vdist,record_vlp_dist,debug_neutv_mod,record_c_detailed_vdist


contains

  subroutine init_velocity_dist(debug_neutv_in,debug_neutv_einmax_in,debug_neutv_nbins_in,crmi)
    implicit none
    integer :: debug_neutv_in,debug_neutv_nbins_in
    real    :: debug_neutv_einmax_in,crmi

    !
    ! Initialize the options which control execution - these start off as unstructured input variables with the same name
    !

    debug_neutv_mod = debug_neutv_in
    debug_neutv_einmax = debug_neutv_einmax_in
    debug_neutv_nbins = debug_neutv_nbins_in

    !
    ! If the option is turned on then allocate the data storage
    !
! slmod begin - bug
! This offends the Sun compiler since DEBUG_NEUTV_MOD (and DEBUG_HEUTV_IN) are 
! declared INTEGER.  -SL 27.03.07
    if (debug_neutv_mod.ne.0) then 
!
!    if (debug_neutv_mod) then 
! slmod end
       !
       ! Assign bin scaling
       !
       vin_scale = 1.38E4 &
       &          * sqrt( debug_neutv_einmax / crmi) &
       &          / debug_neutv_nbins


       !
       ! allocate storage
       !
       call allocate_vdist_data
    endif

  end subroutine init_velocity_dist


  subroutine allocate_vdist_data
    implicit none
    integer :: flag


    allocate(vin_dist(0:debug_neutv_nbins),stat=flag)
    allocate(vlp_dist(-debug_neutv_nbins:debug_neutv_nbins))
    allocate(psi_dist(-maxangs:maxangs),stat=flag)
    allocate(beta_dist(-maxangs:maxangs),stat=flag)
    allocate(anglan_dist(-maxangs:maxangs),stat=flag)
    allocate(tanlan_dist(-maxangs:maxangs),stat=flag)
    allocate(angtot_dist(-maxangs:maxangs),stat=flag)

    !
    ! Detailed
    !
    allocate(vpol_dist(0:debug_neutv_nbins),stat=flag)
    allocate(vtor_dist(-debug_neutv_nbins:debug_neutv_nbins))
    allocate(vr_dist(-debug_neutv_nbins:debug_neutv_nbins))
    allocate(vz_dist(-debug_neutv_nbins:debug_neutv_nbins))



    vin_dist = 0.0
    psi_dist = 0.0
    beta_dist = 0.0
    anglan_dist = 0.0
    tanlan_dist = 0.0
    angtot_dist = 0.0

    vlp_dist = 0.0
    
    vpol_dist = 0.0
    vtor_dist = 0.0
    vr_dist = 0.0
    vz_dist = 0.0


  end subroutine allocate_vdist_data



  subroutine record_vdist(vin,beta,psi,tanlan,anglan,sputy)
    use mod_params
    implicit none
    real :: vin,beta,psi,tanlan,anglan,sputy

    integer :: iv,ia


    !
    ! If neutral velocity debugging is off then return or if the arrays have not been allocated
    ! The arrays are only allocated in neutbatch and are deallocated after the results have been
    ! printed - they are only recorded for one launch at a time and are triggered from neutbatch. 
    ! However, the launch routine can be called from div to launch self-sputtered particles where
    ! the debugging of the initial neutral velocity distribution is not yet supported - to add
    ! this support - a call to the init_velocity_dist routine would need to be added before the 
    ! call to LAUNCH in the div.f module and a call to the print_vdist routine added afterwards. 
    !
    ! In any case it is safer if this code checks to make sure that the data arrays have been 
    ! allocated before trying to write to them. 
    !

    if (debug_neutv_mod.eq.0.or.(.not.allocated(vin_dist))) return


    !
    !       Record detailed velocity debugging information if the 
    !       option is active
    !
    !
    !          Velocity 
    !
    iv = min(int(vin/vin_scale),debug_neutv_nbins)
! slmod begin
!   Not sure what is happening here, but it is a very rare 
!   event.  First observed when I started running a supplimental 
!   neutral launch (i-fwp-0120j) - SL, 01/12/2009
    if (iv.lt.0) then
      write(0,*) 'WARNING record_vdist: IV.LT.0'
      write(6,*) 'WARNING record_vdist: IV.LT.0'
      return
    endif
! slmod end
    vin_dist(iv) = vin_dist(iv) + sputy
    !
    !          Angles
    !
    ia = int(beta*raddeg) 
    beta_dist(ia) = beta_dist(ia) + sputy

    ia = int(psi*raddeg) 
    psi_dist(ia) = psi_dist(ia) + sputy

    ia = int(tanlan*raddeg) 
    tanlan_dist(ia) = tanlan_dist(ia) + sputy

    ia = int(anglan*raddeg) 
    anglan_dist(ia) = anglan_dist(ia) + sputy

    ia = int((anglan+tanlan)*raddeg) 
    angtot_dist(ia) = angtot_dist(ia) + sputy



  end subroutine record_vdist





  subroutine record_vlp_dist(vlp,sputy)
    use mod_params
    implicit none
    real :: vlp,sputy

    integer :: iv


    !
    ! If neutral velocity debugging is off or array not allocated then return
    !

    if (debug_neutv_mod.eq.0.or.(.not.allocated(vlp_dist))) return


    !
    !       Record detailed velocity debugging information if the 
    !       option is active
    !
    !
    !          Velocity 
    !
    iv = min(int(abs(vlp)/vin_scale),debug_neutv_nbins) * sign(1.0,vlp)
    vlp_dist(iv) = vlp_dist(iv) + sputy
    !


  end subroutine record_vlp_dist


  subroutine record_c_detailed_vdist(hc_v,sputy)
    use hc_velocity_type
    implicit none
    type (hc_velocity_type1) :: hc_v
    real :: sputy

    !
    ! local variables
    !
    integer :: iv
    real :: vtor,vpol,vr,vz

    !
    ! If neutral velocity debugging is off then return
    !

    if (debug_neutv_mod.eq.0) return


    vr = hc_v%v(1)
    vz = hc_v%v(2)
    vtor = hc_v%v(3)

    vpol = sqrt(hc_v%v(1)**2 + hc_v%v(2)**2)

       !
       !          Velocity 

    ! toroidal
       !
       iv = min(int(abs(vtor)/vin_scale),debug_neutv_nbins) * sign(1.0,vtor)
       vtor_dist(iv) = vtor_dist(iv) + sputy
       !

       !
       ! poloidal
       !
       iv = min(int(abs(vpol)/vin_scale),debug_neutv_nbins)
       vpol_dist(iv) = vpol_dist(iv) + sputy

       !
       ! Vr
       !
       iv = min(int(abs(vr)/vin_scale),debug_neutv_nbins) * sign(1.0,vr)
       vr_dist(iv) = vr_dist(iv) + sputy

       !
       ! Vz
       !
       iv = min(int(abs(vz)/vin_scale),debug_neutv_nbins) * sign(1.0,vz)
       vz_dist(iv) = vz_dist(iv) + sputy


    
  end subroutine record_c_detailed_vdist




  subroutine print_vdist(crmi)
    implicit none

    real :: crmi

    !
    ! Local variables
    !

    real :: eiv
    integer :: iv,ia

    !
    !     Print out Neutral velocity debug information and 
    !     deallocate the arrays.
    !
    if (debug_neutv_mod.gt.0) then 
       write(6,*) 
       write(6,'(a,f20.5)') 'NEUT: DEBUG VELOCITY: VIN and VPOL',vin_scale
       write(6,*) 
       do iv = 0,debug_neutv_nbins

          eiv = crmi * ((iv*vin_scale+vin_scale/2.0)/1.38e4)**2

          write(6,'(i6,2(1x,g18.6),2(1x,f18.4))') iv, &
          &             iv*vin_scale+vin_scale/2.0,eiv, &
          &             vin_dist(iv),vpol_dist(iv)
       end do
       write(6,*) 
       write(6,*) 'NEUT: DEBUG VELOCITY: ANGLES'
       write(6,*) 
       do ia = -maxangs,maxangs
          write(6,'(i6,1x,f8.1,1x,6(1x,f18.3))') &
          &             ia,ia+0.5,&
          &             beta_dist(ia),psi_dist(ia),&
          &             anglan_dist(ia),tanlan_dist(ia),&
          &             angtot_dist(ia)
       end do
       
       write(6,*) 
       write(6,'(a)') 'LINE PROFILE VELOCITY DISTRIBUTION - LINE PROFILE OPTION MUST BE ACTIVE'
       write(6,*) 
       write(6,'(a,f20.5)') 'NEUT: DEBUG VELOCITY: VLP',vin_scale
       write(6,*) 
       do iv = -debug_neutv_nbins,debug_neutv_nbins

          eiv = crmi * ((iv*vin_scale+vin_scale/2.0)/1.38e4)**2

          write(6,'(i6,2(1x,g18.6),1x,f18.4)') iv, &
          &             iv*vin_scale+vin_scale/2.0,eiv, &
          &             vlp_dist(iv)
       end do

       write(6,*) 
       write(6,'(a)') 'DETAILED INITIAL VELOCITY DISTRIBUTIONS - VR, VZ, VT'
       write(6,*) 
       write(6,'(a,f20.5)') 'NEUT: DEBUG VELOCITY: DETAILED',vin_scale
       write(6,*) 
       do iv = -debug_neutv_nbins,debug_neutv_nbins

          eiv = crmi * ((iv*vin_scale+vin_scale/2.0)/1.38e4)**2

          write(6,'(i6,2(1x,g18.6),3(1x,f18.4))') iv, &
          &             iv*vin_scale+vin_scale/2.0,eiv, &
          &             vr_dist(iv),vz_dist(iv),vtor_dist(iv)
       end do


       !
       !        Deallocate arrays 
       !
       call deallocate_vdist_data

    endif



  end subroutine print_vdist

  subroutine deallocate_vdist_data
    implicit none

    !
    !        Deallocate arrays 
    !

    deallocate(vin_dist)
    deallocate(vlp_dist)
    deallocate(beta_dist)
    deallocate(psi_dist)
    deallocate(anglan_dist)
    deallocate(tanlan_dist)
    deallocate(angtot_dist)

    deallocate(vpol_dist)
    deallocate(vtor_dist)
    deallocate(vr_dist)
    deallocate(vz_dist)


  end subroutine deallocate_vdist_data






end module velocity_dist
