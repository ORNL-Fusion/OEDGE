module mod_promptdep
  use debug_options
  implicit none

  !
  !     common block containing variables related to the prompt
  !     deposition option.
  !
  !     -*-fortran-*-
  ! common /prompt_data/ prompt_depopt,promptdeps,mps_thickness,mps_energy
  !
  ! save /prompt_data/
  !
  real,public,allocatable :: promptdeps(:,:),mps_thickness(:,:),mps_energy(:,:)
  !
  integer,public :: prompt_depopt
  
  ! sazmod - adding generalized prompt redeposition coefficients and
  ! and average charge state near the target option.
  real, public :: prompt_dep_avg_z

  public :: allocate_mod_promptdep,deallocate_mod_promptdep
  
  ! Declaration and variable for the Guterl prompt redepositon scaling.
  public :: w_prob_nonprompt_guterl

contains

  subroutine allocate_mod_promptdep
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_promptdep','ALLOCATE')

    call allocate_array(promptdeps,maxnds,9,'promptdeps',ierr)
    call allocate_array(mps_thickness,maxnrs,2,'mps_thickness',ierr)
    call allocate_array(mps_energy,maxnrs,2,'mps_energy',ierr)

  end subroutine allocate_mod_promptdep


  subroutine deallocate_mod_promptdep
    implicit none

    call pr_trace('mod_promptdep','DEALLOCATE')

    if (allocated(promptdeps)) deallocate(promptdeps)
    if (allocated(mps_thickness)) deallocate(mps_thickness)
    if (allocated(mps_energy)) deallocate(mps_energy)

  end subroutine deallocate_mod_promptdep
  
  real function w_prob_nonprompt_guterl(te, B, crmb, rizb, tau_iz)
    implicit none
    
    ! This implements the prompt redeposition scaling for W from
    ! Guterl NME 2021. The assumptions within the study are made here
    ! as well since the scaling is derived from ERO simulations that
    ! likewise made the same assumptions. This technically puts it at
    ! odds with a few things within DIVIMP, notably:
    !   - Characteristic sheath width assumed to be 5*larmor radius
    !       for deuterium
    !   - Velocity of sputtered W atoms assumed to be sqrt(2*8.68/mz)
    !       despite the actual atom velocity being available within 
    !       DIVIMP.
    ! If the controlling physics could be better elucidated beyond the
    ! lambda_iz_hat parameter used in the paper, then a more detailed
    ! model could be implemented, but nonetheless this is a decent 
    ! start.
  
    real :: te, B, crmb, rizb, tau_iz, lambda_sheath, lambda_iz_hat
    real :: larmor
    real, parameter :: epsilon_sheath = 5.0
    external larmor

    ! lambda_sheath is defined as 5 * the main ion gyroradius in Guterl.
    ! Not completely sure why, but that's what the scaling is built off.
    ! 5.0 has been hardcoded here for now.
    lambda_sheath = epsilon_sheath * larmor(crmb, te, B, rizb)
    
    ! Then, lambda_iz^hat is assumed to be:
    ! vb * tau_iz / lambda_sheath. 
    ! vb is the energy of the sputtered W atom, which Guterl assmes to 
    ! be sqrt(2 * Eb / mz), where Eb is the binding energy (8.68 eV).
    ! tau_iz is the characteristic ionization time for W, calculated
    ! from ADAS (scd50.dat). We could use the actual ion energy for vb, 
    ! but the scaling assumes the above so we stick with it.
    ! Hardcoding vb = sqrt(2*Eb/mz) = 3020.5568 m/s
    lambda_iz_hat = 3020.5568 * tau_iz / lambda_sheath
    
    ! The probability of NOT promptly redepositing is from the scaling:
    !   1 - f_redep = p_nonprompt = exp(-a*x^b)
    !   a = 1.485
    !   b = -0.56 
    w_prob_nonprompt_guterl = exp(-1.485*lambda_iz_hat**(-0.56))
    !write(0,*) 'te, B, crmb, rizb, tau_iz, larmor, lambda_sheath, lambda_iz_hat, prob = ',te, B, crmb, rizb, tau_iz, larmor(crmb, te, B, rizb), lambda_sheath, lambda_iz_hat, w_prob_nonprompt_guterl
    return
    
  
  end function w_prob_nonprompt_guterl

end module mod_promptdep
