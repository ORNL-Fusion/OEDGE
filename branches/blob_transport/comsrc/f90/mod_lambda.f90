module mod_lambda

  !
  ! This module supports redefining the value of lambda used throughout the DIVIMP and LIM codes
  ! There are two inputs - the lambda_opt and the lambda_val for fixed values of lambda
  ! The supplied function returns the lambda value to be used as a function of the option selected and local
  ! density and temperature passed to the routine.
  !
  ! All points in the code using lambda should be modified to consistently use this routine
  ! The default option is 0 which will assign a defaul constant value of 15.0
  !
  ! The default lambda in DIVIMP is option 0 and 15.0
  ! The default lambda in LIM is option 2
  
  integer,public :: lambda_opt=0  ! T47
  real,public :: lambda_val = 15.0  ! T48

  ! LIM currently uses a fixed value of lambda for all cells calculated from the value at the inboard limiter tip.
  ! This isn't very meaningful since there is nothing special about the plasma condition at that location.
  ! However, to maintain backward compatibility a new option is being added where option 0 uses the spatially constant
  ! but possibly varying due to the plasma conditions and option 1 (which is the default in divimp) changes the value
  ! of lambda based on local plasma conditions.
  ! Lambda_vary_opt is not an unstructured input in divimp
  integer,public :: lambda_vary_opt=0
  
  public:: coulomb_lambda,print_lambda_option

contains

  real function coulomb_lambda(n,t)
    implicit none
    real :: n,t

    ! jdemod
    ! Note: Reiser.f uses a separate constant optional input called coulomb_log
    ! This code was not updated because n,t were not always available the way the code was designed and written.
    ! The reiser code only works with constant lambda values without a significant re-write. 
    
    if (n.le.0.0.or.t.le.0.0) then
       ! error conditions on input
       coulomb_lambda = lambda_val
       return
    endif


    if (lambda_opt.eq.0) then
       ! base line constant value of lambda
       coulomb_lambda = lambda_val
    elseif (lambda_opt.eq.1) then 
       ! Used in the GetIonRelTime routine in DIVIMP
       ! called from IonViscosity
       ! Appears to be dead end code added by Steve that isn't
       ! called from anywhere else
       !
       ! Formula used as an option in the hc code  hc_lambda_calc option
       ! Use calculation for Lambda found in Dolan for ion-collision dominated plasma.
       ! Originally in Sivukhin, D.V., Coulomb collisions in a fully ionized plasma in
       ! Review of Plasma Physics (Consultation Bureau, New York, 1966) Vol. 4, p.88.
        coulomb_lambda = 30.0 - 0.5 * LOG(n) + 1.5 * LOG(t)
    elseif (lambda_opt.eq.2) then 
       ! formula for lambda used in LIM - possibly from Wesson? - but that might use 15.2 instead of 17.3
       ! T in this formula is Ti, n is density
       coulomb_lambda = 17.3 - 0.5*LOG(n/1.0E20) + 1.5*LOG(t/1000.0)
    elseif (lambda_opt.eq.3) then 
       ! SOL22 - coulomb_lambda used the electron-ion energy transfer term
       if ((1.5e13 * t**(1.5) / sqrt(n)).le.1.0) then
          coulomb_lambda = lambda_val
       else
          coulomb_lambda = log(1.5e13 * t**(1.5) / sqrt(n))
       endif
    endif

  end function coulomb_lambda

  subroutine print_lambda_option
    implicit none
    
    ! jdemod - print out the lambda options selected to the .dat file

    if (lambda_opt.eq.0) then
       CALL PRR('  COULOMB LOGARITHM OPTION 0: A CONSTANT VALUE OF LAMBDA IS USED =',lambda_val)
    elseif (lambda_opt.eq.1) then
       CALL PRC('  COULOMB LOGARITHM OPTION 1: THE FOLLOWING FORMULA FOR LAMBDA IS USED')
       call prc('                              LAMBDA= 30.0 - 0.5 * LOG(n) + 1.5 * LOG(t)')
       call prc('                              Formula originally by  Sivukhin, D.V.')
       call prc('                              Coulomb collisions in a fully ionized plasma in')
       call prc('                              Review of Plasma Physics (Consultation Bureau, New York, 1966) Vol. 4, p.88.')
       call prc('                              Also from Dolan for ion collision dominated plasma')
    elseif (lambda_opt.eq.2) then 
       CALL PRC('  COULOMB LOGARITHM OPTION 2: LIM BASE OTION: THE FOLLOWING FORMULA FOR LAMBDA IS USED')
       call prc('                              LAMBDA= 17.3 - 0.5*LOG(n/1.0E20) + 1.5*LOG(t/1000.0)')
       call prc('                              IN LIM: LAMBDA is constant at this value using plasma conditions at limiter tip'&
                                          'unless lambda_vary_opt=1')
       call prc('                              IN DIVIMP: Lambda varies with local plasma conditions')        
    elseif (lambda_opt.eq.3) then 
       CALL PRC('  COULOMB LOGARITHM OPTION 3: SOL22 PEI POWER TERM BASE OTION: THE FOLLOWING FORMULA FOR LAMBDA IS USED')
       call prc('                              LAMBDA= log(1.5e13 * t**(1.5) / sqrt(n)) )')
       call prr('                              IF LOG argument is =<1 this uses the constant value lambda =',lambda_val)
    endif
    call prb
    if (lambda_vary_opt.eq.1) then
       call prc(' COULOMB VARIATION OPTION 1: LIM ONLY: LAMBDA values vary spatially based on cell plasma conditions')
    endif
       
    
  end subroutine print_lambda_option


  

end module mod_lambda
